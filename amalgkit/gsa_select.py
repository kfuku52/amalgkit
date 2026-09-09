"""Defer size-dependent selection until public GSA FASTQ sizes are measured."""

import json

import pandas

from amalgkit.gsa_fastq import is_gsa_row

DEFERRED_COLUMN = "gsa_deferred_select"


def reset_deferred_selection(metadata, kind=None):
    """A new selection replaces earlier decisions, including removed rules."""
    df = metadata.df
    if DEFERRED_COLUMN not in df:
        return metadata
    for index, value in df[DEFERRED_COLUMN].items():
        recipes = _recipes(value)
        remaining = [recipe for recipe in recipes if kind is not None and recipe.get("kind") != kind]
        df.at[index, DEFERRED_COLUMN] = json.dumps(remaining, sort_keys=True)
        if "gsa_selection_status" in df:
            df.at[index, "gsa_selection_status"] = "pending_input_counts" if remaining else ""
    return metadata


def validate_selection_ready(df):
    """Explicit and inferred downstream inputs have the same eligibility rules."""
    eligible = pandas.Series(True, index=df.index)
    for column, expected in (("exclusion", "no"), ("is_sampled", "yes")):
        if column in df:
            values = df[column].fillna("").astype(str).str.strip().str.lower()
            if values.ne("").any():
                eligible &= values.eq(expected)
    pending = pandas.Series(False, index=df.index)
    if "gsa_selection_status" in df:
        pending |= df["gsa_selection_status"].fillna("").str.strip().eq("pending_input_counts")
    if DEFERRED_COLUMN in df:
        pending |= df[DEFERRED_COLUMN].map(lambda value: bool(_recipes(value)))
    if (eligible & pending).any():
        raise ValueError("GSA selection is pending input counts; finish getfastq for all selected runs first.")


def _recipes(value):
    if not isinstance(value, str) or not value.strip():
        return []
    recipes = json.loads(value)
    if not isinstance(recipes, list) or any(not isinstance(recipe, dict) for recipe in recipes):
        raise ValueError("Invalid deferred GSA selection rules")
    return recipes


def _append_recipe(df, indices, recipe):
    if len(indices) == 0:
        return
    if DEFERRED_COLUMN not in df:
        df[DEFERRED_COLUMN] = ""
    if "gsa_selection_status" not in df:
        df["gsa_selection_status"] = ""
    for index in indices:
        recipes = _recipes(df.at[index, DEFERRED_COLUMN])
        if recipe not in recipes:
            recipes.append(recipe)
        df.at[index, DEFERRED_COLUMN] = json.dumps(recipes, sort_keys=True)
        df.at[index, "gsa_selection_status"] = "pending_input_counts"


def defer_unknown_count_filter(df, column, threshold, target, outcome, protected):
    unknown = pandas.Series(False, index=df.index)
    if column != "total_spots" or "data_source" not in df:
        return unknown
    numeric = pandas.to_numeric(df[column], errors="coerce")
    unknown = df["data_source"].fillna("").str.lower().eq("gsa") & numeric.isna() & ~protected
    _append_recipe(
        df,
        df.index[unknown],
        {
            "kind": "min_spots",
            "threshold": threshold,
            "target": target,
            "outcome": outcome,
        },
    )
    return unknown


def defer_unknown_count_dedup(df, rule, eligible):
    if rule["target_column"] != "total_spots" or "data_source" not in df:
        return eligible
    numeric = pandas.to_numeric(df["total_spots"], errors="coerce")
    unknown = df["data_source"].fillna("").str.lower().eq("gsa") & numeric.isna() & eligible
    if not unknown.any():
        return eligible
    keys = df[rule["columns"]].fillna("").astype(str).apply(tuple, axis=1)
    pending = eligible & keys.isin(set(keys[unknown]))
    _append_recipe(
        df,
        df.index[pending],
        {
            "kind": "dedup",
            "columns": rule["columns"],
            "outcome": rule["outcome"],
        },
    )
    return eligible & ~pending


def resolve_deferred_selection(metadata):
    """Resolve known counts; leave incomplete array-job groups explicitly pending."""
    df = metadata.df
    if DEFERRED_COLUMN not in df:
        return metadata
    remaining = {index: _recipes(row.get(DEFERRED_COLUMN)) for index, row in df.iterrows()}
    dedups = []
    for index, row in df.iterrows():
        for recipe in list(remaining[index]):
            if recipe.get("kind") == "min_spots":
                if not is_gsa_row(row):
                    raise ValueError("Deferred GSA threshold requires a GSA input")
                if row.get("read_count_status") != "measured":
                    continue
                threshold = recipe.get("threshold")
                if type(threshold) is not int or threshold < 0:
                    raise ValueError("Invalid deferred GSA minimum spots")
                target, outcome = recipe.get("target"), recipe.get("outcome")
                if target not in df or not isinstance(outcome, str):
                    raise ValueError("Invalid deferred GSA threshold outcome")
                if float(row["total_spots"]) < threshold or float(row["total_spots"]) <= 0:
                    df.at[index, target] = outcome
                remaining[index].remove(recipe)
            elif recipe.get("kind") == "dedup":
                if recipe not in dedups:
                    dedups.append(recipe)
            else:
                raise ValueError("Unknown deferred GSA selection rule")
    for recipe in dedups:
        columns = recipe.get("columns")
        if not isinstance(columns, list) or not columns or any(column not in df for column in columns):
            raise ValueError("Invalid deferred GSA deduplication columns")
        if not isinstance(recipe.get("outcome"), str):
            raise ValueError("Invalid deferred GSA deduplication outcome")
        matching = pandas.Series({index: recipe in recipes for index, recipes in remaining.items()})
        selected = df["is_sampled"].ne("no") if "is_sampled" in df else pandas.Series(True, index=df.index)
        eligible = matching & selected & df["exclusion"].eq("no")
        for _, group in df.loc[eligible].groupby(columns, sort=False, dropna=False):
            candidates = group.copy()
            candidates["total_spots"] = pandas.to_numeric(candidates["total_spots"], errors="coerce")
            if candidates["total_spots"].isna().any():
                continue
            candidates = candidates.sort_values(["total_spots", "run"], ascending=[False, True], kind="stable")
            df.loc[candidates.index[1:], "exclusion"] = recipe["outcome"]
            for index in candidates.index:
                remaining[index].remove(recipe)
        for index in df.index[matching & ~eligible]:
            remaining[index].remove(recipe)
    for index, recipes in remaining.items():
        if _recipes(df.at[index, DEFERRED_COLUMN]):
            df.at[index, DEFERRED_COLUMN] = json.dumps(recipes, sort_keys=True)
            df.at[index, "gsa_selection_status"] = "pending_input_counts" if recipes else "resolved"
    # Keep is_sampled stable so --batch indices still identify the same runs.
    # Downstream commands also honor exclusion, including newly resolved rules.
    return metadata
