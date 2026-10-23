"""Selection equivalence for projected groups and batched flag assignment."""

import numpy
import pandas
import pytest

from amalgkit.metadata_utils import Metadata, SELECT_SAMPLING_STRATEGIES


# Captured from 2dbd568 with NumPy 2.5.3 / pandas 3.0.6; these include first-seen groups,
# missing species/group/project labels, sparse exclusions and seeded tie breaks.
EXPECTED = {
    "largest_bioprojects_first:0": ["R4", "R5", "R14", "R16", "R20", "R29", "R33", "R34", "R35", "R42", "R44", "R48"],
    "largest_bioprojects_first:17": ["R4", "R6", "R7", "R8", "R12", "R20", "R28", "R32", "R35", "R40", "R44", "R48"],
    "largest_bioprojects_first:21": ["R0", "R6", "R7", "R12", "R16", "R20", "R28", "R32", "R36", "R40", "R44", "R48"],
    "maximize_bioproject_diversity:0": ["R2", "R4", "R5", "R8", "R20", "R21", "R27", "R34", "R35", "R39", "R42", "R47"],
    "maximize_bioproject_diversity:17": [
        "R4",
        "R6",
        "R12",
        "R19",
        "R21",
        "R23",
        "R28",
        "R29",
        "R35",
        "R37",
        "R41",
        "R48",
    ],
    "maximize_bioproject_diversity:21": [
        "R0",
        "R6",
        "R7",
        "R12",
        "R18",
        "R19",
        "R20",
        "R21",
        "R25",
        "R27",
        "R29",
        "R44",
    ],
    "random:0": ["R8", "R9", "R12", "R14", "R21", "R28", "R34", "R36", "R37", "R41", "R46", "R48"],
    "random:17": ["R0", "R4", "R6", "R7", "R9", "R11", "R12", "R15", "R19", "R20", "R28", "R34"],
    "random:21": ["R0", "R6", "R7", "R13", "R15", "R20", "R22", "R25", "R28", "R37", "R40", "R46"],
    "smallest_bioprojects_first:0": ["R5", "R11", "R14", "R20", "R21", "R22", "R23", "R26", "R27", "R39", "R42", "R48"],
    "smallest_bioprojects_first:17": ["R7", "R9", "R12", "R20", "R21", "R22", "R25", "R26", "R27", "R35", "R37", "R48"],
    "smallest_bioprojects_first:21": ["R0", "R9", "R12", "R20", "R21", "R22", "R25", "R26", "R27", "R28", "R37", "R48"],
}


def sampling_frame():
    ids = numpy.arange(49)
    return pandas.DataFrame(
        {
            "run": ["R" + str(i) for i in ids],
            "scientific_name": numpy.array(["zeta", "alpha", None, "zeta", None, "alpha", "beta"], dtype=object)[
                ids % 7
            ],
            "sample_group": numpy.array(["leaf", "root", None, "leaf", "", "root", " "], dtype=object)[ids % 7],
            "bioproject": numpy.array(["B", "A", None, "C"], dtype=object)[ids % 4],
            "exclusion": numpy.array([" No ", "no", None, "low", "NO", "", "no"], dtype=object)[ids % 7],
            "annotation": ["note" + str(i) for i in ids],
            "counts": pandas.array(ids, dtype="Int64"),
        }
    )


@pytest.mark.parametrize("strategy", SELECT_SAMPLING_STRATEGIES)
@pytest.mark.parametrize("seed", [0, 17, 21])
@pytest.mark.parametrize("index_kind", ["range", "integer", "text"])
def test_sampling_full_metadata_and_seeded_selection(strategy, seed, index_kind):
    metadata = Metadata.from_dataframe(sampling_frame())
    if index_kind == "integer":
        metadata.df.index = numpy.arange(49)[::-1] * 3
    elif index_kind == "text":
        metadata.df.index = ["label" + str(i) for i in range(49)]
    expected = metadata.df.copy()
    expected["sample_group"] = expected["sample_group"].fillna("").astype(str)
    expected["exclusion"] = expected["exclusion"].fillna("no").astype(str).str.strip().str.lower()
    expected["bioproject"] = expected["bioproject"].fillna("unknown").astype(str).values
    selected = EXPECTED[strategy + ":" + str(seed)]
    expected["is_sampled"] = numpy.where(expected["run"].isin(selected), "yes", "no")
    expected["is_qualified"] = numpy.where(
        (expected["exclusion"] == "no") & (expected["sample_group"] != ""), "yes", "no"
    )
    expected["sampling_strategy"] = strategy
    expected["sampling_seed"] = seed
    expected = expected.reset_index(drop=True)

    metadata.label_sampled_data(max_sample=3, sampling_strategy=strategy, random_seed=seed)

    pandas.testing.assert_frame_equal(metadata.df, expected)


@pytest.mark.parametrize("strategy", SELECT_SAMPLING_STRATEGIES)
def test_sampling_empty_table_preserves_schema_and_option(strategy):
    metadata = Metadata.from_dataframe(sampling_frame().iloc[:0])
    with pandas.option_context("mode.chained_assignment", "raise"):
        metadata.label_sampled_data(sampling_strategy=strategy)
        assert pandas.get_option("mode.chained_assignment") == "raise"
    assert metadata.df.empty
    assert metadata.df["sampling_strategy"].dtype == pandas.Series([strategy]).dtype
    assert metadata.df["sampling_seed"].dtype == numpy.dtype("int64")


def test_sampling_completed_groups_retained_after_later_failure(monkeypatch):
    metadata = Metadata.from_dataframe(
        pandas.DataFrame(
            {
                "scientific_name": ["Species A", "Species B"],
                "sample_group": ["leaf", "leaf"],
                "bioproject": ["P", "P"],
                "run": ["R0", "R1"],
                "exclusion": ["no", "no"],
            }
        )
    )
    original = metadata._sample_random
    calls = 0

    def failing_sampler(df, **kwargs):
        nonlocal calls
        calls += 1
        if calls == 2:
            raise RuntimeError("second group")
        return original(df, **kwargs)

    monkeypatch.setattr(metadata, "_sample_random", failing_sampler)
    with pandas.option_context("mode.chained_assignment", "raise"):
        with pytest.raises(RuntimeError, match="second group"):
            metadata.label_sampled_data(sampling_strategy="random")
        assert pandas.get_option("mode.chained_assignment") == "raise"
    assert metadata.df["is_sampled"].tolist() == ["yes", "no"]


@pytest.mark.parametrize("strategy", SELECT_SAMPLING_STRATEGIES)
def test_sampling_duplicate_index_preserves_legacy_assignment_error(strategy):
    metadata = Metadata.from_dataframe(sampling_frame().iloc[:2])
    metadata.df["scientific_name"] = "Species A"
    metadata.df["sample_group"] = "leaf"
    metadata.df.index = [7, 7]
    with pandas.option_context("mode.chained_assignment", "raise"):
        with pytest.raises(ValueError):
            metadata.label_sampled_data(sampling_strategy=strategy)
        assert pandas.get_option("mode.chained_assignment") == "raise"


@pytest.mark.parametrize("strategy", SELECT_SAMPLING_STRATEGIES)
@pytest.mark.parametrize("seed", [-1, None, "invalid"])
def test_sampling_invalid_seed_leaves_metadata_and_option_unchanged(strategy, seed):
    metadata = Metadata.from_dataframe(sampling_frame())
    expected = metadata.df.copy()
    with pandas.option_context("mode.chained_assignment", "raise"):
        with pytest.raises(ValueError, match="non-negative integer"):
            metadata.label_sampled_data(sampling_strategy=strategy, random_seed=seed)
        assert pandas.get_option("mode.chained_assignment") == "raise"
    pandas.testing.assert_frame_equal(metadata.df, expected)
