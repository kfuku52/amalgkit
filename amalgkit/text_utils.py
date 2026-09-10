"""Normalization helpers shared across metadata and workflow modules."""

from __future__ import annotations

from collections.abc import Iterable

import pandas


def normalize_unique_text(values: Iterable[object]) -> list[str]:
    """Return non-empty text values in first-seen order without duplicates."""
    normalized_values = []
    for value in values:
        if value is None or pandas.isna(value):
            continue
        normalized = str(value).strip()
        if normalized == "":
            continue
        normalized_values.append(normalized)
    return list(dict.fromkeys(normalized_values))


def parse_sample_group_argument(sample_group_arg: object) -> list[str]:
    tokens = []
    token = []
    value = str(sample_group_arg)
    index = 0
    while index < len(value):
        char = value[index]
        if char == "\\" and index + 1 < len(value) and value[index + 1] in {"\\", ",", "|"}:
            token.append(value[index + 1])
            index += 2
            continue
        if char in {",", "|"}:
            tokens.append("".join(token))
            token = []
        else:
            token.append(char)
        index += 1
    tokens.append("".join(token))
    return normalize_unique_text(tokens)


def serialize_sample_groups(values: Iterable[object]) -> str:
    groups = normalize_unique_text(values)
    escaped = [value.replace("\\", "\\\\").replace(",", "\\,").replace("|", "\\|") for value in groups]
    return "|".join(escaped)
