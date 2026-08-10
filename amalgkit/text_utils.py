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
