"""TSV readers that keep identifiers lexical and measurements numeric."""

from __future__ import annotations

import os
from collections.abc import Sequence
from typing import Any

import pandas


def read_identifier_tsv(
    path: str | os.PathLike[str],
    *,
    identifier_columns: Sequence[str | int] = (0,),
    index_col: str | int | None = None,
    **kwargs: Any,
) -> pandas.DataFrame:
    """Preserve leading zeroes and literal NA tokens only in identifier columns.

    Set the index *after* parsing: pandas otherwise applies its missing-value
    inference to index columns even when a converter was supplied. Numeric
    columns retain the usual numeric/NA parsing policy.
    """
    frame = pandas.read_csv(path, sep="\t", converters=dict.fromkeys(identifier_columns, str), **kwargs)
    if index_col is not None:
        label = frame.columns[index_col] if isinstance(index_col, int) else index_col
        frame = frame.set_index(label)
    return frame


def read_annotation_tsv(path: str | os.PathLike[str], **kwargs: Any) -> pandas.DataFrame:
    """Read lexical annotations; only empty cells, not literal NA IDs, are missing."""
    return pandas.read_csv(path, sep="\t", dtype=str, keep_default_na=False, na_values=[""], **kwargs)
