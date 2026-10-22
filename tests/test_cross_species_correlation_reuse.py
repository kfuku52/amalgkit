"""Contracts for sharing finite Pearson pairs and writing group scores in batches."""

import numpy
import pandas
import pytest

from amalgkit.cross_species_computation import calculate_correlation_within_group


def _inputs():
    metadata = pandas.DataFrame(
        {
            "species_tag": ["S", "S", "T", "T"],
            "sample_group": ["A", "B", "A", "B"],
            "run": ["a", "b", "a", "b"],
            "exclusion": "no",
        }
    )
    matrix = pandas.DataFrame(
        {
            "S_a": [1.0, 2.0, 4.0, 8.0, 16.0],
            "S_b": [16.0, 9.0, 5.0, 3.0, 2.0],
            "T_a": [2.0, 1.0, 3.0, 7.0, 10.0],
            "T_b": [4.0, 9.0, 1.0, 3.0, 6.0],
        }
    )
    return metadata, matrix


@pytest.mark.parametrize("reference_exclusion", ["run", "species"])
@pytest.mark.parametrize("scale", [1e-200, 1.0, 1e200])
def test_group_scores_keep_scale_invariance_without_mutating_inputs(reference_exclusion, scale):
    metadata, unscaled = _inputs()
    matrix = unscaled * scale
    original_matrix = matrix.copy(deep=True)
    original_metadata = metadata.copy(deep=True)

    actual = calculate_correlation_within_group(metadata, matrix, "test", reference_exclusion)

    within_reference = unscaled["T_a"]
    other_reference = unscaled[["S_b", "T_b"]].median(axis=1) if reference_exclusion == "run" else unscaled["T_b"]
    assert actual.loc[0, "within_group_cor_test"] == pytest.approx(
        numpy.corrcoef(unscaled["S_a"], within_reference)[0, 1], abs=1e-14
    )
    assert actual.loc[0, "max_nongroup_cor_test"] == pytest.approx(
        numpy.corrcoef(unscaled["S_a"], other_reference)[0, 1], abs=1e-14
    )
    assert actual.loc[0, "within_common_genes_test"] == 5
    pandas.testing.assert_frame_equal(matrix, original_matrix)
    pandas.testing.assert_frame_equal(metadata, original_metadata)


@pytest.mark.parametrize("token", ["bad", "1_000"])
def test_text_conversion_keeps_pearson_pairs_separate_from_count_support(token):
    metadata, matrix = _inputs()
    matrix = matrix.astype(object)
    matrix.at[0, "S_a"] = token

    actual = calculate_correlation_within_group(metadata, matrix, "test", "species")

    # NumPy accepts underscores in float strings; pandas.to_numeric does not.
    # The existing scalar Pearson API and support-count policy retain both meanings.
    left = [2.0, 4.0, 8.0, 16.0] if token == "bad" else [1000.0, 2.0, 4.0, 8.0, 16.0]
    right = [1.0, 3.0, 7.0, 10.0] if token == "bad" else [2.0, 1.0, 3.0, 7.0, 10.0]
    assert actual.loc[0, "within_group_cor_test"] == pytest.approx(numpy.corrcoef(left, right)[0, 1], abs=1e-14)
    assert actual.loc[0, "within_common_genes_test"] == 4
    unsupported = calculate_correlation_within_group(metadata, matrix, "test", "species", min_common_genes=5)
    assert numpy.isnan(unsupported.loc[0, "within_group_cor_test"])


@pytest.mark.parametrize("reference_exclusion", ["run", "species"])
def test_duplicate_metadata_labels_keep_last_scored_row_even_for_excluded_rows(reference_exclusion):
    metadata, matrix = _inputs()
    metadata.loc[4] = ["Z", "A", "absent", "old"]
    metadata["within_group_cor"] = 0.123
    metadata["max_nongroup_cor"] = 0.234
    metadata["cs_within_common_genes"] = 17
    metadata["cs_min_nongroup_common_genes"] = 19
    unique = calculate_correlation_within_group(metadata, matrix, "test", reference_exclusion)
    metadata.index = pandas.Index([7, 3, 3, 9, 3], name="row")

    actual = calculate_correlation_within_group(metadata, matrix, "test", reference_exclusion)

    columns = [
        "within_group_cor_test",
        "max_nongroup_cor_test",
        "within_common_genes_test",
        "min_nongroup_common_genes_test",
    ]
    expected = unique.loc[[0, 2, 2, 3, 2], columns].copy()
    expected.index = metadata.index
    pandas.testing.assert_frame_equal(actual[columns], expected, check_exact=True)
    assert actual.iloc[-1]["exclusion"] == "old"


def test_group_scores_do_not_reuse_arrays_after_input_changes():
    metadata, matrix = _inputs()
    initial = calculate_correlation_within_group(metadata, matrix, "test", "species")
    matrix["S_a"] = matrix["T_a"].to_numpy()

    changed = calculate_correlation_within_group(metadata, matrix, "test", "species")

    assert changed.loc[0, "within_group_cor_test"] == pytest.approx(1.0)
    assert changed.loc[0, "within_group_cor_test"] != initial.loc[0, "within_group_cor_test"]
