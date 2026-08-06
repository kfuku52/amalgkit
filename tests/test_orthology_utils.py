import numpy
import pandas
import pytest

from amalgkit.orthology_utils import (
    DEFAULT_SINGLE_COPY_THRESHOLD,
    get_single_copy_orthogroup_mask,
    validate_single_copy_threshold,
)


def test_single_copy_orthogroup_mask_uses_percentage_of_species():
    genecount = pandas.DataFrame(
        {
            'A': [1, 1, 1, 2],
            'B': [1, 1, 1, 2],
            'C': [2, 1, 1, 2],
            'D': [0, 2, 1, 2],
        },
        index=['OG50', 'OG75', 'OG100', 'OG0'],
    )

    default_mask = get_single_copy_orthogroup_mask(genecount, ['A', 'B', 'C', 'D'])
    strict_mask = get_single_copy_orthogroup_mask(genecount, ['A', 'B', 'C', 'D'], 100.0)

    assert DEFAULT_SINGLE_COPY_THRESHOLD == 50.0
    assert default_mask.to_dict() == {'OG50': True, 'OG75': True, 'OG100': True, 'OG0': False}
    assert strict_mask.to_dict() == {'OG50': False, 'OG75': False, 'OG100': True, 'OG0': False}


@pytest.mark.parametrize('value', [0, -1, 100.1, numpy.nan, numpy.inf, 'invalid'])
def test_single_copy_threshold_validation_rejects_invalid_values(value):
    with pytest.raises(ValueError, match='greater than 0 and at most 100'):
        validate_single_copy_threshold(value)


def test_single_copy_orthogroup_mask_requires_species_columns():
    genecount = pandas.DataFrame({'A': [1]}, index=['OG1'])

    with pytest.raises(ValueError, match='At least one species'):
        get_single_copy_orthogroup_mask(genecount, [])
    with pytest.raises(ValueError, match='Species column.*B'):
        get_single_copy_orthogroup_mask(genecount, ['A', 'B'])
