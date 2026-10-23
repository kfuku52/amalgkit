import numpy
import pandas
import pytest

from amalgkit.orthology_utils import (
    DEFAULT_SINGLE_COPY_THRESHOLD,
    generate_multisp_busco_table,
    orthogroup2genecount,
    get_single_copy_orthogroup_mask,
    validate_single_copy_threshold,
)


def test_busco_merge_preserves_quoted_annotation(tmp_path):
    busco_dir = tmp_path / 'busco'
    busco_dir.mkdir()
    description = 'Protein with "quoted" annotation'
    (busco_dir / 'Species_A.tsv').write_text(
        'OG1\tComplete\tgene1\t100\t200\thttp://odb\t' + description + '\n'
    )
    output = tmp_path / 'merged.tsv'

    generate_multisp_busco_table(str(busco_dir), str(output))

    result = pandas.read_csv(output, sep='\t')
    assert result.loc[0, 'description'] == description
    assert result.loc[0, 'Species_A'] == 'gene1'


def test_busco_merge_preserves_full_species_tokens_for_genecount(tmp_path):
    busco_dir = tmp_path / 'busco'
    busco_dir.mkdir()
    row = 'OG1\tComplete\tgene1\t100\t200\thttp://odb\tgene desc\n'
    (busco_dir / 'Homo_sapiens_busco.tsv').write_text(row)
    (busco_dir / 'Homo_sapiens_subsp_busco.tsv').write_text(row)
    orthogroup_path = tmp_path / 'orthogroups.tsv'
    genecount_path = tmp_path / 'genecount.tsv'

    generate_multisp_busco_table(str(busco_dir), str(orthogroup_path))
    merged = pandas.read_csv(orthogroup_path, sep='\t')
    assert {'Homo_sapiens', 'Homo_sapiens_subsp'} <= set(merged.columns)
    orthogroup2genecount(
        str(orthogroup_path), str(genecount_path), ['Homo_sapiens', 'Homo_sapiens_subsp'],
    )
    counts = pandas.read_csv(genecount_path, sep='\t')
    assert counts['Homo_sapiens'].tolist() == [1]
    assert counts['Homo_sapiens_subsp'].tolist() == [1]


def test_orthogroup_gene_counts_preserve_tokens_and_input_species_order(tmp_path):
    orthogroups = tmp_path / 'orthogroups.tsv'
    counts = tmp_path / 'counts.tsv'
    orthogroups.write_text(
        'busco_id\tSpecies_B\tdescription\tSpecies_A\tSpecies_C\n'
        '0001\tNA\tignored, annotation\t0001\t-\n'
        'NA\t\tannotation\tgene-name\tgene1, gene2\n'
        'OG3\t,\tannotation\tgene,\t,gene\n'
        'OG4\t \tannotation\t--\tgene,,other\n'
        'OG5\tnan\tannotation\tNULL\tgène,遺伝子\n'
        'OG6\ta,-\tannotation\t-\t\n'
        'OG7\t' + 'long_gene_' + 'x' * 10000 + '\tannotation\tgene\tgene\n'
    )

    orthogroup2genecount(
        str(orthogroups), str(counts), ['Species_C', 'Species_A', 'Species_B', 'Species_A'],
    )

    assert counts.read_text() == (
        'orthogroup_id\tSpecies_B\tSpecies_A\tSpecies_C\n'
        '0001\t1\t1\t0\n'
        'NA\t0\t1\t2\n'
        'OG3\t2\t2\t2\n'
        'OG4\t1\t1\t3\n'
        'OG5\t1\t1\t2\n'
        'OG6\t2\t0\t0\n'
        'OG7\t1\t1\t1\n'
    )


@pytest.mark.parametrize('rows', ['', '0001\tgene\nNA\t-\n'])
@pytest.mark.parametrize('species', [[], ['Species_A']])
def test_orthogroup_gene_counts_allow_empty_axes(tmp_path, rows, species):
    orthogroups = tmp_path / 'orthogroups.tsv'
    counts = tmp_path / 'counts.tsv'
    orthogroups.write_text('busco_id\tSpecies_A\n' + rows)

    orthogroup2genecount(str(orthogroups), str(counts), species)

    if species:
        expected = 'orthogroup_id\tSpecies_A\n' + ('0001\t1\nNA\t0\n' if rows else '')
    else:
        expected = 'orthogroup_id\n' + ('0001\nNA\n' if rows else '')
    assert counts.read_text() == expected


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
