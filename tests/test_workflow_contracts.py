"""Exercise file contracts through real producers/consumers, stubbing only tools."""

import gzip
import os
from pathlib import Path
from types import SimpleNamespace

import numpy
import pandas
import pytest

from amalgkit import quant
from amalgkit.batch_effect_io import read_expression_matrix_tsv
from amalgkit.cross_species_filter import _select_single_copy_orthogroups
from amalgkit.cstmm_python import _get_df_exp_single_copy_ortholog, _read_est_counts
from amalgkit.getfastq import sequence_extraction_private
from amalgkit.getfastq_stats import read_getfastq_stats_row
from amalgkit.main import build_main_parser
from amalgkit.merge import merge_main, merge_species_quant_tables
from amalgkit.metadata_utils import Metadata
from amalgkit.orthology_utils import orthogroup2genecount, read_busco_species_table
from amalgkit.per_species_finalize_python import _read_expression_tsv


pytestmark = pytest.mark.integration


@pytest.mark.parametrize('run_id', ['0001', 'NA'])
def test_resume_stats_preserve_lexical_run_identifiers(tmp_path, run_id):
    (tmp_path / 'getfastq_stats.tsv').write_text(f'run\tnum_written\n{run_id}\t2\n', encoding='utf-8')
    row = read_getfastq_stats_row(str(tmp_path), run_id)
    assert row is not None
    assert row['run'] == run_id
    assert row['num_written'] == 2


def _metadata(layout='single', read_paths=()):
    return Metadata.from_DataFrame(pandas.DataFrame({
        'run': ['R1'], 'scientific_name': ['Species alpha'],
        'sample_group': ['leaf'], 'bioproject': ['P1'],
        'lib_layout': [layout], 'total_spots': [2],
        'total_bases': [200], 'spot_length': [100],
        'nominal_length': [200], 'nominal_sdev': [20],
        'exclusion': ['no'], 'private_file': ['yes'],
        'read1_path': [str(read_paths[0]) if read_paths else ''],
        'read2_path': [str(read_paths[1]) if len(read_paths) > 1 else ''],
    }))


@pytest.mark.parametrize('compressed', [False, True])
@pytest.mark.parametrize('layout', ['single', 'paired'])
def test_private_fastq_producer_to_quant_default_cleanup(tmp_path, monkeypatch, compressed, layout):
    sources = []
    for mate in range(2 if layout == 'paired' else 1):
        source = tmp_path / ('source' + str(mate) + ('.fastq.gz' if compressed else '.fastq'))
        opener = gzip.open if compressed else open
        with opener(source, 'wt', encoding='utf-8') as handle:
            for read in range(2):
                handle.write('@r{}\n{}\n+\n{}\n'.format(read, 'AC'[mate] * 100, 'I' * 100))
        sources.append(source)
    original_bytes = [path.read_bytes() for path in sources]
    original_inodes = [path.stat().st_ino for path in sources]
    metadata = _metadata(layout, sources)
    run_dir = tmp_path / 'getfastq' / 'R1'
    run_dir.mkdir(parents=True)
    parser = build_main_parser()
    getfastq_args = parser.parse_args([
        'getfastq', '--out_dir', str(tmp_path), '--fastp', 'no', '--rrna_filter', 'no', '--contam_filter', 'no',
    ])
    sra_stat = {'sra_id': 'R1', 'layout': layout, 'getfastq_sra_dir': str(run_dir), 'metadata_idx': 0}
    sequence_extraction_private(metadata, sra_stat, getfastq_args)
    managed = sorted(run_dir.glob('*.amalgkit.fastq.gz'))
    assert len(managed) == len(sources)
    assert all(path.is_symlink() == compressed for path in managed)

    def kallisto(_args, in_files, metadata, sra_stat, output_dir, index):
        assert len(in_files) == len(sources)
        for path in in_files:
            with gzip.open(path, 'rt') as handle:
                assert handle.read().startswith('@r0\n')
        Path(output_dir, 'R1_abundance.tsv').write_text(
            'target_id\tlength\teff_length\test_counts\ttpm\n0001\t100\t90\t2\t1000000\n',
            encoding='utf-8',
        )
        Path(output_dir, 'R1_run_info.json').write_text('{"p_pseudoaligned":100}', encoding='utf-8')

    monkeypatch.setattr(quant, 'call_kallisto', kallisto)
    quant_args = parser.parse_args(['quant', '--out_dir', str(tmp_path)])
    assert quant_args.clean_fastq
    quant.run_quant(quant_args, metadata, 'R1', 'unused.idx', backend='kallisto')

    assert quant.validate_quant_outputs('R1', str(tmp_path / 'quant' / 'R1'))[0]
    assert all(not os.path.lexists(path) for path in managed)
    assert all(Path(str(path) + '.safely_removed').is_file() for path in managed)
    assert [path.read_bytes() for path in sources] == original_bytes
    assert [path.stat().st_ino for path in sources] == original_inodes


@pytest.mark.parametrize('column', ['eff_length', 'est_counts', 'tpm'])
@pytest.mark.parametrize('bad_value', ['-1', 'NaN', 'inf', 'broken', ''])
def test_merge_rejects_invalid_values_without_replacing_previous_output(tmp_path, column, bad_value):
    run_dir = tmp_path / 'quant' / 'R1'
    run_dir.mkdir(parents=True)
    values = dict(target_id='0001', length='100', eff_length='90', est_counts='2', tpm='1000000')
    values[column] = bad_value
    (run_dir / 'R1_abundance.tsv').write_text(
        '\t'.join(values) + '\n' + '\t'.join(values.values()) + '\n', encoding='utf-8',
    )
    prior = tmp_path / 'merge' / 'valuable.tsv'
    prior.parent.mkdir()
    prior.write_text('prior result', encoding='utf-8')
    metadata_path = tmp_path / 'metadata.tsv'
    _metadata().df.to_csv(metadata_path, sep='\t', index=False)
    args = SimpleNamespace(out_dir=str(tmp_path), metadata=str(metadata_path), internal_jobs=1)

    with pytest.raises(ValueError, match='R1.*' + column):
        merge_main(args)

    assert prior.read_text(encoding='utf-8') == 'prior result'
    assert list(prior.parent.iterdir()) == [prior]


def test_merge_rejects_header_only_input(tmp_path):
    run_dir = tmp_path / 'quant' / 'R1'
    run_dir.mkdir(parents=True)
    (run_dir / 'R1_abundance.tsv').write_text('target_id\teff_length\test_counts\ttpm\n', encoding='utf-8')
    with pytest.raises(ValueError, match='did not contain any data rows'):
        merge_species_quant_tables('Species alpha', _metadata(), str(run_dir.parent), str(tmp_path / 'merge'))


def test_lexical_gene_ids_survive_merge_orthology_and_downstream_reads(tmp_path):
    run_dir = tmp_path / 'quant' / 'R1'
    run_dir.mkdir(parents=True)
    ids = ['0001', '1', 'NA']
    pandas.DataFrame({
        'target_id': ids, 'eff_length': [90] * 3, 'est_counts': [2, 3, 4], 'tpm': [100, 200, 300],
    }).to_csv(run_dir / 'R1_abundance.tsv', sep='\t', index=False)
    count_dir = tmp_path / 'merge'
    merge_species_quant_tables('Species alpha', _metadata(), str(run_dir.parent), str(count_dir))
    counts = _read_est_counts(str(count_dir), 'Species_alpha')
    assert counts.index.tolist() == ids
    table_path = count_dir / 'Species_alpha' / 'Species_alpha_est_counts.tsv'
    for reader in (_read_expression_tsv, read_expression_matrix_tsv):
        assert reader(table_path).index.tolist() == ids

    orthogroup_path = tmp_path / 'orthogroups.tsv'
    pandas.DataFrame({'busco_id': ['00001', '1', 'NA'], 'Species_alpha': ids}).to_csv(
        orthogroup_path, sep='\t', index=False,
    )
    genecount_path = tmp_path / 'genecount.tsv'
    orthogroup2genecount(str(orthogroup_path), str(genecount_path), ['Species_alpha'])
    joined = _get_df_exp_single_copy_ortholog(
        str(genecount_path), str(orthogroup_path), str(count_dir), {'Species_alpha': counts},
    )
    numpy.testing.assert_array_equal(joined.iloc[:, 0], [2, 3, 4])
    assert joined.index.tolist() == ['00001', '1', 'NA']
    selected = _select_single_copy_orthogroups(str(orthogroup_path), str(genecount_path), ['Species_alpha'])
    assert selected['Species_alpha'].tolist() == ids

    busco_path = tmp_path / 'busco.tsv'
    busco_path.write_text('00001\tComplete\t0001\t50\t100\t-\t-\n', encoding='utf-8')
    busco = read_busco_species_table(busco_path)
    assert busco.loc[0, 'sequence'] == '0001'
    assert busco.loc[0, 'busco_id'] == '00001'
