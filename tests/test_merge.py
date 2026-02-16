import numpy
import pandas
import pytest
import os
from types import SimpleNamespace

from amalgkit.merge import merge_fastp_stats_into_metadata, merge_species_quant_tables, merge_main
from amalgkit.util import Metadata


class TestMergeFastpStatsIntoMetadata:
    def test_merges_fastp_stats_tsv(self, tmp_path):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001', 'SRR002'],
            'scientific_name': ['sp1', 'sp1'],
            'exclusion': ['no', 'no'],
        }))
        srr001_dir = tmp_path / 'getfastq' / 'SRR001'
        srr001_dir.mkdir(parents=True)
        pandas.DataFrame([{
            'run': 'SRR001',
            'fastp_duplication_rate': 12.5,
            'fastp_insert_size_peak': 250.0,
        }]).to_csv(srr001_dir / 'fastp_stats.tsv', sep='\t', index=False)

        metadata = merge_fastp_stats_into_metadata(metadata, str(tmp_path))

        assert metadata.df.loc[0, 'fastp_duplication_rate'] == 12.5
        assert metadata.df.loc[0, 'fastp_insert_size_peak'] == 250.0
        assert numpy.isnan(metadata.df.loc[1, 'fastp_duplication_rate'])
        assert numpy.isnan(metadata.df.loc[1, 'fastp_insert_size_peak'])

    def test_adds_columns_when_getfastq_dir_missing(self, tmp_path):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'scientific_name': ['sp1'],
            'exclusion': ['no'],
        }))

        metadata = merge_fastp_stats_into_metadata(metadata, str(tmp_path))

        assert 'fastp_duplication_rate' in metadata.df.columns
        assert 'fastp_insert_size_peak' in metadata.df.columns
        assert numpy.isnan(metadata.df.loc[0, 'fastp_duplication_rate'])
        assert numpy.isnan(metadata.df.loc[0, 'fastp_insert_size_peak'])

    def test_keeps_nan_when_column_missing_in_tsv(self, tmp_path):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'scientific_name': ['sp1'],
            'exclusion': ['no'],
        }))
        srr001_dir = tmp_path / 'getfastq' / 'SRR001'
        srr001_dir.mkdir(parents=True)
        pandas.DataFrame([{
            'run': 'SRR001',
            'fastp_duplication_rate': 44.0,
        }]).to_csv(srr001_dir / 'fastp_stats.tsv', sep='\t', index=False)

        metadata = merge_fastp_stats_into_metadata(metadata, str(tmp_path))

        assert metadata.df.loc[0, 'fastp_duplication_rate'] == 44.0
        assert numpy.isnan(metadata.df.loc[0, 'fastp_insert_size_peak'])

    def test_ignores_non_directory_entries_in_getfastq_root(self, tmp_path):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'scientific_name': ['sp1'],
            'exclusion': ['no'],
        }))
        getfastq_root = tmp_path / 'getfastq'
        getfastq_root.mkdir(parents=True)
        (getfastq_root / 'SRR001').write_text('not a directory')

        metadata = merge_fastp_stats_into_metadata(metadata, str(tmp_path))

        assert numpy.isnan(metadata.df.loc[0, 'fastp_duplication_rate'])
        assert numpy.isnan(metadata.df.loc[0, 'fastp_insert_size_peak'])


def test_merge_species_quant_tables_single_pass_reads(tmp_path, monkeypatch):
    quant_dir = tmp_path / 'quant'
    merge_dir = tmp_path / 'merge'
    (quant_dir / 'SRR001').mkdir(parents=True)
    (quant_dir / 'SRR002').mkdir(parents=True)
    for run, base in [('SRR001', 1.0), ('SRR002', 10.0)]:
        pandas.DataFrame({
            'target_id': ['tx1', 'tx2'],
            'eff_length': [base + 0.1, base + 0.2],
            'est_counts': [base + 0.3, base + 0.4],
            'tpm': [base + 0.5, base + 0.6],
        }).to_csv(quant_dir / run / f'{run}_abundance.tsv', sep='\t', index=False)

    metadata = Metadata.from_DataFrame(pandas.DataFrame({
        'run': ['SRR001', 'SRR002'],
        'scientific_name': ['Species A', 'Species A'],
        'exclusion': ['no', 'no'],
    }))

    original_read_csv = pandas.read_csv
    read_calls = {'n': 0}
    usecols_calls = []

    def counting_read_csv(path, *args, **kwargs):
        if str(path).endswith('_abundance.tsv'):
            read_calls['n'] += 1
            usecols_calls.append(tuple(kwargs.get('usecols', [])))
        return original_read_csv(path, *args, **kwargs)

    monkeypatch.setattr('amalgkit.merge.pandas.read_csv', counting_read_csv)
    n = merge_species_quant_tables('Species A', metadata, str(quant_dir), str(merge_dir))
    assert n == 2
    assert read_calls['n'] == 2
    assert usecols_calls[0] == ('target_id', 'eff_length', 'est_counts', 'tpm')
    assert usecols_calls[1] == ('eff_length', 'est_counts', 'tpm')
    assert (merge_dir / 'Species_A' / 'Species_A_eff_length.tsv').exists()
    assert (merge_dir / 'Species_A' / 'Species_A_est_counts.tsv').exists()
    assert (merge_dir / 'Species_A' / 'Species_A_tpm.tsv').exists()
    eff = pandas.read_csv(merge_dir / 'Species_A' / 'Species_A_eff_length.tsv', sep='\t')
    assert list(eff.columns) == ['target_id', 'SRR001', 'SRR002']


def test_merge_main_rejects_nonpositive_species_jobs(tmp_path):
    args = SimpleNamespace(out_dir=str(tmp_path), species_jobs=0, metadata='inferred')
    with pytest.raises(ValueError, match='--species_jobs must be > 0'):
        merge_main(args)


def test_merge_main_parallel_species_jobs(tmp_path, monkeypatch):
    out_dir = tmp_path / 'out'
    out_dir.mkdir()
    metadata = Metadata.from_DataFrame(pandas.DataFrame({
        'run': ['R1', 'R2'],
        'scientific_name': ['Species A', 'Species B'],
        'exclusion': ['no', 'no'],
    }))
    args = SimpleNamespace(out_dir=str(out_dir), species_jobs=2, metadata='inferred')
    processed = []

    def fake_merge_species(sp, _metadata, _quant_dir, _merge_dir):
        processed.append(sp)
        return 1

    monkeypatch.setattr('amalgkit.merge.load_metadata', lambda _args: metadata)
    monkeypatch.setattr('amalgkit.merge.merge_species_quant_tables', fake_merge_species)
    monkeypatch.setattr('amalgkit.merge.merge_fastp_stats_into_metadata', lambda m, _d: m)
    monkeypatch.setattr('amalgkit.merge.write_updated_metadata', lambda _m, _p, _a: None)
    monkeypatch.setattr('amalgkit.merge.subprocess.check_call', lambda _cmd: 0)

    merge_main(args)

    assert set(processed) == {'Species A', 'Species B'}
