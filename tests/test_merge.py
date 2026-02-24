import numpy
import pandas
import pytest
import os
from types import SimpleNamespace

from amalgkit.merge import (
    collect_valid_run_ids,
    merge_fastp_stats_into_metadata,
    merge_species_quant_tables,
    merge_main,
    scan_quant_abundance_paths,
)
from amalgkit.util import Metadata


class TestMergeFastpStatsIntoMetadata:
    def test_collect_valid_run_ids_filters_missing_entries(self):
        run_ids = collect_valid_run_ids([numpy.nan, 'SRR001', ' SRR001 ', '', 'SRR002'])
        assert run_ids == ['SRR001', 'SRR002']

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

    def test_merges_getfastq_stage_stats_tsv(self, tmp_path):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'scientific_name': ['sp1'],
            'exclusion': ['no'],
        }))
        srr001_dir = tmp_path / 'getfastq' / 'SRR001'
        srr001_dir.mkdir(parents=True)
        pandas.DataFrame([{
            'run': 'SRR001',
            'num_dumped': 10,
            'num_rejected': 3,
            'num_written': 7,
            'bp_dumped': 1000,
            'bp_rejected': 300,
            'bp_written': 700,
            'bp_discarded': 500,
            'fastp_duplication_rate': 12.5,
            'fastp_insert_size_peak': 250.0,
        }]).to_csv(srr001_dir / 'getfastq_stats.tsv', sep='\t', index=False)

        metadata = merge_fastp_stats_into_metadata(metadata, str(tmp_path))

        assert metadata.df.loc[0, 'num_dumped'] == 10
        assert metadata.df.loc[0, 'num_rejected'] == 3
        assert metadata.df.loc[0, 'bp_dumped'] == 1000
        assert metadata.df.loc[0, 'bp_rejected'] == 300
        assert metadata.df.loc[0, 'bp_discarded'] == 500
        assert metadata.df.loc[0, 'fastp_duplication_rate'] == 12.5
        assert metadata.df.loc[0, 'fastp_insert_size_peak'] == 250.0

    def test_rejects_metadata_without_run_column(self, tmp_path):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'scientific_name': ['sp1'],
            'exclusion': ['no'],
        }))
        metadata.df = metadata.df.drop(columns=['run'])

        with pytest.raises(ValueError, match='Missing required metadata column\\(s\\) for merge fastp stats: run'):
            merge_fastp_stats_into_metadata(metadata, str(tmp_path))

    def test_avoids_getfastq_root_scandir(self, tmp_path, monkeypatch):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'scientific_name': ['sp1'],
            'exclusion': ['no'],
        }))
        srr001_dir = tmp_path / 'getfastq' / 'SRR001'
        srr001_dir.mkdir(parents=True)
        pandas.DataFrame([{
            'run': 'SRR001',
            'fastp_duplication_rate': 10.0,
            'fastp_insert_size_peak': 111.0,
        }]).to_csv(srr001_dir / 'fastp_stats.tsv', sep='\t', index=False)
        getfastq_root = os.path.realpath(str(tmp_path / 'getfastq'))
        real_scandir = os.scandir

        def fail_on_root_scandir(path):
            if os.path.realpath(path) == getfastq_root:
                raise AssertionError('merge_fastp_stats_into_metadata should not scan getfastq root.')
            return real_scandir(path)

        monkeypatch.setattr('amalgkit.merge.os.scandir', fail_on_root_scandir)

        metadata = merge_fastp_stats_into_metadata(metadata, str(tmp_path))

        assert metadata.df.loc[0, 'fastp_duplication_rate'] == 10.0
        assert metadata.df.loc[0, 'fastp_insert_size_peak'] == 111.0

    def test_adds_columns_when_getfastq_dir_missing(self, tmp_path):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'scientific_name': ['sp1'],
            'exclusion': ['no'],
        }))

        metadata = merge_fastp_stats_into_metadata(metadata, str(tmp_path))

        assert 'fastp_duplication_rate' in metadata.df.columns
        assert 'fastp_insert_size_peak' in metadata.df.columns
        assert 'num_rejected' in metadata.df.columns
        assert 'bp_rejected' in metadata.df.columns
        assert numpy.isnan(metadata.df.loc[0, 'fastp_duplication_rate'])
        assert numpy.isnan(metadata.df.loc[0, 'fastp_insert_size_peak'])

    def test_rejects_getfastq_file_path(self, tmp_path):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'scientific_name': ['sp1'],
            'exclusion': ['no'],
        }))
        (tmp_path / 'getfastq').write_text('not a directory')

        with pytest.raises(NotADirectoryError, match='getfastq path exists but is not a directory'):
            merge_fastp_stats_into_metadata(metadata, str(tmp_path))

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

    def test_ignores_directory_named_fastp_stats_file(self, tmp_path):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'scientific_name': ['sp1'],
            'exclusion': ['no'],
        }))
        srr001_dir = tmp_path / 'getfastq' / 'SRR001'
        srr001_dir.mkdir(parents=True)
        (srr001_dir / 'fastp_stats.tsv').mkdir()

        metadata = merge_fastp_stats_into_metadata(metadata, str(tmp_path))

        assert numpy.isnan(metadata.df.loc[0, 'fastp_duplication_rate'])
        assert numpy.isnan(metadata.df.loc[0, 'fastp_insert_size_peak'])

    def test_ignores_missing_run_ids_while_merging_fastp_stats(self, tmp_path):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': [numpy.nan, 'SRR001'],
            'scientific_name': ['sp1', 'sp1'],
            'exclusion': ['no', 'no'],
        }))
        srr001_dir = tmp_path / 'getfastq' / 'SRR001'
        srr001_dir.mkdir(parents=True)
        pandas.DataFrame([{
            'run': 'SRR001',
            'fastp_duplication_rate': 22.5,
            'fastp_insert_size_peak': 180.0,
        }]).to_csv(srr001_dir / 'fastp_stats.tsv', sep='\t', index=False)

        metadata = merge_fastp_stats_into_metadata(metadata, str(tmp_path))

        assert metadata.df.loc[1, 'fastp_duplication_rate'] == 22.5
        assert metadata.df.loc[1, 'fastp_insert_size_peak'] == 180.0

    def test_merges_fastp_stats_for_duplicate_run_rows(self, tmp_path):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001', 'SRR001'],
            'scientific_name': ['sp1', 'sp1'],
            'exclusion': ['no', 'no'],
        }))
        srr001_dir = tmp_path / 'getfastq' / 'SRR001'
        srr001_dir.mkdir(parents=True)
        pandas.DataFrame([{
            'run': 'SRR001',
            'fastp_duplication_rate': 7.5,
            'fastp_insert_size_peak': 175.0,
        }]).to_csv(srr001_dir / 'fastp_stats.tsv', sep='\t', index=False)

        metadata = merge_fastp_stats_into_metadata(metadata, str(tmp_path))

        assert metadata.df['fastp_duplication_rate'].tolist() == [7.5, 7.5]
        assert metadata.df['fastp_insert_size_peak'].tolist() == [175.0, 175.0]

    def test_merges_fastp_stats_for_whitespace_run_ids(self, tmp_path):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': [' SRR001 '],
            'scientific_name': ['sp1'],
            'exclusion': ['no'],
        }))
        srr001_dir = tmp_path / 'getfastq' / 'SRR001'
        srr001_dir.mkdir(parents=True)
        pandas.DataFrame([{
            'run': 'SRR001',
            'fastp_duplication_rate': 19.0,
            'fastp_insert_size_peak': 210.0,
        }]).to_csv(srr001_dir / 'fastp_stats.tsv', sep='\t', index=False)

        metadata = merge_fastp_stats_into_metadata(metadata, str(tmp_path))
        assert metadata.df.loc[0, 'fastp_duplication_rate'] == 19.0
        assert metadata.df.loc[0, 'fastp_insert_size_peak'] == 210.0

    def test_respects_max_workers_override(self, tmp_path, monkeypatch):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001', 'SRR002'],
            'scientific_name': ['sp1', 'sp1'],
            'exclusion': ['no', 'no'],
        }))
        for run_id in ['SRR001', 'SRR002']:
            run_dir = tmp_path / 'getfastq' / run_id
            run_dir.mkdir(parents=True)
            pandas.DataFrame([{
                'run': run_id,
                'fastp_duplication_rate': 1.0,
                'fastp_insert_size_peak': 100.0,
            }]).to_csv(run_dir / 'fastp_stats.tsv', sep='\t', index=False)
        observed = {}

        def fake_run_tasks(task_items, task_fn, max_workers):
            observed['max_workers'] = max_workers
            return {}, []

        monkeypatch.setattr('amalgkit.merge.run_tasks_with_optional_threads', fake_run_tasks)
        merge_fastp_stats_into_metadata(metadata, str(tmp_path), max_workers=1)
        assert observed['max_workers'] == 1


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
    assert usecols_calls[1] == ('target_id', 'eff_length', 'est_counts', 'tpm')
    assert (merge_dir / 'Species_A' / 'Species_A_eff_length.tsv').exists()
    assert (merge_dir / 'Species_A' / 'Species_A_est_counts.tsv').exists()
    assert (merge_dir / 'Species_A' / 'Species_A_tpm.tsv').exists()
    eff = pandas.read_csv(merge_dir / 'Species_A' / 'Species_A_eff_length.tsv', sep='\t')
    assert list(eff.columns) == ['target_id', 'SRR001', 'SRR002']


def test_merge_species_quant_tables_parallel_reads(tmp_path, monkeypatch):
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
    observed = {'max_workers': None}

    def fake_run_tasks(task_items, task_fn, max_workers=1):
        observed['max_workers'] = max_workers
        results = {}
        failures = []
        for task_item in task_items:
            try:
                results[task_item] = task_fn(task_item)
            except Exception as exc:
                failures.append((task_item, exc))
        return results, failures

    monkeypatch.setattr('amalgkit.merge.run_tasks_with_optional_threads', fake_run_tasks)
    n = merge_species_quant_tables('Species A', metadata, str(quant_dir), str(merge_dir))

    assert n == 2
    assert observed['max_workers'] == 2
    assert (merge_dir / 'Species_A' / 'Species_A_tpm.tsv').exists()


def test_merge_species_quant_tables_rejects_mismatched_target_ids(tmp_path):
    quant_dir = tmp_path / 'quant'
    merge_dir = tmp_path / 'merge'
    (quant_dir / 'SRR001').mkdir(parents=True)
    (quant_dir / 'SRR002').mkdir(parents=True)
    pandas.DataFrame({
        'target_id': ['tx1', 'tx2'],
        'eff_length': [1.1, 1.2],
        'est_counts': [1.3, 1.4],
        'tpm': [1.5, 1.6],
    }).to_csv(quant_dir / 'SRR001' / 'SRR001_abundance.tsv', sep='\t', index=False)
    pandas.DataFrame({
        'target_id': ['tx2', 'tx1'],
        'eff_length': [10.1, 10.2],
        'est_counts': [10.3, 10.4],
        'tpm': [10.5, 10.6],
    }).to_csv(quant_dir / 'SRR002' / 'SRR002_abundance.tsv', sep='\t', index=False)

    metadata = Metadata.from_DataFrame(pandas.DataFrame({
        'run': ['SRR001', 'SRR002'],
        'scientific_name': ['Species A', 'Species A'],
        'exclusion': ['no', 'no'],
    }))

    with pytest.raises(ValueError, match='Mismatched target_id rows across quant files'):
        merge_species_quant_tables('Species A', metadata, str(quant_dir), str(merge_dir))


def test_merge_species_quant_tables_reports_run_when_required_column_missing(tmp_path):
    quant_dir = tmp_path / 'quant'
    merge_dir = tmp_path / 'merge'
    (quant_dir / 'SRR001').mkdir(parents=True)
    pandas.DataFrame({
        'target_id': ['tx1', 'tx2'],
        'eff_length': [1.1, 1.2],
        'est_counts': [1.3, 1.4],
        # tpm column intentionally missing
    }).to_csv(quant_dir / 'SRR001' / 'SRR001_abundance.tsv', sep='\t', index=False)

    metadata = Metadata.from_DataFrame(pandas.DataFrame({
        'run': ['SRR001'],
        'scientific_name': ['Species A'],
        'exclusion': ['no'],
    }))

    with pytest.raises(ValueError, match=r'Failed to read quant output table for run SRR001'):
        merge_species_quant_tables('Species A', metadata, str(quant_dir), str(merge_dir))


def test_merge_species_quant_tables_ignores_missing_run_ids(tmp_path):
    quant_dir = tmp_path / 'quant'
    merge_dir = tmp_path / 'merge'
    (quant_dir / 'SRR001').mkdir(parents=True)
    pandas.DataFrame({
        'target_id': ['tx1', 'tx2'],
        'eff_length': [1.1, 1.2],
        'est_counts': [1.3, 1.4],
        'tpm': [1.5, 1.6],
    }).to_csv(quant_dir / 'SRR001' / 'SRR001_abundance.tsv', sep='\t', index=False)

    metadata = Metadata.from_DataFrame(pandas.DataFrame({
        'run': ['SRR001', numpy.nan, ''],
        'scientific_name': ['Species A', 'Species A', 'Species A'],
        'exclusion': ['no', 'no', 'no'],
    }))

    n = merge_species_quant_tables('Species A', metadata, str(quant_dir), str(merge_dir))

    assert n == 1
    eff = pandas.read_csv(merge_dir / 'Species_A' / 'Species_A_eff_length.tsv', sep='\t')
    assert list(eff.columns) == ['target_id', 'SRR001']


def test_merge_species_quant_tables_handles_whitespace_species_and_exclusion(tmp_path):
    quant_dir = tmp_path / 'quant'
    merge_dir = tmp_path / 'merge'
    (quant_dir / 'SRR001').mkdir(parents=True)
    pandas.DataFrame({
        'target_id': ['tx1', 'tx2'],
        'eff_length': [1.1, 1.2],
        'est_counts': [1.3, 1.4],
        'tpm': [1.5, 1.6],
    }).to_csv(quant_dir / 'SRR001' / 'SRR001_abundance.tsv', sep='\t', index=False)

    metadata = Metadata.from_DataFrame(pandas.DataFrame({
        'run': ['SRR001'],
        'scientific_name': [' Species A '],
        'exclusion': [' no '],
    }))

    n = merge_species_quant_tables('Species A', metadata, str(quant_dir), str(merge_dir))

    assert n == 1
    eff = pandas.read_csv(merge_dir / 'Species_A' / 'Species_A_eff_length.tsv', sep='\t')
    assert list(eff.columns) == ['target_id', 'SRR001']


def test_merge_species_quant_tables_excludes_runs_marked_as_exclusion_yes(tmp_path):
    quant_dir = tmp_path / 'quant'
    merge_dir = tmp_path / 'merge'
    for run_id, base in [('SRR001', 1.0), ('SRR002', 10.0)]:
        run_dir = quant_dir / run_id
        run_dir.mkdir(parents=True)
        pandas.DataFrame({
            'target_id': ['tx1', 'tx2'],
            'eff_length': [base + 0.1, base + 0.2],
            'est_counts': [base + 0.3, base + 0.4],
            'tpm': [base + 0.5, base + 0.6],
        }).to_csv(run_dir / f'{run_id}_abundance.tsv', sep='\t', index=False)

    metadata = Metadata.from_DataFrame(pandas.DataFrame({
        'run': ['SRR001', 'SRR002'],
        'scientific_name': ['Species A', 'Species A'],
        'exclusion': ['no', 'yes'],
    }))

    n = merge_species_quant_tables('Species A', metadata, str(quant_dir), str(merge_dir))

    assert n == 1
    eff = pandas.read_csv(merge_dir / 'Species_A' / 'Species_A_eff_length.tsv', sep='\t')
    assert list(eff.columns) == ['target_id', 'SRR001']


def test_scan_quant_abundance_paths_filters_target_runs(tmp_path):
    quant_dir = tmp_path / 'quant'
    (quant_dir / 'SRR001').mkdir(parents=True)
    (quant_dir / 'SRR002').mkdir(parents=True)
    (quant_dir / 'OTHER').mkdir(parents=True)
    (quant_dir / 'SRR001' / 'SRR001_abundance.tsv').write_text('target_id\teff_length\test_counts\ttpm\n')
    (quant_dir / 'SRR002' / 'SRR002_abundance.tsv').write_text('target_id\teff_length\test_counts\ttpm\n')
    (quant_dir / 'OTHER' / 'OTHER_abundance.tsv').write_text('target_id\teff_length\test_counts\ttpm\n')

    detected = scan_quant_abundance_paths(str(quant_dir), target_runs={'SRR001', 'SRR002'})

    assert set(detected.keys()) == {'SRR001', 'SRR002'}
    assert detected['SRR001'].endswith('SRR001_abundance.tsv')
    assert detected['SRR002'].endswith('SRR002_abundance.tsv')


def test_scan_quant_abundance_paths_ignores_directory_named_abundance_file(tmp_path):
    quant_dir = tmp_path / 'quant'
    run_dir = quant_dir / 'SRR001'
    run_dir.mkdir(parents=True)
    (run_dir / 'SRR001_abundance.tsv').mkdir()

    detected = scan_quant_abundance_paths(str(quant_dir), target_runs={'SRR001'})

    assert detected == {}


def test_scan_quant_abundance_paths_rejects_file_path(tmp_path):
    quant_path = tmp_path / 'quant'
    quant_path.write_text('not a directory')

    with pytest.raises(NotADirectoryError, match='not a directory'):
        scan_quant_abundance_paths(str(quant_path))


def test_merge_main_rejects_nonpositive_species_jobs(tmp_path, monkeypatch):
    args = SimpleNamespace(out_dir=str(tmp_path), internal_jobs=0, metadata='inferred')
    monkeypatch.setattr('amalgkit.merge.check_rscript', lambda: None)
    with pytest.raises(ValueError, match='--internal_jobs must be > 0'):
        merge_main(args)


def test_merge_main_rejects_out_dir_file_path(tmp_path, monkeypatch):
    out_path = tmp_path / 'out_path'
    out_path.write_text('not a directory')
    args = SimpleNamespace(out_dir=str(out_path), internal_jobs=1, metadata='inferred')
    monkeypatch.setattr('amalgkit.merge.check_rscript', lambda: None)

    with pytest.raises(NotADirectoryError, match='Output path exists but is not a directory'):
        merge_main(args)


def test_merge_main_rejects_merge_path_file(tmp_path, monkeypatch):
    out_dir = tmp_path / 'out'
    out_dir.mkdir()
    (out_dir / 'merge').write_text('not a directory')
    args = SimpleNamespace(out_dir=str(out_dir), internal_jobs=1, metadata='inferred')
    monkeypatch.setattr('amalgkit.merge.check_rscript', lambda: None)

    with pytest.raises(NotADirectoryError, match='Merge path exists but is not a directory'):
        merge_main(args)


def test_merge_main_rejects_metadata_without_valid_scientific_name(tmp_path, monkeypatch):
    out_dir = tmp_path / 'out'
    out_dir.mkdir()
    metadata = Metadata.from_DataFrame(pandas.DataFrame({
        'run': ['R1'],
        'scientific_name': ['   '],
        'exclusion': ['no'],
    }))
    args = SimpleNamespace(out_dir=str(out_dir), internal_jobs=1, metadata='inferred')
    monkeypatch.setattr('amalgkit.merge.check_rscript', lambda: None)
    monkeypatch.setattr('amalgkit.merge.load_metadata', lambda _args: metadata)

    with pytest.raises(ValueError, match='No valid scientific_name entries were found'):
        merge_main(args)


def test_merge_main_rejects_metadata_without_run_column(tmp_path, monkeypatch):
    out_dir = tmp_path / 'out'
    out_dir.mkdir()
    metadata = Metadata.from_DataFrame(pandas.DataFrame({
        'scientific_name': ['Species A'],
        'exclusion': ['no'],
    }))
    metadata.df = metadata.df.drop(columns=['run'])
    args = SimpleNamespace(out_dir=str(out_dir), internal_jobs=1, metadata='inferred')
    monkeypatch.setattr('amalgkit.merge.check_rscript', lambda: None)
    monkeypatch.setattr('amalgkit.merge.load_metadata', lambda _args: metadata)

    with pytest.raises(ValueError, match='Missing required metadata column\\(s\\) for merge: run'):
        merge_main(args)


def test_merge_main_raises_when_no_quant_abundance_detected(tmp_path, monkeypatch):
    out_dir = tmp_path / 'out'
    out_dir.mkdir()
    metadata = Metadata.from_DataFrame(pandas.DataFrame({
        'run': ['R1'],
        'scientific_name': ['Species A'],
        'exclusion': ['no'],
    }))
    args = SimpleNamespace(out_dir=str(out_dir), internal_jobs=1, metadata='inferred')

    monkeypatch.setattr('amalgkit.merge.check_rscript', lambda: None)
    monkeypatch.setattr('amalgkit.merge.load_metadata', lambda _args: metadata)
    monkeypatch.setattr('amalgkit.merge.merge_species_quant_tables', lambda **_kwargs: 0)
    monkeypatch.setattr(
        'amalgkit.merge.merge_fastp_stats_into_metadata',
        lambda *_args, **_kwargs: (_ for _ in ()).throw(AssertionError('should not run postprocess when no quant output exists')),
    )

    with pytest.raises(FileNotFoundError, match='No quant abundance file was detected for any species'):
        merge_main(args)


def test_merge_main_parallel_species_jobs(tmp_path, monkeypatch):
    out_dir = tmp_path / 'out'
    out_dir.mkdir()
    metadata = Metadata.from_DataFrame(pandas.DataFrame({
        'run': ['R1', 'R2'],
        'scientific_name': ['Species A', 'Species B'],
        'exclusion': ['no', 'no'],
    }))
    args = SimpleNamespace(out_dir=str(out_dir), internal_jobs=2, metadata='inferred')
    processed = []

    def fake_merge_species(sp, metadata=None, quant_dir=None, merge_dir=None, run_abundance_paths=None):
        processed.append(sp)
        return 1

    monkeypatch.setattr('amalgkit.merge.check_rscript', lambda: None)
    monkeypatch.setattr('amalgkit.merge.load_metadata', lambda _args: metadata)
    monkeypatch.setattr('amalgkit.merge.merge_species_quant_tables', fake_merge_species)
    monkeypatch.setattr('amalgkit.merge.merge_fastp_stats_into_metadata', lambda m, _d, max_workers='auto': m)
    monkeypatch.setattr('amalgkit.merge.write_updated_metadata', lambda _m, _p, _a, max_workers='auto': None)
    monkeypatch.setattr('amalgkit.merge.subprocess.check_call', lambda _cmd: 0)

    merge_main(args)

    assert set(processed) == {'Species A', 'Species B'}


def test_merge_main_cpu_budget_caps_species_jobs_to_serial(tmp_path, monkeypatch):
    out_dir = tmp_path / 'out'
    out_dir.mkdir()
    metadata = Metadata.from_DataFrame(pandas.DataFrame({
        'run': ['R1', 'R2'],
        'scientific_name': ['Species A', 'Species B'],
        'exclusion': ['no', 'no'],
    }))
    args = SimpleNamespace(out_dir=str(out_dir), internal_jobs=4, internal_cpu_budget=1, metadata='inferred')
    processed = []

    def fake_merge_species(sp, metadata=None, quant_dir=None, merge_dir=None, run_abundance_paths=None):
        processed.append(sp)
        return 1

    def fail_if_called(*_args, **_kwargs):
        raise AssertionError('run_tasks_with_optional_threads should not be used when --internal_cpu_budget caps internal_jobs to 1.')

    monkeypatch.setattr('amalgkit.merge.check_rscript', lambda: None)
    monkeypatch.setattr('amalgkit.merge.load_metadata', lambda _args: metadata)
    monkeypatch.setattr('amalgkit.merge.merge_species_quant_tables', fake_merge_species)
    monkeypatch.setattr('amalgkit.merge.run_tasks_with_optional_threads', fail_if_called)
    monkeypatch.setattr('amalgkit.merge.merge_fastp_stats_into_metadata', lambda m, _d, max_workers='auto': m)
    monkeypatch.setattr('amalgkit.merge.write_updated_metadata', lambda _m, _p, _a, max_workers='auto': None)
    monkeypatch.setattr('amalgkit.merge.subprocess.check_call', lambda _cmd: 0)

    merge_main(args)

    assert set(processed) == {'Species A', 'Species B'}


def test_merge_main_calls_check_rscript(tmp_path, monkeypatch):
    out_dir = tmp_path / 'out'
    out_dir.mkdir()
    metadata = Metadata.from_DataFrame(pandas.DataFrame({
        'run': ['R1'],
        'scientific_name': ['Species A'],
        'exclusion': ['no'],
    }))
    args = SimpleNamespace(out_dir=str(out_dir), internal_jobs=1, metadata='inferred')
    called = {'check_rscript': 0}

    def fake_check_rscript():
        called['check_rscript'] += 1

    monkeypatch.setattr('amalgkit.merge.check_rscript', fake_check_rscript)
    monkeypatch.setattr('amalgkit.merge.load_metadata', lambda _args: metadata)
    monkeypatch.setattr('amalgkit.merge.merge_species_quant_tables', lambda **_kwargs: 1)
    monkeypatch.setattr('amalgkit.merge.merge_fastp_stats_into_metadata', lambda m, _d, max_workers='auto': m)
    monkeypatch.setattr('amalgkit.merge.write_updated_metadata', lambda _m, _p, _a, max_workers='auto': None)
    monkeypatch.setattr('amalgkit.merge.subprocess.check_call', lambda _cmd: 0)

    merge_main(args)

    assert called['check_rscript'] == 1
