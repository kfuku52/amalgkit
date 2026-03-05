import shutil
import subprocess
from pathlib import Path

import pandas
import pytest


REQUIRED_R_PACKAGES = [
    'Rtsne',
    'ggplot2',
    'sva',
]


def _has_rscript_and_required_packages():
    if shutil.which('Rscript') is None:
        return False
    quoted = ','.join(['"{}"'.format(p) for p in REQUIRED_R_PACKAGES])
    expr = (
        'pkgs <- c({}); '
        'ok <- all(vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)); '
        'quit(status=ifelse(ok, 0, 1))'
    ).format(quoted)
    out = subprocess.run(['Rscript', '-e', expr], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return out.returncode == 0


@pytest.fixture(scope='module')
def require_finalize_r_runtime():
    if not _has_rscript_and_required_packages():
        pytest.skip('Rscript or required R packages are not available for finalize integration tests.')


def _repo_root():
    return Path(__file__).resolve().parents[1]


def _write_finalize_fixture(tmp_path, sample_groups, bioprojects, species='Finalizus example'):
    runs = ['RUN{:02d}'.format(i + 1) for i in range(len(sample_groups))]
    genes = ['G{:03d}'.format(i + 1) for i in range(50)]
    metadata = pandas.DataFrame({
        'run': runs,
        'scientific_name': [species] * len(runs),
        'sample_group': sample_groups,
        'bioproject': bioprojects,
        'exclusion': ['no'] * len(runs),
        'mapping_rate': [100.0] * len(runs),
        'lib_layout': ['PAIRED'] * len(runs),
        'lib_selection': ['cDNA'] * len(runs),
        'instrument': ['Illumina'] * len(runs),
        'total_spots': [1_000_000 + i * 10_000 for i in range(len(runs))],
        'total_bases': [150_000_000 + i * 1_000_000 for i in range(len(runs))],
    })
    metadata_path = tmp_path / 'metadata.tsv'
    metadata.to_csv(metadata_path, sep='\t', index=False)

    counts = pandas.DataFrame({'target_id': genes})
    for run_id, sample_group in zip(runs, sample_groups):
        if sample_group == 'A':
            values = [100 + i for i in range(25)] + [5 + i for i in range(25)]
        else:
            values = [5 + i for i in range(25)] + [100 + i for i in range(25)]
        counts[run_id] = values
    count_path = tmp_path / 'counts.tsv'
    counts.to_csv(count_path, sep='\t', index=False)

    eff_length = pandas.DataFrame({'target_id': genes})
    for run in runs:
        eff_length[run] = 1000
    eff_length_path = tmp_path / 'eff_length.tsv'
    eff_length.to_csv(eff_length_path, sep='\t', index=False)

    return {
        'species': species,
        'species_tag': species.replace(' ', '_'),
        'metadata_path': metadata_path,
        'count_path': count_path,
        'eff_length_path': eff_length_path,
        'sample_groups': sample_groups,
    }


def _run_finalize_r(tmp_path, fixture, sva_nsv='auto', sva_B='auto', sva_B_auto_max='80', seed='7'):
    repo = _repo_root()
    out_dir = tmp_path / 'out'
    out_dir.mkdir(parents=True, exist_ok=True)
    selected_sample_groups = '|'.join(sorted(set(fixture['sample_groups'])))
    cmd = [
        'Rscript',
        str(repo / 'amalgkit' / 'finalize.r'),
        str(fixture['count_path']),
        str(fixture['metadata_path']),
        str(out_dir),
        str(fixture['eff_length_path']),
        selected_sample_groups,
        'DEFAULT',
        'log2p1-fpkm',
        'sva',
        '1',    # clip_negative
        '1',    # maintain_zero
        str(repo / 'amalgkit' / 'util.r'),
        'auto',  # ruvseq_control_genes
        'auto',  # ruvseq_k
        '5',     # ruvseq_k_max
        '1000',  # ruvseq_control_top_n
        '100',   # ruvseq_min_controls
        str(seed),
        str(sva_nsv),
        str(sva_B),
        str(sva_B_auto_max),
    ]
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    return proc, out_dir


def _read_batch_summary(out_dir, species_tag):
    path = out_dir / 'curate' / species_tag / 'tables' / '{}.sva.batch_effect_summary.tsv'.format(species_tag)
    assert path.exists()
    return pandas.read_csv(path, sep='\t')


def _read_species_metadata(out_dir, species_tag):
    path = out_dir / 'curate' / species_tag / 'tables' / '{}.metadata.tsv'.format(species_tag)
    assert path.exists()
    return pandas.read_csv(path, sep='\t')


def test_finalize_r_handles_confounded_design_and_writes_batch_summary(require_finalize_r_runtime, tmp_path):
    fixture = _write_finalize_fixture(
        tmp_path=tmp_path,
        sample_groups=['A', 'A', 'B', 'B'],
        bioprojects=['BP1', 'BP1', 'BP2', 'BP2'],  # perfectly confounded with sample_group
        species='Confoundus example',
    )
    proc, out_dir = _run_finalize_r(tmp_path=tmp_path, fixture=fixture, sva_nsv='auto', sva_B='auto')
    assert proc.returncode == 0, proc.stdout

    summary_df = _read_batch_summary(out_dir, fixture['species_tag'])
    assert {'resolved_sva_nsv', 'resolved_sva_B', 'skip_reason', 'corrected_run_count'}.issubset(set(summary_df.columns))
    metadata_df = _read_species_metadata(out_dir, fixture['species_tag'])
    assert {'batch_corrected', 'batch_alg_used'}.issubset(set(metadata_df.columns))


def test_finalize_r_single_sample_group_does_not_crash_and_marks_uncorrected(require_finalize_r_runtime, tmp_path):
    fixture = _write_finalize_fixture(
        tmp_path=tmp_path,
        sample_groups=['A', 'A', 'A', 'A'],
        bioprojects=['BP1', 'BP1', 'BP2', 'BP2'],
        species='Singlesamplegroup example',
    )
    proc, out_dir = _run_finalize_r(tmp_path=tmp_path, fixture=fixture, sva_nsv='auto', sva_B='auto')
    assert proc.returncode == 0, proc.stdout

    summary_df = _read_batch_summary(out_dir, fixture['species_tag'])
    assert summary_df.loc[0, 'corrected_run_count'] == 0
    assert str(summary_df.loc[0, 'skip_reason']) != ''
    metadata_df = _read_species_metadata(out_dir, fixture['species_tag'])
    assert set(metadata_df['batch_corrected'].astype(str)) == {'no'}


def test_finalize_r_accepts_explicit_nsv_zero(require_finalize_r_runtime, tmp_path):
    fixture = _write_finalize_fixture(
        tmp_path=tmp_path,
        sample_groups=['A', 'A', 'B', 'B'],
        bioprojects=['BP1', 'BP1', 'BP2', 'BP2'],
        species='Nsvzero example',
    )
    proc, out_dir = _run_finalize_r(tmp_path=tmp_path, fixture=fixture, sva_nsv='0', sva_B='auto')
    assert proc.returncode == 0, proc.stdout

    summary_df = _read_batch_summary(out_dir, fixture['species_tag'])
    assert int(summary_df.loc[0, 'resolved_sva_nsv']) == 0
    assert summary_df.loc[0, 'skip_reason'] == 'sva_nsv_zero'


def test_finalize_r_clamps_large_manual_nsv(require_finalize_r_runtime, tmp_path):
    fixture = _write_finalize_fixture(
        tmp_path=tmp_path,
        sample_groups=['A', 'A', 'B', 'B'],
        bioprojects=['BP1', 'BP1', 'BP2', 'BP2'],
        species='Clampnsv example',
    )
    proc, out_dir = _run_finalize_r(tmp_path=tmp_path, fixture=fixture, sva_nsv='99', sva_B='20')
    assert proc.returncode == 0, proc.stdout

    summary_df = _read_batch_summary(out_dir, fixture['species_tag'])
    assert int(summary_df.loc[0, 'resolved_sva_nsv']) < 99


def test_finalize_r_sva_plots_work_when_optional_metadata_columns_missing(require_finalize_r_runtime, tmp_path):
    fixture = _write_finalize_fixture(
        tmp_path=tmp_path,
        sample_groups=['A', 'A', 'B', 'B'],
        bioprojects=['BP1', 'BP1', 'BP2', 'BP2'],
        species='Missingcolumns example',
    )
    metadata = pandas.read_csv(fixture['metadata_path'], sep='\t')
    metadata = metadata.drop(columns=['lib_layout', 'lib_selection', 'instrument', 'total_spots', 'total_bases'])
    metadata.to_csv(fixture['metadata_path'], sep='\t', index=False)

    proc, out_dir = _run_finalize_r(tmp_path=tmp_path, fixture=fixture, sva_nsv='auto', sva_B='auto')
    assert proc.returncode == 0, proc.stdout

    species_tag = fixture['species_tag']
    plot_dir = out_dir / 'curate' / species_tag / 'plots'
    assert (plot_dir / '{}.batch_compare.sva.pdf'.format(species_tag)).exists()
