import shutil
import subprocess
from types import SimpleNamespace
from pathlib import Path

import pandas
import pytest

from amalgkit.per_species_tables import generate_per_species_tables
from amalgkit.r_config import temporary_r_config


REQUIRED_R_PACKAGES = [
    'Rtsne',
    'ggplot2',
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
def require_r_runtime():
    if not _has_rscript_and_required_packages():
        pytest.skip('Rscript or required R packages are not available for per-species table integration tests.')


def _repo_root():
    return Path(__file__).resolve().parents[1]


def _write_tiny_curate_fixture(tmp_path):
    runs = ['RUN1', 'RUN2', 'RUN3', 'RUN4']
    genes = ['G{:03d}'.format(i) for i in range(1, 21)]
    species = 'Testus example'

    metadata = pandas.DataFrame({
        'run': runs,
        'scientific_name': [species] * 4,
        'sample_group': ['A', 'A', 'B', 'B'],
        'bioproject': ['BP1'] * 4,
        'exclusion': ['no'] * 4,
        'mapping_rate': [100.0] * 4,
    })
    metadata_path = tmp_path / 'metadata.tsv'
    metadata.to_csv(metadata_path, sep='\t', index=False)

    # Construct intentionally discordant profiles so high correlation thresholds can remove all runs.
    counts = pandas.DataFrame({
        'target_id': genes,
        'RUN1': list(range(1, 21)),
        'RUN2': list(range(20, 0, -1)),
        'RUN3': [1 if i % 2 == 0 else 100 for i in range(20)],
        'RUN4': [100 if i % 2 == 0 else 1 for i in range(20)],
    })
    count_path = tmp_path / 'counts.tsv'
    counts.to_csv(count_path, sep='\t', index=False)

    eff_length = pandas.DataFrame({'target_id': genes})
    for run in runs:
        eff_length[run] = 1000
    eff_length_path = tmp_path / 'eff_length.tsv'
    eff_length.to_csv(eff_length_path, sep='\t', index=False)

    return {
        'species': species,
        'metadata_path': metadata_path,
        'count_path': count_path,
        'eff_length_path': eff_length_path,
    }


def _run_curate_r(
    tmp_path,
    fixture,
    one_outlier_per_iteration,
    correlation_threshold,
    selected_sample_groups='A|B',
    sample_group_colors='DEFAULT',
):
    repo = _repo_root()
    out_dir = tmp_path / 'out'
    out_dir.mkdir(parents=True, exist_ok=True)
    config_map = {
        'est_counts_path': str(fixture['count_path']),
        'metadata_path': str(fixture['metadata_path']),
        'out_dir': str(out_dir),
        'eff_length_path': str(fixture['eff_length_path']),
        'dist_method': 'pearson',
        'mapping_rate_cutoff': '0.2',
        'min_dif': '0',
        'plot_intermediate': '0',
        'selected_sample_groups': selected_sample_groups,
        'sample_group_colors': sample_group_colors,
        'transform_method': 'log2p1-fpkm',
        'one_outlier_per_iteration': '1' if one_outlier_per_iteration else '0',
        'correlation_threshold': str(correlation_threshold),
        'batch_effect_alg': 'no',
        'clip_negative': '1',
        'maintain_zero': '1',
        'r_util_path': str(repo / 'amalgkit' / 'util.r'),
        'skip_curation_flag': '0',
        'outlier_method': 'legacy',
        'robust_margin_threshold': '0',
        'robust_z_threshold': '-2.5',
        'disable_auto_outlier_filter_flag': '0',
    }
    with temporary_r_config(config_map, prefix='test_curate_r_') as config_path:
        cmd = ['Rscript', str(repo / 'amalgkit' / 'prepare_tables.r'), config_path]
        return subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True), out_dir


def _write_sample_group_drop_fixture(tmp_path):
    runs = ['ROOT1', 'FLOWER1', 'LEAF1']
    genes = ['G{:03d}'.format(i) for i in range(1, 21)]
    species = 'Dropus samplegroup'

    metadata = pandas.DataFrame({
        'run': runs,
        'scientific_name': [species] * 3,
        'sample_group': ['root', 'flower', 'leaf'],
        'bioproject': ['BP1', 'BP1', 'BP2'],
        'exclusion': ['manual_removal', 'no', 'no'],
        'mapping_rate': [100.0] * 3,
    })
    metadata_path = tmp_path / 'metadata_drop.tsv'
    metadata.to_csv(metadata_path, sep='\t', index=False)

    counts = pandas.DataFrame({
        'target_id': genes,
        'ROOT1': list(range(1, 21)),
        'FLOWER1': list(range(21, 41)),
        'LEAF1': list(range(41, 61)),
    })
    count_path = tmp_path / 'counts_drop.tsv'
    counts.to_csv(count_path, sep='\t', index=False)

    eff_length = pandas.DataFrame({'target_id': genes})
    for run in runs:
        eff_length[run] = 1000
    eff_length_path = tmp_path / 'eff_length_drop.tsv'
    eff_length.to_csv(eff_length_path, sep='\t', index=False)

    return {
        'species': species,
        'metadata_path': metadata_path,
        'count_path': count_path,
        'eff_length_path': eff_length_path,
    }


def test_curate_r_writes_round_and_final_summaries(require_r_runtime, tmp_path):
    fixture = _write_tiny_curate_fixture(tmp_path)
    proc, out_dir = _run_curate_r(
        tmp_path=tmp_path,
        fixture=fixture,
        one_outlier_per_iteration=False,
        correlation_threshold=0.3,
    )
    assert proc.returncode == 0, proc.stdout

    species_tag = fixture['species'].replace(' ', '_')
    tables_dir = out_dir / 'per_species' / species_tag / 'tables'
    round_path = tables_dir / '{}.no.curation_round_summary.tsv'.format(species_tag)
    final_path = tables_dir / '{}.no.curation_final_summary.tsv'.format(species_tag)
    assert round_path.exists(), proc.stdout
    assert final_path.exists(), proc.stdout

    round_df = pandas.read_csv(round_path, sep='\t')
    final_df = pandas.read_csv(final_path, sep='\t')
    assert {'step', 'round', 'num_runs_before', 'num_runs_after', 'num_runs_removed', 'removed_runs'}.issubset(set(round_df.columns))
    assert {'mapping_rate_zero', 'mapping_rate_cutoff'}.issubset(set(round_df['step'].tolist()))
    assert final_df.loc[0, 'num_runs_after_sample_group_filter'] == 4
    assert final_df.loc[0, 'num_runs_final_kept'] + final_df.loc[0, 'num_runs_final_excluded'] == 4
    assert final_df.loc[0, 'total_runtime_sec'] >= 0


def test_curate_r_handles_all_samples_removed_in_iterative_round(require_r_runtime, tmp_path):
    fixture = _write_tiny_curate_fixture(tmp_path)
    proc, out_dir = _run_curate_r(
        tmp_path=tmp_path,
        fixture=fixture,
        one_outlier_per_iteration=False,
        correlation_threshold=1.1,
    )
    assert proc.returncode == 0, proc.stdout
    assert 'Execution halted' not in proc.stdout
    assert 'No sample is available. Outlier removal will be skipped.' in proc.stdout

    species_tag = fixture['species'].replace(' ', '_')
    final_path = out_dir / 'per_species' / species_tag / 'tables' / '{}.no.curation_final_summary.tsv'.format(species_tag)
    final_df = pandas.read_csv(final_path, sep='\t')
    assert final_df.loc[0, 'num_runs_final_kept'] == 0
    assert final_df.loc[0, 'num_runs_final_excluded'] == 4


def test_curate_r_handles_filtered_out_sample_group_with_default_colors(require_r_runtime, tmp_path):
    fixture = _write_sample_group_drop_fixture(tmp_path)
    proc, out_dir = _run_curate_r(
        tmp_path=tmp_path,
        fixture=fixture,
        one_outlier_per_iteration=False,
        correlation_threshold=0.3,
        selected_sample_groups='root|flower|leaf',
        sample_group_colors='DEFAULT',
    )
    assert proc.returncode == 0, proc.stdout
    assert 'Execution halted' not in proc.stdout
    assert 'missing value where TRUE/FALSE needed' not in proc.stdout

    species_tag = fixture['species'].replace(' ', '_')
    final_path = out_dir / 'per_species' / species_tag / 'tables' / '{}.no.curation_final_summary.tsv'.format(species_tag)
    final_df = pandas.read_csv(final_path, sep='\t')
    assert final_df.loc[0, 'num_runs_after_sample_group_filter'] == 3
    assert final_df.loc[0, 'num_runs_final_kept'] == 2


def test_generate_per_species_tables_exits_nonzero_when_species_input_files_are_missing(tmp_path, monkeypatch):
    input_dir = tmp_path / 'input_merge'
    input_dir.mkdir(parents=True, exist_ok=True)
    metadata_path = tmp_path / 'metadata.tsv'
    pandas.DataFrame({
        'run': ['R1'],
        'scientific_name': ['Missing Files Species'],
        'sample_group': ['A'],
        'exclusion': ['no'],
    }).to_csv(metadata_path, sep='\t', index=False)

    out_dir = tmp_path / 'out'
    args = SimpleNamespace(
        input_dir=str(input_dir),
        out_dir=str(out_dir),
        norm='log2p1-fpkm',
        redo=False,
        sample_group=None,
        metadata=str(metadata_path),
        dist_method='pearson',
        mapping_rate=0.2,
        correlation_threshold=0.3,
        plot_intermediate=False,
        sample_group_color='DEFAULT',
        one_outlier_per_iter=False,
        batch_effect_alg='no',
        clip_negative=True,
        maintain_zero=True,
        skip_curation=False,
    )

    monkeypatch.setattr('amalgkit.per_species_tables.check_rscript', lambda: None)
    with pytest.raises(RuntimeError) as e:
        generate_per_species_tables(args)
    assert 'Per-species table generation failed for 1/1 species' in str(e.value)
