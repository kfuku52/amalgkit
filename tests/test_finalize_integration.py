import pandas
import pytest

from amalgkit.command_context import PerSpeciesTableContext
from amalgkit.per_species_tables import generate_per_species_tables
from amalgkit.util import Metadata
from tests.support.per_species import build_per_species_args


pytestmark = pytest.mark.integration


def _write_finalize_fixture(tmp_path, sample_groups, bioprojects, species='Finalizus example', include_optional_columns=True):
    runs = ['RUN{:02d}'.format(i + 1) for i in range(len(sample_groups))]
    genes = ['G{:03d}'.format(i + 1) for i in range(50)]
    species_tag = species.replace(' ', '_')

    input_dir = tmp_path / 'input'
    species_dir = input_dir / species_tag
    species_dir.mkdir(parents=True, exist_ok=True)

    metadata_df = pandas.DataFrame(
        {
            'run': runs,
            'scientific_name': [species] * len(runs),
            'sample_group': sample_groups,
            'bioproject': bioprojects,
            'exclusion': ['no'] * len(runs),
            'mapping_rate': [100.0] * len(runs),
        }
    )
    if include_optional_columns:
        metadata_df['lib_layout'] = ['PAIRED'] * len(runs)
        metadata_df['lib_selection'] = ['cDNA'] * len(runs)
        metadata_df['instrument'] = ['Illumina'] * len(runs)
        metadata_df['total_spots'] = [1_000_000 + i * 10_000 for i in range(len(runs))]
        metadata_df['total_bases'] = [150_000_000 + i * 1_000_000 for i in range(len(runs))]
    metadata = Metadata.from_DataFrame(metadata_df)

    counts_df = pandas.DataFrame({'target_id': genes})
    for run_id, sample_group in zip(runs, sample_groups):
        if sample_group == 'A':
            values = [100 + i for i in range(25)] + [5 + i for i in range(25)]
        else:
            values = [5 + i for i in range(25)] + [100 + i for i in range(25)]
        counts_df[run_id] = values
    eff_length_df = pandas.DataFrame({'target_id': genes})
    for run in runs:
        eff_length_df[run] = 1000

    counts_df.to_csv(species_dir / '{}_est_counts.tsv'.format(species_tag), sep='\t', index=False)
    eff_length_df.to_csv(species_dir / '{}_eff_length.tsv'.format(species_tag), sep='\t', index=False)
    return {
        'metadata': metadata,
        'input_dir': str(input_dir),
        'species': species,
        'species_tag': species_tag,
        'sample_groups': sample_groups,
    }


def _inject_latent_batch_signal(input_dir, species_tag):
    counts_path = input_dir / species_tag / '{}_est_counts.tsv'.format(species_tag)
    counts_df = pandas.read_csv(counts_path, sep='\t')
    run_columns = [column for column in counts_df.columns if column != 'target_id']
    values = {
        'RUN01': [120, 118, 122, 119, 12, 14, 15, 13, 80, 82, 78, 81] + [10 + i for i in range(38)],
        'RUN02': [121, 117, 123, 120, 13, 15, 14, 12, 32, 30, 34, 31] + [10 + i for i in range(38)],
        'RUN03': [14, 13, 12, 11, 131, 129, 133, 130, 79, 83, 77, 80] + [50 + i for i in range(38)],
        'RUN04': [13, 12, 11, 10, 132, 130, 134, 131, 33, 31, 35, 32] + [50 + i for i in range(38)],
    }
    for run_id in run_columns:
        counts_df.loc[:, run_id] = values[run_id]
    counts_df.to_csv(counts_path, sep='\t', index=False)


def _build_finalize_args(tmp_path, **overrides):
    return build_per_species_args(tmp_path, **{
        'batch_effect_alg': 'sva',
        'disable_auto_outlier_filter': True,
        'worker_mode': 'finalize',
        'seed': '7',
        'sva_B_auto_max': 80,
        **overrides,
    })


def _run_finalize_python(tmp_path, fixture, **overrides):
    args = _build_finalize_args(tmp_path, **overrides)
    generate_per_species_tables(
        args,
        context=PerSpeciesTableContext(metadata=fixture['metadata'], input_dir=fixture['input_dir']),
    )
    out_dir = tmp_path / 'out'
    return out_dir


def _read_batch_summary(out_dir, species_tag, batch_effect_alg='sva'):
    path = out_dir / 'per_species' / species_tag / 'tables' / '{}.{}.batch_effect_summary.tsv'.format(species_tag, batch_effect_alg)
    return pandas.read_csv(path, sep='\t')


def _read_species_metadata(out_dir, species_tag):
    path = out_dir / 'per_species' / species_tag / 'tables' / '{}.metadata.tsv'.format(species_tag)
    return pandas.read_csv(path, sep='\t')


def _read_corrected_tc(out_dir, species_tag, batch_effect_alg='sva'):
    path = out_dir / 'per_species' / species_tag / 'tables' / '{}.{}.tc.tsv'.format(species_tag, batch_effect_alg)
    return pandas.read_csv(path, sep='\t', index_col=0)


def test_finalize_python_sva_handles_confounded_design_and_writes_batch_summary(tmp_path):
    fixture = _write_finalize_fixture(
        tmp_path=tmp_path,
        sample_groups=['A', 'A', 'B', 'B'],
        bioprojects=['BP1', 'BP1', 'BP2', 'BP2'],
        species='Confoundus example',
    )
    out_dir = _run_finalize_python(tmp_path=tmp_path, fixture=fixture, sva_nsv='auto', sva_B='auto')

    summary_df = _read_batch_summary(out_dir, fixture['species_tag'])
    metadata_df = _read_species_metadata(out_dir, fixture['species_tag'])
    assert {'resolved_sva_nsv', 'resolved_sva_B', 'skip_reason', 'corrected_run_count'}.issubset(set(summary_df.columns))
    assert {'batch_corrected', 'batch_alg_used'}.issubset(set(metadata_df.columns))


def test_finalize_python_sva_single_sample_group_marks_uncorrected(tmp_path):
    fixture = _write_finalize_fixture(
        tmp_path=tmp_path,
        sample_groups=['A', 'A', 'A', 'A'],
        bioprojects=['BP1', 'BP1', 'BP2', 'BP2'],
        species='Singlesamplegroup example',
    )
    out_dir = _run_finalize_python(tmp_path=tmp_path, fixture=fixture, sva_nsv='auto', sva_B='auto')

    summary_df = _read_batch_summary(out_dir, fixture['species_tag'])
    metadata_df = _read_species_metadata(out_dir, fixture['species_tag'])
    assert summary_df.loc[0, 'corrected_run_count'] == 0
    assert str(summary_df.loc[0, 'skip_reason']) != ''
    assert set(metadata_df['batch_corrected'].astype(str)) == {'no'}


def test_finalize_python_sva_accepts_explicit_nsv_zero(tmp_path):
    fixture = _write_finalize_fixture(
        tmp_path=tmp_path,
        sample_groups=['A', 'A', 'B', 'B'],
        bioprojects=['BP1', 'BP1', 'BP2', 'BP2'],
        species='Nsvzero example',
    )
    out_dir = _run_finalize_python(tmp_path=tmp_path, fixture=fixture, sva_nsv='0', sva_B='5')

    summary_df = _read_batch_summary(out_dir, fixture['species_tag'])
    metadata_df = _read_species_metadata(out_dir, fixture['species_tag'])
    assert int(summary_df.loc[0, 'resolved_sva_nsv']) == 0
    assert summary_df.loc[0, 'skip_reason'] == 'sva_nsv_zero'
    assert set(metadata_df['batch_corrected'].astype(str)) == {'no'}


def test_finalize_python_sva_supports_positive_manual_nsv(tmp_path):
    fixture = _write_finalize_fixture(
        tmp_path=tmp_path,
        sample_groups=['A', 'A', 'B', 'B'],
        bioprojects=['BP1', 'BP1', 'BP2', 'BP2'],
        species='Pythonmanualpositive example',
    )
    out_dir = _run_finalize_python(tmp_path=tmp_path, fixture=fixture, sva_nsv='1', sva_B='5')

    summary_df = _read_batch_summary(out_dir, fixture['species_tag'])
    metadata_df = _read_species_metadata(out_dir, fixture['species_tag'])
    corrected_tc = _read_corrected_tc(out_dir, fixture['species_tag'])
    assert int(summary_df.loc[0, 'resolved_sva_nsv']) == 1
    assert int(summary_df.loc[0, 'resolved_sva_B']) == 5
    assert int(summary_df.loc[0, 'corrected_run_count']) == 4
    assert set(metadata_df['batch_corrected'].astype(str)) == {'yes'}
    assert corrected_tc.shape[1] == 4


@pytest.mark.slow
def test_finalize_python_sva_plots_work_when_optional_metadata_columns_missing(tmp_path):
    fixture = _write_finalize_fixture(
        tmp_path=tmp_path,
        sample_groups=['A', 'A', 'B', 'B'],
        bioprojects=['BP1', 'BP1', 'BP2', 'BP2'],
        species='Missingcolumns example',
        include_optional_columns=False,
    )
    out_dir = _run_finalize_python(tmp_path=tmp_path, fixture=fixture, sva_nsv='1', sva_B='5')

    plot_dir = out_dir / 'per_species' / fixture['species_tag'] / 'plots'
    assert (plot_dir / '{}.before_after.sva.pdf'.format(fixture['species_tag'])).exists()


@pytest.mark.optional_dependency
def test_finalize_python_combatseq_runs_end_to_end(tmp_path):
    pytest.importorskip('inmoose.pycombat')
    fixture = _write_finalize_fixture(
        tmp_path=tmp_path,
        sample_groups=['A', 'B', 'A', 'B'],
        bioprojects=['BP1', 'BP1', 'BP2', 'BP2'],
        species='Combatseq example',
    )
    out_dir = _run_finalize_python(tmp_path=tmp_path, fixture=fixture, batch_effect_alg='combatseq')

    summary_df = _read_batch_summary(out_dir, fixture['species_tag'], batch_effect_alg='combatseq')
    corrected_tc = _read_corrected_tc(out_dir, fixture['species_tag'], batch_effect_alg='combatseq')
    assert int(summary_df.loc[0, 'corrected_run_count']) == 4
    assert corrected_tc.shape[1] == 4


def test_finalize_python_ruvseq_runs_end_to_end(tmp_path):
    fixture = _write_finalize_fixture(
        tmp_path=tmp_path,
        sample_groups=['A', 'A', 'B', 'B'],
        bioprojects=['BP1', 'BP2', 'BP1', 'BP2'],
        species='Pythonruvseq example',
    )
    # Inject a genuine latent batch signal so a real unwanted-variation
    # component exists for RUVseq to estimate. Without this, the design fit
    # leaves only float-noise residuals and the correct resolution is k=0
    # (the regression test in test_batch_effect_ruvseq.py pins that case).
    _inject_latent_batch_signal(tmp_path / 'input', fixture['species_tag'])
    out_dir = _run_finalize_python(
        tmp_path=tmp_path,
        fixture=fixture,
        batch_effect_alg='ruvseq',
        ruvseq_k='1',
    )

    summary_df = _read_batch_summary(out_dir, fixture['species_tag'], batch_effect_alg='ruvseq')
    metadata_df = _read_species_metadata(out_dir, fixture['species_tag'])
    corrected_tc = _read_corrected_tc(out_dir, fixture['species_tag'], batch_effect_alg='ruvseq')
    assert int(summary_df.loc[0, 'resolved_ruv_k']) == 1
    assert int(summary_df.loc[0, 'corrected_run_count']) == 4
    assert set(metadata_df['batch_corrected'].astype(str)) == {'yes'}
    assert corrected_tc.shape[1] == 4


def test_finalize_python_ruvseq_single_group_design_failure(tmp_path):
    fixture = _write_finalize_fixture(
        tmp_path=tmp_path,
        sample_groups=['A', 'A', 'A', 'A'],
        bioprojects=['BP1', 'BP2', 'BP1', 'BP2'],
        species='Ruvseq single-group example',
    )
    out_dir = _run_finalize_python(
        tmp_path=tmp_path,
        fixture=fixture,
        batch_effect_alg='ruvseq',
        ruvseq_k='auto',
    )

    summary_df = _read_batch_summary(out_dir, fixture['species_tag'], batch_effect_alg='ruvseq')
    corrected_tc = _read_corrected_tc(out_dir, fixture['species_tag'], batch_effect_alg='ruvseq')
    assert str(summary_df.loc[0, 'skip_reason']) == 'ruvseq_design_failed'
    assert pandas.isna(summary_df.loc[0, 'resolved_ruv_k'])
    assert corrected_tc.shape[1] == 4


def test_finalize_python_latent_glm_runs_end_to_end_and_stays_nonnegative(tmp_path):
    fixture = _write_finalize_fixture(
        tmp_path=tmp_path,
        sample_groups=['A', 'A', 'B', 'B'],
        bioprojects=['BP1', 'BP2', 'BP1', 'BP2'],
        species='Latentglm example',
    )
    _inject_latent_batch_signal(tmp_path / 'input', fixture['species_tag'])
    out_dir = _run_finalize_python(
        tmp_path=tmp_path,
        fixture=fixture,
        batch_effect_alg='latent_glm',
        latent_family='nb',
        latent_k='1',
        latent_k_max=3,
        latent_max_iter=50,
        latent_tol=1e-6,
    )

    summary_df = _read_batch_summary(out_dir, fixture['species_tag'], batch_effect_alg='latent_glm')
    metadata_df = _read_species_metadata(out_dir, fixture['species_tag'])
    corrected_tc = _read_corrected_tc(out_dir, fixture['species_tag'], batch_effect_alg='latent_glm')
    assert int(summary_df.loc[0, 'resolved_latent_k']) == 1
    assert summary_df.loc[0, 'latent_family'] == 'nb'
    assert int(summary_df.loc[0, 'corrected_run_count']) == 4
    assert set(metadata_df['batch_corrected'].astype(str)) == {'yes'}
    assert (corrected_tc.to_numpy(dtype=float) >= 0.0).all()
