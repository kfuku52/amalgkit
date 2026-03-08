from types import SimpleNamespace

import pandas
import pytest

from amalgkit.command_context import PerSpeciesTableContext
from amalgkit.per_species_tables import generate_per_species_tables
from amalgkit.util import Metadata


def _write_prepare_fixture(tmp_path, species='Testus example'):
    runs = ['RUN1', 'RUN2', 'RUN3', 'RUN4']
    genes = ['G{:03d}'.format(i) for i in range(1, 21)]
    species_tag = species.replace(' ', '_')

    input_dir = tmp_path / 'input'
    species_dir = input_dir / species_tag
    species_dir.mkdir(parents=True, exist_ok=True)

    metadata_df = pandas.DataFrame(
        {
            'run': runs,
            'scientific_name': [species] * 4,
            'sample_group': ['A', 'A', 'B', 'B'],
            'bioproject': ['BP1'] * 4,
            'exclusion': ['no'] * 4,
            'mapping_rate': [100.0, 95.0, 100.0, 95.0],
        }
    )
    metadata = Metadata.from_DataFrame(metadata_df)

    counts_df = pandas.DataFrame(
        {
            'target_id': genes,
            'RUN1': list(range(1, 21)),
            'RUN2': list(range(2, 22)),
            'RUN3': list(range(21, 41)),
            'RUN4': list(range(22, 42)),
        }
    )
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
    }


def _write_sample_group_drop_fixture(tmp_path):
    runs = ['ROOT1', 'FLOWER1', 'LEAF1']
    genes = ['G{:03d}'.format(i) for i in range(1, 21)]
    species = 'Dropus samplegroup'
    species_tag = species.replace(' ', '_')

    input_dir = tmp_path / 'input'
    species_dir = input_dir / species_tag
    species_dir.mkdir(parents=True, exist_ok=True)

    metadata_df = pandas.DataFrame(
        {
            'run': runs,
            'scientific_name': [species] * 3,
            'sample_group': ['root', 'flower', 'leaf'],
            'bioproject': ['BP1', 'BP1', 'BP2'],
            'exclusion': ['manual_removal', 'no', 'no'],
            'mapping_rate': [100.0] * 3,
        }
    )
    metadata = Metadata.from_DataFrame(metadata_df)

    counts_df = pandas.DataFrame(
        {
            'target_id': genes,
            'ROOT1': list(range(1, 21)),
            'FLOWER1': list(range(21, 41)),
            'LEAF1': list(range(41, 61)),
        }
    )
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
    }


def _build_prepare_args(tmp_path):
    return SimpleNamespace(
        out_dir=str(tmp_path / 'out'),
        redo=False,
        metadata='inferred',
        input_dir='inferred',
        sample_group=None,
        sample_group_color='DEFAULT',
        batch=None,
        threads='auto',
        internal_jobs=1,
        internal_cpu_budget='auto',
        dist_method='pearson',
        mapping_rate=0.0,
        correlation_threshold=0.3,
        plot_intermediate=False,
        one_outlier_per_iter=False,
        norm='log2p1-fpkm',
        clip_negative=True,
        maintain_zero=True,
        batch_effect_alg='no',
        skip_curation=False,
        margin_threshold=0.0,
        robust_z_threshold=-2.5,
        worker_mode='prepare_tables',
    )


def test_prepare_tables_python_writes_round_and_final_summaries(tmp_path):
    fixture = _write_prepare_fixture(tmp_path)
    args = _build_prepare_args(tmp_path)

    generate_per_species_tables(
        args,
        context=PerSpeciesTableContext(metadata=fixture['metadata'], input_dir=fixture['input_dir']),
    )

    species_tag = fixture['species_tag']
    tables_dir = tmp_path / 'out' / 'per_species' / species_tag / 'tables'
    round_path = tables_dir / '{}.no.curation_round_summary.tsv'.format(species_tag)
    final_path = tables_dir / '{}.no.curation_final_summary.tsv'.format(species_tag)
    tau_path = tables_dir / '{}.no.tau.tsv'.format(species_tag)
    corr_path = tables_dir / '{}.no.correlation_statistics.tsv'.format(species_tag)
    batch_path = tables_dir / '{}.no.batch_effect_summary.tsv'.format(species_tag)

    assert round_path.exists()
    assert final_path.exists()
    assert tau_path.exists()
    assert corr_path.exists()
    assert batch_path.exists()

    round_df = pandas.read_csv(round_path, sep='\t')
    final_df = pandas.read_csv(final_path, sep='\t')
    assert {'step', 'round', 'num_runs_before', 'num_runs_after', 'num_runs_removed', 'removed_runs'}.issubset(set(round_df.columns))
    assert final_df.loc[0, 'num_runs_after_sample_group_filter'] == 4
    assert final_df.loc[0, 'num_runs_final_kept'] + final_df.loc[0, 'num_runs_final_excluded'] == 4
    assert final_df.loc[0, 'total_runtime_sec'] >= 0


def test_prepare_tables_python_handles_sample_group_drop_with_default_colors(tmp_path):
    fixture = _write_sample_group_drop_fixture(tmp_path)
    args = _build_prepare_args(tmp_path)
    args.sample_group = 'root|flower|leaf'

    generate_per_species_tables(
        args,
        context=PerSpeciesTableContext(metadata=fixture['metadata'], input_dir=fixture['input_dir']),
    )

    species_tag = fixture['species_tag']
    final_path = tmp_path / 'out' / 'per_species' / species_tag / 'tables' / '{}.no.curation_final_summary.tsv'.format(species_tag)
    final_df = pandas.read_csv(final_path, sep='\t')
    assert final_df.loc[0, 'num_runs_after_sample_group_filter'] == 3
    assert final_df.loc[0, 'num_runs_final_kept'] >= 1


def test_generate_per_species_tables_exits_nonzero_when_species_input_files_are_missing(tmp_path):
    input_dir = tmp_path / 'input_merge'
    input_dir.mkdir(parents=True, exist_ok=True)
    metadata = Metadata.from_DataFrame(
        pandas.DataFrame(
            {
                'run': ['R1'],
                'scientific_name': ['Missing Files Species'],
                'sample_group': ['A'],
                'exclusion': ['no'],
            }
        )
    )

    args = _build_prepare_args(tmp_path)
    with pytest.raises(RuntimeError) as exc_info:
        generate_per_species_tables(
            args,
            context=PerSpeciesTableContext(metadata=metadata, input_dir=str(input_dir)),
        )
    assert 'Per-species table generation failed for 1/1 species' in str(exc_info.value)
