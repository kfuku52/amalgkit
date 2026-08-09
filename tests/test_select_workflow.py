import os
import subprocess
import sys
import pandas
import pytest

from types import SimpleNamespace
from pathlib import Path

from amalgkit.select import (
    select_main,
)
from amalgkit.util import Metadata

SELECT_RULE_COLUMNS = [
    'rule_id',
    'enabled',
    'stage',
    'priority',
    'columns',
    'pattern',
    'action',
    'target_column',
    'outcome',
    'scope_column',
    'scope_mode',
    'stop_on_match',
    'parameter_name',
    'parameter_value',
    'note',
]


def write_select_rules(path, rows):
    defaults = {
        'rule_id': '',
        'enabled': 'yes',
        'stage': '',
        'priority': '0',
        'columns': '',
        'pattern': '',
        'action': '',
        'target_column': '',
        'outcome': '',
        'scope_column': '',
        'scope_mode': '',
        'stop_on_match': 'yes',
        'parameter_name': '',
        'parameter_value': '',
        'note': '',
    }
    normalized_rows = []
    for row in rows:
        normalized = defaults.copy()
        normalized.update(row)
        normalized_rows.append(normalized)
    pandas.DataFrame(normalized_rows, columns=SELECT_RULE_COLUMNS).to_csv(path, sep='\t', index=False)


def build_test_validate_rows():
    return [
        {
            'rule_id': 'validate_hint_leaf',
            'stage': 'validate',
            'priority': '10',
            'pattern': r'\b(?:leaf|leaves|leaflet|leaflets|lamina|laminae|leaf[\s_-]?blade|leaf[\s_-]?blades|seedling[\s_-]?leaf|seedling[\s_-]?leaves)\b',
            'action': 'hint_organ',
            'outcome': 'leaf',
        },
        {
            'rule_id': 'validate_hint_root',
            'stage': 'validate',
            'priority': '11',
            'pattern': r'\b(?:root|roots|taproot|taproots|primary root|primary roots|lateral root|lateral roots|radicle|radicles|root[\s_-]?tip|root[\s_-]?tips|seedling[\s_-]?root|seedling[\s_-]?roots|hairy root|hairy roots)\b',
            'action': 'hint_organ',
            'outcome': 'root',
        },
        {
            'rule_id': 'validate_hint_flower',
            'stage': 'validate',
            'priority': '12',
            'pattern': r'\b(?:flower|flowers|floral|inflorescence|inflorescences|spike|spikes|panicle|panicles|catkin|catkins)\b',
            'action': 'hint_organ',
            'outcome': 'flower',
        },
        {
            'rule_id': 'validate_hint_review',
            'stage': 'validate',
            'priority': '13',
            'pattern': r'\b(?:petal|petals|corolla|anther|anthers|ovary|ovaries|pistil|pistils|style|styles|stigma|stigmas|bract|bracts|bud|buds|petiole|petioles|root hair|root hairs|apex|meristem|fruit|fruits|seed|seeds|spine|spines)\b',
            'action': 'hint_review',
            'outcome': 'review',
        },
        {
            'rule_id': 'validate_ignore_safe_metadata',
            'stage': 'validate',
            'priority': '14',
            'pattern': r'^(?:(?:biological|technical)\s+replicates?(?:\s+[A-Za-z0-9._-]+)?|replicates?(?:\s+[A-Za-z0-9._-]+)?|repeats?(?:\s+[A-Za-z0-9._-]+)?|rep(?:licate)?\s*[A-Za-z0-9._-]+|samples?(?:\s+[A-Za-z0-9._-]+)?|libraries?(?:\s+[A-Za-z0-9._-]+)?|libs?(?:\s+[A-Za-z0-9._-]+)?|lanes?(?:\s+[A-Za-z0-9._-]+)?|controls?(?:\s+[A-Za-z0-9._-]+)?|mock(?:\s+[A-Za-z0-9._-]+)?|(?:cold|heat|salt|drought|stress|treated|treatment|cultivar)(?:\s+[A-Za-z0-9._-]+)*)$',
            'action': 'ignore_segment',
            'outcome': 'ignore',
        },
    ]


def build_test_filter_dedup_rows():
    return [
        {
            'rule_id': 'filter_low_nspots',
            'stage': 'filter',
            'priority': '3000',
            'columns': 'total_spots',
            'action': 'exclude_if_lt_parameter',
            'target_column': 'exclusion',
            'outcome': 'low_nspots',
            'parameter_name': 'min_nspots',
        },
        {
            'rule_id': 'filter_missing_taxid',
            'stage': 'filter',
            'priority': '3010',
            'columns': 'taxid_{parameter}',
            'action': 'exclude_if_missing_selected_rank',
            'target_column': 'exclusion',
            'outcome': 'missing_taxid',
            'parameter_name': 'mark_missing_rank',
        },
        {
            'rule_id': 'filter_no_sample_group',
            'stage': 'filter',
            'priority': '3020',
            'columns': 'sample_group',
            'action': 'exclude_if_empty',
            'target_column': 'exclusion',
            'outcome': 'no_sample_group',
        },
        {
            'rule_id': 'dedup_redundant_biosample',
            'stage': 'dedup',
            'priority': '3030',
            'columns': 'bioproject,biosample',
            'action': 'exclude_redundant_by_best_numeric',
            'target_column': 'total_spots',
            'outcome': 'redundant_biosample',
            'parameter_name': 'mark_redundant_biosamples',
        },
    ]


# ---------------------------------------------------------------------------
# write_select_outputs (writes metadata and pivot tables)
# ---------------------------------------------------------------------------


class TestSelectMain:
    def test_rejects_out_dir_file_path(self, tmp_path, monkeypatch):
        out_path = tmp_path / 'out_path'
        out_path.write_text('not a directory')
        args = SimpleNamespace(
            out_dir=str(out_path),
            select_rules_tsv='inferred',
            metadata='inferred',
            sample_group=None,
            min_nspots=0,
            mark_missing_rank='none',
            mark_redundant_biosamples=False,
            max_sample=1,
        )
        monkeypatch.setattr(
            'amalgkit.select.read_select_config',
            lambda *_args, **_kwargs: (_ for _ in ()).throw(AssertionError('read_select_config should not be called')),
        )

        with pytest.raises(NotADirectoryError, match='Output path exists but is not a directory'):
            select_main(args)

    def test_missing_inferred_select_rules_suggests_dataset_commands(self, tmp_path):
        out_dir = tmp_path / 'out'
        args = SimpleNamespace(
            out_dir=str(out_dir),
            select_rules_tsv='inferred',
            metadata='inferred',
            sample_group=None,
            min_nspots=0,
            mark_missing_rank='none',
            mark_redundant_biosamples=False,
            max_sample=1,
        )

        with pytest.raises(FileNotFoundError) as exc:
            select_main(args)

        message = str(exc.value)
        assert 'select rules file not found' in message
        assert 'amalgkit dataset --name init --out_dir' in message
        assert 'amalgkit dataset --rule_set base --out_dir' in message

    def test_uses_runtime_copy_without_mutating_caller_args(self, tmp_path, monkeypatch):
        raw_out_dir = str(tmp_path / 'nested' / '..' / 'out')
        args = SimpleNamespace(
            out_dir=raw_out_dir,
            select_rules_tsv='inferred',
            metadata='inferred',
            sample_group=None,
            min_nspots=0,
            mark_missing_rank='none',
            mark_redundant_biosamples=False,
            max_sample=1,
        )
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'scientific_name': ['Species A'],
            'sample_group': ['brain'],
            'bioproject': ['PRJ1'],
            'biosample': ['SAM1'],
            'total_spots': [100],
            'exclusion': ['no'],
        }))
        observed = {}

        def fake_read_select_config(path_select_rules):
            observed['select_rules_tsv'] = path_select_rules
            return {
                'rules': [{'stage': 'normalize', 'columns': ['sample_group'], 'regex': None}],
                'parameters': {
                    'min_nspots': 1000000,
                    'max_sample': 1,
                    'mark_missing_rank': 'none',
                    'mark_redundant_biosamples': False,
                },
            }

        def fake_load_metadata(runtime_args):
            observed['load_out_dir'] = runtime_args.out_dir
            return metadata

        def fake_prepare_select_metadata(current_metadata, select_rules):
            observed['prepare_rules'] = select_rules
            return current_metadata

        def fake_apply_select_filters(current_metadata, runtime_args, dir_config):
            observed['filter_out_dir'] = runtime_args.out_dir
            observed['filter_rules'] = dir_config
            return current_metadata

        def fake_write_select_outputs(**kwargs):
            observed['write_table'] = kwargs['path_metadata_table']
            observed['preserve_existing_original'] = kwargs['preserve_existing_original']

        monkeypatch.setattr('amalgkit.select.read_select_config', fake_read_select_config)
        monkeypatch.setattr('amalgkit.select.load_metadata', fake_load_metadata)
        monkeypatch.setattr('amalgkit.select.prepare_select_metadata', fake_prepare_select_metadata)
        monkeypatch.setattr('amalgkit.select.apply_select_filters', fake_apply_select_filters)
        monkeypatch.setattr('amalgkit.select.write_select_outputs', fake_write_select_outputs)

        select_main(args)

        normalized_out_dir = os.path.realpath(raw_out_dir)
        assert observed['select_rules_tsv'] == os.path.join(normalized_out_dir, 'select_rules.tsv')
        assert observed['load_out_dir'] == normalized_out_dir
        assert observed['filter_out_dir'] == normalized_out_dir
        assert observed['prepare_rules'] == observed['filter_rules']
        assert observed['write_table'] == os.path.join(normalized_out_dir, 'metadata', 'metadata.tsv')
        assert observed['preserve_existing_original'] is False
        assert args.out_dir == raw_out_dir

    def test_explicit_metadata_refreshes_existing_original(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'out'
        metadata_dir = out_dir / 'metadata'
        metadata_dir.mkdir(parents=True)
        (metadata_dir / 'metadata_original.tsv').write_text('old\n', encoding='utf-8')
        explicit_metadata = tmp_path / 'new.tsv'
        explicit_metadata.write_text('run\nSRR_NEW\n', encoding='utf-8')
        args = SimpleNamespace(
            out_dir=str(out_dir),
            select_rules_tsv='inferred',
            metadata=str(explicit_metadata),
            sample_group=None,
            min_nspots=0,
            mark_missing_rank='none',
            mark_redundant_biosamples=False,
            max_sample=1,
        )
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR_NEW'],
            'scientific_name': ['Species A'],
            'sample_group': ['leaf'],
            'bioproject': ['PRJ1'],
            'biosample': ['SAM1'],
            'total_spots': [100],
            'exclusion': ['no'],
        }))
        observed = {}

        monkeypatch.setattr(
            'amalgkit.select.read_select_config',
            lambda _path: {'rules': [], 'parameters': {}},
        )
        monkeypatch.setattr('amalgkit.select.load_metadata', lambda _args: metadata)
        monkeypatch.setattr('amalgkit.select.prepare_select_metadata', lambda value, _rules: value)
        monkeypatch.setattr('amalgkit.select.apply_select_filters', lambda value, _args, _rules: value)
        monkeypatch.setattr(
            'amalgkit.select.write_select_outputs',
            lambda **kwargs: observed.update(kwargs),
        )

        select_main(args)

        assert observed['preserve_existing_original'] is False
        assert observed['metadata_original_df']['run'].tolist() == ['SRR_NEW']

    def test_explicit_current_select_output_reuses_existing_original(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'out'
        metadata_dir = out_dir / 'metadata'
        metadata_dir.mkdir(parents=True)
        current_metadata = metadata_dir / 'metadata.tsv'
        original_metadata = metadata_dir / 'metadata_original.tsv'
        current_metadata.write_text('run\nSRR_FILTERED\n', encoding='utf-8')
        original_metadata.write_text('run\nSRR_ORIGINAL\n', encoding='utf-8')
        args = SimpleNamespace(
            out_dir=str(out_dir),
            select_rules_tsv='inferred',
            metadata=str(current_metadata),
            sample_group=None,
            min_nspots=0,
            mark_missing_rank='none',
            mark_redundant_biosamples=False,
            max_sample=1,
        )
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR_ORIGINAL'],
            'scientific_name': ['Species A'],
            'sample_group': ['leaf'],
            'total_spots': [100],
            'exclusion': ['no'],
        }))
        observed = {}

        monkeypatch.setattr(
            'amalgkit.select.read_select_config',
            lambda _path: {'rules': [], 'parameters': {}},
        )

        def fake_load_metadata(runtime_args):
            observed['loaded_metadata'] = runtime_args.metadata
            return metadata

        monkeypatch.setattr('amalgkit.select.load_metadata', fake_load_metadata)
        monkeypatch.setattr('amalgkit.select.prepare_select_metadata', lambda value, _rules: value)
        monkeypatch.setattr('amalgkit.select.apply_select_filters', lambda value, _args, _rules: value)
        monkeypatch.setattr(
            'amalgkit.select.write_select_outputs',
            lambda **kwargs: observed.update(kwargs),
        )

        select_main(args)

        assert observed['loaded_metadata'] == str(original_metadata)
        assert observed['preserve_existing_original'] is True


class TestSelectBatchMain:
    def test_batch_mode_writes_summary_queue_and_manifest(self, tmp_path):
        metadata_specieswise_dir = tmp_path / 'metadata_specieswise'
        metadata_specieswise_dir.mkdir()
        species_tsv = tmp_path / 'species.tsv'
        rules_path = tmp_path / 'select_rules.tsv'
        write_select_rules(rules_path, [
            {
                'rule_id': 'param_min_nspots',
                'stage': 'parameter',
                'parameter_name': 'min_nspots',
                'parameter_value': '0',
            },
            {
                'rule_id': 'param_max_sample',
                'stage': 'parameter',
                'parameter_name': 'max_sample',
                'parameter_value': '2',
            },
            {
                'rule_id': 'param_mark_missing_rank',
                'stage': 'parameter',
                'parameter_name': 'mark_missing_rank',
                'parameter_value': 'none',
            },
            {
                'rule_id': 'param_mark_redundant_biosamples',
                'stage': 'parameter',
                'parameter_name': 'mark_redundant_biosamples',
                'parameter_value': 'no',
            },
            {
                'rule_id': 'param_sample_group',
                'stage': 'parameter',
                'parameter_name': 'sample_group',
                'parameter_value': 'flower,leaf,root',
            },
            {
                'rule_id': 'flower_whole',
                'stage': 'normalize',
                'priority': '10',
                'columns': 'sample_attribute_tissue,sample_group',
                'pattern': r'\bflower\b',
                'action': 'assign',
                'outcome': 'flower',
            },
            {
                'rule_id': 'leaf_whole',
                'stage': 'normalize',
                'priority': '20',
                'columns': 'sample_attribute_tissue,sample_group',
                'pattern': r'\bleaf\b',
                'action': 'assign',
                'outcome': 'leaf',
            },
            {
                'rule_id': 'root_whole',
                'stage': 'normalize',
                'priority': '30',
                'columns': 'sample_attribute_tissue,sample_group',
                'pattern': r'\broot\b',
                'action': 'assign',
                'outcome': 'root',
            },
        ])

        species_rows = pandas.DataFrame({
            'scientific_name': ['Species alpha', 'Species beta'],
            'species_token': ['Species_alpha', 'Species_beta'],
        })
        species_rows.to_csv(species_tsv, sep='\t', index=False)

        def write_merged_metadata(species_token, scientific_name, organ_counts):
            species_dir = metadata_specieswise_dir / species_token
            species_dir.mkdir()
            rows = []
            counter = 1
            for organ, count in organ_counts.items():
                for _ in range(count):
                    rows.append({
                        'run': 'SRR{:04d}'.format(counter),
                        'scientific_name': scientific_name,
                        'sample_group': '',
                        'sample_attribute_tissue': organ,
                        'bioproject': 'PRJ{}'.format(counter),
                        'biosample': 'SAM{}'.format(counter),
                        'total_spots': 100,
                        'exclusion': 'no',
                    })
                    counter += 1
            pandas.DataFrame(rows).to_csv(
                species_dir / '{}.metadata.tsv'.format(species_token),
                sep='\t',
                index=False,
            )

        write_merged_metadata(
            species_token='Species_alpha',
            scientific_name='Species alpha',
            organ_counts={'flower': 2, 'leaf': 2, 'root': 2},
        )
        write_merged_metadata(
            species_token='Species_beta',
            scientific_name='Species beta',
            organ_counts={'flower': 1, 'leaf': 2, 'root': 2},
        )

        out_dir = tmp_path / 'select_batch'
        args = SimpleNamespace(
            out_dir=str(out_dir),
            select_rules_tsv=str(rules_path),
            metadata='inferred',
            sample_group=None,
            min_nspots=None,
            mark_missing_rank=None,
            mark_redundant_biosamples=None,
            max_sample=None,
            species_tsv=str(species_tsv),
            metadata_specieswise_dir=str(metadata_specieswise_dir),
            summary_tsv='inferred',
            queue_tsv='inferred',
            manifest_tsv='inferred',
            batch_label='inferred',
            threads='auto',
            internal_jobs='auto',
            internal_cpu_budget='auto',
        )

        select_main(args)

        normalization_summary = tmp_path / 'select_batch' / 'normalization_summary.tsv'
        select_summary = tmp_path / 'select_batch' / 'select_summary.tsv'
        select_queue = tmp_path / 'select_batch' / 'select_queue.tsv'
        manifest = tmp_path / 'select_batch' / 'external_manifest.tsv'
        manifest_all_tissues_ge30 = tmp_path / 'select_batch' / 'external_manifest_all_tissues_ge30.tsv'
        manifest_all_tissues_ge3 = tmp_path / 'select_batch' / 'external_manifest_all_tissues_ge3.tsv'
        manifest_all_tissues_ge1 = tmp_path / 'select_batch' / 'external_manifest_all_tissues_ge1.tsv'
        manifest_any_tissues_ge1 = tmp_path / 'select_batch' / 'external_manifest_any_tissues_ge1.tsv'

        assert normalization_summary.exists()
        assert select_summary.exists()
        assert select_queue.exists()
        assert manifest.exists()
        assert manifest_all_tissues_ge30.exists()
        assert manifest_all_tissues_ge3.exists()
        assert manifest_all_tissues_ge1.exists()
        assert manifest_any_tissues_ge1.exists()

        summary_df = pandas.read_csv(select_summary, sep='\t')
        assert set(summary_df['species_token'].tolist()) == {'Species_alpha', 'Species_beta'}
        queue_by_species = summary_df.set_index('species_token')['queue_tier'].to_dict()
        any_tissues_ge1_by_species = summary_df.set_index('species_token')['any_tissues_ge1'].to_dict()
        all_tissues_ge1_by_species = summary_df.set_index('species_token')['all_tissues_ge1'].to_dict()
        assert queue_by_species['Species_alpha'] == 'all_tissues_ge30'
        assert queue_by_species['Species_beta'] == 'all_tissues_ge1'
        assert bool(any_tissues_ge1_by_species['Species_alpha']) is True
        assert bool(any_tissues_ge1_by_species['Species_beta']) is True
        assert bool(all_tissues_ge1_by_species['Species_alpha']) is True
        assert bool(all_tissues_ge1_by_species['Species_beta']) is True

        manifest_df = pandas.read_csv(manifest, sep='\t')
        assert set(manifest_df['queue_tier'].tolist()) == {'all_tissues_ge30', 'all_tissues_ge1'}
        assert 'selected_metadata_path' in manifest_df.columns

    def test_batch_mode_uses_parallel_worker_path(self, tmp_path, monkeypatch):
        metadata_specieswise_dir = tmp_path / 'metadata_specieswise'
        metadata_specieswise_dir.mkdir()
        species_tsv = tmp_path / 'species.tsv'
        rules_path = tmp_path / 'select_rules.tsv'
        write_select_rules(rules_path, [
            {
                'rule_id': 'param_min_nspots',
                'stage': 'parameter',
                'parameter_name': 'min_nspots',
                'parameter_value': '0',
            },
            {
                'rule_id': 'param_max_sample',
                'stage': 'parameter',
                'parameter_name': 'max_sample',
                'parameter_value': '2',
            },
            {
                'rule_id': 'param_mark_missing_rank',
                'stage': 'parameter',
                'parameter_name': 'mark_missing_rank',
                'parameter_value': 'none',
            },
            {
                'rule_id': 'param_mark_redundant_biosamples',
                'stage': 'parameter',
                'parameter_name': 'mark_redundant_biosamples',
                'parameter_value': 'no',
            },
            {
                'rule_id': 'param_sample_group',
                'stage': 'parameter',
                'parameter_name': 'sample_group',
                'parameter_value': 'flower,leaf,root',
            },
            {
                'rule_id': 'flower_whole',
                'stage': 'normalize',
                'priority': '10',
                'columns': 'sample_attribute_tissue',
                'pattern': r'\bflower\b',
                'action': 'assign',
                'outcome': 'flower',
            },
            {
                'rule_id': 'leaf_whole',
                'stage': 'normalize',
                'priority': '20',
                'columns': 'sample_attribute_tissue',
                'pattern': r'\bleaf\b',
                'action': 'assign',
                'outcome': 'leaf',
            },
            {
                'rule_id': 'root_whole',
                'stage': 'normalize',
                'priority': '30',
                'columns': 'sample_attribute_tissue',
                'pattern': r'\broot\b',
                'action': 'assign',
                'outcome': 'root',
            },
        ])
        species_rows = pandas.DataFrame({
            'scientific_name': ['Species alpha', 'Species beta'],
            'species_token': ['Species_alpha', 'Species_beta'],
        })
        species_rows.to_csv(species_tsv, sep='\t', index=False)
        for species_token, scientific_name in [('Species_alpha', 'Species alpha'), ('Species_beta', 'Species beta')]:
            species_dir = metadata_specieswise_dir / species_token
            species_dir.mkdir()
            pandas.DataFrame([
                {
                    'run': 'SRR1',
                    'scientific_name': scientific_name,
                    'sample_group': '',
                    'sample_attribute_tissue': 'flower',
                    'bioproject': 'PRJ1',
                    'biosample': 'SAM1',
                    'total_spots': 100,
                    'exclusion': 'no',
                },
                {
                    'run': 'SRR2',
                    'scientific_name': scientific_name,
                    'sample_group': '',
                    'sample_attribute_tissue': 'leaf',
                    'bioproject': 'PRJ2',
                    'biosample': 'SAM2',
                    'total_spots': 100,
                    'exclusion': 'no',
                },
                {
                    'run': 'SRR3',
                    'scientific_name': scientific_name,
                    'sample_group': '',
                    'sample_attribute_tissue': 'root',
                    'bioproject': 'PRJ3',
                    'biosample': 'SAM3',
                    'total_spots': 100,
                    'exclusion': 'no',
                },
            ]).to_csv(species_dir / '{}.metadata.tsv'.format(species_token), sep='\t', index=False)

        observed = {}

        def fake_resolve_jobs(args, task_count):
            observed['task_count'] = task_count
            return 2

        def fake_run_tasks(task_items, task_fn, max_workers=1):
            observed['used_parallel'] = True
            observed['max_workers'] = max_workers
            results = {}
            for task in task_items:
                results[task] = task_fn(task)
            return results, []

        monkeypatch.setattr('amalgkit.select._resolve_select_species_jobs', fake_resolve_jobs)
        monkeypatch.setattr('amalgkit.select.run_tasks_with_optional_threads', fake_run_tasks)

        args = SimpleNamespace(
            out_dir=str(tmp_path / 'select_batch'),
            select_rules_tsv=str(rules_path),
            metadata='inferred',
            sample_group=None,
            min_nspots=None,
            mark_missing_rank=None,
            mark_redundant_biosamples=None,
            max_sample=None,
            species_tsv=str(species_tsv),
            metadata_specieswise_dir=str(metadata_specieswise_dir),
            summary_tsv='inferred',
            queue_tsv='inferred',
            manifest_tsv='inferred',
            batch_label='inferred',
            threads='auto',
            internal_jobs='auto',
            internal_cpu_budget='auto',
        )

        select_main(args)

        assert observed['task_count'] == 2
        assert observed['used_parallel'] is True
        assert observed['max_workers'] == 2


class TestCliEntry:
    def test_python_module_entrypoint_supports_help(self):
        repo_root = Path(__file__).resolve().parents[1]
        completed = subprocess.run(
            [sys.executable, '-m', 'amalgkit', 'select', '--help'],
            cwd=str(repo_root),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        assert completed.returncode == 0
        assert 'usage: amalgkit select' in completed.stdout
        assert '--sample_group tissueA,tissueB,tissueC,...' not in completed.stdout
        assert '--max_sample' not in completed.stdout
        assert '--mark_missing_rank' not in completed.stdout
        assert '--mark_redundant_biosamples' not in completed.stdout
