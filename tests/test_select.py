import os
import subprocess
import sys
import numpy
import pandas
import pytest

from types import SimpleNamespace
from pathlib import Path

from amalgkit.select import (
    apply_select_config_parameters,
    apply_select_control_rules,
    apply_select_filters,
    classify_select_text,
    prepare_select_metadata,
    read_select_config,
    read_select_rules,
    write_select_outputs,
    resolve_select_rules_tsv,
    filter_metadata_by_sample_group,
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

class TestWriteSelectOutputs:
    def test_writes_metadata_and_pivots(self, tmp_path, sample_metadata):
        """Writes metadata.tsv, pivot_qualified.tsv, and pivot_selected.tsv."""
        metadata_dir = tmp_path / 'metadata'
        metadata_dir.mkdir()
        path_original = str(tmp_path / 'metadata' / 'metadata_original.tsv')
        path_table = str(tmp_path / 'metadata' / 'metadata.tsv')
        # Write a source metadata table to copy
        sample_metadata.df.to_csv(path_table, sep='\t', index=False)
        sample_metadata.label_sampled_data(max_sample=10)
        write_select_outputs(path_original, path_table, str(metadata_dir), sample_metadata)
        assert os.path.exists(path_original)
        assert os.path.exists(path_table)
        assert os.path.exists(str(metadata_dir / 'pivot_qualified.tsv'))
        assert os.path.exists(str(metadata_dir / 'pivot_selected.tsv'))

    def test_refreshes_original_from_current_input(self, tmp_path, sample_metadata):
        metadata_dir = tmp_path / 'metadata'
        metadata_dir.mkdir()
        path_original = str(tmp_path / 'metadata' / 'metadata_original.tsv')
        path_table = str(tmp_path / 'metadata' / 'metadata.tsv')
        original_df = sample_metadata.df.copy(deep=True)
        with open(path_original, 'w') as f:
            f.write('marker_content\n')
        sample_metadata.df.to_csv(path_table, sep='\t', index=False)
        sample_metadata.label_sampled_data(max_sample=10)
        write_select_outputs(
            path_original,
            path_table,
            str(metadata_dir),
            sample_metadata,
            metadata_original_df=original_df,
        )
        loaded = pandas.read_csv(path_original, sep='\t')
        assert loaded.shape[0] == original_df.shape[0]
        assert set(loaded['run'].tolist()) == set(original_df['run'].tolist())

    def test_writes_original_without_existing_metadata_table(self, tmp_path, sample_metadata):
        metadata_dir = tmp_path / 'metadata'
        path_original = str(tmp_path / 'metadata' / 'metadata_original.tsv')
        path_table = str(tmp_path / 'metadata' / 'metadata.tsv')
        sample_metadata.label_sampled_data(max_sample=10)
        original_df = sample_metadata.df.copy(deep=True)

        write_select_outputs(
            path_metadata_original=path_original,
            path_metadata_table=path_table,
            metadata_dir=str(metadata_dir),
            metadata=sample_metadata,
            metadata_original_df=original_df,
        )

        assert os.path.exists(path_original)
        assert os.path.exists(path_table)
        original_loaded = pandas.read_csv(path_original, sep='\t')
        assert original_loaded.shape[0] == original_df.shape[0]
        assert set(original_loaded['run'].tolist()) == set(original_df['run'].tolist())

    def test_uses_current_input_for_original_even_when_metadata_table_exists(self, tmp_path, sample_metadata):
        metadata_dir = tmp_path / 'metadata'
        metadata_dir.mkdir()
        path_original = str(metadata_dir / 'metadata_original.tsv')
        path_table = str(metadata_dir / 'metadata.tsv')
        existing_filtered = sample_metadata.df.iloc[[0], :].copy()
        existing_filtered.to_csv(path_table, sep='\t', index=False)
        original_df = sample_metadata.df.copy(deep=True)
        sample_metadata.label_sampled_data(max_sample=10)

        write_select_outputs(
            path_metadata_original=path_original,
            path_metadata_table=path_table,
            metadata_dir=str(metadata_dir),
            metadata=sample_metadata,
            metadata_original_df=original_df,
        )

        original_loaded = pandas.read_csv(path_original, sep='\t')
        assert original_loaded.shape[0] == original_df.shape[0]
        assert set(original_loaded['run'].tolist()) == set(original_df['run'].tolist())
        assert original_loaded.shape[0] != existing_filtered.shape[0]

    def test_rejects_metadata_original_directory_path(self, tmp_path, sample_metadata):
        metadata_dir = tmp_path / 'metadata'
        metadata_dir.mkdir()
        path_original = metadata_dir / 'metadata_original.tsv'
        path_original.mkdir()
        path_table = str(metadata_dir / 'metadata.tsv')
        sample_metadata.label_sampled_data(max_sample=10)

        with pytest.raises(NotADirectoryError, match='not a file'):
            write_select_outputs(str(path_original), path_table, str(metadata_dir), sample_metadata)

    def test_rejects_metadata_table_directory_path(self, tmp_path, sample_metadata):
        metadata_dir = tmp_path / 'metadata'
        metadata_dir.mkdir()
        path_original = str(metadata_dir / 'metadata_original.tsv')
        path_table = metadata_dir / 'metadata.tsv'
        path_table.mkdir()
        sample_metadata.label_sampled_data(max_sample=10)

        with pytest.raises(NotADirectoryError, match='not a file'):
            write_select_outputs(path_original, str(path_table), str(metadata_dir), sample_metadata)

    def test_rejects_metadata_dir_file_path(self, tmp_path, sample_metadata):
        metadata_dir = tmp_path / 'metadata'
        metadata_dir.write_text('not a directory')
        path_original = str(tmp_path / 'metadata_original.tsv')
        path_table = str(tmp_path / 'metadata.tsv')
        sample_metadata.label_sampled_data(max_sample=10)

        with pytest.raises(NotADirectoryError, match='not a directory'):
            write_select_outputs(path_original, path_table, str(metadata_dir), sample_metadata)


class TestSelectHelpers:
    def test_resolve_select_rules_tsv_inferred(self):
        args = SimpleNamespace(out_dir='/tmp/out', select_rules_tsv='inferred')
        assert resolve_select_rules_tsv(args) == os.path.realpath('/tmp/out/select_rules.tsv')

    def test_resolve_select_rules_tsv_explicit(self):
        args = SimpleNamespace(out_dir='/tmp/out', select_rules_tsv='/tmp/custom.tsv')
        assert resolve_select_rules_tsv(args) == os.path.realpath('/tmp/custom.tsv')

    def test_filter_metadata_by_sample_group(self):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001', 'SRR002', 'SRR003'],
            'sample_group': ['brain', 'liver', 'brain'],
            'scientific_name': ['sp1', 'sp1', 'sp2'],
            'exclusion': ['no', 'no', 'no'],
        }))
        out = filter_metadata_by_sample_group(metadata, 'brain')
        assert set(out.df['run']) == {'SRR001', 'SRR003'}

    def test_filter_metadata_by_sample_group_strips_tokens(self):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001', 'SRR002', 'SRR003'],
            'sample_group': ['brain', 'liver', 'heart'],
            'scientific_name': ['sp1', 'sp1', 'sp2'],
            'exclusion': ['no', 'no', 'no'],
        }))
        out = filter_metadata_by_sample_group(metadata, 'brain, liver')
        assert set(out.df['run']) == {'SRR001', 'SRR002'}

    def test_filter_metadata_by_sample_group_supports_pipe_separator(self):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001', 'SRR002', 'SRR003'],
            'sample_group': ['brain', 'liver', 'heart'],
            'scientific_name': ['sp1', 'sp1', 'sp2'],
            'exclusion': ['no', 'no', 'no'],
        }))
        out = filter_metadata_by_sample_group(metadata, 'brain|heart')
        assert set(out.df['run']) == {'SRR001', 'SRR003'}

    def test_filter_metadata_by_sample_group_strips_metadata_values(self):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001', 'SRR002'],
            'sample_group': [' brain ', 'liver'],
            'scientific_name': ['sp1', 'sp1'],
            'exclusion': ['no', 'no'],
        }))
        out = filter_metadata_by_sample_group(metadata, 'brain')
        assert set(out.df['run']) == {'SRR001'}

    def test_filter_metadata_by_sample_group_none_keeps_all(self):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001', 'SRR002'],
            'sample_group': ['brain', 'liver'],
            'scientific_name': ['sp1', 'sp1'],
            'exclusion': ['no', 'no'],
        }))
        out = filter_metadata_by_sample_group(metadata, None)
        assert out.df.shape[0] == 2

    def test_filter_metadata_by_sample_group_empty_argument_raises(self):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001', 'SRR002'],
            'sample_group': ['brain', 'liver'],
            'scientific_name': ['sp1', 'sp1'],
            'exclusion': ['no', 'no'],
        }))
        with pytest.raises(ValueError, match='No sample_group was selected'):
            filter_metadata_by_sample_group(metadata, '  ')

    def test_filter_metadata_by_sample_group_requires_column(self):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001', 'SRR002'],
            'scientific_name': ['sp1', 'sp1'],
            'exclusion': ['no', 'no'],
        }))
        metadata.df = metadata.df.drop(columns=['sample_group'])
        with pytest.raises(ValueError, match='Column \"sample_group\" is required'):
            filter_metadata_by_sample_group(metadata, 'brain')

    def test_read_select_config_parses_parameter_rows(self, tmp_path):
        rules_path = tmp_path / 'select_rules.tsv'
        write_select_rules(rules_path, [
            {
                'rule_id': 'param_min_nspots',
                'stage': 'parameter',
                'parameter_name': 'min_nspots',
                'parameter_value': '1000000',
            },
            {
                'rule_id': 'param_sampling_strategy',
                'stage': 'parameter',
                'parameter_name': 'sampling_strategy',
                'parameter_value': 'random',
            },
            {
                'rule_id': 'flower_whole',
                'stage': 'normalize',
                'priority': '10',
                'columns': 'sample_group',
                'pattern': r'\\bflower\\b',
                'action': 'assign',
                'outcome': 'flower',
            },
        ])

        config = read_select_config(str(rules_path))

        assert config['parameters']['min_nspots'] == 1000000
        assert config['parameters']['sampling_strategy'] == 'random'
        assert len(config['rules']) == 1
        assert config['rules'][0]['rule_id'] == 'flower_whole'

    def test_read_select_config_parses_filter_dedup_and_validate_rows(self, tmp_path):
        rules_path = tmp_path / 'select_rules.tsv'
        write_select_rules(rules_path, [
            {
                'rule_id': 'validate_hint_leaf',
                'stage': 'validate',
                'priority': '10',
                'pattern': r'\\bleaf\\b',
                'action': 'hint_organ',
                'outcome': 'leaf',
            },
            {
                'rule_id': 'filter_no_sample_group',
                'stage': 'filter',
                'priority': '20',
                'columns': 'sample_group',
                'action': 'exclude_if_empty',
                'target_column': 'exclusion',
                'outcome': 'no_sample_group',
            },
            {
                'rule_id': 'dedup_redundant_biosample',
                'stage': 'dedup',
                'priority': '30',
                'columns': 'bioproject,biosample',
                'action': 'exclude_redundant_by_best_numeric',
                'target_column': 'total_spots',
                'outcome': 'redundant_biosample',
                'parameter_name': 'mark_redundant_biosamples',
            },
            {
                'rule_id': 'flower_whole',
                'stage': 'normalize',
                'priority': '40',
                'columns': 'sample_group',
                'pattern': r'\\bflower\\b',
                'action': 'assign',
                'outcome': 'flower',
            },
        ])

        config = read_select_config(str(rules_path))

        rule_by_id = {rule['rule_id']: rule for rule in config['rules']}
        assert rule_by_id['validate_hint_leaf']['stage'] == 'validate'
        assert rule_by_id['filter_no_sample_group']['stage'] == 'filter'
        assert rule_by_id['dedup_redundant_biosample']['stage'] == 'dedup'
        assert rule_by_id['dedup_redundant_biosample']['parameter_name'] == 'mark_redundant_biosamples'

    def test_apply_select_config_parameters_uses_config_when_args_missing(self):
        args = SimpleNamespace(
            min_nspots=None,
            max_sample=None,
            mark_missing_rank=None,
            mark_redundant_biosamples=None,
            sample_group=None,
            sampling_strategy=None,
        )

        out = apply_select_config_parameters(
            args,
            {
                'min_nspots': 1000000,
                'max_sample': 30,
                'mark_missing_rank': 'none',
                'mark_redundant_biosamples': False,
                'sample_group': 'flower,leaf,root',
                'sampling_strategy': 'largest_bioprojects_first',
            },
        )

        assert out.min_nspots == 1000000
        assert out.max_sample == 30
        assert out.mark_missing_rank == 'none'
        assert out.mark_redundant_biosamples is False
        assert out.sample_group == 'flower,leaf,root'
        assert out.sampling_strategy == 'largest_bioprojects_first'

    def test_apply_select_config_parameters_keeps_runtime_sample_group_override(self):
        args = SimpleNamespace(
            min_nspots=None,
            max_sample=None,
            mark_missing_rank=None,
            mark_redundant_biosamples=None,
            sample_group='leaf',
            sampling_strategy='random',
        )

        out = apply_select_config_parameters(
            args,
            {
                'min_nspots': 1000000,
                'max_sample': 30,
                'mark_missing_rank': 'none',
                'mark_redundant_biosamples': False,
                'sample_group': 'flower,leaf,root',
                'sampling_strategy': 'largest_bioprojects_first',
            },
        )

        assert out.sample_group == 'leaf'
        assert out.sampling_strategy == 'random'


class TestSamplingStrategies:
    @staticmethod
    def _build_sampling_metadata():
        return Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001', 'SRR002', 'SRR003', 'SRR004', 'SRR005', 'SRR006', 'SRR007'],
            'scientific_name': ['Species A'] * 7,
            'sample_group': ['leaf'] * 7,
            'bioproject': ['PRJ1', 'PRJ1', 'PRJ1', 'PRJ1', 'PRJ2', 'PRJ2', 'PRJ3'],
            'biosample': ['SAM1', 'SAM2', 'SAM3', 'SAM4', 'SAM5', 'SAM6', 'SAM7'],
            'total_spots': [100] * 7,
            'exclusion': ['no'] * 7,
        }))

    def test_maximize_bioproject_diversity_is_current_default(self, monkeypatch):
        metadata = self._build_sampling_metadata()
        monkeypatch.setattr(numpy.random, 'permutation', lambda values: numpy.array(list(values)))

        metadata.label_sampled_data(max_sample=3, sampling_strategy='maximize_bioproject_diversity')

        selected = metadata.df.loc[metadata.df['is_sampled'] == 'yes', 'bioproject'].value_counts().to_dict()
        assert selected == {'PRJ1': 1, 'PRJ2': 1, 'PRJ3': 1}

    def test_largest_bioprojects_first_fills_largest_project_before_others(self, monkeypatch):
        metadata = self._build_sampling_metadata()
        monkeypatch.setattr(numpy.random, 'permutation', lambda values: numpy.array(list(values)))

        metadata.label_sampled_data(max_sample=3, sampling_strategy='largest_bioprojects_first')

        selected = metadata.df.loc[metadata.df['is_sampled'] == 'yes', 'bioproject'].value_counts().to_dict()
        assert selected == {'PRJ1': 3}

    def test_smallest_bioprojects_first_fills_smallest_projects_before_others(self, monkeypatch):
        metadata = self._build_sampling_metadata()
        monkeypatch.setattr(numpy.random, 'permutation', lambda values: numpy.array(list(values)))

        metadata.label_sampled_data(max_sample=3, sampling_strategy='smallest_bioprojects_first')

        selected = metadata.df.loc[metadata.df['is_sampled'] == 'yes', 'bioproject'].value_counts().to_dict()
        assert selected == {'PRJ2': 2, 'PRJ3': 1}

    def test_random_sampling_ignores_bioproject_balance(self, monkeypatch):
        metadata = self._build_sampling_metadata()
        monkeypatch.setattr(numpy.random, 'permutation', lambda values: numpy.array(list(values)))

        metadata.label_sampled_data(max_sample=3, sampling_strategy='random')

        selected_runs = metadata.df.loc[metadata.df['is_sampled'] == 'yes', 'run'].tolist()
        assert selected_runs == ['SRR001', 'SRR002', 'SRR003']

    def test_invalid_sampling_strategy_raises(self):
        metadata = self._build_sampling_metadata()

        with pytest.raises(ValueError, match='Unknown sampling_strategy'):
            metadata.label_sampled_data(max_sample=3, sampling_strategy='unsupported_strategy')

    @pytest.mark.parametrize(
        ('text', 'expected_status', 'expected_organ'),
        [
            ('flower', 'organ', 'flower'),
            ('Flower of Bora mutant(purple-petal) at full-bloom stage rep1', 'organ', 'flower'),
            ('root tips', 'organ', 'root'),
            ('primary root', 'organ', 'root'),
            ('lateral root', 'organ', 'root'),
            ('radicle', 'organ', 'root'),
            ('leaf blade', 'organ', 'leaf'),
            ('seedling leaf', 'organ', 'leaf'),
            ('spikes in seedling stage', 'organ', 'flower'),
            ('seedling roots', 'organ', 'root'),
            ('seedling root', 'organ', 'root'),
            ('Leaves were sampled from fully grown plants with well developed roots and flowers.', 'mixed', ''),
            ('leaf and spine', 'review', ''),
            ('leaf and leaves', 'organ', 'leaf'),
            ('leaf,spike', 'mixed', ''),
            ('flag leaf; young spike', 'mixed', ''),
            ('RNA-Seq of Cannabis Plant 3; Auto-flower; Vegetative; Day 20', 'review', ''),
            ('RNAseq of Panax notoginseng with flower removal1', 'review', ''),
            ('RNAseq of Panax notoginseng with flower kept1', 'review', ''),
            ('flower, biological replicate D1', 'organ', 'flower'),
            ('flower, cold treatment', 'organ', 'flower'),
            ('leaf, cultivar A', 'organ', 'leaf'),
            ('leaf, petal', 'review', ''),
            ('unopened flower buds', 'review', ''),
            ('leaflet', 'organ', 'leaf'),
            ('lamina', 'organ', 'leaf'),
            ('taproot', 'organ', 'root'),
            ('bract', 'review', ''),
            ('inflorescence apex', 'review', ''),
            ('petal', 'review', ''),
            ('corolla', 'review', ''),
            ('anther', 'review', ''),
            ('ovary', 'review', ''),
            ('pistil', 'review', ''),
            ('flower bud', 'review', ''),
            ('petiole', 'review', ''),
            ('root hair', 'review', ''),
            ('hairy root', 'organ', 'root'),
        ],
    )
    def test_classify_select_text_respects_whole_organ_policy(self, tmp_path, text, expected_status, expected_organ):
        rules_path = tmp_path / 'select_rules.tsv'
        write_select_rules(rules_path, [
            {
                'rule_id': 'review_hairy_root_culture',
                'stage': 'normalize',
                'priority': '10',
                'columns': 'sample_group',
                'pattern': r'\bhairy root culture\b',
                'action': 'assign',
                'outcome': 'review',
            },
            {
                'rule_id': 'root_taproot',
                'stage': 'normalize',
                'priority': '15',
                'columns': 'sample_group',
                'pattern': r'\b(?:taproot|taproots)\b',
                'action': 'assign',
                'outcome': 'root',
            },
            {
                'rule_id': 'root_whole_root_tips',
                'stage': 'normalize',
                'priority': '16',
                'columns': 'sample_group',
                'pattern': r'\b(?:root[\s_-]?tip|root[\s_-]?tips)\b',
                'action': 'assign',
                'outcome': 'root',
            },
            {
                'rule_id': 'root_whole_primary_lateral_radicle',
                'stage': 'normalize',
                'priority': '17',
                'columns': 'sample_group',
                'pattern': r'\b(?:primary root|primary roots|lateral root|lateral roots|radicle|radicles)\b',
                'action': 'assign',
                'outcome': 'root',
            },
            {
                'rule_id': 'mixed_leaf_root',
                'stage': 'normalize',
                'priority': '18',
                'columns': 'sample_group',
                'pattern': r'\b(?:leaf|leaves|foliar|foliage|trifoliate|leave|lamina|laminae|leaflet|leaflets|leaf[\s_-]?blades?|seedling[\s_-]?leaves?)\b.*\b(?:root|roots|taproot|taproots|primary root|primary roots|lateral root|lateral roots|radicle|radicles|root[\s_-]?tip|root[\s_-]?tips|seedling[\s_-]?roots?|hairy root|hairy roots)\b|\b(?:root|roots|taproot|taproots|primary root|primary roots|lateral root|lateral roots|radicle|radicles|root[\s_-]?tip|root[\s_-]?tips|seedling[\s_-]?roots?|hairy root|hairy roots)\b.*\b(?:leaf|leaves|foliar|foliage|trifoliate|leave|lamina|laminae|leaflet|leaflets|leaf[\s_-]?blades?|seedling[\s_-]?leaves?)\b',
                'action': 'assign',
                'outcome': 'mixed',
            },
            {
                'rule_id': 'mixed_leaf_flower',
                'stage': 'normalize',
                'priority': '18',
                'columns': 'sample_group',
                'pattern': r'\b(?:leaf|leaves|foliar|foliage|trifoliate|leave|lamina|laminae|leaflet|leaflets|leaf[\s_-]?blades?|seedling[\s_-]?leaves?)\b.*\b(?:flower|flowers|floral|inflorescence|inflorescences|catkin|catkins|spikes?\s+in\s+seedling\s+stage)\b|\b(?:flower|flowers|floral|inflorescence|inflorescences|catkin|catkins|spikes?\s+in\s+seedling\s+stage)\b.*\b(?:leaf|leaves|foliar|foliage|trifoliate|leave|lamina|laminae|leaflet|leaflets|leaf[\s_-]?blades?|seedling[\s_-]?leaves?)\b',
                'action': 'assign',
                'outcome': 'mixed',
            },
            {
                'rule_id': 'mixed_root_flower',
                'stage': 'normalize',
                'priority': '18',
                'columns': 'sample_group',
                'pattern': r'\b(?:root|roots|taproot|taproots|primary root|primary roots|lateral root|lateral roots|radicle|radicles|root[\s_-]?tip|root[\s_-]?tips|seedling[\s_-]?roots?|hairy root|hairy roots)\b.*\b(?:flower|flowers|floral|inflorescence|inflorescences|catkin|catkins|spikes?\s+in\s+seedling\s+stage)\b|\b(?:flower|flowers|floral|inflorescence|inflorescences|catkin|catkins|spikes?\s+in\s+seedling\s+stage)\b.*\b(?:root|roots|taproot|taproots|primary root|primary roots|lateral root|lateral roots|radicle|radicles|root[\s_-]?tip|root[\s_-]?tips|seedling[\s_-]?roots?|hairy root|hairy roots)\b',
                'action': 'assign',
                'outcome': 'mixed',
            },
            {
                'rule_id': 'review_flower_false_positive_phrase',
                'stage': 'normalize',
                'priority': '27',
                'columns': 'sample_group',
                'pattern': r'\bauto[\s_-]?flower(?:ing)?\b|\bflower[\s_-]?(?:removal|kept)\d*\b|\b(?:removal|kept)\d*[\s_-]?flower\b',
                'action': 'assign',
                'outcome': 'review',
            },
            {
                'rule_id': 'flower_whole_full_bloom',
                'stage': 'normalize',
                'priority': '28',
                'columns': 'sample_group',
                'pattern': r'\b(?:flower|flowers)\b(?![\s_-]?buds?)\b.{0,40}\bfull[- ]bloom\b|\bfull[- ]bloom\b.{0,40}\b(?:flower|flowers)\b(?![\s_-]?buds?)',
                'action': 'assign',
                'outcome': 'flower',
            },
            {
                'rule_id': 'review_suborgan',
                'stage': 'normalize',
                'priority': '30',
                'columns': 'sample_group',
                'pattern': r'\b(?:petal|corolla|anther|ovary|pistil|bract|flower[\s_-]?buds?|petiole|root hair)\b',
                'action': 'assign',
                'outcome': 'review',
            },
            {
                'rule_id': 'leaf_whole_leaf_blade',
                'stage': 'normalize',
                'priority': '25',
                'columns': 'sample_group',
                'pattern': r'\bleaf[\s_-]?blades?\b',
                'action': 'assign',
                'outcome': 'leaf',
            },
            {
                'rule_id': 'flower_whole_spikes_seedling_stage',
                'stage': 'normalize',
                'priority': '26',
                'columns': 'sample_group',
                'pattern': r'\bspikes?\s+in\s+seedling\s+stage\b',
                'action': 'assign',
                'outcome': 'flower',
            },
            {
                'rule_id': 'leaf_whole_seedling_leaf',
                'stage': 'normalize',
                'priority': '27',
                'columns': 'sample_group',
                'pattern': r'\bseedling[\s_-]?leaves?\b',
                'action': 'assign',
                'outcome': 'leaf',
            },
            {
                'rule_id': 'root_whole_seedling_root',
                'stage': 'normalize',
                'priority': '28',
                'columns': 'sample_group',
                'pattern': r'\bseedling[\s_-]?roots?\b',
                'action': 'assign',
                'outcome': 'root',
            },
            {
                'rule_id': 'review_structure',
                'stage': 'normalize',
                'priority': '30',
                'columns': 'sample_group',
                'pattern': r'\b(?:inflorescence apex)\b',
                'action': 'assign',
                'outcome': 'review',
            },
            {
                'rule_id': 'root_hairy_root',
                'stage': 'normalize',
                'priority': '100',
                'columns': 'sample_group',
                'pattern': r'\bhairy root\b',
                'action': 'assign',
                'outcome': 'root',
            },
            {
                'rule_id': 'flower_whole',
                'stage': 'normalize',
                'priority': '110',
                'columns': 'sample_group',
                'pattern': r'\b(?:flower|flowers|spike|spikes)\b',
                'action': 'assign',
                'outcome': 'flower',
            },
            {
                'rule_id': 'leaf_whole',
                'stage': 'normalize',
                'priority': '120',
                'columns': 'sample_group',
                'pattern': r'\b(?:leaf|leaves|lamina|laminae|leaflet|leaflets)\b',
                'action': 'assign',
                'outcome': 'leaf',
            },
        ] + build_test_validate_rows())
        parsed_rules = read_select_rules(str(rules_path))
        normalize_rules = [rule for rule in parsed_rules if rule['stage'] == 'normalize']
        validate_rules = [rule for rule in parsed_rules if rule['stage'] == 'validate']
        result = classify_select_text(text, normalize_rules, validate_rules=validate_rules)
        assert result['status'] == expected_status
        assert result['organ'] == expected_organ

    def test_prepare_select_metadata_applies_unknown_segment_validator(self, tmp_path):
        rules_path = tmp_path / 'select_rules.tsv'
        write_select_rules(rules_path, build_test_validate_rows() + [
            {
                'rule_id': 'leaf_whole',
                'stage': 'normalize',
                'priority': '120',
                'columns': 'sample_group',
                'pattern': r'\b(?:leaf|leaves|lamina|laminae|leaflet|leaflets)\b',
                'action': 'assign',
                'outcome': 'leaf',
            },
            {
                'rule_id': 'root_whole',
                'stage': 'normalize',
                'priority': '130',
                'columns': 'sample_group',
                'pattern': r'\b(?:root|roots)\b',
                'action': 'assign',
                'outcome': 'root',
            },
            {
                'rule_id': 'flower_whole',
                'stage': 'normalize',
                'priority': '140',
                'columns': 'sample_group',
                'pattern': r'\b(?:flower|flowers)\b',
                'action': 'assign',
                'outcome': 'flower',
            },
        ])
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001', 'SRR002'],
            'scientific_name': ['sp1', 'sp1'],
            'sample_group': ['leaf and spine', 'leaf and leaves'],
            'exclusion': ['no', 'no'],
        }))
        out = prepare_select_metadata(metadata, read_select_rules(str(rules_path)))
        observed = out.df.set_index('run')[[
            'sample_group',
            'sample_group_normalization_status',
            'sample_group_normalization_rule_id',
        ]]
        assert observed.loc['SRR001', 'sample_group'] == 'review'
        assert observed.loc['SRR001', 'sample_group_normalization_status'] == 'review'
        assert observed.loc['SRR001', 'sample_group_normalization_rule_id'] == 'validate_review_segment'
        assert observed.loc['SRR002', 'sample_group'] == 'leaf'
        assert observed.loc['SRR002', 'sample_group_normalization_status'] == 'organ'

    def test_prepare_select_metadata_applies_list_segment_validator(self, tmp_path):
        rules_path = tmp_path / 'select_rules.tsv'
        write_select_rules(rules_path, build_test_validate_rows() + [
            {
                'rule_id': 'flower_whole',
                'stage': 'normalize',
                'priority': '10',
                'columns': 'sample_group',
                'pattern': r'\b(?:flower|flowers|spike|spikes|inflorescence|inflorescences)\b',
                'action': 'assign',
                'outcome': 'flower',
            },
            {
                'rule_id': 'leaf_whole',
                'stage': 'normalize',
                'priority': '20',
                'columns': 'sample_group',
                'pattern': r'\b(?:leaf|leaves)\b',
                'action': 'assign',
                'outcome': 'leaf',
            },
            {
                'rule_id': 'review_suborgan',
                'stage': 'normalize',
                'priority': '30',
                'columns': 'sample_group',
                'pattern': r'\b(?:petal|petals)\b',
                'action': 'assign',
                'outcome': 'review',
            },
        ])
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001', 'SRR002', 'SRR003', 'SRR004', 'SRR005'],
            'scientific_name': ['sp1', 'sp1', 'sp1', 'sp1', 'sp1'],
            'sample_group': [
                'leaf,spike',
                'flower, biological replicate D1',
                'leaf; petal',
                'flower, cold treatment',
                'leaf, cultivar A',
            ],
            'exclusion': ['no', 'no', 'no', 'no', 'no'],
        }))
        out = prepare_select_metadata(metadata, read_select_rules(str(rules_path)))
        observed = out.df.set_index('run')[[
            'sample_group',
            'sample_group_normalization_status',
            'sample_group_normalization_rule_id',
        ]]
        assert observed.loc['SRR001', 'sample_group'] == 'mixed'
        assert observed.loc['SRR001', 'sample_group_normalization_rule_id'] == 'validate_mixed_segment'
        assert observed.loc['SRR002', 'sample_group'] == 'flower'
        assert observed.loc['SRR002', 'sample_group_normalization_status'] == 'organ'
        assert observed.loc['SRR003', 'sample_group'] == 'review'
        assert observed.loc['SRR003', 'sample_group_normalization_rule_id'] == 'validate_review_segment'
        assert observed.loc['SRR004', 'sample_group'] == 'flower'
        assert observed.loc['SRR004', 'sample_group_normalization_status'] == 'organ'
        assert observed.loc['SRR005', 'sample_group'] == 'leaf'
        assert observed.loc['SRR005', 'sample_group_normalization_status'] == 'organ'

    def test_assign_safe_normalize_rule_skips_validator(self, tmp_path):
        rules_path = tmp_path / 'select_rules.tsv'
        write_select_rules(rules_path, build_test_validate_rows() + [
            {
                'rule_id': 'leaf_safe_floral_induction_context',
                'stage': 'normalize',
                'priority': '10',
                'columns': 'sample_group',
                'pattern': r'\b(?:leaf|leaves)\b.{0,120}\b(?:flower|floral)\s+induction\b|\b(?:flower|floral)\s+induction\b.{0,120}\b(?:leaf|leaves)\b',
                'action': 'assign_safe',
                'outcome': 'leaf',
            },
            {
                'rule_id': 'mixed_leaf_flower',
                'stage': 'normalize',
                'priority': '20',
                'columns': 'sample_group',
                'pattern': r'\b(?:leaf|leaves)\b.*\b(?:flower|flowers|floral)\b|\b(?:flower|flowers|floral)\b.*\b(?:leaf|leaves)\b',
                'action': 'assign',
                'outcome': 'mixed',
            },
            {
                'rule_id': 'flower_whole',
                'stage': 'normalize',
                'priority': '30',
                'columns': 'sample_group',
                'pattern': r'\b(?:flower|flowers)\b',
                'action': 'assign',
                'outcome': 'flower',
            },
            {
                'rule_id': 'leaf_whole',
                'stage': 'normalize',
                'priority': '40',
                'columns': 'sample_group',
                'pattern': r'\b(?:leaf|leaves)\b',
                'action': 'assign',
                'outcome': 'leaf',
            },
        ])
        parsed_rules = read_select_rules(str(rules_path))
        normalize_rules = [rule for rule in parsed_rules if rule['stage'] == 'normalize']
        validate_rules = [rule for rule in parsed_rules if rule['stage'] == 'validate']

        result = classify_select_text(
            'RNA-seq of Cajanus cajan during floral induction: Leaves before induction',
            normalize_rules,
            validate_rules=validate_rules,
        )

        assert result['status'] == 'organ'
        assert result['organ'] == 'leaf'
        assert result['rule_id'] == 'leaf_safe_floral_induction_context'


class TestSelectRuleApplication:
    def test_control_rules_protect_union_of_controls_within_scope(self, tmp_path):
        rules_path = tmp_path / 'select_rules.tsv'
        write_select_rules(rules_path, [
            {
                'rule_id': 'control_wt_title',
                'stage': 'control',
                'priority': '10',
                'columns': 'sample_title',
                'pattern': r'\bWT\b',
                'action': 'mark_non_control',
                'target_column': 'exclusion',
                'outcome': 'non_control',
                'scope_column': 'bioproject',
                'scope_mode': 'mark_other_rows_in_scope',
            },
            {
                'rule_id': 'control_wild_description',
                'stage': 'control',
                'priority': '20',
                'columns': 'sample_description',
                'pattern': r'wild',
                'action': 'mark_non_control',
                'target_column': 'exclusion',
                'outcome': 'non_control',
                'scope_column': 'bioproject',
                'scope_mode': 'mark_other_rows_in_scope',
            },
        ])
        select_rules = read_select_rules(str(rules_path))
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001', 'SRR002', 'SRR003'],
            'scientific_name': ['Species A', 'Species A', 'Species A'],
            'sample_group': ['flower', 'flower', 'flower'],
            'bioproject': ['PRJ1', 'PRJ1', 'PRJ1'],
            'biosample': ['SAM1', 'SAM2', 'SAM3'],
            'total_spots': [100, 100, 100],
            'exclusion': ['no', 'no', 'no'],
            'sample_title': ['WT-Flower-1', 'Flower-2', '397a-Flower-1'],
            'sample_description': ['', 'Flowers of wild plant', 'Transgenic flower'],
        }))

        out = apply_select_control_rules(metadata, select_rules)

        assert out.df.loc[out.df['run'] == 'SRR001', 'exclusion'].iloc[0] == 'no'
        assert out.df.loc[out.df['run'] == 'SRR002', 'exclusion'].iloc[0] == 'no'
        assert out.df.loc[out.df['run'] == 'SRR003', 'exclusion'].iloc[0] == 'non_control'

    def test_control_rules_respect_multi_column_scope(self, tmp_path):
        rules_path = tmp_path / 'select_rules.tsv'
        write_select_rules(rules_path, [
            {
                'rule_id': 'control_mock_treatment',
                'stage': 'control',
                'priority': '10',
                'columns': 'treatment',
                'pattern': r'mock',
                'action': 'mark_non_control',
                'target_column': 'exclusion',
                'outcome': 'non_control',
                'scope_column': 'bioproject,sample_group',
                'scope_mode': 'mark_other_rows_in_scope',
            },
        ])
        select_rules = read_select_rules(str(rules_path))
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001', 'SRR002', 'SRR003', 'SRR004'],
            'scientific_name': ['Species A'] * 4,
            'sample_group': ['leaf', 'leaf', 'flower', 'flower'],
            'bioproject': ['PRJ1'] * 4,
            'biosample': ['SAM1', 'SAM2', 'SAM3', 'SAM4'],
            'total_spots': [100, 100, 100, 100],
            'exclusion': ['no', 'no', 'no', 'no'],
            'treatment': ['mock treatment', 'PepMov infection', '', ''],
        }))

        out = apply_select_control_rules(metadata, select_rules)

        assert out.df.loc[out.df['run'] == 'SRR001', 'exclusion'].iloc[0] == 'no'
        assert out.df.loc[out.df['run'] == 'SRR002', 'exclusion'].iloc[0] == 'non_control'
        assert out.df.loc[out.df['run'] == 'SRR003', 'exclusion'].iloc[0] == 'no'
        assert out.df.loc[out.df['run'] == 'SRR004', 'exclusion'].iloc[0] == 'no'

    def test_prepare_and_filter_metadata_applies_single_file_rules(self, tmp_path):
        rules_path = tmp_path / 'select_rules.tsv'
        write_select_rules(rules_path, [
            {
                'rule_id': 'aggregate_sample_attribute',
                'stage': 'aggregate',
                'priority': '10',
                'columns': 'sample_attribute',
                'action': 'append',
                'target_column': 'sample_group',
            },
            {
                'rule_id': 'exclude_cancer',
                'stage': 'exclude',
                'priority': '20',
                'columns': 'sample_attribute',
                'pattern': 'cancer',
                'action': 'exclude',
                'target_column': 'exclusion',
                'outcome': 'bad_tissue',
            },
        ])
        select_rules = read_select_rules(str(rules_path))
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'scientific_name': ['Species A'],
            'sample_group': [''],
            'sample_attribute': ['cancer'],
            'bioproject': ['PRJ1'],
            'biosample': ['SAM1'],
            'total_spots': [100],
            'exclusion': ['no'],
        }))
        args = SimpleNamespace(
            min_nspots=0,
            mark_missing_rank='none',
            mark_redundant_biosamples=False,
            max_sample=20,
        )

        metadata = prepare_select_metadata(metadata, select_rules)
        out = apply_select_filters(metadata, args, select_rules)

        assert out.df.loc[0, 'exclusion'] == 'bad_tissue'
        assert out.df.loc[0, 'sample_group'] == 'cancer'
        assert 'sample_attribute' in out.df.columns

    def test_prepare_select_metadata_aggregates_sparse_sample_attribute_aliases(self, tmp_path):
        rules_path = tmp_path / 'select_rules.tsv'
        write_select_rules(rules_path, [
            {
                'rule_id': 'aggregate_source_name_alias',
                'stage': 'aggregate',
                'priority': '10',
                'columns': 'sample_attribute_source_name',
                'action': 'append',
                'target_column': 'source_name',
            },
            {
                'rule_id': 'aggregate_sample_title_alias',
                'stage': 'aggregate',
                'priority': '20',
                'columns': 'sample_attribute_sample_title',
                'action': 'append',
                'target_column': 'sample_title',
            },
            {
                'rule_id': 'normalize_root_from_source_name',
                'stage': 'normalize',
                'priority': '30',
                'columns': 'source_name',
                'pattern': r'\broot\b',
                'action': 'assign',
                'outcome': 'root',
            },
            {
                'rule_id': 'normalize_leaf_from_sample_title',
                'stage': 'normalize',
                'priority': '40',
                'columns': 'sample_title',
                'pattern': r'\bleaf\b',
                'action': 'assign',
                'outcome': 'leaf',
            },
        ])
        select_rules = read_select_rules(str(rules_path))
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001', 'SRR002'],
            'scientific_name': ['Species A', 'Species A'],
            'sample_group': ['', ''],
            'source_name': ['', 'existing source'],
            'sample_title': ['existing title', 'existing title'],
            'sample_attribute_source_name': ['primary root', ''],
            'sample_attribute_sample_title': ['', 'leaf replicate 1'],
            'bioproject': ['PRJ1', 'PRJ1'],
            'biosample': ['SAM1', 'SAM2'],
            'total_spots': [100, 100],
            'exclusion': ['no', 'no'],
        }))

        out = prepare_select_metadata(metadata, select_rules)

        assert out.df.loc[0, 'source_name'] == 'primary root'
        assert out.df.loc[0, 'sample_group'] == 'root'
        assert out.df.loc[0, 'sample_group_normalization_source'] == 'source_name'
        assert out.df.loc[1, 'sample_title'] == 'existing title; leaf replicate 1'
        assert out.df.loc[1, 'sample_group'] == 'leaf'
        assert out.df.loc[1, 'sample_group_normalization_source'] == 'sample_title'

    def test_apply_select_filters_prefers_non_excluded_duplicate_biosample(self, tmp_path):
        rules_path = tmp_path / 'select_rules.tsv'
        write_select_rules(rules_path, build_test_filter_dedup_rows())
        select_rules = read_select_rules(str(rules_path))
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001', 'SRR002'],
            'scientific_name': ['Species A', 'Species A'],
            'sample_group': ['flower', 'flower'],
            'bioproject': ['PRJ1', 'PRJ1'],
            'biosample': ['SAM1', 'SAM1'],
            'total_spots': [100, 6000000],
            'exclusion': ['no', 'no'],
        }))
        args = SimpleNamespace(
            min_nspots=1000000,
            mark_missing_rank='none',
            mark_redundant_biosamples=True,
            max_sample=10,
        )

        out = apply_select_filters(metadata, args, select_rules)

        assert out.df.loc[out.df['run'] == 'SRR001', 'exclusion'].iloc[0] == 'low_nspots'
        assert out.df.loc[out.df['run'] == 'SRR002', 'exclusion'].iloc[0] == 'no'
        assert out.df.loc[out.df['run'] == 'SRR002', 'is_sampled'].iloc[0] == 'yes'

    def test_apply_select_filters_keeps_highest_depth_duplicate_biosample(self, tmp_path):
        rules_path = tmp_path / 'select_rules.tsv'
        write_select_rules(rules_path, build_test_filter_dedup_rows())
        select_rules = read_select_rules(str(rules_path))
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001', 'SRR002'],
            'scientific_name': ['Species A', 'Species A'],
            'sample_group': ['flower', 'flower'],
            'bioproject': ['PRJ1', 'PRJ1'],
            'biosample': ['SAM1', 'SAM1'],
            'total_spots': [6000000, 7000000],
            'exclusion': ['no', 'no'],
        }))
        args = SimpleNamespace(
            min_nspots=1000000,
            mark_missing_rank='none',
            mark_redundant_biosamples=True,
            max_sample=10,
        )

        out = apply_select_filters(metadata, args, select_rules)

        assert out.df.loc[out.df['run'] == 'SRR001', 'exclusion'].iloc[0] == 'redundant_biosample'
        assert out.df.loc[out.df['run'] == 'SRR002', 'exclusion'].iloc[0] == 'no'
        assert out.df.loc[out.df['run'] == 'SRR002', 'is_sampled'].iloc[0] == 'yes'

    def test_preserves_taxid_and_extra_columns_after_filtering(self, tmp_path):
        rules_path = tmp_path / 'select_rules.tsv'
        write_select_rules(rules_path, build_test_filter_dedup_rows() + [
            {
                'rule_id': 'leaf_whole',
                'stage': 'normalize',
                'priority': '10',
                'columns': 'sample_attribute_tissue,sample_group',
                'pattern': r'\bleaf\b',
                'action': 'assign',
                'outcome': 'leaf',
            },
        ])
        select_rules = read_select_rules(str(rules_path))
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'scientific_name': ['Species A'],
            'sample_group': ['leaf'],
            'bioproject': ['PRJ1'],
            'biosample': ['SAM1'],
            'total_spots': [100],
            'exclusion': ['no'],
            'taxid_species': [12345],
            'sample_attribute_tissue': ['leaf'],
            'sample_group_normalization_source': ['sample_attribute_tissue'],
        }))
        args = SimpleNamespace(
            min_nspots=0,
            mark_missing_rank='species',
            mark_redundant_biosamples=False,
            max_sample=10,
        )

        metadata = prepare_select_metadata(metadata, select_rules)
        out = apply_select_filters(metadata, args, select_rules)

        assert 'taxid_species' in out.df.columns
        assert out.df.loc[0, 'taxid_species'] == 12345
        assert 'sample_attribute_tissue' in out.df.columns
        assert out.df.loc[0, 'sample_attribute_tissue'] == 'leaf'
        assert 'sample_group_normalization_source' in out.df.columns
        assert out.df.loc[0, 'sample_group_normalization_source'] == 'sample_attribute_tissue'

    def test_apply_select_filters_marks_empty_sample_group_as_no_sample_group(self, tmp_path):
        rules_path = tmp_path / 'select_rules.tsv'
        write_select_rules(rules_path, build_test_filter_dedup_rows())
        select_rules = read_select_rules(str(rules_path))
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'scientific_name': ['Species A'],
            'sample_group': [''],
            'bioproject': ['PRJ1'],
            'biosample': ['SAM1'],
            'total_spots': [100],
            'exclusion': ['no'],
        }))
        args = SimpleNamespace(
            min_nspots=0,
            mark_missing_rank='none',
            mark_redundant_biosamples=False,
            max_sample=20,
        )

        out = apply_select_filters(metadata, args, select_rules)

        assert out.df.loc[0, 'exclusion'] == 'no_sample_group'
        assert out.df.loc[0, 'is_qualified'] == 'no'

    def test_full_bloom_flower_rule_ignores_sample_description_context(self, tmp_path):
        rules_path = tmp_path / 'select_rules.tsv'
        write_select_rules(rules_path, [
            {
                'rule_id': 'flower_whole_full_bloom',
                'stage': 'normalize',
                'priority': '10',
                'columns': 'sample_title,source_name,exp_title',
                'pattern': r'\b(?:flower|flowers)\b(?![\s_-]?buds?)\b.{0,40}\bfull[- ]bloom\b|\bfull[- ]bloom\b.{0,40}\b(?:flower|flowers)\b(?![\s_-]?buds?)',
                'action': 'assign',
                'outcome': 'flower',
            },
            {
                'rule_id': 'fruit_non_target',
                'stage': 'normalize',
                'priority': '20',
                'columns': 'sample_title,source_name,exp_title,sample_description',
                'pattern': r'\bfruit\b',
                'action': 'assign',
                'outcome': 'non_target',
            },
        ])
        select_rules = read_select_rules(str(rules_path))
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'scientific_name': ['Species A'],
            'sample_group': [''],
            'sample_title': ['Plant sample from Species A'],
            'exp_title': ['RNAseq of Species A fruits 5 DAFB'],
            'sample_description': ['Flowers were collected at full bloom before fruit sampling.'],
            'bioproject': ['PRJ1'],
            'biosample': ['SAM1'],
            'total_spots': [100],
            'exclusion': ['no'],
        }))

        out = prepare_select_metadata(metadata, select_rules)

        assert out.df.loc[0, 'sample_group'] == 'non_target'
        assert out.df.loc[0, 'sample_group_normalization_status'] == 'non_target'
        assert out.df.loc[0, 'sample_group_normalization_rule_id'] == 'fruit_non_target'

    def test_default_plantae_pooled_multi_tissue_is_review(self):
        rules_path = Path(__file__).resolve().parents[1] / 'amalgkit' / 'select_rule_sets' / 'plantae' / 'select_rules.tsv'
        select_rules = read_select_rules(str(rules_path))
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR037746'],
            'scientific_name': ['Oryza sativa'],
            'sample_group': [''],
            'tissue': [
                'sample pooled from callus,seedling shoot,seedling root,'
                'tillering leaf,flowering panicle,flowering leaf,filling panicle'
            ],
            'sample_description': ['small RNA from mixture of 8 tissues'],
            'bioproject': ['PRJ1'],
            'biosample': ['SAM1'],
            'total_spots': [1000000],
            'exclusion': ['no'],
        }))

        out = prepare_select_metadata(metadata, select_rules)

        assert out.df.loc[0, 'sample_group'] == 'review'
        assert out.df.loc[0, 'sample_group_normalization_status'] == 'review'
        assert out.df.loc[0, 'sample_group_normalization_rule_id'] == 'normalize_review_pooled_multi_tissue'

    def test_default_plantae_same_organ_pooled_leaf_is_retained(self):
        rules_path = Path(__file__).resolve().parents[1] / 'amalgkit' / 'select_rule_sets' / 'plantae' / 'select_rules.tsv'
        select_rules = read_select_rules(str(rules_path))
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['ERR268475'],
            'scientific_name': ['Quercus robur'],
            'sample_group': [''],
            'sample_description': ['pooled oligo-dt RNA from leaves of two Q. robur clones'],
            'bioproject': ['PRJ1'],
            'biosample': ['SAM1'],
            'total_spots': [1000000],
            'exclusion': ['no'],
        }))

        out = prepare_select_metadata(metadata, select_rules)

        assert out.df.loc[0, 'sample_group'] == 'leaf'
        assert out.df.loc[0, 'sample_group_normalization_status'] == 'organ'
        assert out.df.loc[0, 'sample_group_normalization_rule_id'] == 'normalize_leaf_strict_sample_text'

    def test_default_plantae_delimited_sample_text_is_recovered(self):
        rules_path = Path(__file__).resolve().parents[1] / 'amalgkit' / 'select_rule_sets' / 'plantae' / 'select_rules.tsv'
        select_rules = read_select_rules(str(rules_path))
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRRROOT', 'SRRLEAF', 'SRRFLOWER'],
            'scientific_name': ['Species A', 'Species A', 'Species A'],
            'sample_group': ['', '', ''],
            'sample_description': [
                'TRANSCRIPTOMIC RNA from Castanea dentata Ellis-1, Root JAQP-RNAseq',
                'TRANSCRIPTOMIC RNA from Castanea dentata Ellis-1, Leaf JAQQ-RNAseq',
                'RNA from Species A, Flower JAQF-RNAseq',
            ],
            'bioproject': ['PRJ1', 'PRJ1', 'PRJ1'],
            'biosample': ['SAM1', 'SAM2', 'SAM3'],
            'total_spots': [1000000, 1000000, 1000000],
            'exclusion': ['no', 'no', 'no'],
        }))

        out = prepare_select_metadata(metadata, select_rules)

        groups = dict(zip(out.df['run'], out.df['sample_group']))
        rule_ids = dict(zip(out.df['run'], out.df['sample_group_normalization_rule_id']))
        assert groups == {
            'SRRROOT': 'root',
            'SRRLEAF': 'leaf',
            'SRRFLOWER': 'flower',
        }
        assert rule_ids == {
            'SRRROOT': 'normalize_root_delimited_sample_text',
            'SRRLEAF': 'normalize_leaf_delimited_sample_text',
            'SRRFLOWER': 'normalize_flower_delimited_sample_text',
        }

    def test_default_plantae_sample_specific_tissue_overrides_project_level_group(self):
        rules_path = Path(__file__).resolve().parents[1] / 'amalgkit' / 'select_rule_sets' / 'plantae' / 'select_rules.tsv'
        select_rules = read_select_rules(str(rules_path))
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRRROOT', 'SRRLEAF', 'SRRFLOWER'],
            'scientific_name': ['Species A', 'Species A', 'Species A'],
            'sample_group': ['flower', 'flower', 'leaf'],
            'sample_attribute_tissue': ['Root', 'Leaf', 'Inflorescence'],
            'exp_title': [
                'RNA-seq of Species A flower project',
                'RNA-seq of Species A flower project',
                'RNA-seq of Species A leaf project',
            ],
            'bioproject': ['PRJ1', 'PRJ1', 'PRJ1'],
            'biosample': ['SAM1', 'SAM2', 'SAM3'],
            'total_spots': [1000000, 1000000, 1000000],
            'exclusion': ['no', 'no', 'no'],
        }))

        out = prepare_select_metadata(metadata, select_rules)

        groups = dict(zip(out.df['run'], out.df['sample_group']))
        rule_ids = dict(zip(out.df['run'], out.df['sample_group_normalization_rule_id']))
        assert groups == {
            'SRRROOT': 'root',
            'SRRLEAF': 'leaf',
            'SRRFLOWER': 'flower',
        }
        assert rule_ids == {
            'SRRROOT': 'normalize_root_sample_specific_text',
            'SRRLEAF': 'normalize_leaf_sample_specific_text',
            'SRRFLOWER': 'normalize_flower_sample_specific_text',
        }

    def test_default_plantae_strong_tissue_conflict_is_review(self):
        rules_path = Path(__file__).resolve().parents[1] / 'amalgkit' / 'select_rule_sets' / 'plantae' / 'select_rules.tsv'
        select_rules = read_select_rules(str(rules_path))
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRRCONFLICT'],
            'scientific_name': ['Species A'],
            'sample_group': ['flower'],
            'sample_attribute_tissue': ['Leaf'],
            'lib_name': ['Species_A_RNA_flower_rep1'],
            'exp_title': ['RNA-seq of Species A flower'],
            'bioproject': ['PRJ1'],
            'biosample': ['SAM1'],
            'total_spots': [1000000],
            'exclusion': ['no'],
        }))

        out = prepare_select_metadata(metadata, select_rules)

        assert out.df.loc[0, 'sample_group'] == 'review'
        assert out.df.loc[0, 'sample_group_normalization_status'] == 'review'
        assert out.df.loc[0, 'sample_group_normalization_source'] == 'sample_attribute_tissue'
        assert out.df.loc[0, 'sample_group_normalization_rule_id'] == 'normalize_review_cross_field_tissue_conflict'
        assert 'previous_rule_id=normalize_flower_sample_specific_text' in out.df.loc[0, 'sample_group_normalization_text']

    def test_default_plantae_assign_safe_source_is_not_cross_field_conflict(self):
        rules_path = Path(__file__).resolve().parents[1] / 'amalgkit' / 'select_rule_sets' / 'plantae' / 'select_rules.tsv'
        select_rules = read_select_rules(str(rules_path))
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRRSAFE'],
            'scientific_name': ['Species A'],
            'sample_group': [''],
            'sample_attribute_tissue': ['leaves (covering the inflorescence)'],
            'bioproject': ['PRJ1'],
            'biosample': ['SAM1'],
            'total_spots': [1000000],
            'exclusion': ['no'],
        }))

        out = prepare_select_metadata(metadata, select_rules)

        assert out.df.loc[0, 'sample_group'] == 'leaf'
        assert out.df.loc[0, 'sample_group_normalization_status'] == 'organ'
        assert out.df.loc[0, 'sample_group_normalization_rule_id'] == 'normalize_leaf_safe_covering_inflorescence'

    def test_default_plantae_sample_specific_reproductive_suborgans_are_not_whole_flower(self):
        rules_path = Path(__file__).resolve().parents[1] / 'amalgkit' / 'select_rule_sets' / 'plantae' / 'select_rules.tsv'
        select_rules = read_select_rules(str(rules_path))
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRRPOLLEN', 'SRRSPIKELET'],
            'scientific_name': ['Species A', 'Species A'],
            'sample_group': ['leaf', 'leaf'],
            'sample_attribute_tissue': ['pollen', 'Spikelet'],
            'sample_title': ['Leaf sample collected after heat treatment', 'RNA-seq of rice spikelet'],
            'exp_title': ['RNA-seq of pollen sample', 'RNA-seq of leaf project'],
            'bioproject': ['PRJ1', 'PRJ1'],
            'biosample': ['SAM1', 'SAM2'],
            'total_spots': [1000000, 1000000],
            'exclusion': ['no', 'no'],
        }))

        out = prepare_select_metadata(metadata, select_rules)

        groups = dict(zip(out.df['run'], out.df['sample_group']))
        statuses = dict(zip(out.df['run'], out.df['sample_group_normalization_status']))
        rule_ids = dict(zip(out.df['run'], out.df['sample_group_normalization_rule_id']))
        assert groups == {
            'SRRPOLLEN': 'review',
            'SRRSPIKELET': 'flower',
        }
        assert statuses == {
            'SRRPOLLEN': 'review',
            'SRRSPIKELET': 'organ',
        }
        assert rule_ids == {
            'SRRPOLLEN': 'normalize_review_floral_suborgan_sample_specific_text',
            'SRRSPIKELET': 'normalize_flower_sample_specific_text',
        }

    def test_default_plantae_normalize_rules_use_whitelisted_columns(self):
        rules_path = Path(__file__).resolve().parents[1] / 'amalgkit' / 'select_rule_sets' / 'plantae' / 'select_rules.tsv'
        select_rules = read_select_rules(str(rules_path))
        base_allowed = {
            'sample_attribute_tissue',
            'tissue',
            'sample_group',
            'sample_title',
            'source_name',
            'exp_title',
        }
        sample_specific_allowed = {
            'organism_part',
            'organ',
            'plant_structure',
            'tissue_type',
            'tissue_types',
            'tissues',
            'tissue_source',
            'sample_attribute_source_name',
            'isolation_source',
            'sample_attribute_sample_title',
            'sample_name',
            'lib_name',
        }
        strict_sample_text_allowed = {
            'sample_description',
            'sample_attribute_sample_description',
            'sample_attribute_source_name',
            'sample_attribute_sample_title',
            'isolation_source',
            'sample_name',
            'library_name',
        }
        normalize_rules = [rule for rule in select_rules if rule['stage'] == 'normalize']
        assert len(normalize_rules) > 0
        for rule in normalize_rules:
            allowed = base_allowed
            if rule['rule_id'].endswith('_sample_specific_text'):
                allowed = base_allowed | sample_specific_allowed
            if rule['rule_id'].endswith(('_strict_sample_text', '_delimited_sample_text')):
                allowed = base_allowed | strict_sample_text_allowed
            if rule['rule_id'] == 'normalize_review_pooled_multi_tissue':
                allowed = base_allowed | strict_sample_text_allowed
            assert set(rule['columns']).issubset(allowed), rule['rule_id']

    def test_default_plantae_config_contains_manual_recovery_and_filter_rules(self):
        rules_path = Path(__file__).resolve().parents[1] / 'amalgkit' / 'select_rule_sets' / 'plantae' / 'select_rules.tsv'
        config = read_select_config(str(rules_path))
        select_rules = config['rules']
        rule_ids = {rule['rule_id'] for rule in select_rules}
        validate_rule_ids = {rule['rule_id'] for rule in select_rules if rule['stage'] == 'validate'}
        assert config['parameters']['sampling_strategy'] == 'maximize_bioproject_diversity'
        assert 'normalize_leaf_safe_covering_inflorescence' in rule_ids
        assert 'normalize_leaf_safe_floral_induction_context' in rule_ids
        assert 'normalize_root_safe_leaf_stage_context' in rule_ids
        assert 'normalize_leaf_safe_leaf_blade_and_sheath' in rule_ids
        assert 'normalize_review_pooled_multi_tissue' in rule_ids
        assert 'normalize_flower_delimited_sample_text' in rule_ids
        assert 'normalize_leaf_delimited_sample_text' in rule_ids
        assert 'normalize_root_delimited_sample_text' in rule_ids
        assert 'normalize_flower_sample_specific_text' in rule_ids
        assert 'normalize_leaf_sample_specific_text' in rule_ids
        assert 'normalize_root_sample_specific_text' in rule_ids
        assert 'normalize_review_floral_suborgan_sample_specific_text' in rule_ids
        assert 'filter_low_nspots_1' in rule_ids
        assert 'filter_missing_taxid_2' in rule_ids
        assert 'filter_no_sample_group_3' in rule_ids
        assert 'dedup_redundant_biosample_1' in rule_ids
        assert 'exclude_single_cell_1' in rule_ids
        assert 'exclude_single_cell_1b' in rule_ids
        assert 'exclude_single_nucleus_1a' in rule_ids
        assert 'exclude_single_nucleus_1b' in rule_ids
        assert 'exclude_protoplast_1c' in rule_ids
        assert 'exclude_three_prime_biased_8a' in rule_ids
        assert 'exclude_three_prime_biased_8b' in rule_ids
        assert 'exclude_ribosome_profiling_8c' in rule_ids
        assert 'exclude_bisulfite_or_genome_8d' in rule_ids
        assert 'exclude_non_rna_seq_library_selection_8e' in rule_ids
        assert 'exclude_rnai_4' in rule_ids
        assert 'exclude_cage_8' in rule_ids
        assert 'exclude_cage_8b' in rule_ids
        assert 'exclude_chip_or_rip_3b' in rule_ids
        assert 'control_mock_1' in rule_ids
        assert validate_rule_ids == {
            'validate_hint_flower_1',
            'validate_hint_leaf_2',
            'validate_hint_root_3',
            'validate_hint_review_4',
            'validate_ignore_safe_metadata_5',
        }

    def test_default_plantae_covering_inflorescence_phrase_is_leaf(self):
        rules_path = Path(__file__).resolve().parents[1] / 'amalgkit' / 'select_rule_sets' / 'plantae' / 'select_rules.tsv'
        select_rules = read_select_rules(str(rules_path))
        normalize_rules = [rule for rule in select_rules if rule['stage'] == 'normalize']

        result = classify_select_text(
            'leaves (covering the inflorescence)',
            normalize_rules,
        )

        assert result['status'] == 'organ'
        assert result['organ'] == 'leaf'
        assert result['rule_id'] == 'normalize_leaf_safe_covering_inflorescence'

    def test_default_plantae_excludes_quantseq_and_digital_gene_expression(self):
        rules_path = Path(__file__).resolve().parents[1] / 'amalgkit' / 'select_rule_sets' / 'plantae' / 'select_rules.tsv'
        select_rules = read_select_rules(str(rules_path))
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001', 'SRR002', 'SRR003'],
            'scientific_name': ['Species A'] * 3,
            'sample_group': ['leaf', 'leaf', 'leaf'],
            'sample_title': ['leaf sample 1', 'leaf sample 2', 'leaf sample 3'],
            'sample_description': ['', '', ''],
            'study_title': ['', 'Digital Gene Expression Analysis during Floral Transition', ''],
            'exp_title': ['RNA-seq of leaf with Lexogen QuantSeq 3 mRNA-Seq FWD', '', 'standard bulk RNA-seq'],
            'design': ['', '', ''],
            'lib_name': ['', '', ''],
            'protocol': ['', '3 mRNA-Seq library', 'TruSeq stranded mRNA library'],
            'bioproject': ['PRJ1'] * 3,
            'biosample': ['SAM1', 'SAM2', 'SAM3'],
            'total_spots': [100, 100, 100],
            'exclusion': ['no', 'no', 'no'],
        }))
        args = SimpleNamespace(
            min_nspots=0,
            mark_missing_rank='none',
            mark_redundant_biosamples=False,
            max_sample=20,
        )

        metadata = prepare_select_metadata(metadata, select_rules)
        out = apply_select_filters(metadata, args, select_rules)

        assert out.df.loc[out.df['run'] == 'SRR001', 'exclusion'].iloc[0] == 'three_prime_biased'
        assert out.df.loc[out.df['run'] == 'SRR002', 'exclusion'].iloc[0] == 'three_prime_biased'
        assert out.df.loc[out.df['run'] == 'SRR003', 'exclusion'].iloc[0] == 'no'

    def test_default_plantae_excludes_dge_abbreviation_as_three_prime_biased(self):
        rules_path = Path(__file__).resolve().parents[1] / 'amalgkit' / 'select_rule_sets' / 'plantae' / 'select_rules.tsv'
        select_rules = read_select_rules(str(rules_path))
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRRDGE'],
            'scientific_name': ['Species A'],
            'sample_group': ['leaf'],
            'sample_title': ['DGE analysis of leaf'],
            'sample_description': [''],
            'study_title': [''],
            'exp_title': ['GSM000000: leaf (DGE); Species A; RNA-Seq'],
            'design': [''],
            'lib_name': [''],
            'protocol': [''],
            'bioproject': ['PRJ1'],
            'biosample': ['SAM1'],
            'total_spots': [100],
            'exclusion': ['no'],
        }))
        args = SimpleNamespace(
            min_nspots=0,
            mark_missing_rank='none',
            mark_redundant_biosamples=False,
            max_sample=10,
        )

        metadata = prepare_select_metadata(metadata, select_rules)
        out = apply_select_filters(metadata, args, select_rules)

        assert out.df.loc[out.df['run'] == 'SRRDGE', 'exclusion'].iloc[0] == 'three_prime_biased'

    def test_default_plantae_excludes_single_cell_and_protoplast_without_excluding_lexogen_sense(self):
        rules_path = Path(__file__).resolve().parents[1] / 'amalgkit' / 'select_rule_sets' / 'plantae' / 'select_rules.tsv'
        select_rules = read_select_rules(str(rules_path))
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': [
                'SRRSC',
                'SRRPRO',
                'SRR3P',
                'SRRLEX',
                'SRRSC_SOURCE',
                'SRRPRO_CELL',
                'SRRNUCLEI',
                'SRRCAGE',
                'SRRRNATAG',
                'SRRRIBO',
                'SRRCHIP',
            'SRRWGBS',
            'SRRRACE',
            'SRRSMALL',
        ],
            'scientific_name': ['Species A'] * 14,
            'sample_group': ['leaf'] * 14,
            'sample_title': [
                'leaf sample',
                'leaf protoplast sample',
                'leaf 3 prime sample',
                'leaf bulk sample',
                'leaf sample',
                'leaf sample',
                'sorted root nuclei',
                'leaf sample',
                'leaf sample',
                'leaf sample',
                'leaf sample',
                'leaf sample',
                'leaf sample',
                'leaf small RNA',
            ],
            'sample_attribute_tissue': [
                'leaf',
                'Leaf protoplasts',
                'leaf',
                'leaf',
                'leaf',
                'leaf',
                'root',
                'leaf',
                'leaf',
                'leaf',
                'leaf',
                'leaf',
                'leaf',
                'leaf',
            ],
            'sample_description': [''] * 14,
            'study_title': ['Single-cell transcriptome atlas of leaves', '', '', '', '', '', '', '', '', '', '', '', '', ''],
            'exp_title': [
                '',
                '',
                "3' RNA-seq of leaf",
                'bulk RNA-seq',
                'standard leaf RNA-seq',
                'standard leaf RNA-seq',
                'standard leaf RNA-seq',
                'RAMPAGE of leaf sample',
                '',
                'Ribo-seq of leaf',
                'standard leaf RNA-seq',
                'WGBS sample',
                'standard leaf RNA-seq',
                'Catharanthus roseus leaf small RNA',
            ],
            'design': [
                '',
                '',
                '',
                'Library was prepared with Sense Total RNA-Seq Library Prep Kit, Lexogen',
                '',
                '',
                '',
                '',
                'Libraries were generated using the RNAtag-seq protocol',
                '',
                '',
                '',
                '',
                'Small RNA libraries were constructed using a small RNA sequencing kit',
            ],
            'lib_name': [''] * 14,
            'protocol': [''] * 14,
            'cell_type': ['', '', '', '', '', 'protoplast', '', '', '', '', '', '', '', ''],
            'lib_source': ['', '', '', '', 'TRANSCRIPTOMIC SINGLE CELL', '', '', '', '', '', '', '', '', ''],
            'lib_selection': ['', '', '', '', '', '', '', 'CAGE', '', '', 'ChIP', '', 'RACE', 'size fractionation'],
            'bioproject': ['PRJ1'] * 14,
            'biosample': [f'SAM{i}' for i in range(1, 15)],
            'total_spots': [100] * 14,
            'exclusion': ['no'] * 14,
        }))
        args = SimpleNamespace(
            min_nspots=0,
            mark_missing_rank='none',
            mark_redundant_biosamples=False,
            max_sample=20,
        )

        metadata = prepare_select_metadata(metadata, select_rules)
        out = apply_select_filters(metadata, args, select_rules)

        exclusions = dict(zip(out.df['run'], out.df['exclusion']))
        assert exclusions == {
            'SRRSC': 'single_cell',
            'SRRPRO': 'protoplast',
            'SRR3P': 'three_prime_biased',
            'SRRLEX': 'no',
            'SRRSC_SOURCE': 'single_cell',
            'SRRPRO_CELL': 'protoplast',
            'SRRNUCLEI': 'single_nucleus',
            'SRRCAGE': 'cage',
            'SRRRNATAG': 'three_prime_biased',
            'SRRRIBO': 'ribosome_profiling',
            'SRRCHIP': 'immunoprecipitation',
            'SRRWGBS': 'non_rna_seq_library',
            'SRRRACE': 'non_rna_seq_library',
            'SRRSMALL': 'small_rna',
        }

    def test_review_rule_overrides_original_sample_group(self, tmp_path):
        rules_path = tmp_path / 'select_rules.tsv'
        write_select_rules(rules_path, [
            {
                'rule_id': 'review_structure',
                'stage': 'normalize',
                'priority': '10',
                'columns': 'sample_title,exp_title,sample_group',
                'pattern': r'\bbud\b',
                'action': 'assign',
                'outcome': 'review',
            },
            {
                'rule_id': 'flower_whole',
                'stage': 'normalize',
                'priority': '20',
                'columns': 'sample_title,exp_title,sample_group',
                'pattern': r'\bflower\b',
                'action': 'assign',
                'outcome': 'flower',
            },
        ])
        select_rules = read_select_rules(str(rules_path))
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'scientific_name': ['Species A'],
            'sample_group': ['flower'],
            'exp_title': ['flowers at bud stage'],
            'bioproject': ['PRJ1'],
            'biosample': ['SAM1'],
            'total_spots': [100],
            'exclusion': ['no'],
        }))

        out = prepare_select_metadata(metadata, select_rules)

        assert out.df.loc[0, 'sample_group_original'] == 'flower'
        assert out.df.loc[0, 'sample_group'] == 'review'
        assert out.df.loc[0, 'sample_group_normalization_status'] == 'review'
        assert out.df.loc[0, 'sample_group_normalization_rule_id'] == 'review_structure'


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
        assert args.out_dir == raw_out_dir


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
