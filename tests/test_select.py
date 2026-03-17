import os
import subprocess
import sys
import pandas
import pytest

from types import SimpleNamespace
from pathlib import Path

from amalgkit.select import (
    apply_select_control_rules,
    apply_select_filters,
    classify_select_text,
    prepare_select_metadata,
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
        'note': '',
    }
    normalized_rows = []
    for row in rows:
        normalized = defaults.copy()
        normalized.update(row)
        normalized_rows.append(normalized)
    pandas.DataFrame(normalized_rows, columns=SELECT_RULE_COLUMNS).to_csv(path, sep='\t', index=False)


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

    @pytest.mark.parametrize(
        ('text', 'expected_status', 'expected_organ'),
        [
            ('flower', 'organ', 'flower'),
            ('Flower of Bora mutant(purple-petal) at full-bloom stage rep1', 'organ', 'flower'),
            ('root tips', 'organ', 'root'),
            ('leaf blade', 'organ', 'leaf'),
            ('spikes in seedling stage', 'organ', 'flower'),
            ('seedling roots', 'organ', 'root'),
            ('seedling root', 'organ', 'root'),
            ('unopened flower buds', 'organ', 'flower'),
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
                'rule_id': 'flower_whole_full_bloom',
                'stage': 'normalize',
                'priority': '18',
                'columns': 'sample_group',
                'pattern': r'\b(?:flower|flowers)\b(?![\s_-]?buds?)\b.{0,40}\bfull[- ]bloom\b|\bfull[- ]bloom\b.{0,40}\b(?:flower|flowers)\b(?![\s_-]?buds?)',
                'action': 'assign',
                'outcome': 'flower',
            },
            {
                'rule_id': 'flower_whole_unopened_flower_buds',
                'stage': 'normalize',
                'priority': '19',
                'columns': 'sample_group',
                'pattern': r'\b(?:unopened|closed)\b.{0,20}\bflower[\s_-]?buds?\b|\bflower[\s_-]?buds?\b.{0,20}\b(?:unopened|closed)\b',
                'action': 'assign',
                'outcome': 'flower',
            },
            {
                'rule_id': 'review_suborgan',
                'stage': 'normalize',
                'priority': '20',
                'columns': 'sample_group',
                'pattern': r'\b(?:petal|corolla|anther|ovary|pistil|bract|flower[\s_-]?bud|petiole|root hair)\b',
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
                'rule_id': 'root_whole_seedling_root',
                'stage': 'normalize',
                'priority': '27',
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
                'pattern': r'\bflower\b',
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
        ])
        normalize_rules = [rule for rule in read_select_rules(str(rules_path)) if rule['stage'] == 'normalize']
        result = classify_select_text(text, normalize_rules)
        assert result['status'] == expected_status
        assert result['organ'] == expected_organ


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
            max_sample=10,
        )

        metadata = prepare_select_metadata(metadata, select_rules)
        out = apply_select_filters(metadata, args, select_rules)

        assert out.df.loc[0, 'exclusion'] == 'bad_tissue'
        assert out.df.loc[0, 'sample_group'] == 'cancer'
        assert 'sample_attribute' in out.df.columns

    def test_preserves_taxid_and_extra_columns_after_filtering(self, tmp_path):
        rules_path = tmp_path / 'select_rules.tsv'
        write_select_rules(rules_path, [
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

    def test_default_plantae_normalize_rules_use_whitelisted_columns(self):
        rules_path = Path(__file__).resolve().parents[1] / 'amalgkit' / 'config_dir' / 'plantae' / 'select_rules.tsv'
        select_rules = read_select_rules(str(rules_path))
        allowed = {
            'sample_attribute_tissue',
            'tissue',
            'sample_group',
            'sample_title',
            'source_name',
            'exp_title',
        }
        normalize_rules = [rule for rule in select_rules if rule['stage'] == 'normalize']
        assert len(normalize_rules) > 0
        for rule in normalize_rules:
            assert set(rule['columns']).issubset(allowed), rule['rule_id']

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
            'amalgkit.select.read_select_rules',
            lambda *_args, **_kwargs: (_ for _ in ()).throw(AssertionError('read_select_rules should not be called')),
        )

        with pytest.raises(NotADirectoryError, match='Output path exists but is not a directory'):
            select_main(args)

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

        def fake_read_select_rules(path_select_rules):
            observed['select_rules_tsv'] = path_select_rules
            return [{'stage': 'normalize', 'columns': ['sample_group'], 'regex': None}]

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

        monkeypatch.setattr('amalgkit.select.read_select_rules', fake_read_select_rules)
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
            sample_group='flower,leaf,root',
            min_nspots=0,
            mark_missing_rank='none',
            mark_redundant_biosamples=False,
            max_sample=2,
            species_tsv=str(species_tsv),
            metadata_specieswise_dir=str(metadata_specieswise_dir),
            summary_tsv='inferred',
            queue_tsv='inferred',
            manifest_tsv='inferred',
            batch_label='inferred',
        )

        select_main(args)

        normalization_summary = tmp_path / 'select_batch' / 'normalization_summary.tsv'
        select_summary = tmp_path / 'select_batch' / 'select_summary.tsv'
        select_queue = tmp_path / 'select_batch' / 'select_queue.tsv'
        manifest = tmp_path / 'select_batch' / 'external_manifest.tsv'
        manifest_strict = tmp_path / 'select_batch' / 'external_manifest_strict.tsv'
        manifest_relaxed = tmp_path / 'select_batch' / 'external_manifest_relaxed.tsv'

        assert normalization_summary.exists()
        assert select_summary.exists()
        assert select_queue.exists()
        assert manifest.exists()
        assert manifest_strict.exists()
        assert manifest_relaxed.exists()

        summary_df = pandas.read_csv(select_summary, sep='\t')
        assert set(summary_df['species_token'].tolist()) == {'Species_alpha', 'Species_beta'}
        queue_by_species = summary_df.set_index('species_token')['queue_tier'].to_dict()
        assert queue_by_species['Species_alpha'] == 'strict'
        assert queue_by_species['Species_beta'] == 'defer'

        manifest_df = pandas.read_csv(manifest, sep='\t')
        assert set(manifest_df['queue_tier'].tolist()) == {'strict', 'defer'}
        assert 'selected_metadata_path' in manifest_df.columns


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
