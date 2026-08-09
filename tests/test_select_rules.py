import pandas

from types import SimpleNamespace
from pathlib import Path

from amalgkit.select import (
    apply_select_control_rules,
    apply_select_filters,
    classify_select_text,
    prepare_select_metadata,
    read_select_config,
    read_select_rules,
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

    def test_default_plantae_medium_context_tissue_conflict_is_review(self):
        rules_path = Path(__file__).resolve().parents[1] / 'amalgkit' / 'select_rule_sets' / 'plantae' / 'select_rules.tsv'
        select_rules = read_select_rules(str(rules_path))
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRRROOTLIB'],
            'scientific_name': ['Species A'],
            'sample_group': [''],
            'sample_attribute_tissue': ['Leaf'],
            'lib_name': ['Species_A_RNA_root_rep1'],
            'exp_title': ['RNA-seq of Species A root'],
            'bioproject': ['PRJ1'],
            'biosample': ['SAM1'],
            'total_spots': [1000000],
            'exclusion': ['no'],
        }))

        out = prepare_select_metadata(metadata, select_rules)

        assert out.df.loc[0, 'sample_group'] == 'review'
        assert out.df.loc[0, 'sample_group_normalization_status'] == 'review'
        assert out.df.loc[0, 'sample_group_normalization_source'] == 'lib_name'
        assert out.df.loc[0, 'sample_group_normalization_rule_id'] == 'normalize_review_cross_field_tissue_context_conflict'
        assert 'previous_rule_id=normalize_leaf_sample_specific_text' in out.df.loc[0, 'sample_group_normalization_text']

    def test_default_plantae_medium_context_conflict_without_strong_columns_is_review(self):
        rules_path = Path(__file__).resolve().parents[1] / 'amalgkit' / 'select_rule_sets' / 'plantae' / 'select_rules.tsv'
        select_rules = read_select_rules(str(rules_path))
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRRROOTLIB'],
            'scientific_name': ['Species A'],
            'sample_group': [''],
            'sample_title': ['Leaf sample'],
            'lib_name': ['Species_A_RNA_root_rep1'],
            'bioproject': ['PRJ1'],
            'biosample': ['SAM1'],
            'total_spots': [1000000],
            'exclusion': ['no'],
        }))

        out = prepare_select_metadata(metadata, select_rules)

        assert out.df.loc[0, 'sample_group'] == 'review'
        assert out.df.loc[0, 'sample_group_normalization_source'] == 'lib_name'
        assert out.df.loc[0, 'sample_group_normalization_rule_id'] == 'normalize_review_cross_field_tissue_context_conflict'

    def test_default_plantae_safe_leaf_development_context_is_not_reviewed(self):
        rules_path = Path(__file__).resolve().parents[1] / 'amalgkit' / 'select_rule_sets' / 'plantae' / 'select_rules.tsv'
        select_rules = read_select_rules(str(rules_path))
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRRINFLO', 'SRRPANICLE', 'SRRLAMINA'],
            'scientific_name': ['Species A', 'Species A', 'Species A'],
            'sample_group': ['', '', ''],
            'sample_attribute_tissue': ['Leaf', 'FLAG LEAF', 'leaf'],
            'sample_title': [
                'Leaves from 1 meter inflorescence stage',
                'RNA-seq of Flag leaf at Panicle emergence stage',
                '2391567_lamina_T-4_1st_open_flower',
            ],
            'bioproject': ['PRJ1', 'PRJ1', 'PRJ1'],
            'biosample': ['SAM1', 'SAM2', 'SAM3'],
            'total_spots': [1000000, 1000000, 1000000],
            'exclusion': ['no', 'no', 'no'],
        }))

        out = prepare_select_metadata(metadata, select_rules)

        assert out.df['sample_group'].tolist() == ['leaf', 'leaf', 'leaf']
        assert out.df['sample_group_normalization_status'].tolist() == ['organ', 'organ', 'organ']
        assert out.df['sample_group_normalization_rule_id'].tolist() == [
            'normalize_leaf_sample_specific_text',
            'normalize_leaf_sample_specific_text',
            'normalize_leaf_sample_specific_text',
        ]

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
