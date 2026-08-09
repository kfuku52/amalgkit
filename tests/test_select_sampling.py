import numpy
import pandas
import pytest


from amalgkit.select import (
    classify_select_text,
    prepare_select_metadata,
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
        assert selected_runs == ['SRR003', 'SRR004', 'SRR005']

    def test_random_sampling_is_reproducible_and_records_seed(self):
        first = self._build_sampling_metadata()
        second = self._build_sampling_metadata()
        different = self._build_sampling_metadata()

        first.label_sampled_data(max_sample=3, sampling_strategy='random', random_seed=17)
        second.label_sampled_data(max_sample=3, sampling_strategy='random', random_seed=17)
        different.label_sampled_data(max_sample=3, sampling_strategy='random', random_seed=18)

        first_runs = first.df.loc[first.df['is_sampled'] == 'yes', 'run'].tolist()
        second_runs = second.df.loc[second.df['is_sampled'] == 'yes', 'run'].tolist()
        different_runs = different.df.loc[different.df['is_sampled'] == 'yes', 'run'].tolist()
        assert first_runs == second_runs
        assert first_runs != different_runs
        assert first.df['sampling_seed'].unique().tolist() == [17]
        assert first.df['sampling_strategy'].unique().tolist() == ['random']

    def test_random_sampling_rejects_negative_seed(self):
        metadata = self._build_sampling_metadata()

        with pytest.raises(ValueError, match='non-negative integer'):
            metadata.label_sampled_data(
                max_sample=3,
                sampling_strategy='random',
                random_seed=-1,
            )

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
