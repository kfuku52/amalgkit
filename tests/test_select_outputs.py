import os
import pandas
import pytest

from types import SimpleNamespace

from amalgkit.select import (
    apply_select_aggregate_rules,
    apply_select_config_parameters,
    apply_select_filter_rules,
    read_select_config,
    read_select_rules,
    write_select_outputs,
    resolve_select_rules_tsv,
    filter_metadata_by_sample_group,
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

    def test_rejects_symbolic_link_original_even_when_preserving_it(self, tmp_path, sample_metadata):
        metadata_dir = tmp_path / 'metadata'
        metadata_dir.mkdir()
        outside = tmp_path / 'outside.tsv'
        outside.write_text('do not replace\n', encoding='utf-8')
        path_original = metadata_dir / 'metadata_original.tsv'
        path_original.symlink_to(outside)
        path_table = metadata_dir / 'metadata.tsv'
        sample_metadata.label_sampled_data(max_sample=10)

        with pytest.raises(ValueError, match='symbolic-link output file'):
            write_select_outputs(
                str(path_original),
                str(path_table),
                str(metadata_dir),
                sample_metadata,
                preserve_existing_original=True,
            )

        assert outside.read_text(encoding='utf-8') == 'do not replace\n'


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


def test_aggregate_rules_are_idempotent_and_deduplicate_existing_tokens(tmp_path):
    rules_path = tmp_path / 'select_rules.tsv'
    write_select_rules(
        rules_path,
        [
            {
                'rule_id': 'aggregate_tissue',
                'stage': 'aggregate',
                'priority': '10',
                'columns': 'sample_attribute_tissue',
                'action': 'append',
                'target_column': 'source_name',
            },
        ],
    )
    rules = read_select_rules(str(rules_path))
    frame = pandas.DataFrame(
        {
            'source_name': ['leaf; root'],
            'sample_attribute_tissue': ['root; flower'],
        }
    )

    once = apply_select_aggregate_rules(frame, rules)
    twice = apply_select_aggregate_rules(once, rules)

    assert once.loc[0, 'source_name'] == 'leaf; root; flower'
    pandas.testing.assert_frame_equal(once, twice)


@pytest.mark.parametrize(
    ('first_stop_on_match', 'expected'),
    [
        ('yes', 'low_nspots'),
        ('no', 'no_sample_group'),
    ],
)
def test_filter_rules_honor_stop_on_match(tmp_path, first_stop_on_match, expected):
    rules_path = tmp_path / 'select_rules.tsv'
    write_select_rules(
        rules_path,
        [
            {
                'rule_id': 'low_depth',
                'stage': 'filter',
                'priority': '10',
                'columns': 'total_spots',
                'action': 'exclude_if_lt_parameter',
                'target_column': 'exclusion',
                'outcome': 'low_nspots',
                'parameter_name': 'min_nspots',
                'stop_on_match': first_stop_on_match,
            },
            {
                'rule_id': 'empty_group',
                'stage': 'filter',
                'priority': '20',
                'columns': 'sample_group',
                'action': 'exclude_if_empty',
                'target_column': 'exclusion',
                'outcome': 'no_sample_group',
                'stop_on_match': 'yes',
            },
        ],
    )
    rules = read_select_rules(str(rules_path))
    metadata = Metadata.from_DataFrame(
        pandas.DataFrame(
            {
                'run': ['SRR001'],
                'sample_group': [''],
                'total_spots': [10],
                'exclusion': ['no'],
            }
        )
    )

    out = apply_select_filter_rules(
        metadata,
        SimpleNamespace(min_nspots=100),
        rules,
    )

    assert out.df.loc[0, 'exclusion'] == expected


def test_low_spot_filter_excludes_invalid_nonfinite_and_nonpositive_values(tmp_path):
    rules_path = tmp_path / 'select_rules.tsv'
    write_select_rules(
        rules_path,
        [
            {
                'rule_id': 'low_depth',
                'stage': 'filter',
                'priority': '10',
                'columns': 'total_spots',
                'action': 'exclude_if_lt_parameter',
                'target_column': 'exclusion',
                'outcome': 'low_nspots',
                'parameter_name': 'min_nspots',
                'stop_on_match': 'yes',
            },
        ],
    )
    metadata = Metadata.from_DataFrame(pandas.DataFrame({
        'run': ['missing', 'text', 'zero', 'negative', 'infinite', 'low', 'valid'],
        'total_spots': [None, 'unknown', 0, -1, float('inf'), 99, 100],
        'exclusion': ['no'] * 7,
    }))

    out = apply_select_filter_rules(
        metadata,
        SimpleNamespace(min_nspots=100),
        read_select_rules(str(rules_path)),
    )

    assert out.df['exclusion'].tolist() == [
        'low_nspots',
        'low_nspots',
        'low_nspots',
        'low_nspots',
        'low_nspots',
        'low_nspots',
        'no',
    ]
