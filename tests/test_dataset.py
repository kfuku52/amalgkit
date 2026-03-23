import os

import pytest

from amalgkit import dataset as dataset_module
from amalgkit.dataset import (
    export_rule_set,
    extract_dataset,
    list_available_rule_sets,
    validate_dataset_name,
    validate_rule_set_name,
    resolve_dataset_source_dir,
)
from amalgkit.exceptions import AmalgkitExit
from amalgkit.main import build_main_parser


def _prepare_fake_dataset(tmp_path, monkeypatch):
    dataset_root = tmp_path / 'datasets'
    source_dir = dataset_root / 'toy'
    source_dir.mkdir(parents=True)
    (source_dir / 'toy.fa.gz').write_text('FA')
    (source_dir / 'toy_busco.tsv').write_text('BUSCO')
    (source_dir / 'select_rules.tsv').write_text('rule_id\tenabled\tstage\tpriority\tcolumns\tpattern\taction\ttarget_column\toutcome\tscope_column\tscope_mode\tstop_on_match\tnote\n')
    monkeypatch.setattr(dataset_module, 'DATASETS', {
        'toy': {
            'description': 'toy dataset',
            'species': ['Toy_species'],
            'files': {
                'fasta': ['toy.fa.gz'],
                'busco': ['toy_busco.tsv'],
                'rules': ['select_rules.tsv'],
            },
        },
    })
    monkeypatch.setattr(dataset_module, 'get_dataset_dir', lambda: str(dataset_root))
    return source_dir


def _prepare_fake_rule_sets(tmp_path, monkeypatch):
    rules_root = tmp_path / 'select_rule_sets'
    base_dir = rules_root / 'base'
    base_dir.mkdir(parents=True)
    (base_dir / 'select_rules.tsv').write_text(
        'rule_id\tenabled\tstage\tpriority\tcolumns\tpattern\taction\ttarget_column\toutcome\tscope_column\tscope_mode\tstop_on_match\tnote\n'
    )
    plantae_dir = rules_root / 'plantae'
    plantae_dir.mkdir(parents=True)
    (plantae_dir / 'select_rules.tsv').write_text(
        'rule_id\tenabled\tstage\tpriority\tcolumns\tpattern\taction\ttarget_column\toutcome\tscope_column\tscope_mode\tstop_on_match\tnote\nplantae\tyes\tparameter\t0\t\t\t\t\t\t\t\t\t\n'
    )
    invalid_dir = rules_root / 'invalid'
    invalid_dir.mkdir(parents=True)
    (invalid_dir / 'select_rules.tsv').mkdir()
    monkeypatch.setattr(dataset_module, 'get_rule_set_root', lambda: rules_root)
    return rules_root


class TestDatasetExtract:
    def test_extract_dataset_copies_files(self, tmp_path, monkeypatch):
        _prepare_fake_dataset(tmp_path, monkeypatch)
        out_dir = tmp_path / 'out'

        paths = extract_dataset(name='toy', out_dir=str(out_dir), overwrite=False)

        assert os.path.exists(os.path.join(paths['fasta'], 'toy.fa.gz'))
        assert os.path.exists(os.path.join(paths['busco'], 'toy_busco.tsv'))
        assert os.path.exists(os.path.join(paths['rules'], 'select_rules.tsv'))

    def test_extract_dataset_skips_existing_when_not_overwrite(self, tmp_path, monkeypatch):
        _prepare_fake_dataset(tmp_path, monkeypatch)
        out_dir = tmp_path / 'out'
        existing_fasta = out_dir / 'fasta' / 'toy.fa.gz'
        existing_fasta.parent.mkdir(parents=True)
        existing_fasta.write_text('OLD')

        extract_dataset(name='toy', out_dir=str(out_dir), overwrite=False)

        assert existing_fasta.read_text() == 'OLD'

    def test_extract_dataset_overwrites_when_requested(self, tmp_path, monkeypatch):
        _prepare_fake_dataset(tmp_path, monkeypatch)
        out_dir = tmp_path / 'out'
        existing_fasta = out_dir / 'fasta' / 'toy.fa.gz'
        existing_fasta.parent.mkdir(parents=True)
        existing_fasta.write_text('OLD')

        extract_dataset(name='toy', out_dir=str(out_dir), overwrite=True)

        assert existing_fasta.read_text() == 'FA'

    def test_extract_dataset_rejects_output_path_that_is_file(self, tmp_path, monkeypatch):
        _prepare_fake_dataset(tmp_path, monkeypatch)
        out_file = tmp_path / 'out_file'
        out_file.write_text('x')

        with pytest.raises(NotADirectoryError, match='not a directory'):
            extract_dataset(name='toy', out_dir=str(out_file), overwrite=False)

    def test_extract_dataset_rejects_destination_entry_that_is_directory(self, tmp_path, monkeypatch):
        _prepare_fake_dataset(tmp_path, monkeypatch)
        out_dir = tmp_path / 'out'
        bad_dst = out_dir / 'fasta' / 'toy.fa.gz'
        bad_dst.mkdir(parents=True)

        with pytest.raises(NotADirectoryError, match='not a file'):
            extract_dataset(name='toy', out_dir=str(out_dir), overwrite=False)

    def test_extract_dataset_raises_when_source_file_missing(self, tmp_path, monkeypatch):
        _prepare_fake_dataset(tmp_path, monkeypatch)
        dataset_module.DATASETS['toy']['files']['rules'].append('missing.tsv')
        out_dir = tmp_path / 'out'

        with pytest.raises(FileNotFoundError, match='Dataset source file\\(s\\) not found'):
            extract_dataset(name='toy', out_dir=str(out_dir), overwrite=False)

    def test_extract_dataset_raises_when_source_entry_is_directory(self, tmp_path, monkeypatch):
        source_dir = _prepare_fake_dataset(tmp_path, monkeypatch)
        (source_dir / 'not_a_file').mkdir()
        dataset_module.DATASETS['toy']['files']['rules'].append('not_a_file')
        out_dir = tmp_path / 'out'

        with pytest.raises(FileNotFoundError, match='Dataset source file\\(s\\) not found'):
            extract_dataset(name='toy', out_dir=str(out_dir), overwrite=False)


class TestDatasetValidation:
    def test_validate_dataset_name_unknown_exits(self):
        with pytest.raises(ValueError):
            validate_dataset_name('unknown_dataset')

    def test_resolve_dataset_source_dir_missing_exits(self, tmp_path, monkeypatch):
        monkeypatch.setattr(dataset_module, 'get_dataset_dir', lambda: str(tmp_path / 'missing_root'))
        with pytest.raises(FileNotFoundError):
            resolve_dataset_source_dir('toy')


class TestRuleSetExport:
    def test_export_rule_set_copies_select_rules(self, tmp_path, monkeypatch):
        _prepare_fake_rule_sets(tmp_path, monkeypatch)

        path_rules = export_rule_set('base', str(tmp_path / 'out'), overwrite=False)

        assert path_rules == str(tmp_path / 'out' / 'select_rules.tsv')
        assert os.path.exists(path_rules)

    def test_export_rule_set_skips_existing_when_not_overwrite(self, tmp_path, monkeypatch):
        _prepare_fake_rule_sets(tmp_path, monkeypatch)
        out_dir = tmp_path / 'out'
        out_dir.mkdir()
        existing_rules = out_dir / 'select_rules.tsv'
        existing_rules.write_text('existing')

        with pytest.raises(AmalgkitExit) as exc:
            export_rule_set('base', str(out_dir), overwrite=False)

        assert exc.value.exit_code == 0
        assert existing_rules.read_text() == 'existing'

    def test_export_rule_set_overwrites_when_requested(self, tmp_path, monkeypatch):
        _prepare_fake_rule_sets(tmp_path, monkeypatch)
        out_dir = tmp_path / 'out'
        out_dir.mkdir()
        existing_rules = out_dir / 'select_rules.tsv'
        existing_rules.write_text('existing')

        path_rules = export_rule_set('plantae', str(out_dir), overwrite=True)

        assert path_rules == str(existing_rules)
        assert 'plantae' in existing_rules.read_text()

    def test_export_rule_set_rejects_output_path_that_is_file(self, tmp_path, monkeypatch):
        _prepare_fake_rule_sets(tmp_path, monkeypatch)
        out_file = tmp_path / 'out_file'
        out_file.write_text('x')

        with pytest.raises(NotADirectoryError, match='not a directory'):
            export_rule_set('base', str(out_file), overwrite=False)

    def test_export_rule_set_rejects_destination_entry_that_is_directory(self, tmp_path, monkeypatch):
        _prepare_fake_rule_sets(tmp_path, monkeypatch)
        out_dir = tmp_path / 'out'
        out_dir.mkdir()
        (out_dir / 'select_rules.tsv').mkdir()

        with pytest.raises(IsADirectoryError, match='not a file'):
            export_rule_set('base', str(out_dir), overwrite=True)

    def test_list_available_rule_sets_ignores_invalid_directories(self, tmp_path, monkeypatch):
        _prepare_fake_rule_sets(tmp_path, monkeypatch)

        available = list_available_rule_sets()

        assert 'base' in available
        assert 'plantae' in available
        assert 'invalid' not in available

    def test_validate_rule_set_name_raises_for_unknown(self, tmp_path, monkeypatch):
        _prepare_fake_rule_sets(tmp_path, monkeypatch)

        with pytest.raises(ValueError, match='Unknown rule set'):
            validate_rule_set_name('does_not_exist')


class TestDatasetCli:
    def test_dataset_cli_accepts_rule_set(self):
        parser = build_main_parser()

        args = parser.parse_args(['dataset', '--rule_set', 'base'])

        assert args.rule_set == 'base'

    def test_config_command_is_removed(self):
        parser = build_main_parser()

        with pytest.raises(SystemExit):
            parser.parse_args(['config'])
