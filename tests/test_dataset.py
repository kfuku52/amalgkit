import os
from types import SimpleNamespace

import pytest

from amalgkit import dataset as dataset_module
from amalgkit.dataset import (
    export_rule_set,
    extract_dataset,
    initialize_workspace,
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
        assert paths['select_rules_tsv'] == str(out_dir / 'select_rules.tsv')
        assert os.path.exists(paths['select_rules_tsv'])
        assert not os.path.exists(out_dir / 'rules' / 'select_rules.tsv')
        assert os.path.isdir(paths['private_fastq'])
        assert os.path.isdir(paths['downloads'])
        assert os.path.isdir(paths['metadata'])
        assert os.path.isdir(paths['metadata_specieswise'])

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

    def test_validate_dataset_name_accepts_init(self):
        validate_dataset_name('init')

    def test_resolve_dataset_source_dir_missing_exits(self, tmp_path, monkeypatch):
        monkeypatch.setattr(dataset_module, 'get_dataset_dir', lambda: str(tmp_path / 'missing_root'))
        with pytest.raises(FileNotFoundError):
            resolve_dataset_source_dir('toy')


class TestWorkspaceInit:
    def test_initialize_workspace_creates_templates_and_base_rules(self, tmp_path, monkeypatch):
        _prepare_fake_rule_sets(tmp_path, monkeypatch)
        out_dir = tmp_path / 'out'

        paths = initialize_workspace(str(out_dir), overwrite=False)

        assert os.path.isdir(paths['workspace_dirs']['fasta'])
        assert os.path.isdir(paths['workspace_dirs']['private_fastq'])
        assert os.path.isdir(paths['workspace_dirs']['downloads'])
        assert os.path.isdir(paths['workspace_dirs']['download_locks'])
        assert os.path.isdir(paths['workspace_dirs']['metadata'])
        assert os.path.isdir(paths['workspace_dirs']['metadata_specieswise'])
        assert os.path.exists(paths['workspace_readme'])
        assert os.path.exists(paths['species_tsv'])
        assert os.path.exists(paths['organ_terms_tsv'])
        assert os.path.exists(paths['select_rules_tsv'])
        assert 'species.tsv' in (out_dir / 'WORKSPACE_README.md').read_text()
        assert (out_dir / 'species.tsv').read_text() == 'scientific_name\tspecies_token\n'
        assert (out_dir / 'organ_terms.tsv').read_text() == 'sample_group\ttitle_terms\n'
        assert 'rule_id' in (out_dir / 'select_rules.tsv').read_text()

    def test_initialize_workspace_keeps_existing_templates_when_not_overwrite(self, tmp_path, monkeypatch):
        _prepare_fake_rule_sets(tmp_path, monkeypatch)
        out_dir = tmp_path / 'out'
        out_dir.mkdir()
        workspace_readme = out_dir / 'WORKSPACE_README.md'
        workspace_readme.write_text('keep me\n')
        species_tsv = out_dir / 'species.tsv'
        species_tsv.write_text('existing\n')

        paths = initialize_workspace(str(out_dir), overwrite=False)

        assert paths['workspace_readme'] == str(workspace_readme)
        assert workspace_readme.read_text() == 'keep me\n'
        assert paths['species_tsv'] == str(species_tsv)
        assert species_tsv.read_text() == 'existing\n'


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

    def test_dataset_cli_accepts_init_name(self):
        parser = build_main_parser()

        args = parser.parse_args(['dataset', '--name', 'init'])

        assert args.name == 'init'

    def test_config_command_is_removed(self):
        parser = build_main_parser()

        with pytest.raises(SystemExit):
            parser.parse_args(['config'])


class TestDatasetMain:
    def test_explicit_rule_set_overrides_dataset_rules(self, tmp_path, monkeypatch):
        _prepare_fake_dataset(tmp_path, monkeypatch)
        _prepare_fake_rule_sets(tmp_path, monkeypatch)
        out_dir = tmp_path / 'out'

        dataset_module.dataset_main(SimpleNamespace(
            list=False,
            name='toy',
            rule_set='plantae',
            out_dir=str(out_dir),
            overwrite=False,
        ))

        assert (out_dir / 'fasta' / 'toy.fa.gz').exists()
        assert (out_dir / 'busco' / 'toy_busco.tsv').exists()
        rules_text = (out_dir / 'select_rules.tsv').read_text()
        assert 'plantae' in rules_text

    def test_invalid_rule_set_does_not_create_init_workspace_side_effects(self, tmp_path, monkeypatch):
        _prepare_fake_rule_sets(tmp_path, monkeypatch)
        out_dir = tmp_path / 'out'

        with pytest.raises(ValueError, match='Unknown rule set'):
            dataset_module.dataset_main(SimpleNamespace(
                list=False,
                name='init',
                rule_set='does_not_exist',
                out_dir=str(out_dir),
                overwrite=False,
            ))

        assert not out_dir.exists()

    def test_existing_rule_file_blocks_init_before_templates_are_written(self, tmp_path, monkeypatch):
        _prepare_fake_rule_sets(tmp_path, monkeypatch)
        out_dir = tmp_path / 'out'
        out_dir.mkdir()
        (out_dir / 'select_rules.tsv').write_text('keep-rules\n')

        with pytest.raises(AmalgkitExit) as exc:
            dataset_module.dataset_main(SimpleNamespace(
                list=False,
                name='init',
                rule_set='plantae',
                out_dir=str(out_dir),
                overwrite=False,
            ))

        assert exc.value.exit_code == 0
        assert not (out_dir / 'WORKSPACE_README.md').exists()
        assert not (out_dir / 'species.tsv').exists()
        assert not (out_dir / 'organ_terms.tsv').exists()
        assert (out_dir / 'select_rules.tsv').read_text() == 'keep-rules\n'
