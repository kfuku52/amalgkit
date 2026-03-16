import os
import pytest

from amalgkit.config import check_directory, config_main, list_available_config_sets, validate_config_set
from amalgkit.exceptions import AmalgkitExit


# ---------------------------------------------------------------------------
# check_directory (creates config directory or exits if exists)
# ---------------------------------------------------------------------------

class TestCheckDirectory:
    def test_creates_new_directory(self, tmp_path):
        """Creates select_rules.tsv path when it does not exist."""
        class Args:
            out_dir = str(tmp_path)
            config = 'default'
            overwrite = False
        path_config = check_directory(Args())
        assert path_config == str(tmp_path / 'select_rules.tsv')

    def test_existing_directory_no_overwrite_exits(self, tmp_path):
        """Exits when select_rules.tsv exists and overwrite is False."""
        config_file = tmp_path / 'select_rules.tsv'
        config_file.write_text('existing')
        class Args:
            out_dir = str(tmp_path)
            config = 'default'
            overwrite = False
        with pytest.raises(AmalgkitExit) as exc:
            check_directory(Args())
        assert exc.value.exit_code == 0

    def test_existing_directory_with_overwrite(self, tmp_path):
        """Does not exit when select_rules.tsv exists and overwrite is True."""
        config_file = tmp_path / 'select_rules.tsv'
        config_file.write_text('existing')
        class Args:
            out_dir = str(tmp_path)
            config = 'default'
            overwrite = True
        # Should not raise
        check_directory(Args())

    def test_existing_file_path_raises_not_a_directory(self, tmp_path):
        """If output select rules path exists as a directory, raise clear error."""
        config_path = tmp_path / 'select_rules.tsv'
        config_path.mkdir()

        class Args:
            out_dir = str(tmp_path)
            config = 'default'
            overwrite = True

        with pytest.raises(IsADirectoryError, match='not a file'):
            check_directory(Args())

    def test_out_dir_file_path_raises_not_a_directory(self, tmp_path):
        out_file = tmp_path / 'out_file'
        out_file.write_text('x')

        class Args:
            out_dir = str(out_file)
            config = 'default'
            overwrite = True

        with pytest.raises(NotADirectoryError, match='Output path exists but is not a directory'):
            check_directory(Args())


class TestConfigMain:
    def test_unknown_config_set_raises_clear_error(self, tmp_path):
        class Args:
            out_dir = str(tmp_path)
            config = 'does_not_exist'
            overwrite = True

        with pytest.raises(ValueError, match='Unknown config set'):
            config_main(Args())
        assert not (tmp_path / 'config_does_not_exist').exists()

    def test_lists_available_config_sets(self):
        available = list_available_config_sets()
        assert 'base' in available

    def test_list_available_config_sets_ignores_config_named_directories(self, tmp_path, monkeypatch):
        fake_root = tmp_path / 'config_dir'
        fake_root.mkdir()
        valid_set = fake_root / 'valid'
        valid_set.mkdir()
        (valid_set / 'select_rules.tsv').write_text('rule_id\tenabled\tstage\tpriority\tcolumns\tpattern\taction\ttarget_column\toutcome\tscope_column\tscope_mode\tstop_on_match\tnote\n')
        invalid_set = fake_root / 'invalid'
        invalid_set.mkdir()
        (invalid_set / 'select_rules.tsv').mkdir()

        monkeypatch.setattr('amalgkit.config.ir.files', lambda _pkg: fake_root)
        available = list_available_config_sets()

        assert 'valid' in available
        assert 'invalid' not in available

    def test_validate_config_set_raises_for_unknown(self):
        with pytest.raises(ValueError, match='Unknown config set'):
            validate_config_set('does_not_exist')

    def test_raises_when_destination_config_entry_is_directory(self, tmp_path):
        config_path = tmp_path / 'select_rules.tsv'
        config_path.mkdir()

        class Args:
            out_dir = str(tmp_path)
            config = 'base'
            overwrite = True

        with pytest.raises(IsADirectoryError, match='Config output path exists but is not a file'):
            config_main(Args())
