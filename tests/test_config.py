import os
import pytest

from amalgkit.config import check_directory, config_main, list_available_config_sets, validate_config_set


# ---------------------------------------------------------------------------
# check_directory (creates config directory or exits if exists)
# ---------------------------------------------------------------------------

class TestCheckDirectory:
    def test_creates_new_directory(self, tmp_path):
        """Creates config directory when it does not exist."""
        class Args:
            out_dir = str(tmp_path)
            config = 'default'
            overwrite = False
        check_directory(Args())
        assert os.path.isdir(str(tmp_path / 'config_default'))

    def test_existing_directory_no_overwrite_exits(self, tmp_path):
        """Exits when directory exists and overwrite is False."""
        config_dir = tmp_path / 'config_default'
        config_dir.mkdir()
        class Args:
            out_dir = str(tmp_path)
            config = 'default'
            overwrite = False
        with pytest.raises(SystemExit):
            check_directory(Args())

    def test_existing_directory_with_overwrite(self, tmp_path):
        """Does not exit when directory exists and overwrite is True."""
        config_dir = tmp_path / 'config_default'
        config_dir.mkdir()
        class Args:
            out_dir = str(tmp_path)
            config = 'default'
            overwrite = True
        # Should not raise
        check_directory(Args())

    def test_existing_file_path_raises_not_a_directory(self, tmp_path):
        """If config path exists as a file, raise clear error."""
        config_path = tmp_path / 'config_default'
        config_path.write_text('x')

        class Args:
            out_dir = str(tmp_path)
            config = 'default'
            overwrite = True

        with pytest.raises(NotADirectoryError, match='not a directory'):
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
        (valid_set / 'group_attribute.config').write_text('tissue\tsource_name\n')
        invalid_set = fake_root / 'invalid'
        invalid_set.mkdir()
        (invalid_set / 'group_attribute.config').mkdir()

        monkeypatch.setattr('amalgkit.config.ir.files', lambda _pkg: fake_root)
        available = list_available_config_sets()

        assert 'valid' in available
        assert 'invalid' not in available

    def test_validate_config_set_raises_for_unknown(self):
        with pytest.raises(ValueError, match='Unknown config set'):
            validate_config_set('does_not_exist')

    def test_raises_when_destination_config_entry_is_directory(self, tmp_path):
        config_dir = tmp_path / 'config_base'
        config_dir.mkdir()
        (config_dir / 'group_attribute.config').mkdir()

        class Args:
            out_dir = str(tmp_path)
            config = 'base'
            overwrite = True

        with pytest.raises(IsADirectoryError, match='Config output path exists but is not a file'):
            config_main(Args())
