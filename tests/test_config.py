import os
import pytest

from amalgkit.config import check_directory


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
