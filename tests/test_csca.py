import os
import pytest
import pandas
from types import SimpleNamespace

from amalgkit.csca import get_spp_from_dir, generate_csca_input_symlinks, get_sample_group_string
from amalgkit.util import Metadata


# ---------------------------------------------------------------------------
# get_spp_from_dir (filters hidden and tmp files from directory listing)
# ---------------------------------------------------------------------------

class TestGetSppFromDir:
    def test_basic_listing(self, tmp_path):
        (tmp_path / 'Homo_sapiens').mkdir()
        (tmp_path / 'Mus_musculus').mkdir()
        result = get_spp_from_dir(str(tmp_path))
        assert sorted(result) == ['Homo_sapiens', 'Mus_musculus']

    def test_excludes_hidden_files(self, tmp_path):
        (tmp_path / 'Homo_sapiens').mkdir()
        (tmp_path / '.hidden').mkdir()
        (tmp_path / '.DS_Store').touch()
        result = get_spp_from_dir(str(tmp_path))
        assert result == ['Homo_sapiens']

    def test_excludes_tmp_files(self, tmp_path):
        (tmp_path / 'Homo_sapiens').mkdir()
        (tmp_path / 'tmp.amalgkit.1234').touch()
        result = get_spp_from_dir(str(tmp_path))
        assert result == ['Homo_sapiens']

    def test_empty_directory(self, tmp_path):
        result = get_spp_from_dir(str(tmp_path))
        assert result == []

    def test_mixed_entries(self, tmp_path):
        """Only non-hidden, non-tmp entries are returned."""
        (tmp_path / 'Species_A').mkdir()
        (tmp_path / 'Species_B').mkdir()
        (tmp_path / '.gitkeep').touch()
        (tmp_path / 'tmp.output').touch()
        result = get_spp_from_dir(str(tmp_path))
        assert sorted(result) == ['Species_A', 'Species_B']

    def test_ignores_non_directory_entries(self, tmp_path):
        (tmp_path / 'Species_A').mkdir()
        (tmp_path / 'metadata.tsv').write_text('data')
        result = get_spp_from_dir(str(tmp_path))
        assert result == ['Species_A']


# ---------------------------------------------------------------------------
# generate_csca_input_symlinks (creates symlinks from curate to csca)
# ---------------------------------------------------------------------------

class TestGenerateCscaInputSymlinks:
    def test_creates_symlinks(self, tmp_path):
        """Creates symlinks to curate output tables."""
        dir_curate = tmp_path / 'curate'
        sp_tables = dir_curate / 'Homo_sapiens' / 'tables'
        sp_tables.mkdir(parents=True)
        (sp_tables / 'Homo_sapiens_tpm.tsv').write_text('data')
        (sp_tables / 'Homo_sapiens_est_counts.tsv').write_text('data')
        dir_csca_input = str(tmp_path / 'csca_input')
        generate_csca_input_symlinks(dir_csca_input, str(dir_curate), ['Homo_sapiens'])
        assert os.path.isdir(dir_csca_input)
        assert os.path.islink(os.path.join(dir_csca_input, 'Homo_sapiens_tpm.tsv'))
        assert os.path.islink(os.path.join(dir_csca_input, 'Homo_sapiens_est_counts.tsv'))

    def test_multiple_species(self, tmp_path):
        """Creates symlinks for multiple species."""
        dir_curate = tmp_path / 'curate'
        for sp in ['Homo_sapiens', 'Mus_musculus']:
            tables = dir_curate / sp / 'tables'
            tables.mkdir(parents=True)
            (tables / f'{sp}_tpm.tsv').write_text('data')
        dir_csca_input = str(tmp_path / 'csca_input')
        generate_csca_input_symlinks(dir_csca_input, str(dir_curate), ['Homo_sapiens', 'Mus_musculus'])
        assert os.path.islink(os.path.join(dir_csca_input, 'Homo_sapiens_tpm.tsv'))
        assert os.path.islink(os.path.join(dir_csca_input, 'Mus_musculus_tpm.tsv'))

    def test_recreates_directory(self, tmp_path):
        """Removes and recreates csca input directory."""
        dir_curate = tmp_path / 'curate'
        sp_tables = dir_curate / 'Homo_sapiens' / 'tables'
        sp_tables.mkdir(parents=True)
        (sp_tables / 'Homo_sapiens_tpm.tsv').write_text('data')
        dir_csca_input = tmp_path / 'csca_input'
        dir_csca_input.mkdir()
        (dir_csca_input / 'old_file.txt').write_text('old')
        generate_csca_input_symlinks(str(dir_csca_input), str(dir_curate), ['Homo_sapiens'])
        # Old file should be gone
        assert not os.path.exists(str(dir_csca_input / 'old_file.txt'))
        assert os.path.islink(str(dir_csca_input / 'Homo_sapiens_tpm.tsv'))

    def test_ignores_non_tsv_and_subdirectories(self, tmp_path):
        dir_curate = tmp_path / 'curate'
        sp_tables = dir_curate / 'Homo_sapiens' / 'tables'
        sp_tables.mkdir(parents=True)
        (sp_tables / 'Homo_sapiens_tpm.tsv').write_text('data')
        (sp_tables / 'notes.txt').write_text('ignore')
        (sp_tables / 'nested').mkdir()
        dir_csca_input = str(tmp_path / 'csca_input')

        generate_csca_input_symlinks(dir_csca_input, str(dir_curate), ['Homo_sapiens'])

        assert os.path.islink(os.path.join(dir_csca_input, 'Homo_sapiens_tpm.tsv'))
        assert not os.path.exists(os.path.join(dir_csca_input, 'notes.txt'))
        assert not os.path.exists(os.path.join(dir_csca_input, 'nested'))


class TestGetSampleGroupString:
    def test_uses_cli_sample_group(self):
        args = SimpleNamespace(sample_group='leaf,root')
        assert get_sample_group_string(args) == 'leaf|root'

    def test_cli_sample_group_supports_pipe_and_hyphen(self):
        args = SimpleNamespace(sample_group='non-treated | treated')
        assert get_sample_group_string(args) == 'non-treated|treated'

    def test_reads_sample_group_from_metadata(self, monkeypatch):
        args = SimpleNamespace(sample_group=None)
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['R1', 'R2'],
            'sample_group': ['treated', 'control'],
            'exclusion': ['no', 'no'],
        }))
        monkeypatch.setattr('amalgkit.csca.load_metadata', lambda _args: metadata)
        out = get_sample_group_string(args)
        assert set(out.split('|')) == {'treated', 'control'}

    def test_metadata_sample_groups_drop_blank_and_deduplicate(self, monkeypatch):
        args = SimpleNamespace(sample_group=None)
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['R1', 'R2', 'R3', 'R4'],
            'sample_group': [' treated ', '', 'treated', 'control'],
            'exclusion': ['no', 'no', 'no', 'no'],
        }))
        monkeypatch.setattr('amalgkit.csca.load_metadata', lambda _args: metadata)
        assert get_sample_group_string(args) == 'treated|control'

    def test_metadata_sample_groups_drop_nan(self, monkeypatch):
        args = SimpleNamespace(sample_group=None)
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['R1', 'R2', 'R3'],
            'sample_group': ['treated', float('nan'), 'control'],
            'exclusion': ['no', 'no', 'no'],
        }))
        monkeypatch.setattr('amalgkit.csca.load_metadata', lambda _args: metadata)
        assert get_sample_group_string(args) == 'treated|control'

    def test_raises_when_no_sample_group_selected(self):
        args = SimpleNamespace(sample_group='')
        with pytest.raises(ValueError, match='No sample_group was selected'):
            get_sample_group_string(args)
