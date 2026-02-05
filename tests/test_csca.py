import os

from amalgkit.csca import get_spp_from_dir, generate_csca_input_symlinks


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
