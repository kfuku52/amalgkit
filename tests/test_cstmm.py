import os
import pytest

from amalgkit.cstmm import filepath2spp, get_count_files


class TestFilepath2spp:
    def test_basic_extraction(self):
        paths = [
            '/path/to/merge/Homo_sapiens/Homo_sapiens_est_counts.tsv',
            '/path/to/merge/Mus_musculus/Mus_musculus_est_counts.tsv',
        ]
        result = filepath2spp(paths)
        assert result == ['Homo_sapiens', 'Mus_musculus']

    def test_single_species(self):
        paths = ['/merge/Arabidopsis_thaliana/Arabidopsis_thaliana_est_counts.tsv']
        result = filepath2spp(paths)
        assert result == ['Arabidopsis_thaliana']

    def test_basename_only(self):
        paths = ['Drosophila_melanogaster_est_counts.tsv']
        result = filepath2spp(paths)
        assert result == ['Drosophila_melanogaster']


# ---------------------------------------------------------------------------
# get_count_files (wiki: finds est_counts.tsv in species subdirectories)
# ---------------------------------------------------------------------------

class TestGetCountFiles:
    def test_finds_count_files(self, tmp_path):
        """Finds est_counts.tsv files in species subdirectories."""
        sp_dir = tmp_path / 'Homo_sapiens'
        sp_dir.mkdir()
        count_file = sp_dir / 'Homo_sapiens_est_counts.tsv'
        count_file.write_text('target_id\tSRR001\n')
        result = get_count_files(str(tmp_path))
        assert len(result) == 1
        assert result[0].endswith('Homo_sapiens_est_counts.tsv')

    def test_multiple_species(self, tmp_path):
        """Finds count files across multiple species directories."""
        for sp in ['Homo_sapiens', 'Mus_musculus']:
            sp_dir = tmp_path / sp
            sp_dir.mkdir()
            (sp_dir / f'{sp}_est_counts.tsv').write_text('data')
        result = get_count_files(str(tmp_path))
        assert len(result) == 2

    def test_no_count_files_raises(self, tmp_path):
        """No est_counts.tsv files found should raise Exception."""
        sp_dir = tmp_path / 'Homo_sapiens'
        sp_dir.mkdir()
        (sp_dir / 'other_file.tsv').write_text('data')
        with pytest.raises(Exception, match='No est_counts.tsv'):
            get_count_files(str(tmp_path))

    def test_multiple_count_files_per_species_raises(self, tmp_path):
        """Multiple est_counts.tsv in one species dir should raise Exception."""
        sp_dir = tmp_path / 'Homo_sapiens'
        sp_dir.mkdir()
        (sp_dir / 'Homo_sapiens_est_counts.tsv').write_text('data')
        (sp_dir / 'Homo_sapiens_cstmm_est_counts.tsv').write_text('data')
        with pytest.raises(Exception, match='Multiple est_counts.tsv'):
            get_count_files(str(tmp_path))

    def test_ignores_non_directory_entries(self, tmp_path):
        """Non-directory entries in the count dir should be skipped."""
        sp_dir = tmp_path / 'Homo_sapiens'
        sp_dir.mkdir()
        (sp_dir / 'Homo_sapiens_est_counts.tsv').write_text('data')
        # A regular file at the top level (not a directory)
        (tmp_path / 'metadata.tsv').write_text('data')
        result = get_count_files(str(tmp_path))
        assert len(result) == 1

    def test_empty_directory_raises(self, tmp_path):
        """Empty directory should raise since no count files found."""
        with pytest.raises(Exception, match='No est_counts.tsv'):
            get_count_files(str(tmp_path))
