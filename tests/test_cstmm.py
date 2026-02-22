import os
import pytest
from types import SimpleNamespace

from amalgkit.cstmm import filepath2spp, get_count_files, cstmm_main


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
        """No est_counts.tsv files found should raise FileNotFoundError."""
        sp_dir = tmp_path / 'Homo_sapiens'
        sp_dir.mkdir()
        (sp_dir / 'other_file.tsv').write_text('data')
        with pytest.raises(FileNotFoundError, match='No est_counts.tsv'):
            get_count_files(str(tmp_path))

    def test_multiple_count_files_per_species_raises(self, tmp_path):
        """Multiple est_counts.tsv in one species dir should raise ValueError."""
        sp_dir = tmp_path / 'Homo_sapiens'
        sp_dir.mkdir()
        (sp_dir / 'Homo_sapiens_est_counts.tsv').write_text('data')
        (sp_dir / 'Homo_sapiens_cstmm_est_counts.tsv').write_text('data')
        with pytest.raises(ValueError, match='Multiple est_counts.tsv'):
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

    def test_ignores_hidden_and_tmp_directories(self, tmp_path):
        sp_dir = tmp_path / 'Homo_sapiens'
        sp_dir.mkdir()
        (sp_dir / 'Homo_sapiens_est_counts.tsv').write_text('data')
        hidden_dir = tmp_path / '.cache'
        hidden_dir.mkdir()
        (hidden_dir / 'Hidden_est_counts.tsv').write_text('data')
        tmp_dir = tmp_path / 'tmp.partial'
        tmp_dir.mkdir()
        (tmp_dir / 'Tmp_est_counts.tsv').write_text('data')

        result = get_count_files(str(tmp_path))

        assert len(result) == 1
        assert result[0].endswith('Homo_sapiens_est_counts.tsv')

    def test_empty_directory_raises(self, tmp_path):
        """Empty directory should raise since no count files found."""
        with pytest.raises(FileNotFoundError, match='No est_counts.tsv'):
            get_count_files(str(tmp_path))

    def test_raises_when_any_species_directory_lacks_count_file(self, tmp_path):
        sp_a = tmp_path / 'Species_A'
        sp_b = tmp_path / 'Species_B'
        sp_a.mkdir()
        sp_b.mkdir()
        (sp_a / 'Species_A_est_counts.tsv').write_text('target_id\tR1\nG1\t1\n')
        (sp_b / 'notes.txt').write_text('missing count')

        with pytest.raises(FileNotFoundError, match='No est_counts.tsv file found for species directory\\(ies\\): Species_B'):
            get_count_files(str(tmp_path))

    def test_missing_directory_raises_clear_error(self, tmp_path):
        missing = tmp_path / 'missing'
        with pytest.raises(FileNotFoundError, match='Count directory not found'):
            get_count_files(str(missing))

    def test_file_path_raises_not_a_directory(self, tmp_path):
        path_file = tmp_path / 'counts_path'
        path_file.write_text('x')
        with pytest.raises(NotADirectoryError, match='not a directory'):
            get_count_files(str(path_file))


class TestCstmmMain:
    def test_rejects_out_dir_file_path(self, tmp_path):
        out_path = tmp_path / 'out_path'
        out_path.write_text('not a directory')
        args = SimpleNamespace(
            out_dir=str(out_path),
            dir_count='inferred',
            dir_busco=None,
            orthogroup_table=None,
        )
        with pytest.raises(NotADirectoryError, match='Output path exists but is not a directory'):
            cstmm_main(args)

    def test_rejects_cstmm_path_file(self, tmp_path):
        out_dir = tmp_path / 'out'
        out_dir.mkdir()
        (out_dir / 'cstmm').write_text('not a directory')
        args = SimpleNamespace(
            out_dir=str(out_dir),
            dir_count='inferred',
            dir_busco=None,
            orthogroup_table=None,
        )
        with pytest.raises(NotADirectoryError, match='cstmm path exists but is not a directory'):
            cstmm_main(args)

    def test_single_species_allows_missing_ortholog_inputs(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'out'
        merge_dir = out_dir / 'merge'
        species_dir = merge_dir / 'Species_A'
        species_dir.mkdir(parents=True)
        (species_dir / 'Species_A_est_counts.tsv').write_text('target_id\tR1\nG1\t1\n')
        args = SimpleNamespace(
            out_dir=str(out_dir),
            dir_count='inferred',
            dir_busco=None,
            orthogroup_table=None,
        )
        captured = {}

        monkeypatch.setattr('amalgkit.cstmm.check_rscript', lambda: None)
        monkeypatch.setattr('amalgkit.cstmm.cleanup_tmp_amalgkit_files', lambda work_dir: None)
        monkeypatch.setattr('amalgkit.cstmm.subprocess.check_call', lambda cmd: captured.setdefault('cmd', cmd))
        monkeypatch.setattr(
            'amalgkit.cstmm.orthogroup2genecount',
            lambda **_kwargs: (_ for _ in ()).throw(AssertionError('orthogroup2genecount should not run for single species')),
        )

        cstmm_main(args)

        assert captured['cmd'][6] == 'single_species'
        assert captured['cmd'][3] == ''
        assert captured['cmd'][4] == ''

    def test_multi_species_requires_ortholog_inputs(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'out'
        merge_dir = out_dir / 'merge'
        for species in ['Species_A', 'Species_B']:
            species_dir = merge_dir / species
            species_dir.mkdir(parents=True)
            (species_dir / f'{species}_est_counts.tsv').write_text('target_id\tR1\nG1\t1\n')
        args = SimpleNamespace(
            out_dir=str(out_dir),
            dir_count='inferred',
            dir_busco=None,
            orthogroup_table=None,
        )

        monkeypatch.setattr('amalgkit.cstmm.check_rscript', lambda: None)
        with pytest.raises(ValueError, match='One of --orthogroup_table and --dir_busco should be specified'):
            cstmm_main(args)

    def test_multi_species_rejects_missing_orthogroup_table_path(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'out'
        merge_dir = out_dir / 'merge'
        for species in ['Species_A', 'Species_B']:
            species_dir = merge_dir / species
            species_dir.mkdir(parents=True)
            (species_dir / f'{species}_est_counts.tsv').write_text('target_id\tR1\nG1\t1\n')
        args = SimpleNamespace(
            out_dir=str(out_dir),
            dir_count='inferred',
            dir_busco=None,
            orthogroup_table=str(tmp_path / 'missing.tsv'),
        )

        monkeypatch.setattr('amalgkit.cstmm.check_rscript', lambda: None)
        with pytest.raises(FileNotFoundError, match='Orthogroup table not found'):
            cstmm_main(args)

    def test_multi_species_rejects_directory_orthogroup_table_path(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'out'
        merge_dir = out_dir / 'merge'
        for species in ['Species_A', 'Species_B']:
            species_dir = merge_dir / species
            species_dir.mkdir(parents=True)
            (species_dir / f'{species}_est_counts.tsv').write_text('target_id\tR1\nG1\t1\n')
        og_dir = tmp_path / 'orthogroup_dir'
        og_dir.mkdir()
        args = SimpleNamespace(
            out_dir=str(out_dir),
            dir_count='inferred',
            dir_busco=None,
            orthogroup_table=str(og_dir),
        )

        monkeypatch.setattr('amalgkit.cstmm.check_rscript', lambda: None)
        with pytest.raises(IsADirectoryError, match='Orthogroup table path exists but is not a file'):
            cstmm_main(args)

    def test_rejects_missing_count_file_in_any_species_directory(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'out'
        merge_dir = out_dir / 'merge'
        species_a = merge_dir / 'Species_A'
        species_b = merge_dir / 'Species_B'
        species_a.mkdir(parents=True)
        species_b.mkdir(parents=True)
        (species_a / 'Species_A_est_counts.tsv').write_text('target_id\tR1\nG1\t1\n')
        (species_b / 'notes.txt').write_text('missing count')
        args = SimpleNamespace(
            out_dir=str(out_dir),
            dir_count='inferred',
            dir_busco=None,
            orthogroup_table=None,
        )

        monkeypatch.setattr('amalgkit.cstmm.check_rscript', lambda: None)
        with pytest.raises(FileNotFoundError, match='No est_counts.tsv file found for species directory\\(ies\\): Species_B'):
            cstmm_main(args)
