import os
import pytest
import pandas
from types import SimpleNamespace

from amalgkit.csca import get_spp_from_dir, generate_csca_input_symlinks, get_sample_group_string, csca_main
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

    def test_replaces_symlinked_csca_input_path(self, tmp_path):
        dir_curate = tmp_path / 'curate'
        sp_tables = dir_curate / 'Homo_sapiens' / 'tables'
        sp_tables.mkdir(parents=True)
        (sp_tables / 'Homo_sapiens_tpm.tsv').write_text('data')
        real_target = tmp_path / 'real_target'
        real_target.mkdir()
        (real_target / 'keep.txt').write_text('keep')
        dir_csca_input = tmp_path / 'csca_input'
        os.symlink(real_target, dir_csca_input)

        generate_csca_input_symlinks(str(dir_csca_input), str(dir_curate), ['Homo_sapiens'])

        assert os.path.isdir(str(dir_csca_input))
        assert not os.path.islink(str(dir_csca_input))
        assert os.path.islink(str(dir_csca_input / 'Homo_sapiens_tpm.tsv'))
        assert os.path.exists(str(real_target / 'keep.txt'))

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

    def test_raises_clear_error_when_species_tables_dir_missing(self, tmp_path):
        dir_curate = tmp_path / 'curate'
        (dir_curate / 'Homo_sapiens').mkdir(parents=True)
        dir_csca_input = str(tmp_path / 'csca_input')

        with pytest.raises(FileNotFoundError, match='Curate tables directory not found for species'):
            generate_csca_input_symlinks(dir_csca_input, str(dir_curate), ['Homo_sapiens'])

    def test_raises_when_species_tables_dir_has_no_tsv_files(self, tmp_path):
        dir_curate = tmp_path / 'curate'
        sp_tables = dir_curate / 'Homo_sapiens' / 'tables'
        sp_tables.mkdir(parents=True)
        (sp_tables / 'notes.txt').write_text('not a tsv')
        dir_csca_input = str(tmp_path / 'csca_input')

        with pytest.raises(FileNotFoundError, match='No TSV table file was found in curate tables directory'):
            generate_csca_input_symlinks(dir_csca_input, str(dir_curate), ['Homo_sapiens'])

    def test_raises_on_duplicate_table_filename_across_species(self, tmp_path):
        dir_curate = tmp_path / 'curate'
        hs_tables = dir_curate / 'Homo_sapiens' / 'tables'
        mm_tables = dir_curate / 'Mus_musculus' / 'tables'
        hs_tables.mkdir(parents=True)
        mm_tables.mkdir(parents=True)
        (hs_tables / 'shared.tsv').write_text('human')
        (mm_tables / 'shared.tsv').write_text('mouse')
        dir_csca_input = str(tmp_path / 'csca_input')

        with pytest.raises(ValueError, match='Duplicate table filename across species'):
            generate_csca_input_symlinks(
                dir_csca_input,
                str(dir_curate),
                ['Homo_sapiens', 'Mus_musculus'],
            )

    def test_raises_when_csca_input_path_is_a_file(self, tmp_path):
        dir_curate = tmp_path / 'curate'
        sp_tables = dir_curate / 'Homo_sapiens' / 'tables'
        sp_tables.mkdir(parents=True)
        (sp_tables / 'Homo_sapiens_tpm.tsv').write_text('data')
        dir_csca_input = tmp_path / 'csca_input'
        dir_csca_input.write_text('not a directory')

        with pytest.raises(NotADirectoryError, match='CSCA input path exists but is not a directory'):
            generate_csca_input_symlinks(str(dir_csca_input), str(dir_curate), ['Homo_sapiens'])

    def test_raises_when_symlink_destination_is_directory(self, tmp_path, monkeypatch):
        dir_curate = tmp_path / 'curate'
        sp_tables = dir_curate / 'Homo_sapiens' / 'tables'
        sp_tables.mkdir(parents=True)
        (sp_tables / 'Homo_sapiens_tpm.tsv').write_text('data')
        dir_csca_input = tmp_path / 'csca_input'
        target_dst = os.path.join(str(dir_csca_input), 'Homo_sapiens_tpm.tsv')
        real_lexists = os.path.lexists
        real_isdir = os.path.isdir

        def fake_lexists(path):
            if os.path.realpath(path) == os.path.realpath(target_dst):
                return True
            return real_lexists(path)

        def fake_isdir(path):
            if os.path.realpath(path) == os.path.realpath(target_dst):
                return True
            return real_isdir(path)

        monkeypatch.setattr('amalgkit.csca.os.path.lexists', fake_lexists)
        monkeypatch.setattr('amalgkit.csca.os.path.isdir', fake_isdir)

        with pytest.raises(IsADirectoryError, match='CSCA input destination exists but is a directory'):
            generate_csca_input_symlinks(str(dir_csca_input), str(dir_curate), ['Homo_sapiens'])


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

    def test_raises_when_metadata_missing_sample_group_column(self, monkeypatch):
        args = SimpleNamespace(sample_group=None)
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['R1'],
            'exclusion': ['no'],
        }))
        monkeypatch.setattr('amalgkit.csca.load_metadata', lambda _args: metadata)
        with pytest.raises(ValueError, match='sample_group'):
            get_sample_group_string(args)


class TestCscaMain:
    @staticmethod
    def _base_args(out_dir):
        return SimpleNamespace(
            out_dir=str(out_dir),
            sample_group='leaf',
            sample_group_color='DEFAULT',
            dir_busco=None,
            orthogroup_table='dummy.tsv',
            batch_effect_alg='no',
            missing_strategy='strict',
        )

    def test_raises_clear_error_when_curate_dir_missing(self, tmp_path, monkeypatch):
        args = self._base_args(tmp_path / 'out')
        monkeypatch.setattr('amalgkit.csca.check_rscript', lambda: None)
        monkeypatch.setattr('amalgkit.csca.check_ortholog_parameter_compatibility', lambda _args: None)
        with pytest.raises(FileNotFoundError, match='Curate output directory not found'):
            csca_main(args)

    def test_blank_orthogroup_table_raises_parameter_error(self, tmp_path, monkeypatch):
        args = self._base_args(tmp_path / 'out')
        args.orthogroup_table = '   '
        monkeypatch.setattr('amalgkit.csca.check_rscript', lambda: None)
        with pytest.raises(ValueError, match='One of --orthogroup_table and --dir_busco should be specified'):
            csca_main(args)

    def test_raises_when_out_dir_is_file(self, tmp_path, monkeypatch):
        out_path = tmp_path / 'out_path'
        out_path.write_text('not a directory')
        args = self._base_args(out_path)
        monkeypatch.setattr('amalgkit.csca.check_rscript', lambda: None)
        monkeypatch.setattr('amalgkit.csca.check_ortholog_parameter_compatibility', lambda _args: None)
        with pytest.raises(NotADirectoryError, match='Output path exists but is not a directory'):
            csca_main(args)

    def test_raises_when_curate_path_is_file(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'out'
        out_dir.mkdir()
        (out_dir / 'curate').write_text('not a directory')
        args = self._base_args(out_dir)
        monkeypatch.setattr('amalgkit.csca.check_rscript', lambda: None)
        monkeypatch.setattr('amalgkit.csca.check_ortholog_parameter_compatibility', lambda _args: None)
        with pytest.raises(NotADirectoryError, match='Curate output path exists but is not a directory'):
            csca_main(args)

    def test_raises_when_csca_path_is_file(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'out'
        curate_dir = out_dir / 'curate' / 'Species_A' / 'tables'
        curate_dir.mkdir(parents=True)
        (curate_dir / 'Species_A_est_counts.tsv').write_text('target_id\tR1\nG1\t1\n')
        (out_dir / 'csca').write_text('not a directory')
        args = self._base_args(out_dir)
        args.orthogroup_table = str(tmp_path / 'orthogroup.tsv')
        (tmp_path / 'orthogroup.tsv').write_text('busco_id\tSpecies_A\nOG1\tgene1\n')
        monkeypatch.setattr('amalgkit.csca.check_rscript', lambda: None)
        monkeypatch.setattr('amalgkit.csca.check_ortholog_parameter_compatibility', lambda _args: None)
        with pytest.raises(NotADirectoryError, match='csca path exists but is not a directory'):
            csca_main(args)

    def test_raises_when_curate_dir_has_no_species_subdirs(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'out'
        (out_dir / 'curate').mkdir(parents=True)
        args = self._base_args(out_dir)
        monkeypatch.setattr('amalgkit.csca.check_rscript', lambda: None)
        monkeypatch.setattr('amalgkit.csca.check_ortholog_parameter_compatibility', lambda _args: None)
        with pytest.raises(ValueError, match='No curated species directories were found'):
            csca_main(args)

    def test_raises_when_orthogroup_table_path_missing(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'out'
        tables_dir = out_dir / 'curate' / 'Species_A' / 'tables'
        tables_dir.mkdir(parents=True)
        (tables_dir / 'Species_A_est_counts.tsv').write_text('target_id\tR1\nG1\t1\n')
        args = self._base_args(out_dir)
        args.orthogroup_table = str(tmp_path / 'missing.tsv')
        monkeypatch.setattr('amalgkit.csca.check_rscript', lambda: None)
        monkeypatch.setattr('amalgkit.csca.check_ortholog_parameter_compatibility', lambda _args: None)

        with pytest.raises(FileNotFoundError, match='Orthogroup table not found'):
            csca_main(args)

    def test_raises_when_orthogroup_table_path_is_directory(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'out'
        tables_dir = out_dir / 'curate' / 'Species_A' / 'tables'
        tables_dir.mkdir(parents=True)
        (tables_dir / 'Species_A_est_counts.tsv').write_text('target_id\tR1\nG1\t1\n')
        og_dir = tmp_path / 'orthogroup_dir'
        og_dir.mkdir()
        args = self._base_args(out_dir)
        args.orthogroup_table = str(og_dir)
        monkeypatch.setattr('amalgkit.csca.check_rscript', lambda: None)
        monkeypatch.setattr('amalgkit.csca.check_ortholog_parameter_compatibility', lambda _args: None)

        with pytest.raises(IsADirectoryError, match='Orthogroup table path exists but is not a file'):
            csca_main(args)
