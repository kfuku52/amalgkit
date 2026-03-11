import os
import pytest
import pandas
from types import SimpleNamespace

from amalgkit.cross_species_filter import get_species_from_dir, generate_input_symlinks, get_sample_group_string, run_cross_species_filter
from amalgkit.util import Metadata


# ---------------------------------------------------------------------------
# get_species_from_dir (filters hidden and tmp files from directory listing)
# ---------------------------------------------------------------------------

class TestGetSppFromDir:
    def test_basic_listing(self, tmp_path):
        (tmp_path / 'Homo_sapiens').mkdir()
        (tmp_path / 'Mus_musculus').mkdir()
        result = get_species_from_dir(str(tmp_path))
        assert sorted(result) == ['Homo_sapiens', 'Mus_musculus']

    def test_excludes_hidden_files(self, tmp_path):
        (tmp_path / 'Homo_sapiens').mkdir()
        (tmp_path / '.hidden').mkdir()
        (tmp_path / '.DS_Store').touch()
        result = get_species_from_dir(str(tmp_path))
        assert result == ['Homo_sapiens']

    def test_excludes_tmp_files(self, tmp_path):
        (tmp_path / 'Homo_sapiens').mkdir()
        (tmp_path / 'tmp.amalgkit.1234').touch()
        result = get_species_from_dir(str(tmp_path))
        assert result == ['Homo_sapiens']

    def test_empty_directory(self, tmp_path):
        result = get_species_from_dir(str(tmp_path))
        assert result == []

    def test_mixed_entries(self, tmp_path):
        """Only non-hidden, non-tmp entries are returned."""
        (tmp_path / 'Species_A').mkdir()
        (tmp_path / 'Species_B').mkdir()
        (tmp_path / '.gitkeep').touch()
        (tmp_path / 'tmp.output').touch()
        result = get_species_from_dir(str(tmp_path))
        assert sorted(result) == ['Species_A', 'Species_B']

    def test_ignores_non_directory_entries(self, tmp_path):
        (tmp_path / 'Species_A').mkdir()
        (tmp_path / 'metadata.tsv').write_text('data')
        result = get_species_from_dir(str(tmp_path))
        assert result == ['Species_A']


# ---------------------------------------------------------------------------
# generate_input_symlinks (creates symlinks from curate to csca)
# ---------------------------------------------------------------------------

class TestGenerateCscaInputSymlinks:
    def test_creates_symlinks(self, tmp_path):
        """Creates symlinks to curate output tables."""
        dir_curate = tmp_path / 'per_species'
        sp_tables = dir_curate / 'Homo_sapiens' / 'tables'
        sp_tables.mkdir(parents=True)
        (sp_tables / 'Homo_sapiens_tpm.tsv').write_text('data')
        (sp_tables / 'Homo_sapiens_est_counts.tsv').write_text('data')
        dir_csca_input = str(tmp_path / 'csca_input')
        generate_input_symlinks(dir_csca_input, str(dir_curate), ['Homo_sapiens'])
        assert os.path.isdir(dir_csca_input)
        assert os.path.islink(os.path.join(dir_csca_input, 'Homo_sapiens_tpm.tsv'))
        assert os.path.islink(os.path.join(dir_csca_input, 'Homo_sapiens_est_counts.tsv'))

    def test_multiple_species(self, tmp_path):
        """Creates symlinks for multiple species."""
        dir_curate = tmp_path / 'per_species'
        for sp in ['Homo_sapiens', 'Mus_musculus']:
            tables = dir_curate / sp / 'tables'
            tables.mkdir(parents=True)
            (tables / f'{sp}_tpm.tsv').write_text('data')
        dir_csca_input = str(tmp_path / 'csca_input')
        generate_input_symlinks(dir_csca_input, str(dir_curate), ['Homo_sapiens', 'Mus_musculus'])
        assert os.path.islink(os.path.join(dir_csca_input, 'Homo_sapiens_tpm.tsv'))
        assert os.path.islink(os.path.join(dir_csca_input, 'Mus_musculus_tpm.tsv'))

    def test_recreates_directory(self, tmp_path):
        """Removes and recreates csca input directory."""
        dir_curate = tmp_path / 'per_species'
        sp_tables = dir_curate / 'Homo_sapiens' / 'tables'
        sp_tables.mkdir(parents=True)
        (sp_tables / 'Homo_sapiens_tpm.tsv').write_text('data')
        dir_csca_input = tmp_path / 'csca_input'
        dir_csca_input.mkdir()
        (dir_csca_input / 'old_file.txt').write_text('old')
        generate_input_symlinks(str(dir_csca_input), str(dir_curate), ['Homo_sapiens'])
        # Old file should be gone
        assert not os.path.exists(str(dir_csca_input / 'old_file.txt'))
        assert os.path.islink(str(dir_csca_input / 'Homo_sapiens_tpm.tsv'))

    def test_replaces_symlinked_csca_input_path(self, tmp_path):
        dir_curate = tmp_path / 'per_species'
        sp_tables = dir_curate / 'Homo_sapiens' / 'tables'
        sp_tables.mkdir(parents=True)
        (sp_tables / 'Homo_sapiens_tpm.tsv').write_text('data')
        real_target = tmp_path / 'real_target'
        real_target.mkdir()
        (real_target / 'keep.txt').write_text('keep')
        dir_csca_input = tmp_path / 'csca_input'
        os.symlink(real_target, dir_csca_input)

        generate_input_symlinks(str(dir_csca_input), str(dir_curate), ['Homo_sapiens'])

        assert os.path.isdir(str(dir_csca_input))
        assert not os.path.islink(str(dir_csca_input))
        assert os.path.islink(str(dir_csca_input / 'Homo_sapiens_tpm.tsv'))
        assert os.path.exists(str(real_target / 'keep.txt'))

    def test_ignores_non_tsv_and_subdirectories(self, tmp_path):
        dir_curate = tmp_path / 'per_species'
        sp_tables = dir_curate / 'Homo_sapiens' / 'tables'
        sp_tables.mkdir(parents=True)
        (sp_tables / 'Homo_sapiens_tpm.tsv').write_text('data')
        (sp_tables / 'notes.txt').write_text('ignore')
        (sp_tables / 'nested').mkdir()
        dir_csca_input = str(tmp_path / 'csca_input')

        generate_input_symlinks(dir_csca_input, str(dir_curate), ['Homo_sapiens'])

        assert os.path.islink(os.path.join(dir_csca_input, 'Homo_sapiens_tpm.tsv'))
        assert not os.path.exists(os.path.join(dir_csca_input, 'notes.txt'))
        assert not os.path.exists(os.path.join(dir_csca_input, 'nested'))

    def test_raises_clear_error_when_species_tables_dir_missing(self, tmp_path):
        dir_curate = tmp_path / 'per_species'
        (dir_curate / 'Homo_sapiens').mkdir(parents=True)
        dir_csca_input = str(tmp_path / 'csca_input')

        with pytest.raises(FileNotFoundError, match='Per-species tables directory not found for species'):
            generate_input_symlinks(dir_csca_input, str(dir_curate), ['Homo_sapiens'])

    def test_raises_when_species_tables_dir_has_no_tsv_files(self, tmp_path):
        dir_curate = tmp_path / 'per_species'
        sp_tables = dir_curate / 'Homo_sapiens' / 'tables'
        sp_tables.mkdir(parents=True)
        (sp_tables / 'notes.txt').write_text('not a tsv')
        dir_csca_input = str(tmp_path / 'csca_input')

        with pytest.raises(FileNotFoundError, match='No TSV table file was found in per-species tables directory'):
            generate_input_symlinks(dir_csca_input, str(dir_curate), ['Homo_sapiens'])

    def test_raises_on_duplicate_table_filename_across_species(self, tmp_path):
        dir_curate = tmp_path / 'per_species'
        hs_tables = dir_curate / 'Homo_sapiens' / 'tables'
        mm_tables = dir_curate / 'Mus_musculus' / 'tables'
        hs_tables.mkdir(parents=True)
        mm_tables.mkdir(parents=True)
        (hs_tables / 'shared.tsv').write_text('human')
        (mm_tables / 'shared.tsv').write_text('mouse')
        dir_csca_input = str(tmp_path / 'csca_input')

        with pytest.raises(ValueError, match='Duplicate table filename across species'):
            generate_input_symlinks(
                dir_csca_input,
                str(dir_curate),
                ['Homo_sapiens', 'Mus_musculus'],
            )

    def test_raises_when_csca_input_path_is_a_file(self, tmp_path):
        dir_curate = tmp_path / 'per_species'
        sp_tables = dir_curate / 'Homo_sapiens' / 'tables'
        sp_tables.mkdir(parents=True)
        (sp_tables / 'Homo_sapiens_tpm.tsv').write_text('data')
        dir_csca_input = tmp_path / 'csca_input'
        dir_csca_input.write_text('not a directory')

        with pytest.raises(NotADirectoryError, match='Cross-species input path exists but is not a directory'):
            generate_input_symlinks(str(dir_csca_input), str(dir_curate), ['Homo_sapiens'])

    def test_raises_when_symlink_destination_is_directory(self, tmp_path, monkeypatch):
        dir_curate = tmp_path / 'per_species'
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

        monkeypatch.setattr('amalgkit.cross_species_filter.os.path.lexists', fake_lexists)
        monkeypatch.setattr('amalgkit.cross_species_filter.os.path.isdir', fake_isdir)

        with pytest.raises(IsADirectoryError, match='Cross-species input destination exists but is a directory'):
            generate_input_symlinks(str(dir_csca_input), str(dir_curate), ['Homo_sapiens'])


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
        monkeypatch.setattr('amalgkit.cross_species_filter.load_metadata', lambda _args: metadata)
        out = get_sample_group_string(args)
        assert set(out.split('|')) == {'treated', 'control'}

    def test_metadata_sample_groups_drop_blank_and_deduplicate(self, monkeypatch):
        args = SimpleNamespace(sample_group=None)
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['R1', 'R2', 'R3', 'R4'],
            'sample_group': [' treated ', '', 'treated', 'control'],
            'exclusion': ['no', 'no', 'no', 'no'],
        }))
        monkeypatch.setattr('amalgkit.cross_species_filter.load_metadata', lambda _args: metadata)
        assert get_sample_group_string(args) == 'treated|control'

    def test_metadata_sample_groups_drop_nan(self, monkeypatch):
        args = SimpleNamespace(sample_group=None)
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['R1', 'R2', 'R3'],
            'sample_group': ['treated', float('nan'), 'control'],
            'exclusion': ['no', 'no', 'no'],
        }))
        monkeypatch.setattr('amalgkit.cross_species_filter.load_metadata', lambda _args: metadata)
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
        monkeypatch.setattr('amalgkit.cross_species_filter.load_metadata', lambda _args: metadata)
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
            redo=False,
        )

    @staticmethod
    def _write_cross_species_fixture(out_dir):
        species_specs = [
            {
                'species_tag': 'Species_A',
                'scientific_name': 'Species A',
                'runs': ['A_leaf', 'A_root'],
                'groups': ['leaf', 'root'],
                'uncorrected': pandas.DataFrame(
                    {
                        'A_leaf': [10.0, 9.0, 2.0, 1.0],
                        'A_root': [2.0, 1.0, 10.0, 9.0],
                    },
                    index=['A1', 'A2', 'A3', 'A4'],
                ),
                'corrected': pandas.DataFrame(
                    {
                        'A_leaf': [9.5, 8.5, 2.5, 1.5],
                        'A_root': [2.5, 1.5, 9.5, 8.5],
                    },
                    index=['A1', 'A2', 'A3', 'A4'],
                ),
            },
            {
                'species_tag': 'Species_B',
                'scientific_name': 'Species B',
                'runs': ['B_leaf', 'B_root'],
                'groups': ['leaf', 'root'],
                'uncorrected': pandas.DataFrame(
                    {
                        'B_leaf': [8.0, 7.0, 3.0, 2.0],
                        'B_root': [3.0, 2.0, 8.0, 7.0],
                    },
                    index=['B1', 'B2', 'B3', 'B4'],
                ),
                'corrected': pandas.DataFrame(
                    {
                        'B_leaf': [8.5, 7.5, 2.5, 1.5],
                        'B_root': [2.5, 1.5, 8.5, 7.5],
                    },
                    index=['B1', 'B2', 'B3', 'B4'],
                ),
            },
        ]
        per_species_dir = out_dir / 'per_species'
        for spec in species_specs:
            tables_dir = per_species_dir / spec['species_tag'] / 'tables'
            tables_dir.mkdir(parents=True)
            metadata_df = pandas.DataFrame(
                {
                    'run': spec['runs'],
                    'scientific_name': [spec['scientific_name']] * len(spec['runs']),
                    'sample_group': spec['groups'],
                    'bioproject': ['BP1', 'BP2'],
                    'exclusion': ['no'] * len(spec['runs']),
                }
            )
            metadata_df.to_csv(tables_dir / f"{spec['species_tag']}.metadata.tsv", sep='\t', index=False)
            spec['uncorrected'].to_csv(tables_dir / f"{spec['species_tag']}.uncorrected.tc.tsv", sep='\t', index_label='target_id')
            spec['corrected'].to_csv(tables_dir / f"{spec['species_tag']}.no.tc.tsv", sep='\t', index_label='target_id')
            pandas.DataFrame(
                {
                    'bwbw_n': [1, 1],
                    'bwbw_mean': [0.30, 0.55],
                    'bwbw_median': [0.30, 0.55],
                    'bwbw_variance': [0.0, 0.0],
                    'wibw_n': [1, 1],
                    'wibw_mean': [0.10, 0.20],
                    'wibw_median': [0.10, 0.20],
                    'wibw_variance': [0.0, 0.0],
                    'bwwi_n': [1, 1],
                    'bwwi_mean': [0.40, 0.62],
                    'bwwi_median': [0.40, 0.62],
                    'bwwi_variance': [0.0, 0.0],
                    'wiwi_n': [1, 1],
                    'wiwi_mean': [0.18, 0.22],
                    'wiwi_median': [0.18, 0.22],
                    'wiwi_variance': [0.0, 0.0],
                },
                index=['round_0', 'round_1'],
            ).to_csv(tables_dir / f"{spec['species_tag']}.no.correlation_statistics.tsv", sep='\t')
        orthogroup_path = out_dir.parent / 'orthogroup.tsv'
        pandas.DataFrame(
            {
                'busco_id': ['OG1', 'OG2', 'OG3', 'OG4'],
                'Species_A': ['A1', 'A2', 'A3', 'A4'],
                'Species_B': ['B1', 'B2', 'B3', 'B4'],
            }
        ).to_csv(orthogroup_path, sep='\t', index=False)
        return orthogroup_path

    def test_raises_clear_error_when_curate_dir_missing(self, tmp_path, monkeypatch):
        args = self._base_args(tmp_path / 'out')
        monkeypatch.setattr('amalgkit.cross_species_filter.check_ortholog_parameter_compatibility', lambda _args: None)
        with pytest.raises(FileNotFoundError, match='Per-species table output directory not found'):
            run_cross_species_filter(args)

    def test_blank_orthogroup_table_raises_parameter_error(self, tmp_path, monkeypatch):
        args = self._base_args(tmp_path / 'out')
        args.orthogroup_table = '   '
        with pytest.raises(ValueError, match='One of --orthogroup_table and --dir_busco should be specified'):
            run_cross_species_filter(args)

    def test_raises_when_out_dir_is_file(self, tmp_path, monkeypatch):
        out_path = tmp_path / 'out_path'
        out_path.write_text('not a directory')
        args = self._base_args(out_path)
        monkeypatch.setattr('amalgkit.cross_species_filter.check_ortholog_parameter_compatibility', lambda _args: None)
        with pytest.raises(NotADirectoryError, match='Output path exists but is not a directory'):
            run_cross_species_filter(args)

    def test_raises_when_curate_path_is_file(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'out'
        out_dir.mkdir()
        (out_dir / 'per_species').write_text('not a directory')
        args = self._base_args(out_dir)
        monkeypatch.setattr('amalgkit.cross_species_filter.check_ortholog_parameter_compatibility', lambda _args: None)
        with pytest.raises(NotADirectoryError, match='Per-species output path exists but is not a directory'):
            run_cross_species_filter(args)

    def test_raises_when_csca_path_is_file(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'out'
        curate_dir = out_dir / 'per_species' / 'Species_A' / 'tables'
        curate_dir.mkdir(parents=True)
        (curate_dir / 'Species_A_est_counts.tsv').write_text('target_id\tR1\nG1\t1\n')
        (out_dir / 'cross_species').write_text('not a directory')
        args = self._base_args(out_dir)
        args.orthogroup_table = str(tmp_path / 'orthogroup.tsv')
        (tmp_path / 'orthogroup.tsv').write_text('busco_id\tSpecies_A\nOG1\tgene1\n')
        monkeypatch.setattr('amalgkit.cross_species_filter.check_ortholog_parameter_compatibility', lambda _args: None)
        with pytest.raises(NotADirectoryError, match='Cross-species path exists but is not a directory'):
            run_cross_species_filter(args)

    def test_existing_output_requires_redo(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'out'
        tables_dir = out_dir / 'per_species' / 'Species_A' / 'tables'
        tables_dir.mkdir(parents=True)
        (tables_dir / 'Species_A_est_counts.tsv').write_text('target_id\tR1\nG1\t1\n')
        existing_dir = out_dir / 'cross_species'
        existing_dir.mkdir()
        orthogroup = tmp_path / 'orthogroup.tsv'
        orthogroup.write_text('busco_id\tSpecies_A\nOG1\tgene1\n')
        args = self._base_args(out_dir)
        args.orthogroup_table = str(orthogroup)

        with pytest.raises(FileExistsError, match='Use --redo yes to overwrite'):
            run_cross_species_filter(args)

    def test_raises_when_curate_dir_has_no_species_subdirs(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'out'
        (out_dir / 'per_species').mkdir(parents=True)
        args = self._base_args(out_dir)
        monkeypatch.setattr('amalgkit.cross_species_filter.check_ortholog_parameter_compatibility', lambda _args: None)
        with pytest.raises(ValueError, match='No per-species directories were found'):
            run_cross_species_filter(args)

    def test_raises_when_orthogroup_table_path_missing(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'out'
        tables_dir = out_dir / 'per_species' / 'Species_A' / 'tables'
        tables_dir.mkdir(parents=True)
        (tables_dir / 'Species_A_est_counts.tsv').write_text('target_id\tR1\nG1\t1\n')
        args = self._base_args(out_dir)
        args.orthogroup_table = str(tmp_path / 'missing.tsv')
        monkeypatch.setattr('amalgkit.cross_species_filter.check_ortholog_parameter_compatibility', lambda _args: None)

        with pytest.raises(FileNotFoundError, match='Orthogroup table not found'):
            run_cross_species_filter(args)

    def test_raises_when_orthogroup_table_path_is_directory(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'out'
        tables_dir = out_dir / 'per_species' / 'Species_A' / 'tables'
        tables_dir.mkdir(parents=True)
        (tables_dir / 'Species_A_est_counts.tsv').write_text('target_id\tR1\nG1\t1\n')
        og_dir = tmp_path / 'orthogroup_dir'
        og_dir.mkdir()
        args = self._base_args(out_dir)
        args.orthogroup_table = str(og_dir)
        monkeypatch.setattr('amalgkit.cross_species_filter.check_ortholog_parameter_compatibility', lambda _args: None)

        with pytest.raises(IsADirectoryError, match='Orthogroup table path exists but is not a file'):
            run_cross_species_filter(args)

    def test_restores_existing_output_when_staged_run_fails(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'out'
        tables_dir = out_dir / 'per_species' / 'Species_A' / 'tables'
        tables_dir.mkdir(parents=True)
        (tables_dir / 'Species_A_est_counts.tsv').write_text('target_id\tR1\nG1\t1\n')
        existing_dir = out_dir / 'cross_species'
        existing_dir.mkdir()
        (existing_dir / 'old.txt').write_text('old')
        orthogroup = tmp_path / 'orthogroup.tsv'
        orthogroup.write_text('busco_id\tSpecies_A\nOG1\tgene1\n')
        args = self._base_args(out_dir)
        args.orthogroup_table = str(orthogroup)
        args.redo = True

        monkeypatch.setattr(
            'amalgkit.cross_species_filter.orthogroup2genecount',
            lambda **_kwargs: (_ for _ in ()).throw(RuntimeError('boom')),
        )

        with pytest.raises(RuntimeError, match='boom'):
            run_cross_species_filter(args)

        assert (existing_dir / 'old.txt').read_text() == 'old'

    def test_run_cross_species_filter_writes_restored_plot_outputs(self, tmp_path):
        out_dir = tmp_path / 'out'
        orthogroup_path = self._write_cross_species_fixture(out_dir)
        args = self._base_args(out_dir)
        args.sample_group = 'leaf,root'
        args.orthogroup_table = str(orthogroup_path)
        args.missing_strategy = 'row_mean'

        run_cross_species_filter(args)

        cross_species_dir = out_dir / 'cross_species'
        assert (cross_species_dir / 'metadata.tsv').is_file()
        assert (cross_species_dir / 'cross_species_sample_number_heatmap.pdf').is_file()
        assert (cross_species_dir / 'cross_species_group_cor_scatter.pdf').is_file()
        assert (cross_species_dir / 'cross_species_unaveraged_pca_PC12.pdf').is_file()
        assert (cross_species_dir / 'cross_species_unaveraged_pca_PC12_uncorrected.pdf').is_file()
        assert (cross_species_dir / 'cross_species_unaveraged_pca_PC12_corrected.pdf').is_file()
        assert (cross_species_dir / 'cross_species_unaveraged_pca_PC34_uncorrected.pdf').is_file()
        assert (cross_species_dir / 'cross_species_unaveraged_pca_PC34_corrected.pdf').is_file()
        assert (cross_species_dir / 'cross_species_unaveraged_tsne_uncorrected.pdf').is_file()
        assert (cross_species_dir / 'cross_species_unaveraged_tsne_corrected.pdf').is_file()
        assert (cross_species_dir / 'cross_species_SVA_heatmap.pdf').is_file()
        assert (cross_species_dir / 'cross_species_SVA_dendrogram.pdf').is_file()
        assert (cross_species_dir / 'cross_species_averaged_summary.pdf').is_file()
        assert (cross_species_dir / 'cross_species_boxplot.pdf').is_file()
        assert (cross_species_dir / 'cross_species_averaged_tsne.pdf').is_file()
        assert (cross_species_dir / 'cross_species_delta_pcc_boxplot.pdf').is_file()
