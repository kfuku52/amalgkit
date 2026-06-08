import numpy
import pytest
import pandas
from types import SimpleNamespace

from amalgkit.cstmm import filepath2spp, get_count_files, cstmm_main
from amalgkit.normalization_tmm import run_tmm_rounds_for_cstmm


def _read_dcf(path):
    config = {}
    with open(path, encoding='utf-8') as handle:
        for raw_line in handle:
            line = raw_line.rstrip('\n')
            if line == '':
                continue
            key, value = line.split(':', 1)
            config[key] = value.lstrip()
    return config


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
    @staticmethod
    def _write_multi_species_fixture(base_dir):
        merge_dir = base_dir / 'merge'
        species_a_dir = merge_dir / 'Species_A'
        species_b_dir = merge_dir / 'Species_B'
        species_a_dir.mkdir(parents=True)
        species_b_dir.mkdir(parents=True)
        (merge_dir / 'metadata.tsv').write_text(
            'run\tscientific_name\texclusion\tsample_group\tmapping_rate\n'
            'R1\tSpecies A\tno\tleaf\t100\n'
            'R2\tSpecies A\tno\troot\t100\n'
            'R1\tSpecies B\tno\tleaf\t100\n'
            'R2\tSpecies B\tno\troot\t100\n'
        )
        (species_a_dir / 'Species_A_est_counts.tsv').write_text(
            'target_id\tR1\tR2\n'
            'geneA1\t100\t120\n'
            'geneA2\t50\t60\n'
            'geneA3\t30\t40\n'
            'geneA4\t5\t8\n'
        )
        (species_b_dir / 'Species_B_est_counts.tsv').write_text(
            'target_id\tR1\tR2\n'
            'geneB1\t80\t0\n'
            'geneB2\t60\t0\n'
            'geneB3\t20\t0\n'
            'geneB4\t7\t0\n'
        )
        for species_name, genes in [
            ('Species_A', ['geneA1', 'geneA2', 'geneA3', 'geneA4']),
            ('Species_B', ['geneB1', 'geneB2', 'geneB3', 'geneB4']),
        ]:
            eff_length_rows = ['target_id\tR1\tR2'] + ['{}\t1000\t1000'.format(gene_id) for gene_id in genes]
            (merge_dir / species_name / '{}_eff_length.tsv'.format(species_name)).write_text('\n'.join(eff_length_rows) + '\n')
        orthogroup = base_dir / 'orthogroup.tsv'
        orthogroup.write_text(
            'busco_id\tSpecies_A\tSpecies_B\n'
            'OG1\tgeneA1\tgeneB1\n'
            'OG2\tgeneA2\tgeneB2\n'
            'OG3\tgeneA3\tgeneB3\n'
        )
        return merge_dir, orthogroup

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
        (merge_dir / 'metadata.tsv').write_text(
            'run\tscientific_name\texclusion\tsample_group\tmapping_rate\n'
            'R1\tSpecies A\tno\tleaf\t100\n'
        )
        (species_dir / 'Species_A_est_counts.tsv').write_text('target_id\tR1\nG1\t1\n')
        (species_dir / 'Species_A_eff_length.tsv').write_text('target_id\tR1\nG1\t1000\n')
        args = SimpleNamespace(
            out_dir=str(out_dir),
            dir_count='inferred',
            dir_busco=None,
            orthogroup_table=None,
        )
        monkeypatch.setattr(
            'amalgkit.cstmm.orthogroup2genecount',
            lambda **_kwargs: (_ for _ in ()).throw(AssertionError('orthogroup2genecount should not run for single species')),
        )

        cstmm_main(args)

        assert (out_dir / 'cstmm' / 'Species_A' / 'Species_A_cstmm_counts.tsv').is_file()
        assert (out_dir / 'cstmm' / 'metadata.tsv').is_file()

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

        with pytest.raises(FileNotFoundError, match='No est_counts.tsv file found for species directory\\(ies\\): Species_B'):
            cstmm_main(args)

    def test_existing_output_requires_redo(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'out'
        merge_dir = out_dir / 'merge'
        for species in ['Species_A', 'Species_B']:
            species_dir = merge_dir / species
            species_dir.mkdir(parents=True)
            (species_dir / f'{species}_est_counts.tsv').write_text('target_id\tR1\nG1\t1\n')
        (out_dir / 'cstmm').mkdir()
        orthogroup = tmp_path / 'orthogroup.tsv'
        orthogroup.write_text('busco_id\tSpecies_A\tSpecies_B\nOG1\tgene1\tgene2\n')
        args = SimpleNamespace(
            out_dir=str(out_dir),
            dir_count='inferred',
            dir_busco=None,
            orthogroup_table=str(orthogroup),
            redo=False,
        )

        with pytest.raises(FileExistsError, match='Use --redo yes to overwrite'):
            cstmm_main(args)

    def test_restores_existing_output_when_staged_run_fails(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'out'
        merge_dir = out_dir / 'merge'
        for species in ['Species_A', 'Species_B']:
            species_dir = merge_dir / species
            species_dir.mkdir(parents=True)
            (species_dir / f'{species}_est_counts.tsv').write_text('target_id\tR1\nG1\t1\n')
        existing_dir = out_dir / 'cstmm'
        existing_dir.mkdir()
        (existing_dir / 'old.txt').write_text('old')
        orthogroup = tmp_path / 'orthogroup.tsv'
        orthogroup.write_text('busco_id\tSpecies_A\tSpecies_B\nOG1\tgene1\tgene2\n')
        args = SimpleNamespace(
            out_dir=str(out_dir),
            dir_count='inferred',
            dir_busco=None,
            orthogroup_table=str(orthogroup),
            redo=True,
        )

        monkeypatch.setattr(
            'amalgkit.cstmm.orthogroup2genecount',
            lambda **_kwargs: (_ for _ in ()).throw(RuntimeError('boom')),
        )

        with pytest.raises(RuntimeError, match='boom'):
            cstmm_main(args)

        assert (existing_dir / 'old.txt').read_text() == 'old'

    def test_python_single_species_backend_writes_counts_metadata_and_plots(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'out'
        merge_dir = out_dir / 'merge'
        species_dir = merge_dir / 'Species_A'
        species_dir.mkdir(parents=True)
        (merge_dir / 'metadata.tsv').write_text(
            'run\tscientific_name\texclusion\tsample_group\tmapping_rate\n'
            'R1\tSpecies A\tno\tleaf\t100\n'
            'R2\tSpecies A\tno\troot\t100\n'
        )
        (species_dir / 'Species_A_est_counts.tsv').write_text(
            'target_id\tR1\tR2\n'
            'G1\t100\t90\n'
            'G2\t110\t100\n'
            'G3\t120\t110\n'
            'G4\t130\t120\n'
            'G5\t10\t20\n'
            'G6\t11\t21\n'
        )
        (species_dir / 'Species_A_eff_length.tsv').write_text(
            'target_id\tR1\tR2\n'
            'G1\t1000\t1000\n'
            'G2\t1000\t1000\n'
            'G3\t1000\t1000\n'
            'G4\t1000\t1000\n'
            'G5\t1000\t1000\n'
            'G6\t1000\t1000\n'
        )
        args = SimpleNamespace(
            out_dir=str(out_dir),
            dir_count='inferred',
            dir_busco=None,
            orthogroup_table=None,
            redo=False,
            tmm_backend='python',
        )
        cleaned = {}
        monkeypatch.setattr('amalgkit.cstmm.cleanup_tmp_amalgkit_files', lambda work_dir: cleaned.setdefault('work_dir', work_dir))

        cstmm_main(args)

        counts_path = out_dir / 'cstmm' / 'Species_A' / 'Species_A_cstmm_counts.tsv'
        metadata_path = out_dir / 'cstmm' / 'metadata.tsv'
        eff_length_path = out_dir / 'cstmm' / 'Species_A' / 'Species_A_eff_length.tsv'
        assert counts_path.is_file()
        assert metadata_path.is_file()
        assert eff_length_path.is_file()
        assert (out_dir / 'cstmm' / 'cstmm_normalization_factor_histogram.scientific_name.pdf').is_file()
        assert (out_dir / 'cstmm' / 'cstmm_normalization_factor_histogram.sample_group.pdf').is_file()
        assert (out_dir / 'cstmm' / 'cstmm_normalization_factor_scatter.pdf').is_file()
        assert (out_dir / 'cstmm' / 'cstmm_mean_expression_boxplot.pdf').is_file()
        assert (out_dir / 'cstmm' / 'cstmm_exclusion.pdf').is_file()
        assert cleaned['work_dir'] == '.'

        corrected = pandas.read_csv(counts_path, sep='\t')
        corrected_matrix = corrected.set_index('target_id')
        original = pandas.DataFrame(
            {'Species_A_R1': [100, 110, 120, 130, 10, 11], 'Species_A_R2': [90, 100, 110, 120, 20, 21]},
            index=['G1', 'G2', 'G3', 'G4', 'G5', 'G6'],
        ).astype(float)
        roundtrip = run_tmm_rounds_for_cstmm(original)
        expected = original.copy()
        for sample_id, factor in roundtrip.round2_factors.items():
            expected.loc[:, sample_id] = expected.loc[:, sample_id] / factor
        expected.columns = ['R1', 'R2']
        expected.index.name = 'target_id'
        pandas.testing.assert_frame_equal(corrected_matrix, expected)

        metadata_df = pandas.read_csv(metadata_path, sep='\t')
        assert 'tmm_library_size' in metadata_df.columns
        assert 'tmm_normalization_factor' in metadata_df.columns
        assert set(metadata_df['exclusion']) == {'no'}

    def test_python_multi_species_backend_writes_counts_metadata_and_plots(self, tmp_path, monkeypatch):
        out_dir = tmp_path / 'out'
        merge_dir, orthogroup = self._write_multi_species_fixture(out_dir)
        args = SimpleNamespace(
            out_dir=str(out_dir),
            dir_count='inferred',
            dir_busco=None,
            orthogroup_table=str(orthogroup),
            redo=False,
            tmm_backend='python',
        )
        cleaned = {}
        monkeypatch.setattr('amalgkit.cstmm.cleanup_tmp_amalgkit_files', lambda work_dir: cleaned.setdefault('work_dir', work_dir))

        cstmm_main(args)

        metadata_path = out_dir / 'cstmm' / 'metadata.tsv'
        species_a_counts_path = out_dir / 'cstmm' / 'Species_A' / 'Species_A_cstmm_counts.tsv'
        species_b_counts_path = out_dir / 'cstmm' / 'Species_B' / 'Species_B_cstmm_counts.tsv'
        assert metadata_path.is_file()
        assert species_a_counts_path.is_file()
        assert species_b_counts_path.is_file()
        assert (out_dir / 'cstmm' / 'Species_A' / 'Species_A_eff_length.tsv').is_file()
        assert (out_dir / 'cstmm' / 'Species_B' / 'Species_B_eff_length.tsv').is_file()
        assert (out_dir / 'cstmm' / 'cstmm_normalization_factor_histogram.scientific_name.pdf').is_file()
        assert (out_dir / 'cstmm' / 'cstmm_normalization_factor_histogram.sample_group.pdf').is_file()
        assert (out_dir / 'cstmm' / 'cstmm_normalization_factor_scatter.pdf').is_file()
        assert (out_dir / 'cstmm' / 'cstmm_mean_expression_boxplot.pdf').is_file()
        assert (out_dir / 'cstmm' / 'cstmm_exclusion.pdf').is_file()
        assert cleaned['work_dir'] == '.'

        df_sog = pandas.DataFrame(
            {
                'Species_A_R1': [100.0, 50.0, 30.0],
                'Species_A_R2': [120.0, 60.0, 40.0],
                'Species_B_R1': [80.0, 60.0, 20.0],
            },
            index=['OG1', 'OG2', 'OG3'],
        )
        library_sizes = pandas.Series(
            {
                'Species_A_R1': 185.0,
                'Species_A_R2': 228.0,
                'Species_B_R1': 167.0,
            }
        )
        roundtrip = run_tmm_rounds_for_cstmm(counts=df_sog, lib_size=library_sizes)

        expected_species_a = pandas.DataFrame(
            {
                'R1': [100.0, 50.0, 30.0, 5.0],
                'R2': [120.0, 60.0, 40.0, 8.0],
            },
            index=['geneA1', 'geneA2', 'geneA3', 'geneA4'],
        )
        expected_species_b = pandas.DataFrame(
            {
                'R1': [80.0, 60.0, 20.0, 7.0],
                'R2': [0.0, 0.0, 0.0, 0.0],
            },
            index=['geneB1', 'geneB2', 'geneB3', 'geneB4'],
        )
        expected_species_a.loc[:, 'R1'] = expected_species_a.loc[:, 'R1'] / roundtrip.round2_factors['Species_A_R1']
        expected_species_a.loc[:, 'R2'] = expected_species_a.loc[:, 'R2'] / roundtrip.round2_factors['Species_A_R2']
        expected_species_b.loc[:, 'R1'] = expected_species_b.loc[:, 'R1'] / roundtrip.round2_factors['Species_B_R1']
        expected_species_a.index.name = 'target_id'
        expected_species_b.index.name = 'target_id'

        pandas.testing.assert_frame_equal(
            pandas.read_csv(species_a_counts_path, sep='\t').set_index('target_id'),
            expected_species_a,
        )
        pandas.testing.assert_frame_equal(
            pandas.read_csv(species_b_counts_path, sep='\t').set_index('target_id'),
            expected_species_b,
        )

        metadata_df = pandas.read_csv(metadata_path, sep='\t')
        assert 'tmm_library_size' in metadata_df.columns
        assert 'tmm_normalization_factor' in metadata_df.columns
        retained = metadata_df.loc[metadata_df['scientific_name'] != 'Species B', 'exclusion']
        assert set(retained) == {'no'}
        excluded = metadata_df.loc[
            (metadata_df['scientific_name'] == 'Species B') & (metadata_df['run'] == 'R2'),
            'exclusion',
        ]
        assert excluded.tolist() == ['no_cstmm_output']
