import os
import warnings

import pytest
import pandas
import numpy
import xml.etree.ElementTree as ET

from amalgkit.exceptions import AmalgkitExit
from amalgkit.util import (
    Metadata,
    get_sra_stat,
    check_ortholog_parameter_compatibility,
    orthogroup2genecount,
    load_metadata,
    cleanup_tmp_amalgkit_files,
)


# ---------------------------------------------------------------------------
# strtobool
# ---------------------------------------------------------------------------


class TestMetadataInit:
    def test_empty_metadata(self):
        m = Metadata()
        assert isinstance(m.df, pandas.DataFrame)
        assert m.df.shape[0] == 0
        assert 'scientific_name' in m.df.columns
        assert 'run' in m.df.columns

    def test_column_names_present(self):
        m = Metadata()
        for col in ['tissue', 'sample_group', 'bioproject', 'biosample',
                     'lib_layout', 'total_spots', 'exclusion', 'ENA_SRA_Link', 'DDBJ_SRA_Link']:
            assert col in m.df.columns


class TestMetadataReorder:
    def test_reorder_empty(self):
        m = Metadata()
        result = m.reorder()
        assert result is None  # returns None for empty df

    def test_reorder_empty_normalizes_schema_and_drops_removed_columns(self):
        m = Metadata()
        m.df = pandas.DataFrame(columns=['lab', 'batch', 'misc'])
        result = m.reorder()
        assert result is None
        assert 'sample_group' in m.df.columns
        assert 'exclusion' in m.df.columns
        for col in Metadata.removed_metadata_columns:
            assert col not in m.df.columns

    def test_reorder_sets_exclusion_default(self, sample_metadata_df):
        df = sample_metadata_df.copy()
        df['exclusion'] = ''
        m = Metadata()
        m.df = df
        m.reorder()
        assert (m.df['exclusion'] == 'no').all()

    def test_reorder_sample_group_near_front(self, sample_metadata):
        cols = list(sample_metadata.df.columns)
        assert cols.index('sample_group') == 1


class TestMetadataFromDataFrame:
    def test_roundtrip(self, sample_metadata_df):
        m = Metadata.from_DataFrame(sample_metadata_df)
        assert isinstance(m, Metadata)
        assert m.df.shape[0] == 5
        assert 'scientific_name' in m.df.columns

    def test_drops_removed_legacy_columns(self):
        df = pandas.DataFrame({
            'scientific_name': ['Homo sapiens'],
            'run': ['SRR000001'],
            'exclusion': ['no'],
            'lab': ['Legacy Lab'],
            'biomaterial_provider': ['Legacy Provider'],
            'cell': ['HeLa'],
            'location': ['Tokyo'],
            'antibody': ['H3K4me3'],
            'batch': ['B1'],
            'misc': ['legacy'],
        })
        m = Metadata.from_DataFrame(df)
        for col in Metadata.removed_metadata_columns:
            assert col not in m.df.columns

    def test_exclusion_filled(self, sample_metadata_df):
        df = sample_metadata_df.copy()
        df['exclusion'] = ''
        m = Metadata.from_DataFrame(df)
        assert (m.df['exclusion'] == 'no').all()

    def test_exclusion_nan_filled(self, sample_metadata_df):
        df = sample_metadata_df.copy()
        df['exclusion'] = [numpy.nan] * len(df)
        m = Metadata.from_DataFrame(df)
        assert (m.df['exclusion'] == 'no').all()

    def test_does_not_mutate_input_dataframe(self, sample_metadata_df):
        df = sample_metadata_df.copy()
        df['exclusion'] = ''
        _ = Metadata.from_DataFrame(df)
        assert (df['exclusion'] == '').all()

    def test_empty_input_dataframe_gets_standard_columns(self):
        m = Metadata.from_DataFrame(pandas.DataFrame())
        for col in ['sample_group', 'exclusion', 'is_sampled', 'is_qualified']:
            assert col in m.df.columns


class TestMetadataFromXml:
    def test_parse_minimal_xml(self):
        xml_str = b"""<?xml version="1.0" encoding="UTF-8"?>
        <EXPERIMENT_PACKAGE_SET>
          <EXPERIMENT_PACKAGE>
            <EXPERIMENT alias="test">
              <IDENTIFIERS><PRIMARY_ID>SRX000001</PRIMARY_ID></IDENTIFIERS>
              <TITLE>Test experiment</TITLE>
              <STUDY_REF>
                <IDENTIFIERS><PRIMARY_ID>SRP000001</PRIMARY_ID></IDENTIFIERS>
              </STUDY_REF>
              <DESIGN>
                <DESIGN_DESCRIPTION>test</DESIGN_DESCRIPTION>
                <LIBRARY_DESCRIPTOR>
                  <LIBRARY_NAME>testlib</LIBRARY_NAME>
                  <LIBRARY_STRATEGY>RNA-Seq</LIBRARY_STRATEGY>
                  <LIBRARY_SOURCE>TRANSCRIPTOMIC</LIBRARY_SOURCE>
                  <LIBRARY_SELECTION>cDNA</LIBRARY_SELECTION>
                  <LIBRARY_LAYOUT><PAIRED NOMINAL_LENGTH="200"/></LIBRARY_LAYOUT>
                </LIBRARY_DESCRIPTOR>
              </DESIGN>
              <PLATFORM><ILLUMINA><INSTRUMENT_MODEL>HiSeq 2500</INSTRUMENT_MODEL></ILLUMINA></PLATFORM>
            </EXPERIMENT>
            <SUBMISSION accession="SRA000001" lab_name="TestLab" center_name="TestCenter">
              <IDENTIFIERS>
                <PRIMARY_ID>SRA000001</PRIMARY_ID>
                <SUBMITTER_ID>sub1</SUBMITTER_ID>
              </IDENTIFIERS>
            </SUBMISSION>
            <STUDY>
              <DESCRIPTOR><STUDY_TITLE>Test Study</STUDY_TITLE></DESCRIPTOR>
            </STUDY>
            <SAMPLE>
              <IDENTIFIERS><PRIMARY_ID>SRS000001</PRIMARY_ID></IDENTIFIERS>
              <TITLE>Test sample</TITLE>
              <SAMPLE_NAME>
                <SCIENTIFIC_NAME>Homo sapiens</SCIENTIFIC_NAME>
                <TAXON_ID>9606</TAXON_ID>
              </SAMPLE_NAME>
              <DESCRIPTION>A test sample</DESCRIPTION>
              <SAMPLE_ATTRIBUTES>
                <SAMPLE_ATTRIBUTE><TAG>tissue</TAG><VALUE>brain</VALUE></SAMPLE_ATTRIBUTE>
                <SAMPLE_ATTRIBUTE><TAG>cell</TAG><VALUE>HeLa</VALUE></SAMPLE_ATTRIBUTE>
              </SAMPLE_ATTRIBUTES>
            </SAMPLE>
            <RUN_SET>
              <RUN accession="SRR000001" total_spots="1000000" total_bases="200000000" size="50000000" published="2020-01-01">
                <IDENTIFIERS><PRIMARY_ID>SRR000001</PRIMARY_ID></IDENTIFIERS>
                <SRAFiles>
                  <SRAFile supertype="Primary ETL">
                    <Alternatives org="NCBI" url="https://ncbi.example.com/SRR000001"/>
                  </SRAFile>
                </SRAFiles>
              </RUN>
            </RUN_SET>
            <Pool>
              <Member>
                <EXTERNAL_ID namespace="BioProject">PRJNA000001</EXTERNAL_ID>
                <EXTERNAL_ID namespace="BioSample">SAMN000001</EXTERNAL_ID>
              </Member>
            </Pool>
          </EXPERIMENT_PACKAGE>
        </EXPERIMENT_PACKAGE_SET>"""
        root = ET.fromstring(xml_str)
        tree = ET.ElementTree(root)
        m = Metadata.from_xml(tree)
        assert m.df.shape[0] == 1
        assert m.df.loc[0, 'scientific_name'] == 'Homo sapiens'
        assert m.df.loc[0, 'run'] == 'SRR000001'
        assert m.df.loc[0, 'lib_layout'] == 'paired'
        assert m.df.loc[0, 'bioproject'] == 'PRJNA000001'
        assert m.df.loc[0, 'ENA_SRA_Link'] == ''
        assert m.df.loc[0, 'DDBJ_SRA_Link'] == ''
        for col in Metadata.removed_metadata_columns:
            assert col not in m.df.columns

    def test_parse_drr_xml_derives_ddbj_sra_link(self):
        xml_str = b"""<?xml version="1.0" encoding="UTF-8"?>
        <EXPERIMENT_PACKAGE_SET>
          <EXPERIMENT_PACKAGE>
            <EXPERIMENT>
              <IDENTIFIERS><PRIMARY_ID>DRX000001</PRIMARY_ID></IDENTIFIERS>
              <DESIGN>
                <LIBRARY_DESCRIPTOR>
                  <LIBRARY_LAYOUT><SINGLE/></LIBRARY_LAYOUT>
                </LIBRARY_DESCRIPTOR>
              </DESIGN>
            </EXPERIMENT>
            <SUBMISSION>
              <IDENTIFIERS><PRIMARY_ID>DRA000001</PRIMARY_ID></IDENTIFIERS>
            </SUBMISSION>
            <STUDY>
              <DESCRIPTOR><STUDY_TITLE>Test Study</STUDY_TITLE></DESCRIPTOR>
            </STUDY>
            <SAMPLE>
              <IDENTIFIERS><PRIMARY_ID>DRS000001</PRIMARY_ID></IDENTIFIERS>
              <SAMPLE_NAME>
                <SCIENTIFIC_NAME>Test species</SCIENTIFIC_NAME>
                <TAXON_ID>1234</TAXON_ID>
              </SAMPLE_NAME>
            </SAMPLE>
            <RUN_SET>
              <RUN accession="DRR000001" total_spots="10" total_bases="100" size="1000">
                <IDENTIFIERS><PRIMARY_ID>DRR000001</PRIMARY_ID></IDENTIFIERS>
              </RUN>
            </RUN_SET>
          </EXPERIMENT_PACKAGE>
        </EXPERIMENT_PACKAGE_SET>"""
        root = ET.fromstring(xml_str)
        tree = ET.ElementTree(root)
        m = Metadata.from_xml(tree)

        assert m.df.loc[0, 'ENA_SRA_Link'] == ''
        assert m.df.loc[0, 'DDBJ_SRA_Link'] == (
            'https://ddbj.nig.ac.jp/public/ddbj_database/dra/sra/ByExp/sra/DRX/'
            'DRX000/DRX000001/DRR000001/DRR000001.sra'
        )

    def test_parse_empty_xml(self):
        xml_str = b"""<EXPERIMENT_PACKAGE_SET></EXPERIMENT_PACKAGE_SET>"""
        root = ET.fromstring(xml_str)
        tree = ET.ElementTree(root)
        m = Metadata.from_xml(tree)
        assert m.df.shape[0] == 0


class TestMetadataRemoveSpecialchars:
    def test_removes_special_characters(self):
        data = {
            'scientific_name': ['Homo\rsapiens'],
            'sample_group': ['brain\ntissue'],
            'tissue': ["brain'tissue"],
            'run': ['SRR|001'],
            'exclusion': ['no'],
        }
        df = pandas.DataFrame(data)
        m = Metadata.from_DataFrame(df)
        m.remove_specialchars()
        assert '\r' not in m.df.loc[0, 'scientific_name']
        assert '\n' not in m.df.loc[0, 'sample_group']
        assert "'" not in m.df.loc[0, 'tissue']
        assert '|' not in m.df.loc[0, 'run']


class TestMetadataPivot:
    def test_pivot_shape(self, sample_metadata):
        m = sample_metadata
        pivot = m.pivot(qualified_only=False, sampled_only=False)
        assert isinstance(pivot, pandas.DataFrame)
        # 2 species x 2 sample_groups
        assert pivot.shape[0] == 2
        assert pivot.shape[1] == 2

    def test_pivot_qualified_only(self, sample_metadata):
        m = sample_metadata
        m.df.loc[0, 'is_qualified'] = 'no'
        pivot = m.pivot(qualified_only=True, sampled_only=False)
        assert isinstance(pivot, pandas.DataFrame)


class TestMetadataLabelSampledData:
    def test_labels_samples(self, sample_metadata):
        m = sample_metadata
        m.label_sampled_data(max_sample=2)
        assert 'is_sampled' in m.df.columns
        assert 'is_qualified' in m.df.columns
        sampled = m.df.loc[m.df['is_sampled'] == 'yes']
        assert len(sampled) > 0

    def test_max_sample_respected(self):
        # Create 20 samples in one group
        n = 20
        data = {
            'scientific_name': ['Sp1'] * n,
            'sample_group': ['brain'] * n,
            'bioproject': [f'PRJ{i}' for i in range(n)],
            'biosample': [f'S{i}' for i in range(n)],
            'run': [f'R{i}' for i in range(n)],
            'exclusion': ['no'] * n,
        }
        df = pandas.DataFrame(data)
        m = Metadata.from_DataFrame(df)
        m.label_sampled_data(max_sample=5)
        sampled = m.df.loc[
            (m.df['is_sampled'] == 'yes') & (m.df['exclusion'] == 'no')
        ]
        assert len(sampled) == 5

    def test_normalizes_exclusion_before_qualification(self):
        data = {
            'scientific_name': ['Sp1', 'Sp1'],
            'sample_group': ['brain', 'brain'],
            'bioproject': ['PRJ1', 'PRJ2'],
            'biosample': ['S1', 'S2'],
            'run': ['R1', 'R2'],
            'exclusion': [' NO ', 'No'],
        }
        m = Metadata.from_DataFrame(pandas.DataFrame(data))
        m.label_sampled_data(max_sample=10)
        assert (m.df['is_qualified'] == 'yes').all()
        assert set(m.df['exclusion'].tolist()) == {'no'}


class TestMaximizeBioprojSampling:
    def test_sampling_across_bioprojects(self):
        n = 12
        data = {
            'scientific_name': ['Sp1'] * n,
            'sample_group': ['brain'] * n,
            'bioproject': ['PRJ1'] * 4 + ['PRJ2'] * 4 + ['PRJ3'] * 4,
            'biosample': [f'S{i}' for i in range(n)],
            'run': [f'R{i}' for i in range(n)],
            'exclusion': ['no'] * n,
            'is_sampled': ['no'] * n,
        }
        df = pandas.DataFrame(data)
        m = Metadata()
        result = m._maximize_bioproject_sampling(df, target_n=6)
        sampled = result.loc[result['is_sampled'] == 'yes']
        assert len(sampled) == 6
        # Should have samples from multiple bioprojects
        assert len(sampled['bioproject'].unique()) > 1

    def test_marks_all_eligible_when_eligible_le_target(self):
        data = {
            'scientific_name': ['Sp1'] * 8,
            'sample_group': ['brain'] * 8,
            'bioproject': ['PRJ1'] * 4 + ['PRJ2'] * 4,
            'biosample': [f'S{i}' for i in range(8)],
            'run': [f'R{i}' for i in range(8)],
            'exclusion': ['no', 'no', 'excluded', 'excluded', 'excluded', 'excluded', 'excluded', 'excluded'],
            'is_sampled': ['no'] * 8,
        }
        df = pandas.DataFrame(data)
        m = Metadata()
        result = m._maximize_bioproject_sampling(df, target_n=5)
        sampled = result.loc[(result['is_sampled'] == 'yes') & (result['exclusion'] == 'no')]
        assert len(sampled) == 2

    def test_respects_existing_selected_rows(self):
        data = {
            'scientific_name': ['Sp1'] * 6,
            'sample_group': ['brain'] * 6,
            'bioproject': ['PRJ1', 'PRJ1', 'PRJ2', 'PRJ2', 'PRJ3', 'PRJ3'],
            'biosample': [f'S{i}' for i in range(6)],
            'run': [f'R{i}' for i in range(6)],
            'exclusion': ['no'] * 6,
            'is_sampled': ['yes', 'yes', 'no', 'no', 'no', 'no'],
        }
        df = pandas.DataFrame(data)
        m = Metadata()
        result = m._maximize_bioproject_sampling(df, target_n=3)
        sampled = result.loc[(result['is_sampled'] == 'yes') & (result['exclusion'] == 'no')]
        assert len(sampled) == 3

    def test_handles_case_and_whitespace_in_exclusion_flags(self):
        data = {
            'scientific_name': ['Sp1'] * 4,
            'sample_group': ['brain'] * 4,
            'bioproject': ['PRJ1', 'PRJ1', 'PRJ2', 'PRJ2'],
            'biosample': [f'S{i}' for i in range(4)],
            'run': [f'R{i}' for i in range(4)],
            'exclusion': [' NO ', 'No', 'excluded', ' low_nspots '],
            'is_sampled': ['no'] * 4,
        }
        df = pandas.DataFrame(data)
        m = Metadata()
        result = m._maximize_bioproject_sampling(df, target_n=10)
        sampled = result.loc[result['is_sampled'] == 'yes', 'run'].tolist()
        assert sampled == ['R0', 'R1']

    def test_raises_when_exclusion_column_is_missing(self):
        df = pandas.DataFrame({
            'scientific_name': ['Sp1'],
            'sample_group': ['brain'],
            'bioproject': ['PRJ1'],
            'biosample': ['S1'],
            'run': ['R1'],
            'is_sampled': ['no'],
        })
        m = Metadata()
        with pytest.raises(ValueError, match='exclusion'):
            m._maximize_bioproject_sampling(df, target_n=1)


# ---------------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------------

class TestCleanupTmpAmalgkitFiles:
    def test_removes_matching_files_and_directories(self, tmp_path):
        (tmp_path / 'tmp.amalgkit.file1').write_text('x')
        tmp_dir = tmp_path / 'tmp.amalgkit.dir1'
        tmp_dir.mkdir()
        (tmp_dir / 'inner.txt').write_text('y')
        (tmp_path / 'keep.txt').write_text('z')

        cleanup_tmp_amalgkit_files(work_dir=str(tmp_path))

        assert not (tmp_path / 'tmp.amalgkit.file1').exists()
        assert not (tmp_path / 'tmp.amalgkit.dir1').exists()
        assert (tmp_path / 'keep.txt').exists()

    def test_ignores_file_not_found_during_cleanup(self, tmp_path, monkeypatch):
        disappearing = tmp_path / 'tmp.amalgkit.disappearing'
        disappearing.write_text('x')
        stable = tmp_path / 'tmp.amalgkit.stable'
        stable.write_text('y')
        real_remove = os.remove

        def flaky_remove(path):
            if os.path.realpath(path) == os.path.realpath(str(disappearing)):
                if os.path.exists(path):
                    real_remove(path)
                raise FileNotFoundError(path)
            return real_remove(path)

        monkeypatch.setattr('amalgkit.util.os.remove', flaky_remove)

        cleanup_tmp_amalgkit_files(work_dir=str(tmp_path))

        assert not disappearing.exists()
        assert not stable.exists()

class TestGetSraStat:
    def test_basic_stats(self, sample_metadata):
        stat = get_sra_stat('SRR001', sample_metadata)
        assert stat['sra_id'] == 'SRR001'
        assert stat['layout'] == 'paired'
        assert stat['total_spot'] == 10000000
        assert stat['spot_length'] == 200

    def test_inferred_spot_length(self, sample_metadata):
        m = sample_metadata
        m.df.loc[m.df['run'] == 'SRR001', 'spot_length'] = 0
        stat = get_sra_stat('SRR001', m)
        # total_bases / total_spots = 2000000000 / 10000000 = 200
        assert stat['spot_length'] == 200

    def test_num_read_per_sra(self, sample_metadata):
        stat = get_sra_stat('SRR001', sample_metadata, num_bp_per_sra=1000000)
        # num_bp_per_sra / spot_length = 1000000 / 200 = 5000
        assert stat['num_read_per_sra'] == 5000

    def test_prefers_layout_amalgkit_when_available(self, sample_metadata):
        m = sample_metadata
        m.df.loc[m.df['run'] == 'SRR001', 'layout_amalgkit'] = 'single'
        stat = get_sra_stat('SRR001', m)
        assert stat['layout'] == 'single'

    def test_normalizes_layout_case_and_whitespace(self, sample_metadata):
        m = sample_metadata
        m.df.loc[m.df['run'] == 'SRR001', 'lib_layout'] = '  PAIRED  '
        stat = get_sra_stat('SRR001', m)
        assert stat['layout'] == 'paired'

    def test_raises_when_layout_is_invalid(self):
        m = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'lib_layout': ['unknown'],
            'total_spots': [10],
            'spot_length': [100],
            'total_bases': [1000],
            'exclusion': ['no'],
        }))
        with pytest.raises(ValueError, match='Unsupported lib_layout'):
            get_sra_stat('SRR001', m)

    def test_duplicate_run_raises(self):
        m = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001', 'SRR001'],
            'lib_layout': ['single', 'single'],
            'total_spots': [10, 10],
            'spot_length': [100, 100],
            'total_bases': [1000, 1000],
            'exclusion': ['no', 'no'],
        }))
        with pytest.raises(AssertionError, match='multiple metadata rows'):
            get_sra_stat('SRR001', m)

    def test_missing_run_raises(self):
        m = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'lib_layout': ['single'],
            'total_spots': [10],
            'spot_length': [100],
            'total_bases': [1000],
            'exclusion': ['no'],
        }))
        with pytest.raises(AssertionError, match='SRA ID not found'):
            get_sra_stat('SRR999', m)

    def test_requires_run_column(self):
        m = Metadata.from_DataFrame(pandas.DataFrame({
            'lib_layout': ['single'],
            'total_spots': [10],
            'spot_length': [100],
            'total_bases': [1000],
            'exclusion': ['no'],
        }))
        m.df = m.df.drop(columns=['run'])
        with pytest.raises(ValueError, match='Missing required metadata column\\(s\\) for get_sra_stat: run'):
            get_sra_stat('SRR001', m)

    def test_requires_lib_layout_column(self):
        m = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'total_spots': [10],
            'spot_length': [100],
            'total_bases': [1000],
            'exclusion': ['no'],
        }))
        m.df = m.df.drop(columns=['lib_layout'])

        with pytest.raises(ValueError, match='Missing required metadata column\\(s\\) for get_sra_stat: lib_layout'):
            get_sra_stat('SRR001', m)

    def test_raises_when_total_spots_is_zero(self):
        m = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'lib_layout': ['single'],
            'total_spots': [0],
            'spot_length': [100],
            'total_bases': [1000],
            'exclusion': ['no'],
        }))
        with pytest.raises(ValueError, match='total_spots must be > 0'):
            get_sra_stat('SRR001', m)

    def test_raises_when_spot_length_cannot_be_inferred(self):
        m = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'lib_layout': ['single'],
            'total_spots': [10],
            'spot_length': [0],
            'total_bases': [0],
            'exclusion': ['no'],
        }))
        with pytest.raises(ValueError, match='spot_length cannot be inferred'):
            get_sra_stat('SRR001', m)


class TestCheckOrthologParameterCompatibility:
    def test_both_none_raises(self):
        class Args:
            orthogroup_table = None
            dir_busco = None
        with pytest.raises(ValueError, match="One of"):
            check_ortholog_parameter_compatibility(Args())

    def test_both_set_raises(self):
        class Args:
            orthogroup_table = 'table.tsv'
            dir_busco = '/path/to/busco'
        with pytest.raises(ValueError, match="Only one"):
            check_ortholog_parameter_compatibility(Args())

    def test_only_orthogroup_ok(self):
        class Args:
            orthogroup_table = 'table.tsv'
            dir_busco = None
        assert check_ortholog_parameter_compatibility(Args()) == ('table.tsv', None)

    def test_only_busco_ok(self):
        class Args:
            orthogroup_table = None
            dir_busco = '/path/to/busco'
        assert check_ortholog_parameter_compatibility(Args()) == (None, '/path/to/busco')

    def test_blank_orthogroup_table_is_treated_as_none(self):
        class Args:
            orthogroup_table = '   '
            dir_busco = None
        with pytest.raises(ValueError, match='One of'):
            check_ortholog_parameter_compatibility(Args())

    def test_does_not_mutate_args(self):
        class Args:
            orthogroup_table = ' table.tsv '
            dir_busco = None

        args = Args()
        normalized = check_ortholog_parameter_compatibility(args)

        assert normalized == ('table.tsv', None)
        assert args.orthogroup_table == ' table.tsv '
        assert args.dir_busco is None


class TestOrthogroup2Genecount:
    def test_basic_transformation(self, tmp_path):
        ortho_file = tmp_path / 'orthogroups.tsv'
        ortho_file.write_text(
            'busco_id\tSpecies_A\tSpecies_B\n'
            'OG0001\tgene1\tgene2,gene3\n'
            'OG0002\tgene4,gene5,gene6\t-\n'
            'OG0003\t\tgene7\n'
        )
        out_file = tmp_path / 'genecount.tsv'
        orthogroup2genecount(
            str(ortho_file), str(out_file),
            spp=['Species_A', 'Species_B']
        )
        result = pandas.read_csv(str(out_file), sep='\t')
        assert list(result.columns) == ['orthogroup_id', 'Species_A', 'Species_B']
        assert result.loc[0, 'Species_A'] == 1  # gene1
        assert result.loc[0, 'Species_B'] == 2  # gene2,gene3
        assert result.loc[1, 'Species_A'] == 3  # gene4,gene5,gene6
        assert result.loc[1, 'Species_B'] == 0  # '-' -> empty -> 0
        assert result.loc[2, 'Species_A'] == 0  # empty
        assert result.loc[2, 'Species_B'] == 1  # gene7

    def test_handles_blank_hyphen_and_multi_comma_values(self, tmp_path):
        ortho_file = tmp_path / 'orthogroups.tsv'
        ortho_file.write_text(
            'busco_id\tSpecies_A\tSpecies_B\n'
            'OG1000\t-\tgene1\n'
            'OG1001\t\tgene2,gene3,gene4\n'
            'OG1002\tgene5,gene6\t-\n'
        )
        out_file = tmp_path / 'genecount.tsv'
        orthogroup2genecount(
            str(ortho_file), str(out_file),
            spp=['Species_A', 'Species_B']
        )
        result = pandas.read_csv(str(out_file), sep='\t')
        assert result.loc[0, 'Species_A'] == 0
        assert result.loc[0, 'Species_B'] == 1
        assert result.loc[1, 'Species_A'] == 0
        assert result.loc[1, 'Species_B'] == 3
        assert result.loc[2, 'Species_A'] == 2
        assert result.loc[2, 'Species_B'] == 0

    def test_raises_when_busco_id_column_missing(self, tmp_path):
        ortho_file = tmp_path / 'orthogroups.tsv'
        ortho_file.write_text(
            'orthogroup\tSpecies_A\n'
            'OG0001\tgene1\n'
        )
        out_file = tmp_path / 'genecount.tsv'

        with pytest.raises(ValueError, match='Column \"busco_id\" is required'):
            orthogroup2genecount(
                str(ortho_file),
                str(out_file),
                spp=['Species_A'],
            )

    def test_raises_when_species_column_missing(self, tmp_path):
        ortho_file = tmp_path / 'orthogroups.tsv'
        ortho_file.write_text(
            'busco_id\tSpecies_A\n'
            'OG0001\tgene1\n'
        )
        out_file = tmp_path / 'genecount.tsv'

        with pytest.raises(ValueError, match='Species column\\(s\\) not found in orthogroup table'):
            orthogroup2genecount(
                str(ortho_file),
                str(out_file),
                spp=['Species_A', 'Species_B'],
            )


# ---------------------------------------------------------------------------
# label_sampled_data: empty sample_group handling (wiki: select)
# ---------------------------------------------------------------------------

class TestMetadataLabelSampledDataEdgeCases:
    def test_empty_sample_group_marked_unqualified(self):
        """Empty sample_group rows remain unqualified and are not auto-assigned an exclusion label."""
        data = {
            'scientific_name': ['Sp1', 'Sp1'],
            'sample_group': ['brain', ''],
            'bioproject': ['PRJ1', 'PRJ2'],
            'biosample': ['S1', 'S2'],
            'run': ['R1', 'R2'],
            'exclusion': ['no', 'no'],
        }
        m = Metadata.from_DataFrame(pandas.DataFrame(data))
        m.label_sampled_data(max_sample=10)
        empty_sg = m.df.loc[m.df['run'] == 'R2']
        assert empty_sg['exclusion'].values[0] == 'no'
        assert empty_sg['is_qualified'].values[0] == 'no'

    def test_no_futurewarning_when_is_qualified_starts_float(self):
        """Assigning qualification flags should be dtype-safe even from float columns."""
        data = {
            'scientific_name': ['Sp1', 'Sp1'],
            'sample_group': ['brain', ''],
            'bioproject': ['PRJ1', 'PRJ2'],
            'biosample': ['S1', 'S2'],
            'run': ['R1', 'R2'],
            'exclusion': ['no', 'no'],
            # Simulate pandas float inference from all-missing text column.
            'is_qualified': [numpy.nan, numpy.nan],
        }
        m = Metadata.from_DataFrame(pandas.DataFrame(data))
        with warnings.catch_warnings(record=True) as captured:
            warnings.simplefilter('always')
            m.label_sampled_data(max_sample=10)
        future_warnings = [w for w in captured if issubclass(w.category, FutureWarning)]
        assert len(future_warnings) == 0

    def test_preserves_chained_assignment_option(self):
        data = {
            'scientific_name': ['Sp1', 'Sp1'],
            'sample_group': ['brain', 'brain'],
            'bioproject': ['PRJ1', 'PRJ2'],
            'biosample': ['S1', 'S2'],
            'run': ['R1', 'R2'],
            'exclusion': ['no', 'no'],
        }
        m = Metadata.from_DataFrame(pandas.DataFrame(data))
        previous_mode = pandas.get_option('mode.chained_assignment')
        pandas.set_option('mode.chained_assignment', 'raise')
        try:
            m.label_sampled_data(max_sample=10)
            assert pandas.get_option('mode.chained_assignment') == 'raise'
        finally:
            pandas.set_option('mode.chained_assignment', previous_mode)

    def test_restores_chained_assignment_option_after_error(self, monkeypatch):
        data = {
            'scientific_name': ['Sp1', 'Sp1'],
            'sample_group': ['brain', 'brain'],
            'bioproject': ['PRJ1', 'PRJ2'],
            'biosample': ['S1', 'S2'],
            'run': ['R1', 'R2'],
            'exclusion': ['no', 'no'],
        }
        m = Metadata.from_DataFrame(pandas.DataFrame(data))
        previous_mode = pandas.get_option('mode.chained_assignment')
        pandas.set_option('mode.chained_assignment', 'raise')
        monkeypatch.setattr(
            m,
            '_maximize_bioproject_sampling',
            lambda *args, **kwargs: (_ for _ in ()).throw(RuntimeError('forced failure')),
        )
        try:
            with pytest.raises(RuntimeError, match='forced failure'):
                m.label_sampled_data(max_sample=10)
            assert pandas.get_option('mode.chained_assignment') == 'raise'
        finally:
            pandas.set_option('mode.chained_assignment', previous_mode)


# ---------------------------------------------------------------------------
# Metadata.from_xml: SAMPLE_ATTRIBUTES extraction
# ---------------------------------------------------------------------------

class TestMetadataFromXmlAttributes:
    def test_sample_attributes_extracted(self):
        """Wiki: XML SAMPLE_ATTRIBUTES are parsed into extra columns."""
        xml_str = b"""<?xml version="1.0" encoding="UTF-8"?>
        <EXPERIMENT_PACKAGE_SET>
          <EXPERIMENT_PACKAGE>
            <EXPERIMENT alias="test">
              <IDENTIFIERS><PRIMARY_ID>SRX000001</PRIMARY_ID></IDENTIFIERS>
              <TITLE>Test</TITLE>
              <STUDY_REF>
                <IDENTIFIERS><PRIMARY_ID>SRP000001</PRIMARY_ID></IDENTIFIERS>
              </STUDY_REF>
              <DESIGN>
                <DESIGN_DESCRIPTION/>
                <LIBRARY_DESCRIPTOR>
                  <LIBRARY_NAME/>
                  <LIBRARY_STRATEGY>RNA-Seq</LIBRARY_STRATEGY>
                  <LIBRARY_SOURCE>TRANSCRIPTOMIC</LIBRARY_SOURCE>
                  <LIBRARY_SELECTION>cDNA</LIBRARY_SELECTION>
                  <LIBRARY_LAYOUT><SINGLE/></LIBRARY_LAYOUT>
                </LIBRARY_DESCRIPTOR>
              </DESIGN>
              <PLATFORM><ILLUMINA><INSTRUMENT_MODEL>HiSeq 2500</INSTRUMENT_MODEL></ILLUMINA></PLATFORM>
            </EXPERIMENT>
            <SUBMISSION accession="SRA1" lab_name="" center_name="">
              <IDENTIFIERS><PRIMARY_ID>SRA1</PRIMARY_ID><SUBMITTER_ID/></IDENTIFIERS>
            </SUBMISSION>
            <STUDY><DESCRIPTOR><STUDY_TITLE>Study</STUDY_TITLE></DESCRIPTOR></STUDY>
            <SAMPLE>
              <IDENTIFIERS><PRIMARY_ID>SRS1</PRIMARY_ID></IDENTIFIERS>
              <TITLE/>
              <SAMPLE_NAME>
                <SCIENTIFIC_NAME>Mus musculus</SCIENTIFIC_NAME>
                <TAXON_ID>10090</TAXON_ID>
              </SAMPLE_NAME>
              <DESCRIPTION/>
              <SAMPLE_ATTRIBUTES>
                <SAMPLE_ATTRIBUTE><TAG>tissue</TAG><VALUE>liver</VALUE></SAMPLE_ATTRIBUTE>
                <SAMPLE_ATTRIBUTE><TAG>sex</TAG><VALUE>female</VALUE></SAMPLE_ATTRIBUTE>
                <SAMPLE_ATTRIBUTE><TAG>age</TAG><VALUE>8 weeks</VALUE></SAMPLE_ATTRIBUTE>
              </SAMPLE_ATTRIBUTES>
            </SAMPLE>
            <RUN_SET>
              <RUN accession="SRR1" total_spots="5000000" total_bases="500000000" size="100000000" published="2021-01-01">
                <IDENTIFIERS><PRIMARY_ID>SRR1</PRIMARY_ID></IDENTIFIERS>
              </RUN>
            </RUN_SET>
            <Pool>
              <Member>
                <EXTERNAL_ID namespace="BioProject">PRJNA1</EXTERNAL_ID>
                <EXTERNAL_ID namespace="BioSample">SAMN1</EXTERNAL_ID>
              </Member>
            </Pool>
          </EXPERIMENT_PACKAGE>
        </EXPERIMENT_PACKAGE_SET>"""
        root = ET.fromstring(xml_str)
        tree = ET.ElementTree(root)
        m = Metadata.from_xml(tree)
        assert m.df.loc[0, 'tissue'] == 'liver'
        assert m.df.loc[0, 'lib_layout'] == 'single'
        assert m.df.loc[0, 'scientific_name'] == 'Mus musculus'


# ---------------------------------------------------------------------------
# Metadata.reorder: extra columns preserved
# ---------------------------------------------------------------------------

class TestMetadataReorderExtraCols:
    def test_extra_columns_preserved(self):
        """Extra columns not in column_names should be kept after reorder."""
        data = {
            'scientific_name': ['Sp1'],
            'run': ['R1'],
            'exclusion': ['no'],
            'my_custom_column': ['custom_value'],
        }
        m = Metadata.from_DataFrame(pandas.DataFrame(data))
        assert 'my_custom_column' in m.df.columns

    def test_reorder_preserves_data(self, sample_metadata):
        """Reorder should not lose any rows."""
        original_rows = sample_metadata.df.shape[0]
        sample_metadata.reorder()
        assert sample_metadata.df.shape[0] == original_rows


# ---------------------------------------------------------------------------
# Metadata.pivot: sampled_only filter
# ---------------------------------------------------------------------------

class TestMetadataPivotSampledOnly:
    def test_pivot_sampled_only(self, sample_metadata):
        m = sample_metadata
        m.label_sampled_data(max_sample=2)
        pivot = m.pivot(qualified_only=True, sampled_only=True)
        assert isinstance(pivot, pandas.DataFrame)

    def test_pivot_n_sp_cutoff(self, sample_metadata):
        """n_sp_cutoff filters columns with fewer species than cutoff."""
        m = sample_metadata
        pivot = m.pivot(n_sp_cutoff=3, qualified_only=False, sampled_only=False)
        # With cutoff=3, columns where fewer than 3 species appear are dropped
        assert isinstance(pivot, pandas.DataFrame)


# ---------------------------------------------------------------------------
# load_metadata (loads metadata from file)
# ---------------------------------------------------------------------------

class TestLoadMetadata:
    def test_load_from_explicit_path(self, tmp_path, sample_metadata):
        """Loads metadata from an explicit file path."""
        path = tmp_path / 'metadata.tsv'
        sample_metadata.df.to_csv(str(path), sep='\t', index=False)
        class Args:
            metadata = str(path)
            out_dir = str(tmp_path)
        m = load_metadata(Args())
        assert isinstance(m, Metadata)
        assert m.df.shape[0] == 5

    def test_load_metadata_requests_utf8_encoding(self, tmp_path, sample_metadata, monkeypatch):
        path = tmp_path / 'metadata.tsv'
        sample_metadata.df.to_csv(str(path), sep='\t', index=False, encoding='utf-8')
        observed = {}
        original_read_csv = pandas.read_csv

        def capture_read_csv(*args, **kwargs):
            observed['encoding'] = kwargs.get('encoding')
            return original_read_csv(*args, **kwargs)

        monkeypatch.setattr(pandas, 'read_csv', capture_read_csv)

        class Args:
            metadata = str(path)
            out_dir = str(tmp_path)

        load_metadata(Args())

        assert observed['encoding'] == 'utf-8'

    def test_load_from_inferred_path(self, tmp_path, sample_metadata):
        """When metadata='inferred', loads from out_dir/metadata/metadata.tsv."""
        meta_dir = tmp_path / 'metadata'
        meta_dir.mkdir()
        path = meta_dir / 'metadata.tsv'
        sample_metadata.df.to_csv(str(path), sep='\t', index=False)
        class Args:
            metadata = 'inferred'
            out_dir = str(tmp_path)
        m = load_metadata(Args())
        assert isinstance(m, Metadata)
        assert m.df.shape[0] == 5

    def test_raises_when_metadata_file_missing(self, tmp_path):
        class Args:
            metadata = str(tmp_path / 'missing.tsv')
            out_dir = str(tmp_path)
        with pytest.raises(FileNotFoundError, match='Metadata file not found'):
            load_metadata(Args())

    def test_raises_when_metadata_path_is_directory(self, tmp_path):
        metadata_dir = tmp_path / 'metadata_path'
        metadata_dir.mkdir()

        class Args:
            metadata = str(metadata_dir)
            out_dir = str(tmp_path)

        with pytest.raises(IsADirectoryError, match='Metadata path exists but is not a file'):
            load_metadata(Args())

    def test_batch_mode_treats_missing_is_sampled_as_no(self, tmp_path):
        path = tmp_path / 'metadata.tsv'
        pandas.DataFrame({
            'run': ['R1', 'R2'],
            'scientific_name': ['Sp1', 'Sp1'],
            'is_sampled': ['yes', None],
            'exclusion': ['no', 'no'],
        }).to_csv(str(path), sep='\t', index=False)

        class Args:
            metadata = str(path)
            out_dir = str(tmp_path)
            batch = 1

        m = load_metadata(Args())
        assert m.df.shape[0] == 1
        assert m.df.iloc[0]['run'] == 'R1'

    def test_batch_mode_raises_for_invalid_is_sampled_value(self, tmp_path):
        path = tmp_path / 'metadata.tsv'
        pandas.DataFrame({
            'run': ['R1', 'R2'],
            'scientific_name': ['Sp1', 'Sp1'],
            'is_sampled': ['yes', 'maybe'],
            'exclusion': ['no', 'no'],
        }).to_csv(str(path), sep='\t', index=False)

        class Args:
            metadata = str(path)
            out_dir = str(tmp_path)
            batch = 1

        with pytest.raises(ValueError, match='is_sampled'):
            load_metadata(Args())

    def test_batch_mode_rejects_nonpositive_batch(self, tmp_path):
        path = tmp_path / 'metadata.tsv'
        pandas.DataFrame({
            'run': ['R1'],
            'scientific_name': ['Sp1'],
            'is_sampled': ['yes'],
            'exclusion': ['no'],
        }).to_csv(str(path), sep='\t', index=False)

        class Args:
            metadata = str(path)
            out_dir = str(tmp_path)
            batch = 0

        with pytest.raises(ValueError, match='--batch must be >= 1'):
            load_metadata(Args())

    def test_batch_mode_exits_nonzero_when_no_sampled_rows(self, tmp_path):
        path = tmp_path / 'metadata.tsv'
        pandas.DataFrame({
            'run': ['R1', 'R2'],
            'scientific_name': ['Sp1', 'Sp1'],
            'is_sampled': ['no', 'no'],
            'exclusion': ['no', 'no'],
        }).to_csv(str(path), sep='\t', index=False)

        class Args:
            metadata = str(path)
            out_dir = str(tmp_path)
            batch = 1

        with pytest.raises(ValueError, match='No sample is "sampled"'):
            load_metadata(Args())

    def test_batch_mode_uses_explicit_run_scope(self, tmp_path):
        path = tmp_path / 'metadata.tsv'
        pandas.DataFrame({
            'run': ['R1', 'R2'],
            'scientific_name': ['Sp1', 'Sp1'],
            'is_sampled': ['yes', 'no'],
            'exclusion': ['no', 'no'],
        }).to_csv(str(path), sep='\t', index=False)

        class Args:
            metadata = str(path)
            out_dir = str(tmp_path)
            batch = 1

        m = load_metadata(Args(), batch_scope='run')
        assert m.df.shape[0] == 1
        assert m.df.iloc[0]['run'] == 'R1'

    def test_species_batch_ignores_missing_species_names(self, tmp_path):
        path = tmp_path / 'metadata.tsv'
        pandas.DataFrame({
            'run': ['R1', 'R2', 'R3'],
            'scientific_name': ['Sp2', '', 'Sp1'],
            'is_sampled': ['yes', 'yes', 'yes'],
            'exclusion': ['no', 'no', 'no'],
        }).to_csv(str(path), sep='\t', index=False)

        class Args:
            metadata = str(path)
            out_dir = str(tmp_path)
            batch = 1

        m = load_metadata(Args(), batch_scope='species')
        assert m.df.shape[0] == 1
        assert m.df.iloc[0]['scientific_name'] == 'Sp1'

    def test_species_batch_raises_when_no_valid_species(self, tmp_path):
        path = tmp_path / 'metadata.tsv'
        pandas.DataFrame({
            'run': ['R1', 'R2'],
            'scientific_name': ['', None],
            'is_sampled': ['yes', 'yes'],
            'exclusion': ['no', 'no'],
        }).to_csv(str(path), sep='\t', index=False)

        class Args:
            metadata = str(path)
            out_dir = str(tmp_path)
            batch = 1

        with pytest.raises(ValueError, match='No valid scientific_name'):
            load_metadata(Args(), batch_scope='species')

    def test_batch_mode_raises_amalgkit_exit_when_batch_too_large(self, tmp_path):
        path = tmp_path / 'metadata.tsv'
        pandas.DataFrame({
            'run': ['R1'],
            'scientific_name': ['Sp1'],
            'is_sampled': ['yes'],
            'exclusion': ['no'],
        }).to_csv(str(path), sep='\t', index=False)

        class Args:
            metadata = str(path)
            out_dir = str(tmp_path)
            batch = 2

        with pytest.raises(AmalgkitExit) as exc:
            load_metadata(Args())
        assert exc.value.exit_code == 0


# ---------------------------------------------------------------------------
# detect_layout_from_file (corrects layout based on actual files)
# ---------------------------------------------------------------------------
