import json
import os
import warnings
import warnings
import pytest
import pandas
import numpy
import xml.etree.ElementTree as ET
from types import SimpleNamespace

from amalgkit.util import (
    strtobool,
    parse_bool_flags,
    Metadata,
    read_config_file,
    get_sra_stat,
    check_ortholog_parameter_compatibility,
    orthogroup2genecount,
    check_config_dir,
    load_metadata,
    detect_layout_from_file,
    is_there_unpaired_file,
    get_newest_intermediate_file_extension,
    get_mapping_rate,
    get_getfastq_run_dir,
    generate_multisp_busco_table,
    cleanup_tmp_amalgkit_files,
    check_rscript,
    run_tasks_with_optional_threads,
    validate_positive_int_option,
    resolve_cpu_budget,
    resolve_thread_worker_allocation,
    resolve_worker_allocation,
    find_prefixed_entries,
    find_species_prefixed_entries,
    find_run_prefixed_entries,
)


# ---------------------------------------------------------------------------
# strtobool
# ---------------------------------------------------------------------------

class TestStrtobool:
    @pytest.mark.parametrize("val", ["y", "yes", "Yes", "YES", "t", "true", "True", "on", "1"])
    def test_true_values(self, val):
        assert strtobool(val) is True

    @pytest.mark.parametrize("val", ["n", "no", "No", "NO", "f", "false", "False", "off", "0"])
    def test_false_values(self, val):
        assert strtobool(val) is False

    def test_invalid_value(self):
        with pytest.raises(ValueError):
            strtobool("maybe")


class TestParseBoolFlags:
    def test_fills_missing_values_with_default(self):
        result = parse_bool_flags(['yes', None, ''], column_name='is_sampled', default='no')
        assert result.tolist() == [True, False, False]

    def test_raises_for_invalid_values(self):
        with pytest.raises(ValueError, match='invalid boolean flag'):
            parse_bool_flags(['yes', 'maybe'], column_name='is_sampled', default='no')

    def test_accepts_generator_input(self):
        values = (v for v in ['yes', 'no', None])
        result = parse_bool_flags(values, column_name='is_sampled', default='no')
        assert result.tolist() == [True, False, False]

class TestRunTasksWithOptionalThreads:
    def test_empty_tasks(self):
        results, failures = run_tasks_with_optional_threads([], lambda x: x, max_workers=4)
        assert results == {}
        assert failures == []

    def test_serial_collects_successes_and_failures(self):
        def worker(x):
            if x == 2:
                raise RuntimeError('boom')
            return x * 10

        results, failures = run_tasks_with_optional_threads([1, 2, 3], worker, max_workers=1)

        assert results == {1: 10, 3: 30}
        assert len(failures) == 1
        assert failures[0][0] == 2
        assert isinstance(failures[0][1], RuntimeError)

    def test_parallel_collects_successes_and_failures(self):
        def worker(x):
            if x == 3:
                raise ValueError('bad')
            return x + 1

        results, failures = run_tasks_with_optional_threads([1, 2, 3, 4], worker, max_workers=3)

        assert results == {1: 2, 2: 3, 4: 5}
        assert len(failures) == 1
        assert failures[0][0] == 3
        assert isinstance(failures[0][1], ValueError)

    def test_parallel_converts_system_exit_to_failure(self):
        def worker(x):
            if x == 2:
                raise SystemExit(2)
            return x

        results, failures = run_tasks_with_optional_threads([1, 2], worker, max_workers=2)

        assert results == {1: 1}
        assert len(failures) == 1
        assert failures[0][0] == 2
        assert isinstance(failures[0][1], RuntimeError)
        assert 'Task requested exit with code 2.' in str(failures[0][1])

    def test_accepts_string_max_workers(self):
        results, failures = run_tasks_with_optional_threads([1, 2], lambda x: x * 2, max_workers='2')
        assert results == {1: 2, 2: 4}
        assert failures == []

    def test_auto_max_workers_falls_back_to_serial(self):
        results, failures = run_tasks_with_optional_threads([1, 2], lambda x: x + 1, max_workers='auto')
        assert results == {1: 2, 2: 3}
        assert failures == []

    def test_invalid_max_workers_raises_clear_error(self):
        with pytest.raises(ValueError, match='max_workers must be an integer'):
            run_tasks_with_optional_threads([1], lambda x: x, max_workers='two')

class TestValidatePositiveIntOption:
    def test_accepts_positive(self):
        assert validate_positive_int_option(3, 'internal_jobs') == 3
        assert validate_positive_int_option('2', 'internal_jobs') == 2

    def test_rejects_nonpositive(self):
        with pytest.raises(ValueError, match='--internal_jobs must be > 0'):
            validate_positive_int_option(0, 'internal_jobs')


class TestCpuBudgetHelpers:
    def test_resolve_cpu_budget_auto_uses_os_cpu_count(self, monkeypatch):
        monkeypatch.setattr('amalgkit.util.os.cpu_count', lambda: 12)
        assert resolve_cpu_budget(0) == 12

    def test_resolve_cpu_budget_auto_fallbacks_to_one(self, monkeypatch):
        monkeypatch.setattr('amalgkit.util.os.cpu_count', lambda: None)
        assert resolve_cpu_budget(0) == 1

    def test_resolve_cpu_budget_rejects_negative(self):
        with pytest.raises(ValueError, match='--internal_cpu_budget must be >= 0'):
            resolve_cpu_budget(-1)

    def test_resolve_thread_worker_allocation_caps_workers(self):
        threads, workers, budget = resolve_thread_worker_allocation(
            requested_threads=4,
            requested_workers=4,
            internal_cpu_budget=8,
            worker_option_name='internal_jobs',
        )
        assert threads == 1
        assert workers == 4
        assert budget == 4

    def test_resolve_thread_worker_allocation_caps_threads(self):
        threads, workers, budget = resolve_thread_worker_allocation(
            requested_threads=16,
            requested_workers=2,
            internal_cpu_budget=8,
            worker_option_name='internal_jobs',
        )
        assert threads == 4
        assert workers == 2
        assert budget == 8

    def test_resolve_thread_worker_allocation_respects_total_core_budget(self):
        threads, workers, budget = resolve_thread_worker_allocation(
            requested_threads=8,
            requested_workers='auto',
            internal_cpu_budget=64,
            worker_option_name='internal_jobs',
        )
        assert threads == 1
        assert workers == 8
        assert budget == 8

    def test_resolve_worker_allocation_caps_workers(self):
        workers, budget = resolve_worker_allocation(
            requested_workers=10,
            internal_cpu_budget=3,
            worker_option_name='internal_jobs',
        )
        assert workers == 3
        assert budget == 3
        with pytest.raises(ValueError, match='--internal_jobs must be > 0'):
            validate_positive_int_option(-1, 'internal_jobs')

class TestFindPrefixedEntries:
    def test_sorted_list_entries(self):
        entries = ['Homo_sapiens.fa', 'Homo_sapiens.idx', 'Mus_musculus.fa']
        out = find_prefixed_entries(entries, 'Homo_sapiens')
        assert out == ['Homo_sapiens.fa', 'Homo_sapiens.idx']

    def test_set_entries_returns_sorted_output(self):
        entries = {'Homo_sapiens_b.idx', 'Mus_musculus.idx', 'Homo_sapiens_a.idx'}
        out = find_prefixed_entries(entries, 'Homo_sapiens')
        assert out == ['Homo_sapiens_a.idx', 'Homo_sapiens_b.idx']

    def test_no_match(self):
        entries = ['Arabidopsis_thaliana.fa']
        out = find_prefixed_entries(entries, 'Homo_sapiens')
        assert out == []

    def test_unsorted_list_entries(self):
        entries = ['Aardvark.fa', 'Homo_sapiens.fa', 'Mus_musculus.fa', 'Homo_sapiens.idx']
        out = find_prefixed_entries(entries, 'Homo_sapiens')
        assert out == ['Homo_sapiens.fa', 'Homo_sapiens.idx']


class TestFindSpeciesPrefixedEntries:
    def test_rejects_similar_species_prefix(self):
        entries = ['Homo_sapiens.fa', 'Homo_sapiens2.fa', 'Homo_sapiens_k31.idx']
        out = find_species_prefixed_entries(entries, 'Homo_sapiens')
        assert out == ['Homo_sapiens.fa', 'Homo_sapiens_k31.idx']


class TestFindRunPrefixedEntries:
    def test_rejects_similar_run_prefix(self):
        entries = ['SRR001.fastq.gz', 'SRR0010.fastq.gz', 'SRR001_1.fastq.gz']
        out = find_run_prefixed_entries(entries, 'SRR001')
        assert out == ['SRR001.fastq.gz', 'SRR001_1.fastq.gz']

    def test_rejects_hyphen_suffix_variants(self):
        entries = ['SRR001.fastq.gz', 'SRR001-legacy.fastq.gz', 'SRR001_1.fastq.gz']
        out = find_run_prefixed_entries(entries, 'SRR001')
        assert out == ['SRR001.fastq.gz', 'SRR001_1.fastq.gz']


# ---------------------------------------------------------------------------
# Metadata class
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
                     'lib_layout', 'total_spots', 'exclusion']:
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
        for col in Metadata.removed_metadata_columns:
            assert col not in m.df.columns

    def test_parse_empty_xml(self):
        xml_str = b"""<EXPERIMENT_PACKAGE_SET></EXPERIMENT_PACKAGE_SET>"""
        root = ET.fromstring(xml_str)
        tree = ET.ElementTree(root)
        m = Metadata.from_xml(tree)
        assert m.df.shape[0] == 0


class TestMetadataNspotCutoff:
    def test_marks_low_spots(self, sample_metadata):
        m = sample_metadata
        m.nspot_cutoff(1000000)
        # SRR003 has 100 spots, SRR005 has 200 spots - both below cutoff
        low = m.df.loc[m.df['run'].isin(['SRR003', 'SRR005']), 'exclusion']
        assert (low == 'low_nspots').all()
        # SRR001 has 10M spots - should remain 'no'
        high = m.df.loc[m.df['run'] == 'SRR001', 'exclusion'].values[0]
        assert high == 'no'


class TestMetadataMarkExcludeKeywords:
    def test_marks_matching_keyword(self, sample_metadata, tmp_config_dir):
        m = sample_metadata
        m.df.loc[0, 'sample_description'] = 'cancer tissue sample'
        m.mark_exclude_keywords(tmp_config_dir)
        assert m.df.loc[0, 'exclusion'] == 'disease'

    def test_no_false_positives(self, sample_metadata, tmp_config_dir):
        m = sample_metadata
        m.mark_exclude_keywords(tmp_config_dir)
        assert (m.df['exclusion'] == 'no').all()

    def test_strips_config_column_tokens(self, sample_metadata, tmp_path):
        (tmp_path / 'exclude_keyword.config').write_text('sample_description, tissue\tdisease\tcancer\n')
        m = sample_metadata
        m.df.loc[0, 'tissue'] = 'cancer tissue'
        m.mark_exclude_keywords(str(tmp_path))
        assert m.df.loc[0, 'exclusion'] == 'disease'

    def test_raises_for_unknown_config_column(self, sample_metadata, tmp_path):
        (tmp_path / 'exclude_keyword.config').write_text('unknown_column\tdisease\tcancer\n')
        m = sample_metadata
        with pytest.raises(ValueError, match='exclude_keyword.config were not found in metadata'):
            m.mark_exclude_keywords(str(tmp_path))

    def test_raises_for_invalid_keyword_regex(self, sample_metadata, tmp_path):
        (tmp_path / 'exclude_keyword.config').write_text('sample_description\tdisease\t[invalid\n')
        m = sample_metadata
        with pytest.raises(ValueError, match='Invalid regex pattern in exclude_keyword.config row 1'):
            m.mark_exclude_keywords(str(tmp_path))

    def test_ignores_empty_keyword_row_without_mass_exclusion(self, sample_metadata, tmp_path):
        (tmp_path / 'exclude_keyword.config').write_text('sample_description\tdisease\t\n')
        m = sample_metadata
        m.mark_exclude_keywords(str(tmp_path))
        assert (m.df['exclusion'] == 'no').all()

    def test_ignores_empty_reason_row_without_mass_exclusion(self, sample_metadata, tmp_path):
        (tmp_path / 'exclude_keyword.config').write_text('sample_description\t\tcancer\n')
        m = sample_metadata
        m.df.loc[0, 'sample_description'] = 'cancer tissue sample'
        m.mark_exclude_keywords(str(tmp_path))
        assert (m.df['exclusion'] == 'no').all()

    def test_raises_when_exclude_keyword_config_path_is_directory(self, sample_metadata, tmp_path):
        (tmp_path / 'exclude_keyword.config').mkdir()
        m = sample_metadata
        with pytest.raises(IsADirectoryError, match='Config path exists but is not a file'):
            m.mark_exclude_keywords(str(tmp_path))

    def test_raises_for_malformed_config_with_too_few_columns(self, sample_metadata, tmp_path):
        (tmp_path / 'exclude_keyword.config').write_text('sample_description_only\n')
        m = sample_metadata
        with pytest.raises(ValueError, match='exclude_keyword.config must contain at least 3'):
            m.mark_exclude_keywords(str(tmp_path))


class TestMetadataMarkTreatmentTerms:
    def test_marks_non_control(self, tmp_config_dir):
        data = {
            'scientific_name': ['Sp1'] * 4,
            'sample_group': ['leaf'] * 4,
            'treatment': ['wild type', 'wild type', 'drought', 'heat'],
            'bioproject': ['PRJ1'] * 4,
            'biosample': ['S1', 'S2', 'S3', 'S4'],
            'run': ['R1', 'R2', 'R3', 'R4'],
            'exclusion': ['no'] * 4,
        }
        df = pandas.DataFrame(data)
        m = Metadata.from_DataFrame(df)
        m.mark_treatment_terms(tmp_config_dir)
        assert m.df.loc[m.df['run'] == 'R1', 'exclusion'].values[0] == 'no'
        assert m.df.loc[m.df['run'] == 'R2', 'exclusion'].values[0] == 'no'
        assert m.df.loc[m.df['run'] == 'R3', 'exclusion'].values[0] == 'non_control'
        assert m.df.loc[m.df['run'] == 'R4', 'exclusion'].values[0] == 'non_control'

    def test_raises_for_malformed_config_with_too_few_columns(self, sample_metadata, tmp_path):
        (tmp_path / 'control_term.config').write_text('treatment_only\n')
        m = sample_metadata
        with pytest.raises(ValueError, match='control_term.config must contain at least 2'):
            m.mark_treatment_terms(str(tmp_path))

    def test_strips_control_term_config_column_tokens(self, tmp_path):
        (tmp_path / 'control_term.config').write_text('treatment, source_name\twild.type\n')
        data = {
            'scientific_name': ['Sp1', 'Sp1', 'Sp1'],
            'sample_group': ['leaf', 'leaf', 'leaf'],
            'treatment': ['drought', 'heat', 'drought'],
            'source_name': ['wild type', '', ''],
            'bioproject': ['PRJ1', 'PRJ1', 'PRJ2'],
            'biosample': ['S1', 'S2', 'S3'],
            'run': ['R1', 'R2', 'R3'],
            'exclusion': ['no', 'no', 'no'],
        }
        m = Metadata.from_DataFrame(pandas.DataFrame(data))
        m.mark_treatment_terms(str(tmp_path))
        assert m.df.loc[m.df['run'] == 'R1', 'exclusion'].values[0] == 'no'
        assert m.df.loc[m.df['run'] == 'R2', 'exclusion'].values[0] == 'non_control'
        assert m.df.loc[m.df['run'] == 'R3', 'exclusion'].values[0] == 'no'

    def test_raises_for_unknown_control_term_config_column(self, tmp_path):
        (tmp_path / 'control_term.config').write_text('unknown_column\twild.type\n')
        data = {
            'scientific_name': ['Sp1'],
            'sample_group': ['leaf'],
            'treatment': ['wild type'],
            'bioproject': ['PRJ1'],
            'biosample': ['S1'],
            'run': ['R1'],
            'exclusion': ['no'],
        }
        m = Metadata.from_DataFrame(pandas.DataFrame(data))
        with pytest.raises(ValueError, match='control_term.config were not found in metadata'):
            m.mark_treatment_terms(str(tmp_path))

    def test_raises_for_invalid_control_term_regex(self, tmp_path):
        (tmp_path / 'control_term.config').write_text('treatment\t[invalid\n')
        data = {
            'scientific_name': ['Sp1'],
            'sample_group': ['leaf'],
            'treatment': ['wild type'],
            'bioproject': ['PRJ1'],
            'biosample': ['S1'],
            'run': ['R1'],
            'exclusion': ['no'],
        }
        m = Metadata.from_DataFrame(pandas.DataFrame(data))
        with pytest.raises(ValueError, match='Invalid regex pattern in control_term.config row 1'):
            m.mark_treatment_terms(str(tmp_path))


class TestMetadataMarkRedundantBiosample:
    def test_marks_duplicates(self, sample_metadata_df):
        df = sample_metadata_df.copy()
        # Add a duplicate biosample within same bioproject
        new_row = df.iloc[0].copy()
        new_row['run'] = 'SRR006'
        new_row['experiment'] = 'SRX6'
        # same bioproject PRJNA1, same biosample SAMN1
        df = pandas.concat([df, pandas.DataFrame([new_row])], ignore_index=True)
        m = Metadata.from_DataFrame(df)
        m.mark_redundant_biosample(True)
        dup = m.df.loc[m.df['run'] == 'SRR006', 'exclusion'].values[0]
        assert dup == 'redundant_biosample'

    def test_no_action_when_disabled(self, sample_metadata):
        m = sample_metadata
        m.mark_redundant_biosample(False)
        assert (m.df['exclusion'] == 'no').all()


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

class TestReadConfigFile:
    def test_reads_config(self, tmp_path):
        config = tmp_path / 'test.config'
        config.write_text('col1\tcol2\tcol3\n')
        df = read_config_file('test.config', str(tmp_path))
        assert df.shape == (1, 3)

    def test_missing_file_returns_empty(self, tmp_path):
        df = read_config_file('nonexistent.config', str(tmp_path))
        assert isinstance(df, pandas.DataFrame)
        assert df.empty

    def test_single_column_returns_series(self, tmp_path):
        config = tmp_path / 'single.config'
        config.write_text('value1\nvalue2\nvalue3\n')
        result = read_config_file('single.config', str(tmp_path))
        assert isinstance(result, pandas.Series)
        assert len(result) == 3

    def test_malformed_config_raises_parser_error(self, tmp_path):
        config = tmp_path / 'bad.config'
        config.write_text('"unterminated\nvalue2\n')
        with pytest.raises(pandas.errors.ParserError):
            read_config_file('bad.config', str(tmp_path))

    def test_raises_when_config_path_is_directory(self, tmp_path):
        (tmp_path / 'dir.config').mkdir()
        with pytest.raises(IsADirectoryError, match='Config path exists but is not a file'):
            read_config_file('dir.config', str(tmp_path))


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


class TestCheckRscript:
    def test_exits_when_rscript_missing(self, monkeypatch):
        monkeypatch.setattr(
            'amalgkit.util.subprocess.run',
            lambda *_args, **_kwargs: (_ for _ in ()).throw(FileNotFoundError('Rscript')),
        )
        with pytest.raises(SystemExit) as exc:
            check_rscript()
        assert exc.value.code == 1

    def test_exits_when_rscript_probe_returns_nonzero(self, monkeypatch):
        monkeypatch.setattr(
            'amalgkit.util.subprocess.run',
            lambda *_args, **_kwargs: SimpleNamespace(returncode=127, stdout=b'', stderr=b''),
        )
        with pytest.raises(SystemExit) as exc:
            check_rscript()
        assert exc.value.code == 1

    def test_passes_when_rscript_probe_returns_zero(self, monkeypatch):
        monkeypatch.setattr(
            'amalgkit.util.subprocess.run',
            lambda *_args, **_kwargs: SimpleNamespace(returncode=0, stdout=b'', stderr=b''),
        )
        check_rscript()


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
        check_ortholog_parameter_compatibility(Args())  # should not raise

    def test_only_busco_ok(self):
        class Args:
            orthogroup_table = None
            dir_busco = '/path/to/busco'
        check_ortholog_parameter_compatibility(Args())  # should not raise

    def test_blank_orthogroup_table_is_treated_as_none(self):
        class Args:
            orthogroup_table = '   '
            dir_busco = None
        with pytest.raises(ValueError, match='One of'):
            check_ortholog_parameter_compatibility(Args())


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
# group_attributes (wiki: column merging via group_attribute.config)
# ---------------------------------------------------------------------------

class TestMetadataGroupAttributes:
    def test_aggregates_source_to_target(self, tmp_path):
        """Wiki: group_attribute.config merges heterogeneous SRA attribute columns."""
        ga = tmp_path / 'group_attribute.config'
        ga.write_text('tissue\tsource_name\n')
        data = {
            'scientific_name': ['Sp1'],
            'sample_group': [''],
            'tissue': [''],
            'source_name': ['brain cortex'],
            'run': ['R1'],
            'exclusion': ['no'],
        }
        m = Metadata.from_DataFrame(pandas.DataFrame(data))
        m.group_attributes(str(tmp_path))
        assert 'brain cortex[source_name]' in m.df.loc[0, 'tissue']

    def test_aggregates_appends_to_nonempty_target(self, tmp_path):
        """When target already has a value, append with semicolon separator."""
        ga = tmp_path / 'group_attribute.config'
        ga.write_text('tissue\tsource_name\n')
        data = {
            'scientific_name': ['Sp1'],
            'sample_group': [''],
            'tissue': ['brain'],
            'source_name': ['frontal lobe'],
            'run': ['R1'],
            'exclusion': ['no'],
        }
        m = Metadata.from_DataFrame(pandas.DataFrame(data))
        m.group_attributes(str(tmp_path))
        assert 'brain' in m.df.loc[0, 'tissue']
        assert 'frontal lobe[source_name]' in m.df.loc[0, 'tissue']

    def test_missing_config_no_error(self, tmp_path):
        """Issue #108: Missing config should not raise error."""
        m = Metadata.from_DataFrame(pandas.DataFrame({
            'scientific_name': ['Sp1'], 'run': ['R1'], 'exclusion': ['no'],
        }))
        m.group_attributes(str(tmp_path))  # no config file present

    def test_no_futurewarning_when_target_column_is_float(self, tmp_path):
        """String aggregation into float target must not emit pandas FutureWarning."""
        ga = tmp_path / 'group_attribute.config'
        ga.write_text('treatment\tgrowth_condition\n')
        data = {
            'scientific_name': ['Sp1', 'Sp1'],
            'sample_group': ['g1', 'g2'],
            'run': ['R1', 'R2'],
            'exclusion': ['no', 'no'],
            'growth_condition': ['sd medium', 'sd medium'],
            # Simulate dtype inferred from all-missing values at file-load time.
            'treatment': [numpy.nan, numpy.nan],
        }
        m = Metadata.from_DataFrame(pandas.DataFrame(data))
        with warnings.catch_warnings(record=True) as captured:
            warnings.simplefilter('always')
            m.group_attributes(str(tmp_path))
        future_warnings = [w for w in captured if issubclass(w.category, FutureWarning)]
        assert len(future_warnings) == 0

    def test_raises_for_malformed_config_with_too_few_columns(self, tmp_path):
        ga = tmp_path / 'group_attribute.config'
        ga.write_text('tissue_only\n')
        m = Metadata.from_DataFrame(pandas.DataFrame({
            'scientific_name': ['Sp1'],
            'run': ['R1'],
            'exclusion': ['no'],
        }))
        with pytest.raises(ValueError, match='group_attribute.config must contain at least 2'):
            m.group_attributes(str(tmp_path))


# ---------------------------------------------------------------------------
# mark_missing_rank
# ---------------------------------------------------------------------------

class TestMetadataMarkMissingRank:
    def test_marks_missing_taxid(self):
        data = {
            'scientific_name': ['Sp1', 'Sp2'],
            'run': ['R1', 'R2'],
            'exclusion': ['no', 'no'],
            'taxid_species': [9606, pandas.NA],
        }
        df = pandas.DataFrame(data)
        df['taxid_species'] = df['taxid_species'].astype('Int64')
        m = Metadata.from_DataFrame(df)
        m.mark_missing_rank('species')
        assert m.df.loc[m.df['run'] == 'R1', 'exclusion'].values[0] == 'no'
        assert m.df.loc[m.df['run'] == 'R2', 'exclusion'].values[0] == 'missing_taxid'

    def test_none_rank_skips(self, sample_metadata):
        """rank_name='none' should do nothing."""
        m = sample_metadata
        m.mark_missing_rank('none')
        assert (m.df['exclusion'] == 'no').all()

    def test_raises_when_required_rank_column_missing(self, sample_metadata):
        m = sample_metadata
        with pytest.raises(ValueError, match='Column \"taxid_species\" is required'):
            m.mark_missing_rank('species')


# ---------------------------------------------------------------------------
# label_sampled_data: empty sample_group handling (wiki: select)
# ---------------------------------------------------------------------------

class TestMetadataLabelSampledDataEdgeCases:
    def test_empty_sample_group_marked_unqualified(self):
        """Wiki/select: samples with empty sample_group get exclusion=no_tissue_label."""
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
        assert empty_sg['exclusion'].values[0] == 'no_tissue_label'
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
# check_config_dir (issue #108: missing configs should warn, not crash)
# ---------------------------------------------------------------------------

class TestCheckConfigDir:
    def test_all_configs_present(self, tmp_config_dir):
        """No error when all config files are present."""
        check_config_dir(tmp_config_dir, mode='select')

    def test_missing_config_warns(self, tmp_path):
        """Issue #108: Missing config files should print warning, not raise."""
        # Create only one of the three required files
        ga = tmp_path / 'group_attribute.config'
        ga.write_text('tissue\tsource_name\n')
        # Should not raise an exception
        check_config_dir(str(tmp_path), mode='select')

    def test_invalid_mode_raises(self, tmp_path):
        with pytest.raises(ValueError, match='Unsupported config check mode'):
            check_config_dir(str(tmp_path), mode='unknown')

    def test_missing_config_dir_raises(self, tmp_path):
        missing_dir = tmp_path / 'missing'
        with pytest.raises(FileNotFoundError, match='Config directory not found'):
            check_config_dir(str(missing_dir), mode='select')

    def test_config_dir_file_path_raises(self, tmp_path):
        file_path = tmp_path / 'config_path'
        file_path.write_text('not a directory')
        with pytest.raises(NotADirectoryError, match='not a directory'):
            check_config_dir(str(file_path), mode='select')

    def test_config_entry_directory_is_treated_as_missing(self, tmp_path, capsys):
        (tmp_path / 'group_attribute.config').write_text('tissue\tsource_name\n')
        (tmp_path / 'exclude_keyword.config').mkdir()
        check_config_dir(str(tmp_path), mode='select')
        captured = capsys.readouterr()
        assert 'Config entry exists but is not a file: exclude_keyword.config' in captured.err


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
# Metadata.nspot_cutoff edge cases (issue #96, #110)
# ---------------------------------------------------------------------------

class TestMetadataNspotCutoffEdgeCases:
    def test_zero_spots_not_marked(self):
        """Rows with total_spots=0 should NOT be marked low_nspots (they are unknown)."""
        data = {
            'scientific_name': ['Sp1'],
            'run': ['R1'],
            'exclusion': ['no'],
            'total_spots': [0],
        }
        m = Metadata.from_DataFrame(pandas.DataFrame(data))
        m.nspot_cutoff(1000000)
        # Zero spots: the negation -(0==0) is False, so row is NOT marked
        assert m.df.loc[0, 'exclusion'] == 'no'

    def test_empty_string_spots(self):
        """Empty string total_spots should be handled gracefully."""
        data = {
            'scientific_name': ['Sp1'],
            'run': ['R1'],
            'exclusion': ['no'],
            'total_spots': [''],
        }
        m = Metadata.from_DataFrame(pandas.DataFrame(data))
        m.nspot_cutoff(1000000)
        # Empty string -> converted to 0, not marked
        assert m.df.loc[0, 'exclusion'] == 'no'


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

        with pytest.raises(SystemExit) as exc:
            load_metadata(Args())
        assert exc.value.code == 1

    def test_batch_mode_handles_missing_caller_module(self, tmp_path, monkeypatch):
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

        monkeypatch.setattr('amalgkit.util.inspect.getmodule', lambda *_args, **_kwargs: None)
        m = load_metadata(Args())
        assert m.df.shape[0] == 1
        assert m.df.iloc[0]['run'] == 'R1'

    def test_curate_batch_ignores_missing_species_names(self, tmp_path, monkeypatch):
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

        monkeypatch.setattr(
            'amalgkit.util.inspect.getmodule',
            lambda *_args, **_kwargs: type('DummyModule', (), {'__name__': 'amalgkit.curate'})(),
        )
        m = load_metadata(Args())
        assert m.df.shape[0] == 1
        assert m.df.iloc[0]['scientific_name'] == 'Sp1'

    def test_curate_batch_raises_when_no_valid_species(self, tmp_path, monkeypatch):
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

        monkeypatch.setattr(
            'amalgkit.util.inspect.getmodule',
            lambda *_args, **_kwargs: type('DummyModule', (), {'__name__': 'amalgkit.curate'})(),
        )
        with pytest.raises(ValueError, match='No valid scientific_name'):
            load_metadata(Args())


# ---------------------------------------------------------------------------
# detect_layout_from_file (corrects layout based on actual files)
# ---------------------------------------------------------------------------

class TestDetectLayoutFromFile:
    def test_paired_files_detected(self, tmp_path):
        """Paired-end files detected -> layout corrected to 'paired'."""
        sra_stat = {
            'sra_id': 'SRR001',
            'layout': 'single',
            'getfastq_sra_dir': str(tmp_path),
        }
        (tmp_path / 'SRR001_1.amalgkit.fastq.gz').write_text('data')
        (tmp_path / 'SRR001_2.amalgkit.fastq.gz').write_text('data')
        result = detect_layout_from_file(sra_stat)
        assert result['layout'] == 'paired'

    def test_single_files_detected(self, tmp_path):
        """Single-end file detected -> layout corrected to 'single'."""
        sra_stat = {
            'sra_id': 'SRR001',
            'layout': 'paired',
            'getfastq_sra_dir': str(tmp_path),
        }
        (tmp_path / 'SRR001.amalgkit.fastq.gz').write_text('data')
        result = detect_layout_from_file(sra_stat)
        assert result['layout'] == 'single'

    def test_paired_plain_fastq_files_detected(self, tmp_path):
        sra_stat = {
            'sra_id': 'SRR001',
            'layout': 'single',
            'getfastq_sra_dir': str(tmp_path),
        }
        (tmp_path / 'SRR001_1.fastq').write_text('data')
        (tmp_path / 'SRR001_2.fastq').write_text('data')
        result = detect_layout_from_file(sra_stat)
        assert result['layout'] == 'paired'

    def test_no_files_keeps_layout(self, tmp_path):
        """No fastq files -> layout unchanged."""
        sra_stat = {
            'sra_id': 'SRR001',
            'layout': 'paired',
            'getfastq_sra_dir': str(tmp_path),
        }
        result = detect_layout_from_file(sra_stat)
        assert result['layout'] == 'paired'

    def test_paired_safely_removed_files_detected(self, tmp_path):
        sra_stat = {
            'sra_id': 'SRR001',
            'layout': 'single',
            'getfastq_sra_dir': str(tmp_path),
        }
        (tmp_path / 'SRR001_1.fastq.gz.safely_removed').write_text('')
        (tmp_path / 'SRR001_2.fastq.gz.safely_removed').write_text('')
        result = detect_layout_from_file(sra_stat)
        assert result['layout'] == 'paired'

    def test_single_safely_removed_file_detected(self, tmp_path):
        sra_stat = {
            'sra_id': 'SRR001',
            'layout': 'paired',
            'getfastq_sra_dir': str(tmp_path),
        }
        (tmp_path / 'SRR001.fastq.gz.safely_removed').write_text('')
        result = detect_layout_from_file(sra_stat)
        assert result['layout'] == 'single'


# ---------------------------------------------------------------------------
# is_there_unpaired_file
# ---------------------------------------------------------------------------

class TestIsThereUnpairedFile:
    def test_unpaired_file_present(self, tmp_path):
        sra_stat = {
            'sra_id': 'SRR001',
            'getfastq_sra_dir': str(tmp_path),
        }
        (tmp_path / 'SRR001.amalgkit.fastq.gz').write_text('data')
        assert is_there_unpaired_file(sra_stat, ['.amalgkit.fastq.gz']) is True

    def test_no_unpaired_file(self, tmp_path):
        sra_stat = {
            'sra_id': 'SRR001',
            'getfastq_sra_dir': str(tmp_path),
        }
        assert is_there_unpaired_file(sra_stat, ['.amalgkit.fastq.gz']) is False


# ---------------------------------------------------------------------------
# get_newest_intermediate_file_extension
# ---------------------------------------------------------------------------

class TestGetNewestIntermediateFileExtension:
    def test_finds_amalgkit_extension(self, tmp_path):
        """Finds .amalgkit.fastq.gz as the newest extension."""
        sra_stat = {
            'sra_id': 'SRR001',
            'layout': 'paired',
            'getfastq_sra_dir': str(tmp_path),
        }
        (tmp_path / 'SRR001_1.amalgkit.fastq.gz').write_text('data')
        (tmp_path / 'SRR001_2.amalgkit.fastq.gz').write_text('data')
        (tmp_path / 'SRR001_1.fastq.gz').write_text('data')
        (tmp_path / 'SRR001_2.fastq.gz').write_text('data')
        ext = get_newest_intermediate_file_extension(sra_stat, str(tmp_path))
        assert ext == '.amalgkit.fastq.gz'

    def test_finds_fastq_extension(self, tmp_path):
        """Falls back to .fastq.gz when no downstream extensions exist."""
        sra_stat = {
            'sra_id': 'SRR001',
            'layout': 'single',
            'getfastq_sra_dir': str(tmp_path),
        }
        (tmp_path / 'SRR001.fastq.gz').write_text('data')
        ext = get_newest_intermediate_file_extension(sra_stat, str(tmp_path))
        assert ext == '.fastq.gz'

    def test_finds_plain_fastq_extension(self, tmp_path):
        sra_stat = {
            'sra_id': 'SRR001',
            'layout': 'single',
            'getfastq_sra_dir': str(tmp_path),
        }
        (tmp_path / 'SRR001.fastq').write_text('data')
        ext = get_newest_intermediate_file_extension(sra_stat, str(tmp_path))
        assert ext == '.fastq'

    def test_safely_removed(self, tmp_path):
        """Detects .safely_removed flag."""
        sra_stat = {
            'sra_id': 'SRR001',
            'layout': 'paired',
            'getfastq_sra_dir': str(tmp_path),
        }
        (tmp_path / 'SRR001_1.amalgkit.fastq.gz.safely_removed').write_text('')
        (tmp_path / 'SRR001_2.amalgkit.fastq.gz.safely_removed').write_text('')
        ext = get_newest_intermediate_file_extension(sra_stat, str(tmp_path))
        assert ext == '.safely_removed'

    def test_safely_removed_ignores_similar_run_prefix(self, tmp_path):
        sra_stat = {
            'sra_id': 'SRR001',
            'layout': 'paired',
            'getfastq_sra_dir': str(tmp_path),
        }
        (tmp_path / 'SRR0010_1.amalgkit.fastq.gz.safely_removed').write_text('')
        ext = get_newest_intermediate_file_extension(sra_stat, str(tmp_path))
        assert ext == 'no_extension_found'

    def test_safely_removed_requires_both_mates_for_paired_layout(self, tmp_path):
        sra_stat = {
            'sra_id': 'SRR001',
            'layout': 'paired',
            'getfastq_sra_dir': str(tmp_path),
        }
        (tmp_path / 'SRR001_1.amalgkit.fastq.gz.safely_removed').write_text('')
        ext = get_newest_intermediate_file_extension(sra_stat, str(tmp_path))
        assert ext == 'no_extension_found'

    def test_safely_removed_detected_despite_layout_mismatch_single_to_paired(self, tmp_path):
        sra_stat = {
            'sra_id': 'SRR001',
            'layout': 'single',
            'getfastq_sra_dir': str(tmp_path),
        }
        (tmp_path / 'SRR001_1.fastq.gz.safely_removed').write_text('')
        (tmp_path / 'SRR001_2.fastq.gz.safely_removed').write_text('')
        ext = get_newest_intermediate_file_extension(sra_stat, str(tmp_path))
        assert ext == '.safely_removed'

    def test_safely_removed_detected_despite_layout_mismatch_paired_to_single(self, tmp_path):
        sra_stat = {
            'sra_id': 'SRR001',
            'layout': 'paired',
            'getfastq_sra_dir': str(tmp_path),
        }
        (tmp_path / 'SRR001.fastq.gz.safely_removed').write_text('')
        ext = get_newest_intermediate_file_extension(sra_stat, str(tmp_path))
        assert ext == '.safely_removed'

    def test_paired_extension_detection_requires_both_fastq_mates(self, tmp_path):
        sra_stat = {
            'sra_id': 'SRR001',
            'layout': 'paired',
            'getfastq_sra_dir': str(tmp_path),
        }
        (tmp_path / 'SRR001_1.fastq.gz').write_text('data')
        ext = get_newest_intermediate_file_extension(sra_stat, str(tmp_path))
        assert ext == 'no_extension_found'

    def test_no_extension_found(self, tmp_path):
        """No matching files -> 'no_extension_found'."""
        sra_stat = {
            'sra_id': 'SRR001',
            'layout': 'single',
            'getfastq_sra_dir': str(tmp_path),
        }
        ext = get_newest_intermediate_file_extension(sra_stat, str(tmp_path))
        assert ext == 'no_extension_found'

    def test_ignores_directory_named_as_fastq_file(self, tmp_path):
        sra_stat = {
            'sra_id': 'SRR001',
            'layout': 'single',
            'getfastq_sra_dir': str(tmp_path),
        }
        (tmp_path / 'SRR001.fastq.gz').mkdir()
        ext = get_newest_intermediate_file_extension(sra_stat, str(tmp_path))
        assert ext == 'no_extension_found'


# ---------------------------------------------------------------------------
# get_mapping_rate (extracts p_pseudoaligned from quant run_info.json)
# ---------------------------------------------------------------------------

class TestGetMappingRate:
    def test_extracts_mapping_rate(self, tmp_path, sample_metadata):
        """Reads p_pseudoaligned from run_info.json into mapping_rate column."""
        quant_dir = tmp_path / 'quant'
        sra_dir = quant_dir / 'SRR001'
        sra_dir.mkdir(parents=True)
        run_info = {'p_pseudoaligned': 85.5}
        (sra_dir / 'SRR001_run_info.json').write_text(json.dumps(run_info))
        m = get_mapping_rate(sample_metadata, str(quant_dir))
        assert m.df.loc[m.df['run'] == 'SRR001', 'mapping_rate'].values[0] == 85.5

    def test_avoids_quant_root_listdir_scan(self, tmp_path, sample_metadata, monkeypatch):
        quant_dir = tmp_path / 'quant'
        sra_dir = quant_dir / 'SRR001'
        sra_dir.mkdir(parents=True)
        (sra_dir / 'SRR001_run_info.json').write_text(json.dumps({'p_pseudoaligned': 75.0}))

        def fail_if_listdir_called(_path):
            raise AssertionError('get_mapping_rate should not scan quant root with os.listdir.')

        monkeypatch.setattr('amalgkit.util.os.listdir', fail_if_listdir_called)

        m = get_mapping_rate(sample_metadata, str(quant_dir))

        assert m.df.loc[m.df['run'] == 'SRR001', 'mapping_rate'].values[0] == 75.0

    def test_skips_nan_run_ids_when_scanning_quant_subdirs(self, tmp_path):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001', numpy.nan],
            'scientific_name': ['sp1', 'sp1'],
            'exclusion': ['no', 'no'],
        }))
        quant_dir = tmp_path / 'quant'
        sra_dir = quant_dir / 'SRR001'
        sra_dir.mkdir(parents=True)
        (sra_dir / 'SRR001_run_info.json').write_text(json.dumps({'p_pseudoaligned': 66.6}))

        m = get_mapping_rate(metadata, str(quant_dir))

        assert m.df.loc[m.df['run'] == 'SRR001', 'mapping_rate'].values[0] == 66.6

    def test_missing_quant_dir(self, sample_metadata, tmp_path):
        """Missing quant directory does not raise."""
        m = get_mapping_rate(sample_metadata, str(tmp_path / 'nonexistent'))
        assert isinstance(m, Metadata)

    def test_missing_run_info(self, tmp_path, sample_metadata):
        """Missing run_info.json for an SRA -> mapping_rate stays NaN."""
        quant_dir = tmp_path / 'quant'
        sra_dir = quant_dir / 'SRR001'
        sra_dir.mkdir(parents=True)
        # No run_info.json
        m = get_mapping_rate(sample_metadata, str(quant_dir))
        assert numpy.isnan(m.df.loc[m.df['run'] == 'SRR001', 'mapping_rate'].values[0])

    def test_quant_path_is_file_does_not_crash(self, tmp_path, sample_metadata):
        """quant path that exists as a file should be treated as missing directory."""
        quant_file = tmp_path / 'quant'
        quant_file.write_text('not a directory')
        m = get_mapping_rate(sample_metadata, str(quant_file))
        assert isinstance(m, Metadata)

    def test_raises_when_run_column_missing(self, tmp_path):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'scientific_name': ['sp1'],
            'exclusion': ['no'],
        }))
        metadata.df = metadata.df.drop(columns=['run'])

        with pytest.raises(ValueError, match='Column \"run\" is required in metadata to compute mapping_rate'):
            get_mapping_rate(metadata, str(tmp_path / 'quant'))

    def test_duplicate_run_ids_are_all_updated(self, tmp_path):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001', 'SRR001'],
            'scientific_name': ['sp1', 'sp1'],
            'exclusion': ['no', 'no'],
        }))
        quant_dir = tmp_path / 'quant'
        sra_dir = quant_dir / 'SRR001'
        sra_dir.mkdir(parents=True)
        (sra_dir / 'SRR001_run_info.json').write_text(json.dumps({'p_pseudoaligned': 42.0}))

        m = get_mapping_rate(metadata, str(quant_dir))

        assert m.df['mapping_rate'].tolist() == [42.0, 42.0]

    def test_metadata_run_whitespace_is_trimmed_for_matching(self, tmp_path):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': [' SRR001 '],
            'scientific_name': ['sp1'],
            'exclusion': ['no'],
        }))
        quant_dir = tmp_path / 'quant'
        sra_dir = quant_dir / 'SRR001'
        sra_dir.mkdir(parents=True)
        (sra_dir / 'SRR001_run_info.json').write_text(json.dumps({'p_pseudoaligned': 33.3}))

        m = get_mapping_rate(metadata, str(quant_dir))

        assert m.df.loc[0, 'mapping_rate'] == 33.3

    def test_invalid_pseudoaligned_value_is_skipped(self, tmp_path, sample_metadata):
        quant_dir = tmp_path / 'quant'
        sra_dir = quant_dir / 'SRR001'
        sra_dir.mkdir(parents=True)
        (sra_dir / 'SRR001_run_info.json').write_text(json.dumps({'p_pseudoaligned': 'not-a-number'}))

        m = get_mapping_rate(sample_metadata, str(quant_dir))

        assert numpy.isnan(m.df.loc[m.df['run'] == 'SRR001', 'mapping_rate'].values[0])

    def test_numeric_string_pseudoaligned_value_is_parsed(self, tmp_path, sample_metadata):
        quant_dir = tmp_path / 'quant'
        sra_dir = quant_dir / 'SRR001'
        sra_dir.mkdir(parents=True)
        (sra_dir / 'SRR001_run_info.json').write_text(json.dumps({'p_pseudoaligned': '12.5'}))

        m = get_mapping_rate(sample_metadata, str(quant_dir))

        assert m.df.loc[m.df['run'] == 'SRR001', 'mapping_rate'].values[0] == 12.5

    def test_respects_max_workers_override(self, tmp_path, sample_metadata, monkeypatch):
        quant_dir = tmp_path / 'quant'
        for sra_id, value in [('SRR001', 12.5), ('SRR002', 33.3)]:
            sra_dir = quant_dir / sra_id
            sra_dir.mkdir(parents=True)
            (sra_dir / f'{sra_id}_run_info.json').write_text(json.dumps({'p_pseudoaligned': value}))
        observed = {}

        def fake_run_tasks(task_items, task_fn, max_workers):
            observed['max_workers'] = max_workers
            return {}, []

        monkeypatch.setattr('amalgkit.util.run_tasks_with_optional_threads', fake_run_tasks)
        get_mapping_rate(sample_metadata, str(quant_dir), max_workers=1)
        assert observed['max_workers'] == 1


# ---------------------------------------------------------------------------
# get_getfastq_run_dir (creates SRA-specific output directory)
# ---------------------------------------------------------------------------

class TestGetGetfastqRunDir:
    def test_creates_directory(self, tmp_path):
        """Creates getfastq/SRR001 directory and returns path."""
        class Args:
            out_dir = str(tmp_path)
        result = get_getfastq_run_dir(Args(), 'SRR001')
        assert os.path.isdir(result)
        assert result.endswith(os.path.join('getfastq', 'SRR001'))

    def test_existing_directory(self, tmp_path):
        """Returns existing directory without error."""
        class Args:
            out_dir = str(tmp_path)
        gf_dir = tmp_path / 'getfastq' / 'SRR001'
        gf_dir.mkdir(parents=True)
        result = get_getfastq_run_dir(Args(), 'SRR001')
        assert os.path.isdir(result)

    def test_rejects_out_dir_file_path(self, tmp_path):
        out_file = tmp_path / 'out_file'
        out_file.write_text('not a directory')

        class Args:
            out_dir = str(out_file)

        with pytest.raises(NotADirectoryError, match='Output path exists but is not a directory'):
            get_getfastq_run_dir(Args(), 'SRR001')

    def test_rejects_getfastq_root_file_path(self, tmp_path):
        (tmp_path / 'getfastq').write_text('not a directory')

        class Args:
            out_dir = str(tmp_path)

        with pytest.raises(NotADirectoryError, match='getfastq path exists but is not a directory'):
            get_getfastq_run_dir(Args(), 'SRR001')


# ---------------------------------------------------------------------------
# generate_multisp_busco_table (merges BUSCO full_table.tsv files)
# ---------------------------------------------------------------------------

class TestGenerateMultispBuscoTable:
    def test_raises_when_busco_dir_missing(self, tmp_path):
        outfile = tmp_path / 'merged.tsv'
        missing_dir = tmp_path / 'busco_missing'

        with pytest.raises(FileNotFoundError, match='BUSCO directory not found'):
            generate_multisp_busco_table(str(missing_dir), str(outfile))

    def test_raises_when_busco_path_is_file(self, tmp_path):
        outfile = tmp_path / 'merged.tsv'
        busco_file = tmp_path / 'busco.tsv'
        busco_file.write_text('not a directory')

        with pytest.raises(NotADirectoryError, match='BUSCO path exists but is not a directory'):
            generate_multisp_busco_table(str(busco_file), str(outfile))

    def test_raises_when_no_busco_table_files_detected(self, tmp_path):
        outfile = tmp_path / 'merged.tsv'
        busco_dir = tmp_path / 'busco'
        busco_dir.mkdir()

        with pytest.raises(FileNotFoundError, match='No BUSCO full table file'):
            generate_multisp_busco_table(str(busco_dir), str(outfile))

    def test_merges_two_species(self, tmp_path):
        """Merges BUSCO tables from two species into one output file."""
        busco_dir = tmp_path / 'busco'
        busco_dir.mkdir()
        content_a = (
            '# comment line\n'
            'OG0001\tComplete\tgene1\t100\t200\thttp://odb\tgene desc\n'
            'OG0002\tComplete\tgene2\t90\t150\thttp://odb2\tgene desc2\n'
        )
        content_b = (
            '# comment line\n'
            'OG0001\tComplete\tgeneA\t95\t180\t-\t-\n'
            'OG0002\tMissing\t-\t0\t0\t-\t-\n'
        )
        (busco_dir / 'Species_A.tsv').write_text(content_a)
        (busco_dir / 'Species_B.tsv').write_text(content_b)
        outfile = tmp_path / 'merged.tsv'
        generate_multisp_busco_table(str(busco_dir), str(outfile))
        result = pandas.read_csv(str(outfile), sep='\t')
        assert 'Species_A' in result.columns
        assert 'Species_B' in result.columns
        assert result.shape[0] == 2

    def test_raises_when_all_busco_tables_fail_to_parse(self, tmp_path, monkeypatch):
        busco_dir = tmp_path / 'busco'
        busco_dir.mkdir()
        (busco_dir / 'Species_A.tsv').write_text('x')
        outfile = tmp_path / 'merged.tsv'

        def fake_run_tasks(task_items, task_fn, max_workers):
            return {}, [(task_items[0], RuntimeError('bad table'))]

        monkeypatch.setattr('amalgkit.util.run_tasks_with_optional_threads', fake_run_tasks)

        with pytest.warns(UserWarning, match='Failed to parse BUSCO table'):
            with pytest.raises(ValueError, match='Failed to parse any BUSCO table'):
                generate_multisp_busco_table(str(busco_dir), str(outfile))

    def test_raises_on_duplicate_species_label_after_filename_parsing(self, tmp_path):
        busco_dir = tmp_path / 'busco'
        busco_dir.mkdir()
        content = (
            '# comment line\n'
            'OG0001\tComplete\tgene1\t100\t200\thttp://odb\tgene desc\n'
        )
        (busco_dir / 'Homo_sapiens_strain1.tsv').write_text(content)
        (busco_dir / 'Homo_sapiens_strain2.tsv').write_text(content)
        outfile = tmp_path / 'merged.tsv'

        with pytest.raises(ValueError, match='Duplicate species label was detected across BUSCO tables'):
            generate_multisp_busco_table(str(busco_dir), str(outfile))

    def test_ignores_uncommented_busco_header_row(self, tmp_path):
        busco_dir = tmp_path / 'busco'
        busco_dir.mkdir()
        content = (
            'Busco id\tStatus\tSequence\tScore\tLength\tOrthoDB url\tDescription\n'
            'OG0001\tComplete\tgene1\t100\t200\thttp://odb\tgene desc\n'
        )
        (busco_dir / 'Species_A.tsv').write_text(content)
        (busco_dir / 'Species_B.tsv').write_text(content)
        outfile = tmp_path / 'merged.tsv'

        generate_multisp_busco_table(str(busco_dir), str(outfile))

        result = pandas.read_csv(str(outfile), sep='\t')
        assert result['busco_id'].tolist() == ['OG0001']

    def test_accepts_gzipped_tsv_inputs(self, tmp_path):
        import gzip
        busco_dir = tmp_path / 'busco'
        busco_dir.mkdir()
        content = (
            '# comment line\n'
            'OG0001\tComplete\tgene1\t100\t200\thttp://odb\tgene desc\n'
        )
        with gzip.open(busco_dir / 'Species_A.tsv.gz', 'wt') as f:
            f.write(content)
        outfile = tmp_path / 'merged.tsv'

        generate_multisp_busco_table(str(busco_dir), str(outfile))

        result = pandas.read_csv(str(outfile), sep='\t')
        assert 'Species_A' in result.columns
        assert result.shape[0] == 1

    def test_ignores_directory_named_like_busco_table(self, tmp_path):
        busco_dir = tmp_path / 'busco'
        busco_dir.mkdir()
        content = (
            '# comment line\n'
            'OG0001\tComplete\tgene1\t100\t200\thttp://odb\tgene desc\n'
        )
        (busco_dir / 'Species_A.tsv').write_text(content)
        (busco_dir / 'Species_B.tsv').mkdir()
        outfile = tmp_path / 'merged.tsv'

        with warnings.catch_warnings(record=True) as records:
            warnings.simplefilter('always')
            generate_multisp_busco_table(str(busco_dir), str(outfile))

        assert len(records) == 0
        result = pandas.read_csv(str(outfile), sep='\t')
        assert 'Species_A' in result.columns
        assert 'Species_B' not in result.columns


class TestMetadataTaxidValidation:
    def test_add_standard_rank_taxids_requires_nullable_int64_taxid(self, monkeypatch):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'scientific_name': ['Homo sapiens'],
            'exclusion': ['no'],
            'taxid': ['9606'],
        }))

        def fail_if_called():
            raise AssertionError('NCBITaxa should not be initialized for invalid taxid dtype.')

        monkeypatch.setattr('amalgkit.util.ete4.NCBITaxa', fail_if_called)

        with pytest.raises(TypeError, match='taxid column must be Int64 dtype'):
            metadata.add_standard_rank_taxids()

    def test_resolve_scientific_names_requires_nullable_int64_taxid(self, monkeypatch):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'scientific_name': ['H. sapiens'],
            'exclusion': ['no'],
            'taxid': ['9606'],
        }))

        def fail_if_called():
            raise AssertionError('NCBITaxa should not be initialized for invalid taxid dtype.')

        monkeypatch.setattr('amalgkit.util.ete4.NCBITaxa', fail_if_called)

        with pytest.raises(TypeError, match='taxid column must be Int64 dtype'):
            metadata.resolve_scientific_names()

    def test_add_standard_rank_taxids_does_not_swallow_keyboard_interrupt(self, monkeypatch):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'scientific_name': ['Homo sapiens'],
            'exclusion': ['no'],
            'taxid': [9606],
        }))
        metadata.df['taxid'] = metadata.df['taxid'].astype('Int64')

        class InterruptingNcbi:
            def get_lineage(self, _taxid):
                raise KeyboardInterrupt()

            def get_rank(self, _lineage):
                return {}

        monkeypatch.setattr('amalgkit.util.ete4.NCBITaxa', lambda: InterruptingNcbi())

        with pytest.raises(KeyboardInterrupt):
            metadata.add_standard_rank_taxids()

    def test_add_standard_rank_taxids_warns_on_lineage_lookup_error(self, monkeypatch):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'scientific_name': ['Homo sapiens'],
            'exclusion': ['no'],
            'taxid': [9606],
        }))
        metadata.df['taxid'] = metadata.df['taxid'].astype('Int64')

        class FailingNcbi:
            def get_lineage(self, _taxid):
                raise RuntimeError('boom')

            def get_rank(self, _lineage):
                return {}

        monkeypatch.setattr('amalgkit.util.ete4.NCBITaxa', lambda: FailingNcbi())

        with pytest.warns(UserWarning, match='Failed to resolve NCBI lineage'):
            metadata.add_standard_rank_taxids()
        assert 'taxid_species' in metadata.df.columns
        assert metadata.df['taxid_species'].isna().all()
