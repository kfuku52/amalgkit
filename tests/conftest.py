import sys
import types
import os
import pytest
import pandas
import numpy

# Mock ete4 if not installed, so amalgkit.util can be imported in test environments
if 'ete4' not in sys.modules:
    try:
        import ete4
    except ModuleNotFoundError:
        ete4_mock = types.ModuleType('ete4')
        class _DummyNcbiTaxa:
            def __init__(self, *args, **kwargs):
                pass
        ete4_mock.NCBITaxa = _DummyNcbiTaxa
        sys.modules['ete4'] = ete4_mock

from amalgkit.util import Metadata


@pytest.fixture
def sample_metadata_df():
    """A small realistic metadata DataFrame."""
    data = {
        'scientific_name': ['Homo sapiens', 'Homo sapiens', 'Homo sapiens',
                            'Mus musculus', 'Mus musculus'],
        'sample_group': ['brain', 'brain', 'liver', 'brain', 'liver'],
        'tissue': ['brain', 'brain', 'liver', 'brain', 'liver'],
        'genotype': ['', '', '', '', ''],
        'sex': ['male', 'female', 'male', 'female', 'male'],
        'age': ['', '', '', '', ''],
        'treatment': ['', '', '', '', ''],
        'source_name': ['', '', '', '', ''],
        'is_sampled': ['yes', 'yes', 'yes', 'yes', 'yes'],
        'is_qualified': ['yes', 'yes', 'yes', 'yes', 'yes'],
        'exclusion': ['no', 'no', 'no', 'no', 'no'],
        'protocol': ['', '', '', '', ''],
        'bioproject': ['PRJNA1', 'PRJNA1', 'PRJNA2', 'PRJNA3', 'PRJNA3'],
        'biosample': ['SAMN1', 'SAMN2', 'SAMN3', 'SAMN4', 'SAMN5'],
        'experiment': ['SRX1', 'SRX2', 'SRX3', 'SRX4', 'SRX5'],
        'run': ['SRR001', 'SRR002', 'SRR003', 'SRR004', 'SRR005'],
        'sra_primary': ['SRA1', 'SRA2', 'SRA3', 'SRA4', 'SRA5'],
        'sra_sample': ['SRS1', 'SRS2', 'SRS3', 'SRS4', 'SRS5'],
        'sra_study': ['SRP1', 'SRP1', 'SRP2', 'SRP3', 'SRP3'],
        'study_title': ['Study1', 'Study1', 'Study2', 'Study3', 'Study3'],
        'exp_title': ['', '', '', '', ''],
        'design': ['', '', '', '', ''],
        'sample_title': ['', '', '', '', ''],
        'sample_description': ['', '', '', '', ''],
        'lib_name': ['', '', '', '', ''],
        'lib_layout': ['paired', 'paired', 'single', 'paired', 'single'],
        'lib_strategy': ['RNA-Seq', 'RNA-Seq', 'RNA-Seq', 'RNA-Seq', 'RNA-Seq'],
        'lib_source': ['TRANSCRIPTOMIC', 'TRANSCRIPTOMIC', 'TRANSCRIPTOMIC',
                        'TRANSCRIPTOMIC', 'TRANSCRIPTOMIC'],
        'lib_selection': ['cDNA', 'cDNA', 'cDNA', 'cDNA', 'cDNA'],
        'instrument': ['Illumina HiSeq 2500'] * 5,
        'total_spots': [10000000, 5000000, 100, 8000000, 200],
        'total_bases': [2000000000, 1000000000, 20000, 1600000000, 40000],
        'size': [500000000, 250000000, 5000, 400000000, 10000],
        'nominal_length': [200, 200, 0, 200, 0],
        'nominal_sdev': [0, 0, 0, 0, 0],
        'spot_length': [200, 200, 100, 200, 100],
        'read_index': ['', '', '', '', ''],
        'read_class': ['', '', '', '', ''],
        'read_type': ['', '', '', '', ''],
        'base_coord': ['', '', '', '', ''],
        'center': ['', '', '', '', ''],
        'submitter_id': ['', '', '', '', ''],
        'pubmed_id': ['', '', '', '', ''],
        'taxid': [9606, 9606, 9606, 10090, 10090],
        'published_date': ['2020-01-01', '2020-02-01', '2020-03-01',
                           '2020-04-01', '2020-05-01'],
        'NCBI_Link': ['', '', '', '', ''],
        'AWS_Link': ['', '', '', '', ''],
        'GCP_Link': ['', '', '', '', ''],
    }
    return pandas.DataFrame(data)


@pytest.fixture
def sample_metadata(sample_metadata_df):
    """A Metadata object built from sample data."""
    return Metadata.from_DataFrame(sample_metadata_df)


@pytest.fixture
def tmp_config_dir(tmp_path):
    """Temporary directory with mock config files."""
    # group_attribute.config: aggregate 'source_name' into 'tissue'
    ga = tmp_path / 'group_attribute.config'
    ga.write_text('tissue\tsource_name\n')

    # exclude_keyword.config: exclude samples with 'cancer' in sample_description
    ek = tmp_path / 'exclude_keyword.config'
    ek.write_text('sample_description\tdisease\tcancer\n')

    # control_term.config: mark 'wild type' as control in treatment column
    ct = tmp_path / 'control_term.config'
    ct.write_text('treatment\twild.type\n')

    return str(tmp_path)
