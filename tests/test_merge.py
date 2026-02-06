import numpy
import pandas

from amalgkit.merge import merge_fastp_stats_into_metadata
from amalgkit.util import Metadata


class TestMergeFastpStatsIntoMetadata:
    def test_merges_fastp_stats_tsv(self, tmp_path):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001', 'SRR002'],
            'scientific_name': ['sp1', 'sp1'],
            'exclusion': ['no', 'no'],
        }))
        srr001_dir = tmp_path / 'getfastq' / 'SRR001'
        srr001_dir.mkdir(parents=True)
        pandas.DataFrame([{
            'run': 'SRR001',
            'fastp_duplication_rate': 12.5,
            'fastp_insert_size_peak': 250.0,
        }]).to_csv(srr001_dir / 'fastp_stats.tsv', sep='\t', index=False)

        metadata = merge_fastp_stats_into_metadata(metadata, str(tmp_path))

        assert metadata.df.loc[0, 'fastp_duplication_rate'] == 12.5
        assert metadata.df.loc[0, 'fastp_insert_size_peak'] == 250.0
        assert numpy.isnan(metadata.df.loc[1, 'fastp_duplication_rate'])
        assert numpy.isnan(metadata.df.loc[1, 'fastp_insert_size_peak'])

    def test_adds_columns_when_getfastq_dir_missing(self, tmp_path):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'scientific_name': ['sp1'],
            'exclusion': ['no'],
        }))

        metadata = merge_fastp_stats_into_metadata(metadata, str(tmp_path))

        assert 'fastp_duplication_rate' in metadata.df.columns
        assert 'fastp_insert_size_peak' in metadata.df.columns
        assert numpy.isnan(metadata.df.loc[0, 'fastp_duplication_rate'])
        assert numpy.isnan(metadata.df.loc[0, 'fastp_insert_size_peak'])

    def test_keeps_nan_when_column_missing_in_tsv(self, tmp_path):
        metadata = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['SRR001'],
            'scientific_name': ['sp1'],
            'exclusion': ['no'],
        }))
        srr001_dir = tmp_path / 'getfastq' / 'SRR001'
        srr001_dir.mkdir(parents=True)
        pandas.DataFrame([{
            'run': 'SRR001',
            'fastp_duplication_rate': 44.0,
        }]).to_csv(srr001_dir / 'fastp_stats.tsv', sep='\t', index=False)

        metadata = merge_fastp_stats_into_metadata(metadata, str(tmp_path))

        assert metadata.df.loc[0, 'fastp_duplication_rate'] == 44.0
        assert numpy.isnan(metadata.df.loc[0, 'fastp_insert_size_peak'])
