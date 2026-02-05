import re
import pytest
import pandas

from amalgkit.curate import get_sample_group
from amalgkit.util import Metadata


# ---------------------------------------------------------------------------
# get_sample_group (wiki: curate extracts sample groups from args or metadata)
# ---------------------------------------------------------------------------

class TestGetSampleGroup:
    def test_from_args(self):
        """When --sample_group is specified, parse it."""
        class Args:
            sample_group = 'brain,liver,heart'
            metadata = 'metadata.tsv'
        m = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['R1'], 'exclusion': ['no'],
            'sample_group': ['brain'],
        }))
        result = get_sample_group(Args(), m)
        assert result == 'brain|liver|heart'

    def test_from_metadata(self):
        """When --sample_group is None, extract from metadata."""
        class Args:
            sample_group = None
            metadata = 'metadata.tsv'
        m = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['R1', 'R2', 'R3'],
            'exclusion': ['no', 'no', 'no'],
            'sample_group': ['brain', 'liver', 'brain'],
        }))
        result = get_sample_group(Args(), m)
        # Should contain brain and liver separated by pipe
        groups = result.split('|')
        assert 'brain' in groups
        assert 'liver' in groups

    def test_empty_sample_group_exits(self):
        """Wiki: curate exits with error if sample_group is empty."""
        class Args:
            sample_group = None
            metadata = 'metadata.tsv'
        m = Metadata.from_DataFrame(pandas.DataFrame({
            'run': ['R1'],
            'exclusion': ['no'],
            'sample_group': [float('nan')],
        }))
        with pytest.raises(SystemExit):
            get_sample_group(Args(), m)
