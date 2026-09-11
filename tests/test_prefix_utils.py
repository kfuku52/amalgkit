import pytest
from amalgkit.prefix_utils import find_prefixed_entries, find_run_prefixed_entries, find_species_prefixed_entries


@pytest.mark.parametrize('sorted_input', [False, True])
def test_prefix_matching_preserves_delimiter_boundaries(sorted_input):
    entries = ['SRR10', 'SRR1_1.fastq', 'SRR1', 'SRR1-v1', 'SRR1.fastq', 'other']
    if sorted_input:
        entries.sort()
    assert find_run_prefixed_entries(entries, 'SRR1', sorted_input) == ['SRR1', 'SRR1.fastq', 'SRR1_1.fastq']
    assert find_species_prefixed_entries(entries, 'SRR1', sorted_input) == ['SRR1', 'SRR1-v1', 'SRR1.fastq', 'SRR1_1.fastq']
    assert find_prefixed_entries(entries, 'absent', sorted_input) == []
    assert find_prefixed_entries([], '', sorted_input) == []
