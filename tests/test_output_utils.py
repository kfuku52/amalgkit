import pandas
import pytest

from amalgkit import output_utils
from amalgkit.table_io import read_identifier_tsv


@pytest.mark.parametrize('suffix', ['.tsv', '.tsv.gz', '.tsv.bz2', '.tsv.xz', '.tar.gz'])
def test_atomic_dataframe_roundtrip_preserves_compression_and_identifiers(tmp_path, suffix):
    output = tmp_path / ('metadata' + suffix)
    frame = pandas.DataFrame({'run': ['0001', '1', 'NA'], 'count': [1, 2, 3]})
    output_utils.atomic_write_dataframe(frame, output, sep='\t', index=False)
    pandas.testing.assert_frame_equal(read_identifier_tsv(output), frame)


def test_atomic_output_cleans_temporary_file_when_mode_probe_fails(tmp_path, monkeypatch):
    output = tmp_path / 'metadata.tsv'
    output.write_text('previous data\n')

    def fail_mode_probe(*args, **kwargs):
        raise OSError('mode probe failed')

    monkeypatch.setattr(output_utils, 'get_default_creation_mode', fail_mode_probe)
    with pytest.raises(OSError, match='mode probe failed'):
        with output_utils.atomic_output_path(output):
            pytest.fail('writer must not start after preparation failed')
    assert output.read_text() == 'previous data\n'
    assert list(tmp_path.iterdir()) == [output]
