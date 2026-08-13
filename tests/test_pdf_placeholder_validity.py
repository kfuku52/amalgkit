"""Regression test: the integration-test PDF placeholders must be parseable.

Previously tests/conftest.py stubbed plot helpers by writing a two-line
``b'%PDF-1.4\\n% amalgkit test placeholder\\n'`` blob that is *not* a valid
PDF (it has no object table and no ``%%EOF`` trailer).  Integration tests
that only asserted ``.is_file()`` were therefore false-green: a file existed
but could never be opened by a PDF reader.  The stub now writes a genuinely
parseable minimal PDF and self-asserts its own output.
"""

from io import BytesIO

from tests.conftest import _assert_parseable_pdf, _valid_pdf_bytes, _write_valid_pdf


def test_placeholder_bytes_are_parseable_pdf():
    data = _valid_pdf_bytes()
    _assert_parseable_pdf(data)


def test_placeholder_written_to_path_is_parseable(tmp_path):
    out_path = tmp_path / 'sample.pdf'
    returned = _write_valid_pdf(out_path)
    assert returned == str(out_path)
    assert out_path.is_file()
    _assert_parseable_pdf(out_path.read_bytes())


def test_placeholder_written_to_file_object_is_parseable():
    buffer = BytesIO()
    _write_valid_pdf(buffer)
    _assert_parseable_pdf(buffer.getvalue())
