import pytest

from io import BytesIO

from amalgkit.sra_sources import (
    build_ena_sra_download_url,
    fetch_ena_run_file_report,
    is_allowed_download_url,
    normalize_sra_download_url,
    parse_ena_run_file_report,
    read_bounded_response,
)


def test_normalize_sra_download_url_upgrades_exact_ena_ftp_host():
    source_url = 'ftp://ftp.sra.ebi.ac.uk/vol1/srr/SRR123/SRR123456/SRR123456.sra?download=1#archive'

    normalized = normalize_sra_download_url('ENA', source_url)

    assert normalized == (
        'https://ftp.sra.ebi.ac.uk/vol1/srr/SRR123/SRR123456/SRR123456.sra?download=1#archive'
    )


def test_normalize_sra_download_url_upgrades_case_insensitive_ena_ftp_host():
    source_url = 'FTP://FTP.SRA.EBI.AC.UK/vol1/srr/SRR123/SRR123456/SRR123456.sra'

    normalized = normalize_sra_download_url('ENA', source_url)

    assert normalized == 'https://ftp.sra.ebi.ac.uk/vol1/srr/SRR123/SRR123456/SRR123456.sra'


def test_normalize_sra_download_url_does_not_upgrade_ena_host_lookalikes():
    lookalike_urls = (
        'ftp://ftp.sra.ebi.ac.uk.attacker.example/archive.sra',
        'ftp://attacker.example/ftp.sra.ebi.ac.uk/archive.sra',
        'ftp://ftp.sra.ebi.ac.uk@attacker.example/archive.sra',
    )

    assert [normalize_sra_download_url('ENA', url) for url in lookalike_urls] == list(lookalike_urls)


def test_normalize_sra_download_url_preserves_malformed_url():
    source_url = 'ftp://[invalid-ipv6-host/archive.sra'

    assert normalize_sra_download_url('ENA', source_url) == source_url


def test_ena_sra_url_is_not_guessed_from_accession():
    assert build_ena_sra_download_url('SRR000001') == ''
    assert normalize_sra_download_url('ENA', '', run_accession='SRR000001') == ''


def test_ena_reported_url_is_not_extended_with_a_guessed_sra_filename():
    reported_url = 'ftp.sra.ebi.ac.uk/vol1/srr/SRR000/SRR000001'

    assert normalize_sra_download_url(
        'ENA',
        reported_url,
        run_accession='SRR000001',
    ) == 'https://ftp.sra.ebi.ac.uk/vol1/srr/SRR000/SRR000001'


def test_parse_ena_report_uses_reported_fastq_paths_when_sra_is_absent():
    report = (
        'run_accession\tsra_ftp\tfastq_ftp\n'
        'ERR6090701\t\t'
        'ftp.sra.ebi.ac.uk/vol1/fastq/ERR609/001/ERR6090701/ERR6090701_1.fastq.gz;'
        'ftp.sra.ebi.ac.uk/vol1/fastq/ERR609/001/ERR6090701/ERR6090701_2.fastq.gz\n'
    )

    parsed = parse_ena_run_file_report(report, 'ERR6090701')

    assert parsed['sra_urls'] == []
    assert parsed['fastq_urls'] == [
        'https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR609/001/ERR6090701/ERR6090701_1.fastq.gz',
        'https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR609/001/ERR6090701/ERR6090701_2.fastq.gz',
    ]


def test_fetch_ena_report_uses_official_filereport_endpoint():
    observed = {}
    response = (
        'run_accession\tsra_ftp\tfastq_ftp\n'
        'SRR000001\tftp.sra.ebi.ac.uk/path/run.sra\t\n'
    ).encode('utf-8')

    def fake_urlopen(url, timeout):
        observed['url'] = url
        observed['timeout'] = timeout
        return BytesIO(response)

    parsed = fetch_ena_run_file_report(
        'SRR000001',
        timeout=12,
        urlopen_fn=fake_urlopen,
    )

    assert 'www.ebi.ac.uk/ena/portal/api/filereport?' in observed['url']
    assert 'accession=SRR000001' in observed['url']
    assert 'fastq_md5' in observed['url'] and 'fastq_bytes' in observed['url']
    assert observed['timeout'] == 12.0
    assert parsed['sra_urls'] == ['https://ftp.sra.ebi.ac.uk/path/run.sra']


def test_ena_integrity_is_bound_to_each_reported_file():
    report = ('run_accession\tsra_ftp\tfastq_ftp\tfastq_md5\tfastq_bytes\n'
              'SRR1\t\tftp.sra.ebi.ac.uk/a.gz;ftp.sra.ebi.ac.uk/b.gz\t'
              + 'a' * 32 + ';' + 'b' * 32 + '\t100;200\n')
    result = parse_ena_run_file_report(report, 'SRR1')
    assert result['fastq_integrity'] == {
        'https://ftp.sra.ebi.ac.uk/a.gz': {'expected_md5': 'a' * 32, 'expected_bytes': 100},
        'https://ftp.sra.ebi.ac.uk/b.gz': {'expected_md5': 'b' * 32, 'expected_bytes': 200},
    }
    for changed in (report.replace('100;200', '100'), report.replace('100;200', '100;0'),
                    report.replace('b' * 32, 'invalid'), report.replace('/b.gz', '/a.gz')):
        with pytest.raises(ValueError, match='ENA FASTQ'):
            parse_ena_run_file_report(changed, 'SRR1')

def test_is_allowed_download_url_accepts_only_exact_https_endpoints():
    assert is_allowed_download_url('https://sra-pub-run-odp.s3.amazonaws.com/x.sra')
    assert is_allowed_download_url('https://ftp.sra.ebi.ac.uk/x.sra')
    assert is_allowed_download_url('https://storage.googleapis.com/b/x')
    assert is_allowed_download_url('https://ddbj.nig.ac.jp/x.sra')
    assert is_allowed_download_url('https://sra-downloadb.be-md.ncbi.nlm.nih.gov/x.sra')
    for bad in [
        'http://169.254.169.254/latest/meta-data/',
        'http://127.0.0.1/x',
        'http://192.168.1.5/x',
        'https://8.8.8.8/x',
        'ftp://ftp.sra.ebi.ac.uk/x.sra',
        'https://user:secret@ftp.sra.ebi.ac.uk/x.sra',
        'https://ftp.sra.ebi.ac.uk:444/x.sra',
        'http://evil.example/run.sra',
        'https://s3.amazonaws.com/foo.sra',
        'https://attacker.s3.amazonaws.com/x.sra',
        'https://sra-pub-run-attacker.s3.amazonaws.com/x.sra',
        '',
        'not a url',
    ]:
        assert not is_allowed_download_url(bad)

def test_fetch_ena_report_rejects_oversized_body():
    big = (b'run_accession\tsra_ftp\tfastq_ftp\nSRR1\t\t\n' + b'x' * (9 * 1024 * 1024))
    with pytest.raises(ValueError, match='byte metadata read limit'):
        fetch_ena_run_file_report(
            'SRR1', timeout=30,
            urlopen_fn=lambda url, timeout: BytesIO(big),
        )

def test_read_bounded_response_caps_bytes():
    from io import BytesIO
    with pytest.raises(ValueError, match='byte metadata read limit'):
        read_bounded_response(BytesIO(b'z' * 100), max_bytes=10, timeout=30)
    assert read_bounded_response(BytesIO(b'hello'), max_bytes=100, timeout=30) == b'hello'
    assert read_bounded_response(BytesIO(b'z' * 10), max_bytes=10, timeout=30) == b'z' * 10
