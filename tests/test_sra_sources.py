from io import BytesIO

from amalgkit.sra_sources import (
    build_ena_sra_download_url,
    fetch_ena_run_file_report,
    normalize_sra_download_url,
    parse_ena_run_file_report,
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
    assert observed['timeout'] == 12.0
    assert parsed['sra_urls'] == ['https://ftp.sra.ebi.ac.uk/path/run.sra']
