from amalgkit.sra_sources import normalize_sra_download_url


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
