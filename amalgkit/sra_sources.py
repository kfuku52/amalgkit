import re
from urllib.parse import urlsplit, urlunsplit


ENA_SRA_LINK_COLUMN = 'ENA_SRA_Link'
DDBJ_SRA_LINK_COLUMN = 'DDBJ_SRA_Link'

_RUN_ACCESSION_PATTERN = re.compile(r'^(?:SRR|ERR|DRR)\d+$', re.IGNORECASE)
_DDBJ_RUN_ACCESSION_PATTERN = re.compile(r'^DRR\d+$', re.IGNORECASE)
_DDBJ_EXPERIMENT_ACCESSION_PATTERN = re.compile(r'^DRX\d+$', re.IGNORECASE)

_ENA_SRA_HOST = 'ftp.sra.ebi.ac.uk'
_ENA_SRA_ROOT = 'https://{}/vol1'.format(_ENA_SRA_HOST)
_DDBJ_DRA_PUBLIC_ROOT = 'https://ddbj.nig.ac.jp/public/ddbj_database/dra'


def normalize_accession_text(value):
    if value is None:
        return ''
    return str(value).strip().upper()


def build_ena_sra_download_url(run_accession):
    run_accession = normalize_accession_text(run_accession)
    if (_RUN_ACCESSION_PATTERN.match(run_accession) is None) or (len(run_accession) < 6):
        return ''
    run_prefix = run_accession[:3].lower()
    accession_prefix = run_accession[:6]
    return '{}/{}/{}/{}/{}.sra'.format(
        _ENA_SRA_ROOT,
        run_prefix,
        accession_prefix,
        run_accession,
        run_accession,
    )


def build_ddbj_sra_download_url(run_accession, experiment_accession):
    run_accession = normalize_accession_text(run_accession)
    experiment_accession = normalize_accession_text(experiment_accession)
    if _DDBJ_RUN_ACCESSION_PATTERN.match(run_accession) is None:
        return ''
    if (_DDBJ_EXPERIMENT_ACCESSION_PATTERN.match(experiment_accession) is None) or (len(experiment_accession) < 6):
        return ''
    experiment_prefix = experiment_accession[:6]
    return '{}/sra/ByExp/sra/DRX/{}/{}/{}/{}.sra'.format(
        _DDBJ_DRA_PUBLIC_ROOT,
        experiment_prefix,
        experiment_accession,
        run_accession,
        run_accession,
    )


def normalize_sra_download_url(source_name, source_url, run_accession='', experiment_accession=''):
    source_name = str(source_name or '').strip().upper()
    source_url = str(source_url or '').strip()
    run_accession = normalize_accession_text(run_accession)
    experiment_accession = normalize_accession_text(experiment_accession)
    if source_name == 'ENA':
        if source_url == '':
            return build_ena_sra_download_url(run_accession)
        if '://' not in source_url:
            source_url = 'https://' + source_url.lstrip('/')
        else:
            try:
                parsed_url = urlsplit(source_url)
            except ValueError:
                parsed_url = None
            if (
                parsed_url is not None
                and parsed_url.scheme.casefold() == 'ftp'
                and parsed_url.netloc.casefold() == _ENA_SRA_HOST
            ):
                source_url = urlunsplit(
                    ('https', _ENA_SRA_HOST, parsed_url.path, parsed_url.query, parsed_url.fragment)
                )
        source_url = source_url.rstrip('/')
        if (run_accession != '') and (not source_url.lower().endswith('.sra')):
            last_component = source_url.rsplit('/', 1)[-1]
            if last_component.upper() == run_accession:
                source_url = source_url + '/' + run_accession + '.sra'
        return source_url
    if source_name == 'DDBJ':
        if source_url == '':
            return build_ddbj_sra_download_url(
                run_accession=run_accession,
                experiment_accession=experiment_accession,
            )
        if '://' not in source_url:
            source_url = 'https://' + source_url.lstrip('/')
        return source_url
    return source_url
