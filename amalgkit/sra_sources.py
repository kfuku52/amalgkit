import csv
import io
import re
import time
import urllib.parse
import urllib.request
from urllib.parse import urlsplit, urlunsplit


ENA_SRA_LINK_COLUMN = 'ENA_SRA_Link'
DDBJ_SRA_LINK_COLUMN = 'DDBJ_SRA_Link'

_RUN_ACCESSION_PATTERN = re.compile(r'^(?:SRR|ERR|DRR)\d+$', re.IGNORECASE)
_DDBJ_RUN_ACCESSION_PATTERN = re.compile(r'^DRR\d+$', re.IGNORECASE)
_DDBJ_EXPERIMENT_ACCESSION_PATTERN = re.compile(r'^DRX\d+$', re.IGNORECASE)

_ENA_SRA_HOST = 'ftp.sra.ebi.ac.uk'
_ENA_FILEREPORT_URL = 'https://www.ebi.ac.uk/ena/portal/api/filereport'
_DDBJ_DRA_PUBLIC_ROOT = 'https://ddbj.nig.ac.jp/public/ddbj_database/dra'

# Hosts amalgkit is allowed to contact for SRA/ENA/cloud downloads and for its
# own fixed metadata endpoints. Anything else (including private/link-local IP
# literals) is rejected as an SSRF vector.
ALLOWED_DOWNLOAD_HOSTS = frozenset({
    'ftp.sra.ebi.ac.uk',
    'ddbj.nig.ac.jp',
    'storage.googleapis.com',
    'sra-pub-run-odp.s3.amazonaws.com',
    'sra-downloadb.be-md.ncbi.nlm.nih.gov',
    'sra-download.be-md.ncbi.nlm.nih.gov',
    'www.ebi.ac.uk',            # ENA filereport endpoint (fixed)
    'trace.ncbi.nlm.nih.gov',   # NCBI trace XML endpoint (fixed)
})

def is_allowed_download_url(url):
    try:
        parsed = urlsplit(str(url or ''))
        host = parsed.hostname
        port = parsed.port
    except ValueError:
        return False
    if parsed.scheme.casefold() != 'https' or not host:
        return False
    host = host.lower().rstrip('.')
    return (
        parsed.username is None
        and parsed.password is None
        and port in (None, 443)
        and host in ALLOWED_DOWNLOAD_HOSTS
    )

def assert_allowed_download_url(url):
    if not is_allowed_download_url(url):
        raise ValueError('URL scheme or host is not an allowed SRA/ENA/cloud download endpoint.')

class _AllowedHostRedirectHandler(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        if not is_allowed_download_url(newurl):
            raise ValueError('Redirect to a non-allowed SRA/ENA/cloud download endpoint.')
        return super().redirect_request(req, fp, code, msg, headers, newurl)

def build_allowed_host_opener():
    """Return an opener that refuses redirects to non-allow-listed hosts."""
    return urllib.request.build_opener(_AllowedHostRedirectHandler)

MAX_METADATA_READ_BYTES = 8 * 1024 * 1024  # 8 MiB fail-closed cap for metadata bodies
_METADATA_READ_CHUNK = 64 * 1024

def read_bounded_response(response, max_bytes=MAX_METADATA_READ_BYTES, timeout=30):
    """Read a response in chunks, capped by byte count and a wall-clock deadline.

    Fails closed (raises) if either limit is exceeded so an unbounded remote
    body can never be fully slurped into memory.
    """
    deadline = float(timeout)
    started_at = time.monotonic()
    total = 0
    chunks = []
    max_bytes = int(max_bytes)
    if max_bytes < 0:
        raise ValueError('Metadata read byte limit must not be negative.')
    while True:
        if (time.monotonic() - started_at) >= deadline:
            raise TimeoutError(
                'Metadata response exceeded the {:.0f} sec read deadline.'.format(deadline)
            )
        chunk = response.read(min(_METADATA_READ_CHUNK, max_bytes - total + 1))
        if not chunk:
            break
        chunks.append(chunk)
        total += len(chunk)
        if total > max_bytes:
            raise ValueError(
                'Response exceeded the {:,} byte metadata read limit.'.format(max_bytes)
            )
        if (time.monotonic() - started_at) >= deadline:
            raise TimeoutError(
                'Metadata response exceeded the {:.0f} sec read deadline.'.format(deadline)
            )
    return b''.join(chunks)


def normalize_accession_text(value):
    if value is None:
        return ''
    return str(value).strip().upper()


def build_ena_sra_download_url(run_accession):
    """Return no URL: ENA file locations are not derivable from an accession.

    Keep this compatibility function because metadata readers imported it in
    older releases. Callers should use ``fetch_ena_run_file_report`` instead.
    """
    return ''


def normalize_ena_file_url(raw_url):
    raw_url = str(raw_url or '').strip()
    if raw_url == '':
        return ''
    if '://' not in raw_url:
        raw_url = 'https://' + raw_url.lstrip('/')
    try:
        parsed_url = urlsplit(raw_url)
    except ValueError:
        return raw_url
    if (
        parsed_url.scheme.casefold() == 'ftp'
        and parsed_url.netloc.casefold() == _ENA_SRA_HOST
    ):
        return urlunsplit(
            ('https', _ENA_SRA_HOST, parsed_url.path, parsed_url.query, parsed_url.fragment)
        )
    return raw_url


def _split_ena_file_urls(raw_value):
    urls = []
    seen = set()
    for raw_url in str(raw_value or '').split(';'):
        normalized = normalize_ena_file_url(raw_url)
        if (normalized == '') or (normalized in seen):
            continue
        urls.append(normalized)
        seen.add(normalized)
    return urls


def parse_ena_run_file_report(report_text, run_accession):
    run_accession = normalize_accession_text(run_accession)
    if _RUN_ACCESSION_PATTERN.fullmatch(run_accession) is None:
        raise ValueError('Invalid ENA run accession: {}'.format(run_accession))
    report_text = str(report_text or '')
    if report_text.strip() == '':
        return {'sra_urls': [], 'fastq_urls': []}
    reader = csv.DictReader(io.StringIO(report_text), delimiter='\t')
    required_fields = {'run_accession', 'sra_ftp', 'fastq_ftp'}
    if (reader.fieldnames is None) or (not required_fields.issubset(set(reader.fieldnames))):
        raise ValueError('ENA filereport response is missing required columns.')
    for row in reader:
        if normalize_accession_text(row.get('run_accession')) != run_accession:
            continue
        return {
            'sra_urls': _split_ena_file_urls(row.get('sra_ftp')),
            'fastq_urls': _split_ena_file_urls(row.get('fastq_ftp')),
        }
    return {'sra_urls': [], 'fastq_urls': []}


def fetch_ena_run_file_report(run_accession, timeout=30, urlopen_fn=None):
    run_accession = normalize_accession_text(run_accession)
    if _RUN_ACCESSION_PATTERN.fullmatch(run_accession) is None:
        raise ValueError('Invalid ENA run accession: {}'.format(run_accession))
    if urlopen_fn is None:
        urlopen_fn = urllib.request.urlopen
    query = urllib.parse.urlencode({
        'accession': run_accession,
        'result': 'read_run',
        'fields': 'run_accession,sra_ftp,fastq_ftp',
        'format': 'tsv',
    })
    report_url = '{}?{}'.format(_ENA_FILEREPORT_URL, query)
    with urlopen_fn(report_url, timeout=float(timeout)) as response:  # noqa: S310 - fixed ENA HTTPS endpoint
        report_bytes = read_bounded_response(response, timeout=float(timeout))
    report_text = report_bytes.decode('utf-8', errors='strict')
    return parse_ena_run_file_report(report_text=report_text, run_accession=run_accession)


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
            return ''
        source_url = normalize_ena_file_url(source_url)
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
