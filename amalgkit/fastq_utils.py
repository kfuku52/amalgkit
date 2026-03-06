import gzip
import os


DEFAULT_FASTQ_CHUNK_SIZE = 16 * 1024 * 1024


def is_gzip_fastq_path(path_fastq):
    path_fastq_lower = str(path_fastq).lower()
    return path_fastq_lower.endswith(('.fq.gz', '.fastq.gz'))


def open_fastq_text(path_fastq):
    if is_gzip_fastq_path(path_fastq):
        return gzip.open(path_fastq, 'rt')
    return open(path_fastq, 'rt')


def open_fastq_binary(path_fastq):
    if is_gzip_fastq_path(path_fastq):
        return gzip.open(path_fastq, 'rb')
    return open(path_fastq, 'rb')


def count_fastq_lines(path_fastq, chunk_size=DEFAULT_FASTQ_CHUNK_SIZE):
    num_lines = 0
    with open_fastq_binary(path_fastq) as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            num_lines += chunk.count(b'\n')
    return num_lines


def sequence_line_length(line_bytes):
    seq_len = len(line_bytes)
    if seq_len == 0:
        return 0
    if line_bytes.endswith(b'\n'):
        seq_len -= 1
        if (seq_len > 0) and line_bytes.endswith(b'\r\n'):
            seq_len -= 1
    return seq_len


def sample_fastq_reads(path_fastq, max_reads=None):
    num_reads = 0
    total_bases = 0
    total_record_chars = 0
    reached_eof = True
    with open_fastq_binary(path_fastq) as handle:
        while True:
            line1 = handle.readline()
            if line1 == b'':
                break
            line2 = handle.readline()
            line3 = handle.readline()
            line4 = handle.readline()
            if (line2 == b'') or (line3 == b'') or (line4 == b''):
                raise ValueError('Malformed FASTQ (record truncated): {}'.format(path_fastq))
            num_reads += 1
            total_bases += sequence_line_length(line2)
            total_record_chars += (len(line1) + len(line2) + len(line3) + len(line4))
            if (max_reads is not None) and (num_reads >= max_reads):
                reached_eof = (handle.readline() == b'')
                break
    avg_len = int(total_bases / num_reads) if num_reads > 0 else 0
    return num_reads, avg_len, total_record_chars, reached_eof


def count_fastq_records(path_fastq, chunk_size=DEFAULT_FASTQ_CHUNK_SIZE, warning_writer=None):
    num_lines = count_fastq_lines(path_fastq, chunk_size=chunk_size)
    if ((num_lines % 4) != 0) and callable(warning_writer):
        warning_writer('FASTQ line count is not divisible by 4 and may be truncated: {}\n'.format(path_fastq))
    return num_lines // 4


def count_fastq_records_and_bases(path_fastq, warning_writer=None):
    num_records = 0
    num_bases = 0
    with open_fastq_binary(path_fastq) as handle:
        while True:
            line1 = handle.readline()
            if line1 == b'':
                break
            line2 = handle.readline()
            line3 = handle.readline()
            line4 = handle.readline()
            if (line2 == b'') or (line3 == b'') or (line4 == b''):
                if callable(warning_writer):
                    warning_writer('FASTQ record seems truncated while counting bases: {}\n'.format(path_fastq))
                break
            num_records += 1
            num_bases += len(line2.rstrip(b'\r\n'))
    return num_records, num_bases


def parse_seqkit_stats_rows(stdout_txt, empty_data_message='seqkit stats output did not include data rows.'):
    lines = [line.rstrip('\n') for line in (stdout_txt or '').splitlines() if line.strip() != '']
    if len(lines) < 2:
        raise RuntimeError('seqkit stats output was empty.')
    header_fields = lines[0].split('\t')
    rows = []
    for line in lines[1:]:
        value_fields = line.split('\t')
        row = {}
        for idx, field_name in enumerate(header_fields):
            row[field_name] = value_fields[idx] if idx < len(value_fields) else ''
        rows.append(row)
    if len(rows) == 0:
        raise RuntimeError(empty_data_message)
    return rows


def parse_seqkit_stats_numeric(row, key, path_fastq):
    if key not in row:
        raise RuntimeError('seqkit stats output is missing "{}" for {}'.format(key, path_fastq))
    raw = str(row.get(key, '')).strip().replace(',', '')
    if raw == '':
        raise RuntimeError('seqkit stats output has an empty "{}" value for {}'.format(key, path_fastq))
    try:
        return int(float(raw))
    except (TypeError, ValueError) as exc:
        raise RuntimeError('Failed to parse seqkit stats "{}" for {}'.format(key, path_fastq)) from exc


def parse_seqkit_stats_row_num_reads_avg_len(row, path_fastq):
    if ('num_seqs' not in row) or ('avg_len' not in row):
        raise RuntimeError('seqkit stats output was missing required fields for {}'.format(path_fastq))
    try:
        num_reads = int(float(row['num_seqs']))
        avg_len = int(round(float(row['avg_len'])))
    except (TypeError, ValueError) as exc:
        raise RuntimeError('Failed to parse seqkit stats numeric fields for {}.'.format(path_fastq)) from exc
    return num_reads, avg_len


def parse_seqkit_stats_row_records_and_bases(row, path_fastq):
    num_records = parse_seqkit_stats_numeric(row=row, key='num_seqs', path_fastq=path_fastq)
    if ('sum_len' in row) and (str(row.get('sum_len', '')).strip() not in ['', 'NA', 'nan', 'NaN']):
        num_bases = parse_seqkit_stats_numeric(row=row, key='sum_len', path_fastq=path_fastq)
    else:
        avg_len = parse_seqkit_stats_numeric(row=row, key='avg_len', path_fastq=path_fastq)
        num_bases = int(num_records * avg_len)
    return num_records, num_bases


def map_seqkit_stats_rows(stdout_txt, requested_paths, row_parser, missing_message):
    deduplicated_paths = []
    seen = set()
    for path_fastq in requested_paths:
        if path_fastq is None:
            continue
        abs_path = os.path.abspath(path_fastq)
        if abs_path in seen:
            continue
        seen.add(abs_path)
        deduplicated_paths.append(path_fastq)
    if len(deduplicated_paths) == 0:
        return {}
    rows = parse_seqkit_stats_rows(stdout_txt=stdout_txt)
    stats_by_path = {}
    fallback_idx = 0
    fallback_paths = [os.path.abspath(path_fastq) for path_fastq in deduplicated_paths]
    for row in rows:
        row_path = None
        for key in ['file', 'path', 'filename']:
            raw_path = str(row.get(key, '')).strip()
            if raw_path != '':
                row_path = os.path.abspath(raw_path)
                break
        if (row_path is None) and (fallback_idx < len(fallback_paths)):
            row_path = fallback_paths[fallback_idx]
            fallback_idx += 1
        if row_path is None:
            continue
        stats_by_path[row_path] = row_parser(row=row, path_fastq=row_path)
    missing_paths = [path_fastq for path_fastq in fallback_paths if path_fastq not in stats_by_path]
    if len(missing_paths) > 0:
        raise RuntimeError(missing_message.format(missing_paths))
    return stats_by_path
