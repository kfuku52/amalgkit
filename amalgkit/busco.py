import gzip
import os
import re
import shlex
import shutil
import subprocess
import sys

import pandas

from amalgkit.util import (
    acquire_exclusive_lock,
    find_species_prefixed_entries,
    load_metadata,
    run_tasks_with_optional_threads,
    resolve_thread_worker_allocation,
)


REQUIRED_COLUMNS = [
    'busco_id',
    'status',
    'sequence',
    'score',
    'length',
    'orthodb_url',
    'description',
]


FASTA_SUFFIXES = ('.fa', '.fasta', '.fa.gz', '.fasta.gz')


def validate_busco_metadata_columns(metadata, required_columns):
    missing = [col for col in required_columns if col not in metadata.df.columns]
    if len(missing) > 0:
        raise ValueError(
            'Missing required metadata column(s) for busco: {}'.format(', '.join(missing))
        )


def normalize_busco_columns(df):
    def keyify(name):
        return re.sub(r'[^a-z0-9]', '', str(name).lower())

    lookup = {keyify(col): col for col in df.columns}
    required_keys = {
        'buscoid': 'busco_id',
        'status': 'status',
        'sequence': 'sequence',
        'score': 'score',
        'length': 'length',
        'orthodburl': 'orthodb_url',
        'description': 'description',
    }
    missing = [target for key, target in required_keys.items() if key not in lookup]
    if missing:
        raise ValueError('Missing required BUSCO columns: {}'.format(', '.join(missing)))
    renamed = {lookup[key]: target for key, target in required_keys.items()}
    df = df.rename(columns=renamed)
    return df


def _is_header_like_busco_row(values):
    normalized_values = {
        re.sub(r'[^a-z0-9]', '', str(value).lower())
        for value in values
    }
    required_keys = {
        'buscoid',
        'status',
        'sequence',
        'score',
        'length',
        'orthodburl',
        'description',
    }
    return required_keys.issubset(normalized_values)


def normalize_busco_table(src_path, dest_path):
    header = None
    src_path_lower = src_path.lower()
    open_func = gzip.open if src_path_lower.endswith('.gz') else open
    compression = 'gzip' if src_path_lower.endswith('.gz') else None
    with open_func(src_path, 'rt') as f:
        for line in f:
            if not line.startswith('#'):
                continue
            candidate = line.lstrip('#').strip()
            candidate_fields = candidate.split('\t')
            if _is_header_like_busco_row(candidate_fields):
                header = candidate_fields
    if header:
        df = pandas.read_table(
            src_path,
            sep='\t',
            header=None,
            comment='#',
            names=header,
            dtype=str,
            low_memory=False,
            compression=compression,
        )
    else:
        df = pandas.read_table(
            src_path,
            sep='\t',
            header=None,
            comment='#',
            dtype=str,
            low_memory=False,
            compression=compression,
        )
        if df.shape[0] > 0 and _is_header_like_busco_row(df.iloc[0].tolist()):
            inferred_header = [str(value) for value in df.iloc[0].tolist()]
            df = df.iloc[1:, :].reset_index(drop=True)
            if df.shape[0] == 0:
                raise ValueError('BUSCO table is empty: {}'.format(src_path))
            df.columns = inferred_header
    if df.shape[0] == 0:
        raise ValueError('BUSCO table is empty: {}'.format(src_path))
    if not header and df.shape[1] == len(REQUIRED_COLUMNS):
        df.columns = REQUIRED_COLUMNS
    df = normalize_busco_columns(df)
    df = df.loc[:, REQUIRED_COLUMNS]
    with open(dest_path, 'w') as f:
        f.write('# Busco id\tStatus\tSequence\tScore\tLength\tOrthoDB url\tDescription\n')
        df.to_csv(f, sep='\t', index=False, header=False)


def find_full_table(output_dir):
    found = None
    for root, _, files in os.walk(output_dir):
        for filename in files:
            filename_lower = filename.lower()
            if ('full_table' not in filename_lower):
                continue
            if not (filename_lower.endswith('.tsv') or filename_lower.endswith('.tsv.gz')):
                continue
            path = os.path.join(root, filename)
            if found is None:
                found = path
                continue
            if path != found:
                raise ValueError('Multiple BUSCO full_table files detected: {}, {}'.format(found, path))
    if found is None:
        raise FileNotFoundError('No full_table.tsv found under: {}'.format(output_dir))
    return found


def resolve_fasta_dir(args):
    if args.fasta_dir == 'inferred':
        return os.path.join(args.out_dir, 'fasta')
    return os.path.realpath(args.fasta_dir)


def list_fasta_filenames(fasta_dir):
    files = []
    with os.scandir(fasta_dir) as entries:
        for entry in entries:
            if entry.is_file():
                files.append(entry.name)
    return files


def _species_prefix_candidates(sci_name):
    raw = str(sci_name).strip()
    if raw == '':
        return []
    normalized = re.sub(r'\s+', '_', raw)
    normalized = re.sub(r'_+', '_', normalized)
    candidates = [normalized]
    raw_underscored = raw.replace(' ', '_')
    if raw_underscored != normalized:
        candidates.append(raw_underscored)
    normalized_no_dot = normalized.replace('.', '')
    if normalized_no_dot != normalized:
        candidates.append(normalized_no_dot)
    normalized_parts = [part for part in normalized.split('_') if part != '']
    if len(normalized_parts) > 2:
        fallback = '_'.join(normalized_parts[0:2])
        if fallback not in candidates:
            candidates.append(fallback)
        fallback_no_dot = fallback.replace('.', '')
        if fallback_no_dot not in candidates:
            candidates.append(fallback_no_dot)
    return candidates


def resolve_species_fasta(sci_name, fasta_dir, fasta_filenames=None):
    sci_name = str(sci_name).strip()
    if fasta_filenames is None:
        fasta_filenames = list_fasta_filenames(fasta_dir)
    candidates = _species_prefix_candidates(sci_name)
    primary = candidates[0] if len(candidates) > 0 else sci_name.replace(' ', '_')
    for prefix in candidates:
        fasta_files = []
        for filename in find_species_prefixed_entries(fasta_filenames, prefix):
            if not filename.lower().endswith(FASTA_SUFFIXES):
                continue
            fasta_files.append(os.path.join(fasta_dir, filename))
        fasta_files = sorted(set(fasta_files))
        if len(fasta_files) > 1:
            raise ValueError('Found multiple reference fasta files for {}: {}'.format(prefix, ', '.join(fasta_files)))
        if len(fasta_files) == 1:
            if prefix != primary:
                print(
                    "Reference fasta fallback prefix '{}' was used for species '{}'.".format(prefix, sci_name),
                    flush=True,
                )
            return fasta_files[0]
    raise FileNotFoundError('Could not find reference fasta file for {} in: {}'.format(primary, fasta_dir))


def select_tool(args):
    tool = args.tool
    if tool == 'auto':
        if shutil.which(args.compleasm_exe):
            return 'compleasm'
        if shutil.which(args.busco_exe):
            return 'busco'
        raise FileNotFoundError('Neither compleasm nor busco was found on PATH.')
    if tool == 'compleasm':
        if not shutil.which(args.compleasm_exe):
            raise FileNotFoundError('compleasm executable not found: {}'.format(args.compleasm_exe))
        return 'compleasm'
    if tool == 'busco':
        if not shutil.which(args.busco_exe):
            raise FileNotFoundError('busco executable not found: {}'.format(args.busco_exe))
        return 'busco'
    raise ValueError('Unknown tool: {}'.format(tool))


def ensure_clean_dir(path, redo):
    if os.path.lexists(path):
        if os.path.islink(path):
            is_dir_link = os.path.isdir(path)
            if redo:
                os.remove(path)
            elif not is_dir_link:
                raise NotADirectoryError('BUSCO output path exists but is not a directory: {}'.format(path))
        else:
            if not os.path.isdir(path):
                raise NotADirectoryError('BUSCO output path exists but is not a directory: {}'.format(path))
            if redo:
                shutil.rmtree(path)
    os.makedirs(path, exist_ok=True)


def run_command(cmd):
    print('Command: {}'.format(' '.join(cmd)), flush=True)
    out = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    print(out.stdout.decode('utf8', errors='replace'))
    print(out.stderr.decode('utf8', errors='replace'))
    if out.returncode != 0:
        raise RuntimeError('Command failed with exit code {}: {}'.format(out.returncode, ' '.join(cmd)))

def resolve_download_dir(args):
    raw_dir = getattr(args, 'download_dir', 'inferred')
    if raw_dir is None:
        return os.path.join(os.path.realpath(args.out_dir), 'downloads')
    normalized = str(raw_dir).strip()
    if normalized.lower() in ['', 'inferred']:
        return os.path.join(os.path.realpath(args.out_dir), 'downloads')
    return os.path.realpath(normalized)


def _sanitize_lock_suffix(text):
    normalized = re.sub(r'[^A-Za-z0-9._-]+', '_', str(text).strip())
    if normalized == '':
        return 'lineage'
    return normalized


def _dir_has_entries(path_dir):
    if not os.path.isdir(path_dir):
        return False
    try:
        with os.scandir(path_dir) as entries:
            for _ in entries:
                return True
    except OSError:
        return False
    return False


def has_busco_lineage_cache(download_path, lineage):
    lineage = str(lineage).strip()
    if lineage == '':
        return False
    candidates = [
        os.path.join(download_path, lineage),
        os.path.join(download_path, 'lineages', lineage),
    ]
    for candidate in candidates:
        if _dir_has_entries(candidate):
            return True
        if os.path.isfile(candidate):
            return True
        if os.path.isfile(candidate + '.tar.gz'):
            return True
        if os.path.isfile(candidate + '.tgz'):
            return True
    return False


def run_busco(fasta_path, sci_name, output_root, args, extra_args):
    out_name = sci_name.replace(' ', '_')
    download_path = resolve_download_dir(args)
    if os.path.exists(download_path) and (not os.path.isdir(download_path)):
        raise NotADirectoryError(
            'BUSCO download path exists but is not a directory: {}'.format(download_path)
        )
    os.makedirs(download_path, exist_ok=True)
    cmd = [
        args.busco_exe,
        '-i', fasta_path,
        '-o', out_name,
        '-l', args.lineage,
        '-m', 'transcriptome',
        '--out_path', output_root,
        '--download_path', download_path,
        '--cpu', str(args.threads),
    ]
    if args.redo:
        cmd.append('--force')
    cmd.extend(extra_args)
    if has_busco_lineage_cache(download_path=download_path, lineage=args.lineage):
        run_command(cmd)
        return os.path.join(output_root, out_name)

    lock_filename = '.busco_{}.download.lock'.format(_sanitize_lock_suffix(args.lineage))
    lock_path = os.path.join(download_path, lock_filename)
    cache_ready_from_other_process = False
    with acquire_exclusive_lock(lock_path=lock_path, lock_label='BUSCO lineage download'):
        if has_busco_lineage_cache(download_path=download_path, lineage=args.lineage):
            cache_ready_from_other_process = True
        else:
            run_command(cmd)
            return os.path.join(output_root, out_name)
    if cache_ready_from_other_process:
        run_command(cmd)
    return os.path.join(output_root, out_name)


def run_compleasm(fasta_path, sci_name, output_root, args, extra_args):
    out_dir = os.path.join(output_root, sci_name.replace(' ', '_'))
    ensure_clean_dir(out_dir, args.redo)
    cmd = [
        args.compleasm_exe,
        'run',
        '-i', fasta_path,
        '-o', out_dir,
        '-l', args.lineage,
        '-t', str(args.threads),
        '--mode', 'transcriptome',
    ]
    cmd.extend(extra_args)
    run_command(cmd)
    return out_dir


def collect_species(args, metadata):
    if args.fasta is not None:
        if args.species is None:
            raise ValueError('--species is required when --fasta is provided.')
        species_name = str(args.species).strip()
        if species_name == '':
            raise ValueError('--species must not be empty when --fasta is provided.')
        fasta_path = os.path.realpath(args.fasta)
        if not os.path.exists(fasta_path):
            raise FileNotFoundError('FASTA file not found: {}'.format(fasta_path))
        if not os.path.isfile(fasta_path):
            raise IsADirectoryError('FASTA path exists but is not a file: {}'.format(fasta_path))
        return [species_name], {species_name: fasta_path}
    fasta_dir = resolve_fasta_dir(args)
    if not os.path.exists(fasta_dir):
        raise FileNotFoundError('FASTA directory not found: {}'.format(fasta_dir))
    if not os.path.isdir(fasta_dir):
        raise NotADirectoryError('FASTA path exists but is not a directory: {}'.format(fasta_dir))
    validate_busco_metadata_columns(
        metadata=metadata,
        required_columns=['scientific_name'],
    )
    species = []
    for species_name in metadata.df.loc[:, 'scientific_name'].tolist():
        if pandas.isna(species_name):
            continue
        normalized = str(species_name).strip()
        if normalized == '':
            continue
        species.append(normalized)
    # Preserve order while de-duplicating.
    species = list(dict.fromkeys(species))
    if len(species) == 0:
        raise ValueError('No valid scientific_name entries were found in metadata.')
    fasta_filenames = list_fasta_filenames(fasta_dir)
    fasta_map = {}
    for sp in species:
        fasta_map[sp] = resolve_species_fasta(sp, fasta_dir, fasta_filenames=fasta_filenames)
    return species, fasta_map


def process_species_busco(sp, fasta_path, busco_dir, tool, args, extra_args):
    print('Processing species: {}'.format(sp), flush=True)
    if tool == 'busco':
        output_root = busco_dir
        ensure_clean_dir(os.path.join(output_root, sp.replace(' ', '_')), args.redo)
        tool_out_dir = run_busco(fasta_path, sp, output_root, args, extra_args)
    else:
        tool_out_dir = run_compleasm(fasta_path, sp, busco_dir, args, extra_args)
    full_table = find_full_table(tool_out_dir)
    out_table = os.path.join(busco_dir, sp.replace(' ', '_') + '_busco.tsv')
    normalize_busco_table(full_table, out_table)
    print('BUSCO table written: {}'.format(out_table), flush=True)
    return out_table


def load_metadata_if_needed_for_busco(args):
    if args.fasta is not None:
        return None
    return load_metadata(args)


def run_busco_species_jobs(species, fasta_map, busco_dir, tool, args, extra_args, species_jobs):
    if (species_jobs == 1) or (len(species) <= 1):
        for sp in species:
            process_species_busco(
                sp=sp,
                fasta_path=fasta_map[sp],
                busco_dir=busco_dir,
                tool=tool,
                args=args,
                extra_args=extra_args,
            )
        return

    max_workers = min(species_jobs, len(species))
    print('Running BUSCO for {:,} species with {:,} parallel jobs.'.format(len(species), max_workers), flush=True)
    _, failures = run_tasks_with_optional_threads(
        task_items=species,
        task_fn=lambda sp: process_species_busco(
            sp=sp,
            fasta_path=fasta_map[sp],
            busco_dir=busco_dir,
            tool=tool,
            args=args,
            extra_args=extra_args,
        ),
        max_workers=max_workers,
    )
    if failures:
        details = '; '.join(['{}: {}'.format(sp, err) for sp, err in failures])
        raise RuntimeError('BUSCO failed for {}/{} species. {}'.format(len(failures), len(species), details))


def busco_main(args):
    lineage = str(getattr(args, 'lineage', '')).strip()
    if lineage == '':
        raise ValueError('--lineage is required.')
    args.lineage = lineage
    threads, species_jobs, _ = resolve_thread_worker_allocation(
        requested_threads=getattr(args, 'threads', 'auto'),
        requested_workers=getattr(args, 'internal_jobs', 'auto'),
        internal_cpu_budget=getattr(args, 'internal_cpu_budget', 'auto'),
        worker_option_name='internal_jobs',
        context='busco:',
    )
    args.threads = threads
    args.internal_jobs = species_jobs
    out_dir = os.path.realpath(args.out_dir)
    if os.path.exists(out_dir) and (not os.path.isdir(out_dir)):
        raise NotADirectoryError('Output path exists but is not a directory: {}'.format(out_dir))
    busco_dir = os.path.join(out_dir, 'busco')
    if os.path.exists(busco_dir) and (not os.path.isdir(busco_dir)):
        raise NotADirectoryError('BUSCO path exists but is not a directory: {}'.format(busco_dir))
    tool = select_tool(args)
    extra_args = shlex.split(args.tool_args) if args.tool_args else []
    metadata = load_metadata_if_needed_for_busco(args)
    species, fasta_map = collect_species(args, metadata)
    os.makedirs(busco_dir, exist_ok=True)
    run_busco_species_jobs(
        species=species,
        fasta_map=fasta_map,
        busco_dir=busco_dir,
        tool=tool,
        args=args,
        extra_args=extra_args,
        species_jobs=species_jobs,
    )
