import os
import re
import subprocess
import sys
import time

from amalgkit.util import *

INDEX_BUILD_LOCK_POLL_SECONDS = 5
INDEX_BUILD_LOCK_TIMEOUT_SECONDS = 3600

INDEX_FASTA_SUFFIXES = ('.fa', '.fasta', '.fa.gz', '.fasta.gz')

def purge_existing_quant_outputs(sra_id, output_dir):
    stale_names = [
        sra_id + '_abundance.tsv',
        sra_id + '_run_info.json',
        sra_id + '_abundance.h5',
        'abundance.tsv',
        'run_info.json',
        'abundance.h5',
    ]
    for stale_name in stale_names:
        stale_path = os.path.join(output_dir, stale_name)
        if not os.path.exists(stale_path):
            continue
        if not os.path.isfile(stale_path):
            raise IsADirectoryError('Quant output path exists but is not a file: {}'.format(stale_path))
        os.remove(stale_path)

def quant_output_exists(sra_id, output_dir):
    abundance_path = os.path.join(output_dir, sra_id + '_abundance.tsv')
    run_info_path = os.path.join(output_dir, sra_id + '_run_info.json')
    has_abundance = os.path.isfile(abundance_path)
    has_run_info = os.path.isfile(run_info_path)
    if has_abundance and has_run_info:
        print('Output files detected: {}, {}'.format(abundance_path, run_info_path))
        return True
    if not has_abundance:
        print('Output file was not detected: {}'.format(abundance_path))
    if not has_run_info:
        print('Output file was not detected: {}'.format(run_info_path))
    return False

def resolve_nominal_length_for_kallisto(metadata, sra_id, sra_stat):
    nominal_length = sra_stat.get('nominal_length', numpy.nan)
    # Backward-compatible fallback for direct call_kallisto callers.
    if numpy.isnan(pandas.to_numeric(nominal_length, errors='coerce')):
        try:
            idx = get_metadata_row_index_by_run(metadata, sra_id)
            nominal_length = metadata.df.at[idx, 'nominal_length']
        except AssertionError:
            nominal_length = numpy.nan
    nominal_length = pandas.to_numeric(nominal_length, errors='coerce')
    if numpy.isnan(nominal_length) or nominal_length <= 0:
        print("Could not find nominal length in metadata. Assuming fragment length.")
        nominal_length = 200
    elif nominal_length < 200:
        print('Nominal length in metadata is unusually small ({}). Setting it to 200.'.format(nominal_length))
        nominal_length = 200
    print("Fragment length set to: {}".format(nominal_length))
    fragment_sd = nominal_length / 10
    print("Fragment length standard deviation set to: {}".format(fragment_sd))
    return nominal_length, fragment_sd


def build_kallisto_quant_command(args, in_files, lib_layout, output_dir, index, nominal_length=None, fragment_sd=None):
    if lib_layout == 'single':
        if len(in_files) != 1:
            txt = "Library layout: {} and expected 1 input file. " \
                  "Received {} input file[s]. Please check your inputs and metadata."
            raise ValueError(txt.format(lib_layout, len(in_files)))
        return [
            "kallisto", "quant", "--threads", str(args.threads), "--index", index, "-o", output_dir,
            "--single", "-l", str(nominal_length), "-s", str(fragment_sd), in_files[0],
        ]
    if lib_layout == 'paired':
        if len(in_files) != 2:
            txt = "Library layout: {} and expected 2 input files. " \
                  "Received {} input file[s]. Please check your inputs and metadata."
            raise ValueError(txt.format(lib_layout, len(in_files)))
        return [
            "kallisto", "quant", "--threads", str(args.threads), "-i", index, "-o", output_dir,
            in_files[0], in_files[1],
        ]
    raise ValueError("Unsupported library layout: {}. Expected 'single' or 'paired'.".format(lib_layout))


def rename_kallisto_outputs(output_dir, sra_id):
    output_specs = [
        ('run_info.json', sra_id + '_run_info.json', True),
        ('abundance.tsv', sra_id + '_abundance.tsv', True),
        ('abundance.h5', sra_id + '_abundance.h5', False),
    ]
    for src_name, dst_name, is_required in output_specs:
        src_path = os.path.join(output_dir, src_name)
        dst_path = os.path.join(output_dir, dst_name)
        if not os.path.exists(src_path):
            if is_required:
                raise FileNotFoundError('kallisto output file was not generated: {}'.format(src_path))
            continue
        if not os.path.isfile(src_path):
            raise IsADirectoryError('kallisto output path exists but is not a file: {}'.format(src_path))
        if os.path.exists(dst_path) and (not os.path.isfile(dst_path)):
            raise IsADirectoryError('kallisto renamed output path exists but is not a file: {}'.format(dst_path))
        os.replace(src_path, dst_path)


def call_kallisto(args, in_files, metadata, sra_stat, output_dir, index):
    sra_id = sra_stat['sra_id']
    lib_layout = sra_stat['layout']
    if lib_layout == 'single':
        print("Single end reads detected. Proceeding in single mode")
        nominal_length, fragment_sd = resolve_nominal_length_for_kallisto(metadata, sra_id, sra_stat)
        kallisto_cmd = build_kallisto_quant_command(
            args=args,
            in_files=in_files,
            lib_layout=lib_layout,
            output_dir=output_dir,
            index=index,
            nominal_length=nominal_length,
            fragment_sd=fragment_sd,
        )
    else:
        if lib_layout == 'paired':
            print("Paired-end reads detected. Running in paired read mode.")
        kallisto_cmd = build_kallisto_quant_command(
            args=args,
            in_files=in_files,
            lib_layout=lib_layout,
            output_dir=output_dir,
            index=index,
        )

    print('Command: {}'.format(' '.join(kallisto_cmd)))
    kallisto_out = subprocess.run(kallisto_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout_txt = kallisto_out.stdout.decode('utf8', errors='replace')
    stderr_txt = kallisto_out.stderr.decode('utf8', errors='replace')
    print('kallisto quant stdout:')
    print(stdout_txt)
    print('kallisto quant stderr:')
    print(stderr_txt)
    if kallisto_out.returncode != 0:
        sys.stderr.write("kallisto did not finish safely.\n")
        if 'Zero reads pseudoaligned' in stderr_txt:
            sys.stderr.write('No reads are mapped to the reference. This sample will be removed by `amalgkit curate`.')
        raise RuntimeError(
            'kallisto quant failed with exit code {} for {}.'.format(kallisto_out.returncode, sra_id)
        )

    # move output to results with unique name
    rename_kallisto_outputs(output_dir=output_dir, sra_id=sra_id)

    return kallisto_out


def check_kallisto_dependency():
    try:
        probe = subprocess.run(['kallisto', 'version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except FileNotFoundError as exc:
        raise FileNotFoundError('kallisto executable not found: kallisto') from exc
    if probe.returncode != 0:
        raise RuntimeError('kallisto dependency probe failed with exit code {}: kallisto version'.format(probe.returncode))


def list_getfastq_run_files(output_dir):
    return list_dir_entries(output_dir)


def list_dir_entries(path_dir):
    try:
        with os.scandir(path_dir) as entries:
            return {
                entry.name
                for entry in entries
                if entry.is_file()
            }
    except FileNotFoundError:
        return set()


def find_species_index_files(index_dir, sci_name, entries=None):
    if entries is None:
        entries = list_dir_entries(index_dir)
    matched = find_species_prefixed_entries(entries, sci_name)
    return [
        os.path.join(index_dir, entry)
        for entry in matched
        if os.path.isfile(os.path.join(index_dir, entry))
    ]


def find_species_fasta_files(path_fasta_dir, sci_name, entries=None):
    if os.path.exists(path_fasta_dir) and (not os.path.isdir(path_fasta_dir)):
        raise NotADirectoryError('Fasta path exists but is not a directory: {}'.format(path_fasta_dir))
    if entries is None:
        entries = list_dir_entries(path_fasta_dir)
    matched = [
        entry for entry in find_species_prefixed_entries(entries, sci_name)
        if entry.lower().endswith(INDEX_FASTA_SUFFIXES)
        and os.path.isfile(os.path.join(path_fasta_dir, entry))
    ]
    return [os.path.join(path_fasta_dir, entry) for entry in matched]


def check_layout_mismatch(sra_stat, output_dir, files=None):
    if files is None:
        files = list_getfastq_run_files(output_dir)
    if sra_stat['layout'] == 'paired':
        fastq_files = [
            f for f in find_run_prefixed_entries(files, sra_stat['sra_id'])
            if '.fastq' in f
        ]
        if len(fastq_files) == 1:
            sys.stderr.write('Single-end fastq was detected even though layout = {}\n'.format(sra_stat['layout']))
            sys.stderr.write('This sample will be treated as single-end sequencing.\n')
            sra_stat['layout'] = 'single'
    return sra_stat


def resolve_input_fastq_files(sra_stat, output_dir_getfastq, ext, files=None):
    if files is None:
        files = list_getfastq_run_files(output_dir_getfastq)
    sra_id = sra_stat['sra_id']
    if sra_stat['layout'] == 'paired':
        pair1 = sra_id + '_1' + ext
        pair2 = sra_id + '_2' + ext
        if (pair1 in files) and (pair2 in files):
            return [
                os.path.join(output_dir_getfastq, pair1),
                os.path.join(output_dir_getfastq, pair2),
            ]
    elif sra_stat['layout'] == 'single':
        single = sra_id + ext
        if single in files:
            return [os.path.join(output_dir_getfastq, single)]

    matched = sorted([
        f for f in find_run_prefixed_entries(files, sra_id)
        if f.endswith(ext)
    ])
    return [os.path.join(output_dir_getfastq, f) for f in matched]


def run_quant(args, metadata, sra_id, index):
    # make results directory, if not already there
    output_dir = os.path.join(args.out_dir, 'quant', sra_id)
    if os.path.exists(output_dir) and (not os.path.isdir(output_dir)):
        raise NotADirectoryError('Quant run output path exists but is not a directory: {}'.format(output_dir))
    os.makedirs(output_dir, exist_ok=True)
    if args.redo:
        purge_existing_quant_outputs(sra_id=sra_id, output_dir=output_dir)
    is_quant_output_available = quant_output_exists(sra_id, output_dir)
    if is_quant_output_available:
        if args.redo:
            print('The output will be overwritten. Set "--redo no" to not overwrite results.')
        else:
            print('Continued. The output will not be overwritten. If you want to overwrite the results, set "--redo yes".')
            return
    output_dir_getfastq = os.path.join(args.out_dir, 'getfastq', sra_id)
    if os.path.exists(output_dir_getfastq) and (not os.path.isdir(output_dir_getfastq)):
        raise NotADirectoryError(
            'getfastq run path exists but is not a directory: {}'.format(output_dir_getfastq)
        )
    sra_stat = get_sra_stat(sra_id, metadata, num_bp_per_sra=None)
    prefetched_run_files_map = getattr(args, '_prefetched_getfastq_run_files', None)
    if isinstance(prefetched_run_files_map, dict):
        run_files = prefetched_run_files_map.get(sra_id, None)
    else:
        run_files = None
    if run_files is None:
        run_files = list_getfastq_run_files(output_dir_getfastq)
    sra_stat = check_layout_mismatch(sra_stat, output_dir_getfastq, files=run_files)
    sra_stat['getfastq_sra_dir'] = output_dir_getfastq
    ext = get_newest_intermediate_file_extension(sra_stat, work_dir=output_dir_getfastq, files=run_files)
    if ext == '.safely_removed':
        print('These files have been safe-deleted. If you wish to re-obtain the .fastq file(s), run: getfastq --id ', sra_id, ' -w ', args.out_dir)
        print('Skipping.')
        return
    if ext == 'no_extension_found':
        sys.stderr.write('getfastq output not found in: {}, layout = {}\n'.format(sra_stat['getfastq_sra_dir'], sra_stat['layout']))
        txt = 'Exiting. If you wish to obtain the .fastq file(s), run: getfastq --id {}\n'
        sys.stderr.write(txt.format(sra_stat['sra_id']))
        sys.exit(1)
    in_files = resolve_input_fastq_files(sra_stat, output_dir_getfastq, ext, files=run_files)
    if len(in_files) == 0:
        # Refresh once in case files changed after initial snapshot.
        run_files = list_getfastq_run_files(output_dir_getfastq)
        in_files = resolve_input_fastq_files(sra_stat, output_dir_getfastq, ext, files=run_files)
    if not in_files:
        raise FileNotFoundError('{}: Fastq file not found. Check {}'.format(sra_id, output_dir_getfastq))
    print('Input fastq detected:', ', '.join(in_files))
    call_kallisto(args, in_files, metadata, sra_stat, output_dir, index)
    if (args.clean_fastq and quant_output_exists(sra_id, output_dir)):
        print('Safe-deleting getfastq files.', flush=True)
        for in_file in in_files:
            print('Output file detected. Safely removing fastq:', in_file)
            os.remove(in_file)
            with open(in_file + '.safely_removed', "w") as f:
                f.write("This fastq file was safely removed after `amalgkit quant`.")
    else:
        print('Skipping the deletion of getfastq files.', flush=True)

def _resolve_index_lock_options(args):
    lock_poll_seconds = int(getattr(args, 'index_lock_poll', INDEX_BUILD_LOCK_POLL_SECONDS))
    lock_timeout_seconds = int(getattr(args, 'index_lock_timeout', INDEX_BUILD_LOCK_TIMEOUT_SECONDS))
    if lock_poll_seconds <= 0:
        raise ValueError('--index_lock_poll must be > 0 (seconds).')
    if lock_timeout_seconds <= 0:
        raise ValueError('--index_lock_timeout must be > 0 (seconds).')
    return lock_poll_seconds, lock_timeout_seconds

def _resolve_index_dir(args):
    if args.index_dir is not None:
        index_dir = args.index_dir
    else:
        index_dir = os.path.join(args.out_dir, 'index')
    if (not os.path.exists(index_dir)) and args.build_index:
        os.makedirs(index_dir, exist_ok=True)
    if not os.path.exists(index_dir):
        raise FileNotFoundError("Could not find index folder at: {}".format(index_dir))
    if not os.path.isdir(index_dir):
        raise NotADirectoryError("Index path exists but is not a directory: {}".format(index_dir))
    return index_dir

def _assert_lock_path_is_regular_file(lock_path):
    if not os.path.lexists(lock_path):
        return
    if os.path.islink(lock_path) or (not os.path.isfile(lock_path)):
        raise IsADirectoryError('Index lock path exists but is not a file: {}'.format(lock_path))

def _acquire_index_lock(lock_path):
    _assert_lock_path_is_regular_file(lock_path)
    try:
        fd = os.open(lock_path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError:
        return False
    with os.fdopen(fd, 'w') as lock_handle:
        lock_handle.write('{}\n'.format(os.getpid()))
    return True

def _wait_for_existing_builder(index_path, lock_path, sci_name, lock_poll_seconds, lock_timeout_seconds):
    _assert_lock_path_is_regular_file(lock_path)
    if not os.path.lexists(lock_path):
        return False
    print('Another process is building index for {}. Waiting for lock release: {}'.format(sci_name, lock_path), flush=True)
    wait_start = time.time()
    while os.path.lexists(lock_path):
        _assert_lock_path_is_regular_file(lock_path)
        elapsed = time.time() - wait_start
        if elapsed > lock_timeout_seconds:
            txt = 'Timed out after {:,} sec while waiting for index lock: {}\n'
            txt += 'Remove stale lock if no builder is running and rerun quant.'
            raise TimeoutError(txt.format(lock_timeout_seconds, lock_path))
        time.sleep(lock_poll_seconds)
    if os.path.exists(index_path):
        print('Detected completed index after waiting: {}'.format(index_path), flush=True)
        return True
    return False

def _normalize_species_identifier(text):
    normalized = str(text).strip()
    if normalized == '':
        return ''
    normalized = re.sub(r'\s+', '_', normalized)
    normalized = re.sub(r'_+', '_', normalized)
    return normalized

def _find_single_index_file(index_dir, sci_name, entries=None):
    candidates = []
    raw = str(sci_name).strip()
    normalized = _normalize_species_identifier(raw)
    for prefix in [raw, normalized]:
        if (prefix != '') and (prefix not in candidates):
            candidates.append(prefix)
        prefix_no_dot = prefix.replace('.', '')
        if (prefix_no_dot != '') and (prefix_no_dot != prefix) and (prefix_no_dot not in candidates):
            candidates.append(prefix_no_dot)
    normalized_parts = [part for part in normalized.split('_') if part != '']
    if len(normalized_parts) > 2:
        fallback = '_'.join(normalized_parts[:2])
        if fallback not in candidates:
            candidates.append(fallback)
        fallback_no_dot = fallback.replace('.', '')
        if (fallback_no_dot != fallback) and (fallback_no_dot not in candidates):
            candidates.append(fallback_no_dot)
    for prefix in candidates:
        index_files = find_species_index_files(index_dir=index_dir, sci_name=prefix, entries=entries)
        if len(index_files) > 1:
            raise ValueError(
                "Found multiple index files for species. Please make sure there is only one index file for this species.")
        if len(index_files) == 1:
            index_file = index_files[0]
            if prefix != sci_name:
                print(
                    "Kallisto index fallback prefix '{}' was used for species '{}'.".format(prefix, sci_name),
                    flush=True,
                )
            print("Kallisto index file found: {}".format(index_file), flush=True)
            return index_file
    return None

def _resolve_prefetched_index_entries(args, index_dir):
    prefetched_index_entries = getattr(args, '_prefetched_index_entries', None)
    prefetched_index_entries_sorted = getattr(args, '_prefetched_index_entries_sorted', None)
    prefetched_index_dir = getattr(args, '_prefetched_index_dir', None)
    is_matching_prefetch_dir = isinstance(prefetched_index_dir, str) and (os.path.realpath(index_dir) == prefetched_index_dir)
    if (not is_matching_prefetch_dir):
        return None
    if isinstance(prefetched_index_entries_sorted, (list, tuple)):
        return prefetched_index_entries_sorted
    if isinstance(prefetched_index_entries, (set, list, tuple)):
        return prefetched_index_entries
    return None

def _resolve_single_fasta_file(args, sci_name):
    if args.fasta_dir == 'inferred':
        path_fasta_dir = os.path.join(args.out_dir, 'fasta')
    else:
        path_fasta_dir = args.fasta_dir
    prefetched_fasta_entries = getattr(args, '_prefetched_fasta_entries', None)
    prefetched_fasta_entries_sorted = getattr(args, '_prefetched_fasta_entries_sorted', None)
    prefetched_fasta_dir = getattr(args, '_prefetched_fasta_dir', None)
    prefetched_entries = None
    if (
        isinstance(prefetched_fasta_entries_sorted, (list, tuple))
        and isinstance(prefetched_fasta_dir, str)
        and (os.path.realpath(path_fasta_dir) == prefetched_fasta_dir)
    ):
        prefetched_entries = prefetched_fasta_entries_sorted
    elif (
        isinstance(prefetched_fasta_entries, (set, list, tuple))
        and isinstance(prefetched_fasta_dir, str)
        and (os.path.realpath(path_fasta_dir) == prefetched_fasta_dir)
    ):
        prefetched_entries = prefetched_fasta_entries
    fasta_files = find_species_fasta_files(
        path_fasta_dir=path_fasta_dir,
        sci_name=sci_name,
        entries=prefetched_entries,
    )
    if len(fasta_files) > 1:
        txt = "Found multiple reference fasta files for this species: {}\n"
        txt += "Please make sure there is only one index file for this species.\n{}"
        raise ValueError(txt.format(sci_name, ', '.join(fasta_files)))
    if len(fasta_files) == 0:
        txt = "Could not find reference fasta file for this species: {}\n".format(sci_name)
        txt += 'If the reference fasta file is correctly placed, the column "scientific_name" of the --metadata file may need to be edited.'
        raise FileNotFoundError(txt)
    return fasta_files[0]

def _build_kallisto_index(index_path, fasta_file, sci_name):
    print('Reference fasta file found: {}'.format(fasta_file), flush=True)
    print('Building index: {}'.format(index_path), flush=True)
    kallisto_build_cmd = ["kallisto", "index", "-i", index_path, fasta_file]
    print('Command: {}'.format(' '.join(kallisto_build_cmd)), flush=True)
    index_out = subprocess.run(kallisto_build_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout_txt = index_out.stdout.decode('utf8', errors='replace')
    stderr_txt = index_out.stderr.decode('utf8', errors='replace')
    print('kallisto index stdout:')
    print(stdout_txt)
    print('kallisto index stderr:')
    print(stderr_txt)
    if index_out.returncode != 0:
        raise RuntimeError("kallisto index failed for {}.".format(sci_name))
    if not os.path.isfile(index_path):
        raise RuntimeError("Index file was not generated: {}".format(index_path))

def _store_prefetched_entries(args, attr_prefix, entries, path_dir):
    setattr(args, '_{}_entries'.format(attr_prefix), entries)
    if isinstance(entries, set):
        setattr(args, '_{}_entries_sorted'.format(attr_prefix), sorted(entries))
    else:
        setattr(args, '_{}_entries_sorted'.format(attr_prefix), None)
    if isinstance(path_dir, str):
        setattr(args, '_{}_dir'.format(attr_prefix), os.path.realpath(path_dir))
    else:
        setattr(args, '_{}_dir'.format(attr_prefix), None)

def get_index(args, sci_name):
    lock_poll_seconds, lock_timeout_seconds = _resolve_index_lock_options(args)
    index_dir = _resolve_index_dir(args)
    index_path = os.path.join(index_dir, sci_name + '.idx')
    lock_path = os.path.join(index_dir, '.{}.idx.lock'.format(sci_name))
    prefetched_index_lookup_entries = _resolve_prefetched_index_entries(args, index_dir)
    use_prefetched_index = prefetched_index_lookup_entries is not None

    while True:
        if use_prefetched_index:
            index = _find_single_index_file(
                index_dir=index_dir,
                sci_name=sci_name,
                entries=prefetched_index_lookup_entries,
            )
            use_prefetched_index = False
        else:
            index = _find_single_index_file(index_dir=index_dir, sci_name=sci_name)
        if index is not None:
            return index

        if not args.build_index:
            sys.stderr.write('No index file was found in: {}\n'.format(index_dir))
            sys.stderr.write('Try --fasta_dir PATH and --build_index yes\n')
            raise FileNotFoundError("Could not find index file.")

        if _wait_for_existing_builder(
            index_path=index_path,
            lock_path=lock_path,
            sci_name=sci_name,
            lock_poll_seconds=lock_poll_seconds,
            lock_timeout_seconds=lock_timeout_seconds,
        ):
            continue

        if not _acquire_index_lock(lock_path):
            continue

        try:
            # Another process may have finished between lock release and acquisition.
            index = _find_single_index_file(index_dir=index_dir, sci_name=sci_name)
            if index is not None:
                return index

            print("--build_index set. Building index for {}".format(sci_name))
            fasta_file = _resolve_single_fasta_file(args, sci_name)
            _build_kallisto_index(index_path=index_path, fasta_file=fasta_file, sci_name=sci_name)
            print("Kallisto index file found: {}".format(index_path), flush=True)
            return index_path
        finally:
            if os.path.lexists(lock_path):
                _assert_lock_path_is_regular_file(lock_path)
                os.remove(lock_path)


def run_quant_for_sra(args, metadata, sra_id, sci_name):
    print('')
    print('Species: {}'.format(sci_name))
    print('SRA Run ID: {}'.format(sra_id))
    sci_name = _normalize_species_identifier(sci_name)
    index_cache = getattr(args, '_resolved_index_cache', None)
    if isinstance(index_cache, dict) and (sci_name in index_cache):
        index = index_cache[sci_name]
        print('Using pre-resolved index: {}'.format(index))
    else:
        print('Looking for index folder in ', args.out_dir)
        index = get_index(args, sci_name)
    run_quant(args, metadata, sra_id, index)

def pre_resolve_species_indices(args, tasks):
    # Resolve one index per species once to avoid redundant index lookups across SRA runs.
    species_list = sorted(set([_normalize_species_identifier(sci_name) for _, sci_name in tasks]))
    if len(species_list) == 0:
        return {}
    prefetched_fasta_entries = None
    prefetched_fasta_dir = None
    if getattr(args, 'build_index', False):
        if args.fasta_dir == 'inferred':
            path_fasta_dir = os.path.join(args.out_dir, 'fasta')
        else:
            path_fasta_dir = args.fasta_dir
        if os.path.exists(path_fasta_dir) and (not os.path.isdir(path_fasta_dir)):
            raise NotADirectoryError('Fasta path exists but is not a directory: {}'.format(path_fasta_dir))
        prefetched_fasta_entries = list_dir_entries(path_fasta_dir)
        prefetched_fasta_dir = path_fasta_dir
    _store_prefetched_entries(
        args=args,
        attr_prefix='prefetched_fasta',
        entries=prefetched_fasta_entries,
        path_dir=prefetched_fasta_dir,
    )
    if args.index_dir is not None:
        index_dir = args.index_dir
    else:
        index_dir = os.path.join(args.out_dir, 'index')
    if os.path.exists(index_dir) and (not os.path.isdir(index_dir)):
        raise NotADirectoryError('Index path exists but is not a directory: {}'.format(index_dir))
    if (not os.path.exists(index_dir)) and args.build_index:
        os.makedirs(index_dir, exist_ok=True)
    prefetched_index_entries = list_dir_entries(index_dir)
    _store_prefetched_entries(
        args=args,
        attr_prefix='prefetched_index',
        entries=prefetched_index_entries,
        path_dir=index_dir,
    )
    print('Resolving kallisto index for {:,} species.'.format(len(species_list)), flush=True)
    requested_jobs = getattr(args, 'internal_jobs', 'auto')
    if is_auto_parallel_option(requested_jobs):
        requested_jobs = resolve_detected_cpu_count()
    try:
        max_workers = min(max(1, int(requested_jobs)), len(species_list))
    except (TypeError, ValueError):
        max_workers = 1

    if max_workers <= 1:
        resolved = dict()
        for sci_name in species_list:
            print('Pre-resolving index for species: {}'.format(sci_name), flush=True)
            resolved[sci_name] = get_index(args, sci_name)
        return resolved

    print('Pre-resolving indices with {:,} parallel jobs.'.format(max_workers), flush=True)
    resolved_by_species, failures = run_tasks_with_optional_threads(
        task_items=species_list,
        task_fn=lambda sci_name: get_index(args, sci_name),
        max_workers=max_workers,
    )
    if failures:
        details = '; '.join(['{}: {}'.format(sci_name, err) for sci_name, err in failures])
        raise RuntimeError(
            'Failed to pre-resolve index for {}/{} species. {}'.format(
                len(failures),
                len(species_list),
                details,
            )
        )
    resolved = dict()
    for sci_name in species_list:
        resolved[sci_name] = resolved_by_species[sci_name]
    return resolved


def prefetch_getfastq_run_files(args, tasks):
    run_ids = sorted(set([sra_id for sra_id, _ in tasks]))
    if len(run_ids) == 0:
        return {}
    getfastq_root = os.path.join(args.out_dir, 'getfastq')
    prefetched = {}
    for run_id in run_ids:
        run_dir = os.path.join(getfastq_root, run_id)
        try:
            with os.scandir(run_dir) as run_entries:
                prefetched[run_id] = {
                    run_entry.name
                    for run_entry in run_entries
                    if run_entry.is_file()
                }
        except (FileNotFoundError, NotADirectoryError):
            continue
    return prefetched


def build_quant_tasks(metadata):
    required_columns = ['run', 'scientific_name']
    missing_columns = [col for col in required_columns if col not in metadata.df.columns]
    if len(missing_columns) > 0:
        raise ValueError('Missing required metadata column(s) for quant: {}'.format(', '.join(missing_columns)))
    runs = metadata.df['run'].fillna('').astype(str).str.strip()
    species = metadata.df['scientific_name'].fillna('').astype(str).str.strip()
    metadata.df['run'] = runs
    metadata.df['scientific_name'] = species
    missing_species_runs = []
    missing_run_count = 0
    duplicate_runs = []
    seen_runs = set()
    tasks = []
    for run_id, sci_name in zip(runs.tolist(), species.tolist()):
        if run_id == '':
            missing_run_count += 1
            continue
        if sci_name == '':
            missing_species_runs.append(run_id)
            continue
        if run_id in seen_runs:
            duplicate_runs.append(run_id)
            continue
        seen_runs.add(run_id)
        tasks.append((run_id, sci_name))
    if missing_run_count > 0:
        raise ValueError('Missing run ID in metadata for {:,} row(s).'.format(missing_run_count))
    if len(missing_species_runs) > 0:
        raise ValueError('Missing scientific_name in metadata for run(s): {}'.format(', '.join(missing_species_runs)))
    if len(duplicate_runs) > 0:
        duplicate_runs = list(dict.fromkeys(duplicate_runs))
        raise ValueError('Duplicate run ID in metadata for run(s): {}'.format(', '.join(duplicate_runs)))
    if len(tasks) == 0:
        raise ValueError('No valid run/scientific_name entries were found in metadata.')
    return tasks


def quant_main(args):
    threads, jobs, _ = resolve_thread_worker_allocation(
        requested_threads=getattr(args, 'threads', 'auto'),
        requested_workers=getattr(args, 'internal_jobs', 'auto'),
        internal_cpu_budget=getattr(args, 'internal_cpu_budget', 'auto'),
        worker_option_name='internal_jobs',
        context='quant:',
        disable_workers=(getattr(args, 'batch', None) is not None),
    )
    args.threads = threads
    args.internal_jobs = jobs
    out_dir = os.path.realpath(args.out_dir)
    if os.path.exists(out_dir) and (not os.path.isdir(out_dir)):
        raise NotADirectoryError('Output path exists but is not a directory: {}'.format(out_dir))
    quant_dir = os.path.join(out_dir, 'quant')
    if os.path.exists(quant_dir) and (not os.path.isdir(quant_dir)):
        raise NotADirectoryError('Quant path exists but is not a directory: {}'.format(quant_dir))
    args.out_dir = out_dir
    check_kallisto_dependency()
    metadata = load_metadata(args)
    tasks = build_quant_tasks(metadata)
    setattr(args, '_prefetched_getfastq_run_files', prefetch_getfastq_run_files(args, tasks))
    setattr(args, '_resolved_index_cache', pre_resolve_species_indices(args, tasks))
    if (jobs == 1) or (len(tasks) <= 1):
        for sra_id, sci_name in tasks:
            run_quant_for_sra(args, metadata, sra_id, sci_name)
        return

    max_workers = min(jobs, len(tasks))
    print('Running quant for {:,} SRA runs with {:,} parallel jobs.'.format(len(tasks), max_workers), flush=True)
    _, failures = run_tasks_with_optional_threads(
        task_items=tasks,
        task_fn=lambda task: run_quant_for_sra(args, metadata, task[0], task[1]),
        max_workers=max_workers,
    )
    if failures:
        details = '; '.join(['{}: {}'.format(task[0], err) for task, err in failures])
        raise RuntimeError('quant failed for {}/{} SRA runs. {}'.format(len(failures), len(tasks), details))
