import os
import re
import subprocess
import sys
import time
from bisect import bisect_left
from concurrent.futures import ThreadPoolExecutor, as_completed

from amalgkit.util import *

INDEX_BUILD_LOCK_POLL_SECONDS = 5
INDEX_BUILD_LOCK_TIMEOUT_SECONDS = 3600

INDEX_FASTA_SUFFIXES = ('.fa', '.fasta', '.fa.gz', '.fasta.gz')

def quant_output_exists(sra_id, output_dir):
    out_path = os.path.join(output_dir, sra_id + '_abundance.tsv')
    is_output_present = os.path.exists(out_path)
    if is_output_present:
        print('Output file detected: {}'.format(out_path))
        return True
    else:
        print('Output file was not detected: {}'.format(out_path))
        return False

def call_kallisto(args, in_files, metadata, sra_stat, output_dir, index):
    sra_id = sra_stat['sra_id']
    lib_layout = sra_stat['layout']
    kallisto_cmd = ''
    if lib_layout == 'single':
        print("Single end reads detected. Proceeding in single mode")
        if len(in_files) != 1:
            txt = "Library layout: {} and expected 1 input file. " \
                  "Received {} input file[s]. Please check your inputs and metadata."
            raise ValueError(txt.format(lib_layout, len(in_files)))
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
        kallisto_cmd = ["kallisto", "quant", "--threads", str(args.threads), "--index", index, "-o", output_dir,
                        "--single", "-l", str(nominal_length), "-s", str(fragment_sd), in_files[0]]
    elif lib_layout == 'paired':
        if len(in_files) != 2:
            txt = "Library layout: {} and expected 2 input files. " \
                  "Received {} input file[s]. Please check your inputs and metadata."
            raise ValueError(txt.format(lib_layout, len(in_files)))
        print("Paired-end reads detected. Running in paired read mode.")
        kallisto_cmd = ["kallisto", "quant", "--threads", str(args.threads), "-i", index, "-o",
                        output_dir, in_files[0], in_files[1]]
    else:
        raise ValueError("Unsupported library layout: {}. Expected 'single' or 'paired'.".format(lib_layout))

    print('Command: {}'.format(' '.join(kallisto_cmd)))
    kallisto_out = subprocess.run(kallisto_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    print('kallisto quant stdout:')
    print(kallisto_out.stdout.decode('utf8'))
    print('kallisto quant stderr:')
    print(kallisto_out.stderr.decode('utf8'))
    if kallisto_out.returncode != 0:
        sys.stderr.write("kallisto did not finish safely.\n")
        if 'Zero reads pseudoaligned' in kallisto_out.stderr.decode('utf8'):
            sys.stderr.write('No reads are mapped to the reference. This sample will be removed by `amalgkit curate`.')

    # move output to results with unique name
    try:
        os.rename(os.path.join(output_dir, "run_info.json"), os.path.join(output_dir, sra_id + "_run_info.json"))
        os.rename(os.path.join(output_dir, "abundance.tsv"), os.path.join(output_dir, sra_id + "_abundance.tsv"))
        os.rename(os.path.join(output_dir, "abundance.h5"), os.path.join(output_dir, sra_id + "_abundance.h5"))
    except FileNotFoundError:
        pass

    return kallisto_out


def check_kallisto_dependency():
    try:
        subprocess.run(['kallisto', '-h'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except FileNotFoundError:
        raise Exception("kallisto is not installed.")


def list_getfastq_run_files(output_dir):
    try:
        return set(os.listdir(output_dir))
    except FileNotFoundError:
        return set()


def list_dir_entries(path_dir):
    try:
        return set(os.listdir(path_dir))
    except FileNotFoundError:
        return set()


def find_prefixed_entries(entries, prefix):
    if isinstance(entries, (list, tuple)):
        left = bisect_left(entries, prefix)
        right = bisect_left(entries, prefix + '\uffff')
        matched = []
        for i in range(left, right):
            entry = entries[i]
            if entry.startswith(prefix):
                matched.append(entry)
        return matched
    return sorted([
        entry for entry in entries
        if entry.startswith(prefix)
    ])


def find_species_index_files(index_dir, sci_name, entries=None):
    if entries is None:
        entries = list_dir_entries(index_dir)
    matched = find_prefixed_entries(entries, sci_name)
    return [os.path.join(index_dir, entry) for entry in matched]


def find_species_fasta_files(path_fasta_dir, sci_name, entries=None):
    if entries is None:
        entries = list_dir_entries(path_fasta_dir)
    matched = [
        entry for entry in find_prefixed_entries(entries, sci_name)
        if entry.endswith(INDEX_FASTA_SUFFIXES)
    ]
    return [os.path.join(path_fasta_dir, entry) for entry in matched]


def check_layout_mismatch(sra_stat, output_dir, files=None):
    if files is None:
        files = list_getfastq_run_files(output_dir)
    if sra_stat['layout'] == 'paired':
        fastq_files = [
            f for f in files
            if f.startswith(sra_stat['sra_id']) and ('.fastq' in f)
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
        f for f in files
        if f.startswith(sra_id) and f.endswith(ext)
    ])
    return [os.path.join(output_dir_getfastq, f) for f in matched]


def run_quant(args, metadata, sra_id, index):
    # make results directory, if not already there
    output_dir = os.path.join(args.out_dir, 'quant', sra_id)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    is_quant_output_available = quant_output_exists(sra_id, output_dir)
    if is_quant_output_available:
        if args.redo:
            print('The output will be overwritten. Set "--redo no" to not overwrite results.')
        else:
            print('Continued. The output will not be overwritten. If you want to overwrite the results, set "--redo yes".')
            return
    output_dir_getfastq = os.path.join(args.out_dir, 'getfastq', sra_id)
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
    assert in_files, '{}: Fastq file not found. Check {}'.format(sra_id, output_dir)
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

def get_index(args, sci_name):
    lock_poll_seconds = int(getattr(args, 'index_lock_poll', INDEX_BUILD_LOCK_POLL_SECONDS))
    lock_timeout_seconds = int(getattr(args, 'index_lock_timeout', INDEX_BUILD_LOCK_TIMEOUT_SECONDS))
    if lock_poll_seconds <= 0:
        raise ValueError('--index_lock_poll must be > 0 (seconds).')
    if lock_timeout_seconds <= 0:
        raise ValueError('--index_lock_timeout must be > 0 (seconds).')

    def find_index_files(index_dir, sci_name):
        return find_species_index_files(index_dir=index_dir, sci_name=sci_name)

    def acquire_index_lock(lock_path):
        try:
            fd = os.open(lock_path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        except FileExistsError:
            return False
        with os.fdopen(fd, 'w') as lock_handle:
            lock_handle.write('{}\n'.format(os.getpid()))
        return True

    def wait_for_existing_builder(index_path, lock_path):
        if not os.path.exists(lock_path):
            return False
        print('Another process is building index for {}. Waiting for lock release: {}'.format(sci_name, lock_path), flush=True)
        wait_start = time.time()
        while os.path.exists(lock_path):
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

    if args.index_dir is not None:
        index_dir = args.index_dir
    else:
        index_dir = os.path.join(args.out_dir, 'index')
    if not os.path.exists(index_dir) and args.build_index:
            os.makedirs(index_dir, exist_ok=True)
    if not os.path.exists(index_dir):
        raise FileNotFoundError("Could not find index folder at: {}".format(index_dir))

    index_path = os.path.join(index_dir, sci_name + '.idx')
    lock_path = os.path.join(index_dir, '.{}.idx.lock'.format(sci_name))
    prefetched_index_entries = getattr(args, '_prefetched_index_entries', None)
    prefetched_index_entries_sorted = getattr(args, '_prefetched_index_entries_sorted', None)
    prefetched_index_dir = getattr(args, '_prefetched_index_dir', None)
    use_prefetched_index = isinstance(prefetched_index_dir, str) and (os.path.realpath(index_dir) == prefetched_index_dir)
    if use_prefetched_index and isinstance(prefetched_index_entries_sorted, (list, tuple)):
        prefetched_index_lookup_entries = prefetched_index_entries_sorted
    else:
        prefetched_index_lookup_entries = prefetched_index_entries

    while True:
        if use_prefetched_index and isinstance(prefetched_index_lookup_entries, (set, list, tuple)):
            index = find_species_index_files(index_dir=index_dir, sci_name=sci_name, entries=prefetched_index_lookup_entries)
            use_prefetched_index = False
        else:
            index = find_index_files(index_dir, sci_name)
        if len(index) > 1:
            raise ValueError(
                "Found multiple index files for species. Please make sure there is only one index file for this species.")
        elif len(index) == 1:
            index = index[0]
            print("Kallisto index file found: {}".format(index), flush=True)
            return index

        if not args.build_index:
            sys.stderr.write('No index file was found in: {}\n'.format(index_dir))
            sys.stderr.write('Try --fasta_dir PATH and --build_index yes\n')
            raise FileNotFoundError("Could not find index file.")

        if wait_for_existing_builder(index_path=index_path, lock_path=lock_path):
            continue

        if not acquire_index_lock(lock_path):
            continue

        try:
            # Another process may have finished between lock release and acquisition.
            index = find_index_files(index_dir, sci_name)
            if len(index) > 1:
                raise ValueError(
                    "Found multiple index files for species. Please make sure there is only one index file for this species.")
            if len(index) == 1:
                index = index[0]
                print("Kallisto index file found: {}".format(index), flush=True)
                return index

            print("--build_index set. Building index for {}".format(sci_name))
            if (args.fasta_dir=='inferred'):
                path_fasta_dir = os.path.join(args.out_dir, 'fasta')
            else:
                path_fasta_dir = args.fasta_dir
            prefetched_fasta_entries = getattr(args, '_prefetched_fasta_entries', None)
            prefetched_fasta_entries_sorted = getattr(args, '_prefetched_fasta_entries_sorted', None)
            prefetched_fasta_dir = getattr(args, '_prefetched_fasta_dir', None)
            if (
                isinstance(prefetched_fasta_entries_sorted, (list, tuple))
                and isinstance(prefetched_fasta_dir, str)
                and (os.path.realpath(path_fasta_dir) == prefetched_fasta_dir)
            ):
                fasta_files = find_species_fasta_files(
                    path_fasta_dir=path_fasta_dir,
                    sci_name=sci_name,
                    entries=prefetched_fasta_entries_sorted,
                )
            elif (
                isinstance(prefetched_fasta_entries, (set, list, tuple))
                and isinstance(prefetched_fasta_dir, str)
                and (os.path.realpath(path_fasta_dir) == prefetched_fasta_dir)
            ):
                fasta_files = find_species_fasta_files(
                    path_fasta_dir=path_fasta_dir,
                    sci_name=sci_name,
                    entries=prefetched_fasta_entries,
                )
            else:
                fasta_files = find_species_fasta_files(path_fasta_dir=path_fasta_dir, sci_name=sci_name)
            if len(fasta_files) > 1:
                txt = "Found multiple reference fasta files for this species: {}\n"
                txt += "Please make sure there is only one index file for this species.\n{}"
                raise ValueError(txt.format(sci_name, ', '.join(fasta_files)))
            elif len(fasta_files) == 0:
                txt = "Could not find reference fasta file for this species: {}\n".format(sci_name)
                txt += 'If the reference fasta file is correctly placed, the column "scientific_name" of the --metadata file may need to be edited.'
                raise FileNotFoundError(txt)
            fasta_file = fasta_files[0]
            print('Reference fasta file found: {}'.format(fasta_file), flush=True)
            print('Building index: {}'.format(index_path), flush=True)
            kallisto_build_cmd = ["kallisto", "index", "-i", index_path, fasta_file]
            print('Command: {}'.format(' '.join(kallisto_build_cmd)), flush=True)
            index_out = subprocess.run(kallisto_build_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            print('kallisto index stdout:')
            print(index_out.stdout.decode('utf8'))
            print('kallisto index stderr:')
            print(index_out.stderr.decode('utf8'))
            if index_out.returncode != 0:
                raise RuntimeError("kallisto index failed for {}.".format(sci_name))
            if not os.path.exists(index_path):
                raise RuntimeError("Index file was not generated: {}".format(index_path))
            print("Kallisto index file found: {}".format(index_path), flush=True)
            return index_path
        finally:
            if os.path.exists(lock_path):
                os.remove(lock_path)


def run_quant_for_sra(args, metadata, sra_id, sci_name):
    print('')
    print('Species: {}'.format(sci_name))
    print('SRA Run ID: {}'.format(sra_id))
    sci_name = sci_name.replace(" ", "_")
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
    species_list = sorted(set([sci_name.replace(' ', '_') for _, sci_name in tasks]))
    if len(species_list) == 0:
        return {}
    prefetched_fasta_entries = None
    prefetched_fasta_dir = None
    if getattr(args, 'build_index', False):
        if args.fasta_dir == 'inferred':
            path_fasta_dir = os.path.join(args.out_dir, 'fasta')
        else:
            path_fasta_dir = args.fasta_dir
        prefetched_fasta_entries = list_dir_entries(path_fasta_dir)
        prefetched_fasta_dir = os.path.realpath(path_fasta_dir)
    setattr(args, '_prefetched_fasta_entries', prefetched_fasta_entries)
    if isinstance(prefetched_fasta_entries, set):
        setattr(args, '_prefetched_fasta_entries_sorted', sorted(prefetched_fasta_entries))
    else:
        setattr(args, '_prefetched_fasta_entries_sorted', None)
    setattr(args, '_prefetched_fasta_dir', prefetched_fasta_dir)
    if args.index_dir is not None:
        index_dir = args.index_dir
    else:
        index_dir = os.path.join(args.out_dir, 'index')
    prefetched_index_entries = list_dir_entries(index_dir)
    setattr(args, '_prefetched_index_entries', prefetched_index_entries)
    if isinstance(prefetched_index_entries, set):
        setattr(args, '_prefetched_index_entries_sorted', sorted(prefetched_index_entries))
    else:
        setattr(args, '_prefetched_index_entries_sorted', None)
    setattr(args, '_prefetched_index_dir', os.path.realpath(index_dir))
    print('Resolving kallisto index for {:,} species.'.format(len(species_list)), flush=True)
    resolved = dict()
    for sci_name in species_list:
        print('Pre-resolving index for species: {}'.format(sci_name), flush=True)
        resolved[sci_name] = get_index(args, sci_name)
    return resolved


def prefetch_getfastq_run_files(args, tasks):
    run_ids = set([sra_id for sra_id, _ in tasks])
    if len(run_ids) == 0:
        return {}
    getfastq_root = os.path.join(args.out_dir, 'getfastq')
    prefetched = {}
    try:
        with os.scandir(getfastq_root) as entries:
            for entry in entries:
                if (not entry.is_dir()) or (entry.name not in run_ids):
                    continue
                try:
                    with os.scandir(entry.path) as run_entries:
                        prefetched[entry.name] = {run_entry.name for run_entry in run_entries}
                except (FileNotFoundError, NotADirectoryError):
                    continue
    except FileNotFoundError:
        return {}
    return prefetched


def quant_main(args):
    jobs = int(getattr(args, 'jobs', 1))
    if jobs <= 0:
        raise ValueError('--jobs must be > 0.')
    check_kallisto_dependency()
    metadata = load_metadata(args)
    tasks = list(zip(metadata.df['run'].tolist(), metadata.df['scientific_name'].tolist()))
    setattr(args, '_prefetched_getfastq_run_files', prefetch_getfastq_run_files(args, tasks))
    setattr(args, '_resolved_index_cache', pre_resolve_species_indices(args, tasks))
    if (jobs == 1) or (len(tasks) <= 1):
        for sra_id, sci_name in tasks:
            run_quant_for_sra(args, metadata, sra_id, sci_name)
        return

    max_workers = min(jobs, len(tasks))
    print('Running quant for {:,} SRA runs with {:,} parallel jobs.'.format(len(tasks), max_workers), flush=True)
    failures = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(run_quant_for_sra, args, metadata, sra_id, sci_name): sra_id
            for sra_id, sci_name in tasks
        }
        for future in as_completed(futures):
            sra_id = futures[future]
            try:
                future.result()
            except Exception as exc:
                failures.append((sra_id, exc))
    if failures:
        details = '; '.join(['{}: {}'.format(sra_id, err) for sra_id, err in failures])
        raise RuntimeError('quant failed for {}/{} SRA runs. {}'.format(len(failures), len(tasks), details))
