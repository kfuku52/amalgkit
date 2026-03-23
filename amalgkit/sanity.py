import datetime
import json
import numpy as np
import os
import pandas
import re

from amalgkit.exceptions import AmalgkitExit
from amalgkit.fastq_utils import open_fastq_binary
from amalgkit.metadata_utils import (
    get_newest_intermediate_file_extension,
    get_sra_stat,
    load_metadata,
)
from amalgkit.parallel_utils import (
    is_auto_parallel_option,
    resolve_detected_cpu_count,
    run_tasks_with_optional_threads,
    validate_positive_int_option,
)
from amalgkit.prefix_utils import find_run_prefixed_entries, find_species_prefixed_entries


def list_duplicates(seq):
    seen = set()
    seen_add = seen.add
    # adds all elements it doesn't know yet to seen and all other to seen_twice
    seen_twice = set(x for x in seq if x in seen or seen_add(x))
    # turn the set into a list (as requested)
    return list(seen_twice)


def _scan_target_run_dirs(root_path, target_runs):
    run_files_map = {}
    non_dir_runs = set()
    try:
        with os.scandir(root_path) as entries:
            for entry in entries:
                run_id = entry.name
                if run_id not in target_runs:
                    continue
                if not entry.is_dir():
                    non_dir_runs.add(run_id)
                    continue
                try:
                    with os.scandir(entry.path) as run_entries:
                        run_files_map[run_id] = {
                            run_entry.name
                            for run_entry in run_entries
                            if run_entry.is_file()
                        }
                except (FileNotFoundError, NotADirectoryError):
                    pass
    except FileNotFoundError:
        pass
    return run_files_map, non_dir_runs

def _normalize_species_prefix(species):
    normalized = re.sub(r'\s+', '_', str(species).strip())
    normalized = re.sub(r'_+', '_', normalized)
    return normalized


def _get_species_prefix_candidates(species_or_prefix):
    prefix = str(species_or_prefix).strip()
    if prefix == '':
        return []
    normalized = _normalize_species_prefix(prefix)
    candidates = [normalized]
    no_dot = normalized.replace(".", "")
    if no_dot != normalized:
        candidates.append(no_dot)
    # Preserve order while de-duplicating.
    return list(dict.fromkeys(candidates))

def _get_species_fallback_prefix(species):
    parts = str(species).strip().split()
    if len(parts) <= 2:
        return None
    return _normalize_species_prefix(' '.join(parts[0:2]))

def _should_log_per_item(args, num_items):
    if bool(getattr(args, 'quiet', False)):
        return False
    verbose_runs = int(getattr(args, 'verbose_runs', 20))
    if verbose_runs < 0:
        return True
    return num_items <= verbose_runs

def _print_log_mode(args, num_items, singular_label='run', plural_label=None):
    if _should_log_per_item(args, num_items):
        return
    verbose_runs = int(getattr(args, 'verbose_runs', 20))
    if plural_label is None:
        plural_label = singular_label + 's'
    print('Per-{} logs suppressed for {:,} {} (--verbose_runs={}, --quiet={}).'.format(
        singular_label,
        num_items,
        plural_label,
        verbose_runs,
        bool(getattr(args, 'quiet', False)),
    ))

def _write_unavailable_items(output_dir, filename, item_ids, label):
    if not item_ids:
        return ''
    if os.path.exists(output_dir) and (not os.path.isdir(output_dir)):
        raise NotADirectoryError('Sanity output path exists but is not a directory: {}'.format(output_dir))
    os.makedirs(output_dir, exist_ok=True)
    outpath = os.path.join(output_dir, filename)
    print("writing {} to: {}".format(label, outpath))
    with open(outpath, "w") as f:
        for item_id in item_ids:
            f.write(item_id + "\n")
    return outpath


def _normalize_sra_ids(sra_ids):
    normalized = []
    seen = set()
    num_dropped = 0
    for sra_id in list(sra_ids):
        if sra_id is None or pandas.isna(sra_id):
            num_dropped += 1
            continue
        normalized_id = str(sra_id).strip()
        if normalized_id == '':
            num_dropped += 1
            continue
        if normalized_id in seen:
            continue
        seen.add(normalized_id)
        normalized.append(normalized_id)
    return normalized, num_dropped

def _normalize_metadata_run_column(metadata):
    if 'run' not in metadata.df.columns:
        return
    metadata.df.loc[:, 'run'] = (
        metadata.df.loc[:, 'run']
        .fillna('')
        .astype(str)
        .str.strip()
    )

def _prepare_run_output_scan(args, root_path, sra_ids, found_msg, missing_msg):
    if not os.path.exists(root_path):
        print(missing_msg)
        return None, None, None
    if not os.path.isdir(root_path):
        print(missing_msg)
        print('Path exists but is not a directory: {}'.format(root_path))
        return None, None, None
    print(found_msg)
    target_runs = set(sra_ids)
    verbose_run_logs = _should_log_per_item(args, len(target_runs))
    _print_log_mode(args, len(target_runs), singular_label='run')
    run_files_map, non_dir_runs = _scan_target_run_dirs(root_path, target_runs)
    return verbose_run_logs, run_files_map, non_dir_runs

def _resolve_sanity_workers(args, num_tasks, verbose_run_logs):
    if num_tasks <= 1:
        return 1
    if verbose_run_logs:
        return 1
    requested = getattr(args, 'threads', 'auto')
    if is_auto_parallel_option(requested):
        worker_cap = min(8, resolve_detected_cpu_count())
    else:
        worker_cap = validate_positive_int_option(requested, 'threads')
    return max(1, min(worker_cap, num_tasks))

def _run_sanity_tasks(task_items, task_fn, max_workers):
    if (max_workers <= 1) or (len(task_items) <= 1):
        results = dict()
        failures = list()
        for task_item in task_items:
            try:
                results[task_item] = task_fn(task_item)
            except Exception as exc:
                failures.append((task_item, exc))
        return results, failures
    return run_tasks_with_optional_threads(
        task_items=task_items,
        task_fn=task_fn,
        max_workers=max_workers,
    )

def _build_getfastq_sra_stat_cache(sra_ids, metadata):
    cache = dict()
    for sra_id in sra_ids:
        try:
            cache[sra_id] = get_sra_stat(sra_id, metadata)
        except AssertionError:
            cache[sra_id] = None
    return cache


def _print_getfastq_missing_output_message(args, sra_id):
    print("Could not find getfastq output for: ", sra_id, "\n")
    print(
        "Suggested command for rerun: getfastq -e email@adress.com --id ",
        sra_id,
        " -w ",
        args.out_dir,
        "--redo yes --gcp yes --aws yes --ncbi yes",
    )


def _check_single_getfastq_run(
    args,
    sra_id,
    metadata,
    getfastq_path,
    run_files_map,
    non_dir_runs,
    verbose_run_logs,
    sra_stat_cache=None,
):
    sra_path = os.path.join(getfastq_path, sra_id)
    run_files = run_files_map.get(sra_id)
    if (run_files is None) or (sra_id in non_dir_runs):
        if verbose_run_logs:
            _print_getfastq_missing_output_message(args, sra_id)
        return False

    if isinstance(sra_stat_cache, dict):
        sra_stat = sra_stat_cache.get(sra_id, None)
        if sra_stat is None:
            if verbose_run_logs:
                print('Skipping {} due to metadata inconsistency.'.format(sra_id))
            return False
    else:
        try:
            sra_stat = get_sra_stat(sra_id, metadata)
        except AssertionError as exc:
            if verbose_run_logs:
                print('Skipping {} due to metadata inconsistency: {}'.format(sra_id, exc))
            return False
    try:
        ext = get_newest_intermediate_file_extension(sra_stat, sra_path, files=run_files)
    except FileNotFoundError:
        if verbose_run_logs:
            print("could not find any fastq files for ", sra_id, "Please make sure amalgkit getfastq ran properly")
        return False
    if ext == 'no_extension_found':
        if verbose_run_logs:
            print("could not find any fastq files for ", sra_id, "Please make sure amalgkit getfastq ran properly")
        return False

    if verbose_run_logs and (ext != '.safely_removed'):
        files = sorted([
            os.path.join(sra_path, f)
            for f in find_run_prefixed_entries(run_files, sra_id)
            if f.endswith(ext)
        ])
        print("Found:", files)
    return True


def parse_metadata(args, metadata):
    print("Checking essential entries from metadata file.")
    required_columns = ['scientific_name', 'run']
    missing_columns = [col for col in required_columns if col not in metadata.df.columns]
    if len(missing_columns) > 0:
        raise ValueError(
            'Missing required metadata column(s): {}'.format(', '.join(missing_columns))
        )
    species = metadata.df.loc[:, 'scientific_name'].fillna('').astype(str).str.strip()
    species = species.loc[species != '']
    uni_species = np.unique(species.to_numpy(dtype=str))
    if len(uni_species):
        print(len(uni_species), " species detected:")
        print(uni_species)
    else:
        txt = "{} species detected. Please check if --metadata ({}) has a 'scientific_name' column."
        raise ValueError(txt.format(len(uni_species), args.metadata))

    sra_ids = metadata.df.loc[:, 'run'].fillna('').astype(str).str.strip()
    is_missing_run = (sra_ids == '')
    if is_missing_run.any():
        raise ValueError(
            'Found {} metadata entr{} without Run ID in --metadata ({}).'.format(
                int(is_missing_run.sum()),
                'y' if int(is_missing_run.sum()) == 1 else 'ies',
                args.metadata,
            )
        )
    uni_sra_ids = np.unique(sra_ids.to_numpy(dtype=str))
    if len(sra_ids):
        print(len(uni_sra_ids), " SRA runs detected:")
        print(uni_sra_ids)
        # check for duplicate runs
        if len(sra_ids) > len(uni_sra_ids):
            dupes = list_duplicates(sra_ids.tolist())
            raise ValueError(
                "Duplicate SRA IDs detected, where IDs should be unique. Please check these entries: {}".format(dupes)
            )
    else:
        txt = "{} SRA runs detected. Please check if --metadata ({}) has a 'run' column."
        raise ValueError(txt.format(len(uni_sra_ids), args.metadata))
    return uni_species, sra_ids


def check_getfastq_outputs(args, sra_ids, metadata, output_dir):
    print("checking for getfastq outputs: ")
    _normalize_metadata_run_column(metadata)
    sra_ids, num_dropped = _normalize_sra_ids(sra_ids)
    if num_dropped > 0:
        print('Ignored {:,} metadata run entr{} with missing Run ID.'.format(
            num_dropped,
            'y' if num_dropped == 1 else 'ies',
        ))
    if len(sra_ids) == 0:
        raise ValueError('No valid Run IDs were found while checking getfastq outputs.')
    if args.getfastq_dir:
        getfastq_path = args.getfastq_dir
    else:
        getfastq_path = os.path.join(args.out_dir, "getfastq")
    data_available = []
    data_unavailable = []
    verbose_run_logs, run_files_map, non_dir_runs = _prepare_run_output_scan(
        args=args,
        root_path=getfastq_path,
        sra_ids=sra_ids,
        found_msg="amalgkit getfastq output folder detected. Checking presence of output files.",
        missing_msg="Could not find getfastq output folder {}. Have you run getfastq yet?".format(getfastq_path),
    )
    if verbose_run_logs is not None:
        if verbose_run_logs:
            for sra_id in sra_ids:
                print("\n")
                print("Looking for {}".format(sra_id))
                is_available = _check_single_getfastq_run(
                    args=args,
                    sra_id=sra_id,
                    metadata=metadata,
                    getfastq_path=getfastq_path,
                    run_files_map=run_files_map,
                    non_dir_runs=non_dir_runs,
                    verbose_run_logs=verbose_run_logs,
                )
                if is_available:
                    data_available.append(sra_id)
                else:
                    data_unavailable.append(sra_id)
        else:
            max_workers = _resolve_sanity_workers(args=args, num_tasks=len(sra_ids), verbose_run_logs=verbose_run_logs)
            sra_stat_cache = _build_getfastq_sra_stat_cache(sra_ids=sra_ids, metadata=metadata)
            if max_workers > 1:
                print(
                    'Checking getfastq outputs for {:,} runs with {:,} parallel worker(s).'.format(
                        len(sra_ids),
                        max_workers,
                    ),
                    flush=True,
                )
            results_by_sra, failures = _run_sanity_tasks(
                task_items=sra_ids,
                task_fn=lambda sra_id: _check_single_getfastq_run(
                    args=args,
                    sra_id=sra_id,
                    metadata=metadata,
                    getfastq_path=getfastq_path,
                    run_files_map=run_files_map,
                    non_dir_runs=non_dir_runs,
                    verbose_run_logs=False,
                    sra_stat_cache=sra_stat_cache,
                ),
                max_workers=max_workers,
            )
            for sra_id, _exc in failures:
                results_by_sra[sra_id] = False
            for sra_id in sra_ids:
                is_available = bool(results_by_sra.get(sra_id, False))
                if is_available:
                    data_available.append(sra_id)
                else:
                    data_unavailable.append(sra_id)

    else:
        data_unavailable = list(sra_ids)

    if data_unavailable:
        _write_unavailable_items(
            output_dir=output_dir,
            filename="SRA_IDs_without_fastq.txt",
            item_ids=data_unavailable,
            label="SRA IDs without getfastq output",
        )
    else:
        txt = "The getfastq output files for all SRA IDs in --metadata ({}) were found."
        print(txt.format(args.metadata))

    return data_available, data_unavailable

def _find_index_files(index_entries, index_dir_path, prefix):
    matched = find_species_prefixed_entries(index_entries, prefix, entries_sorted=True)
    return [
        os.path.join(index_dir_path, entry)
        for entry in matched
        if os.path.isfile(os.path.join(index_dir_path, entry))
    ]

def _resolve_species_index_files(species, index_entries, index_dir_path, verbose_run_logs=True):
    candidates = _get_species_prefix_candidates(species)
    sci_name = candidates[0]
    index_path = os.path.join(index_dir_path, sci_name + "*")
    if verbose_run_logs:
        print("\n")
        print("Looking for index file {} for species {}".format(index_path, species))
    for candidate in candidates:
        index_files = _find_index_files(index_entries, index_dir_path, candidate)
        if index_files:
            if verbose_run_logs:
                print("Found ", index_files, "!")
            return index_files

    if verbose_run_logs:
        print("could not find anything in", index_path)
    # Deprecate subspecies or variants and look again
    # I.e. if Gorilla_gorilla_gorilla.idx was not found, we look for Gorilla_gorilla.idx instead.
    fallback_prefix = _get_species_fallback_prefix(species)
    if fallback_prefix is None:
        if verbose_run_logs:
            print("Could not find any index files for ", species)
        return []
    if verbose_run_logs:
        print("Ignoring subspecies.")
    fallback_candidates = _get_species_prefix_candidates(fallback_prefix)
    index_path = os.path.join(index_dir_path, fallback_candidates[0] + "*")
    if verbose_run_logs:
        print("Looking for {}".format(index_path))
    for candidate in fallback_candidates:
        index_files = _find_index_files(index_entries, index_dir_path, candidate)
        if index_files:
            if verbose_run_logs:
                print("Found ", index_files, "!")
            return index_files
    if verbose_run_logs:
        print("Could not find any index files for ", species)
    return []


def _classify_index_files(species, index_files, index_available, index_unavailable, verbose_run_logs):
    if len(index_files) == 1:
        index_available.append(species)
        return False
    if len(index_files) > 1:
        if verbose_run_logs:
            print(
                "Multiple possible index files detected for ",
                species,
                ": ",
                index_files,
                ". Please keep only one index file per species.",
            )
        index_unavailable.append(species)
        return True
    index_unavailable.append(species)
    return False


def _check_single_quant_run(sra_id, quant_path, quant_run_files_map, non_dir_runs, verbose_run_logs):
    sra_path = os.path.join(quant_path, sra_id)
    quant_run_files = quant_run_files_map.get(sra_id)
    if (quant_run_files is None) or (sra_id in non_dir_runs):
        if verbose_run_logs:
            print("Could not find output folder ", sra_path, " for ", sra_id)
        return False

    if verbose_run_logs:
        print("Found output folder ", sra_path, " for ", sra_id)
        print("Checking for output files.")
    abundance_file = sra_id + "_abundance.tsv"
    run_info_file = sra_id + "_run_info.json"
    has_abundance = abundance_file in quant_run_files
    has_run_info = run_info_file in quant_run_files
    if has_abundance and has_run_info:
        if verbose_run_logs:
            print("All quant output files present for", sra_id, "!")
        return True

    if verbose_run_logs and (not has_abundance):
        print(os.path.join(sra_path, abundance_file), " is missing! Please check if quant ran correctly")
    if verbose_run_logs and (not has_run_info):
        print(os.path.join(sra_path, run_info_file), " is missing! Please check if quant ran correctly")
    return False


def check_quant_index(args, uni_species, output_dir):
    if args.index_dir:
        index_dir_path = args.index_dir
    else:
        index_dir_path = os.path.join(args.out_dir, "index")
    index_unavailable = []
    index_available = []
    if os.path.exists(index_dir_path):
        if not os.path.isdir(index_dir_path):
            print("Could not find index directory ", index_dir_path, " . Did you provide the correct Path?")
            print('Path exists but is not a directory: {}'.format(index_dir_path))
            index_unavailable = list(uni_species)
            _write_unavailable_items(
                output_dir=output_dir,
                filename="species_without_index.txt",
                item_ids=index_unavailable,
                label="species without index",
            )
            return index_available, index_unavailable
        index_entries = sorted(os.listdir(index_dir_path))
        species_items = list(uni_species)
        verbose_run_logs = _should_log_per_item(args, len(species_items))
        _print_log_mode(args, len(species_items), singular_label='species', plural_label='species')
        ambiguous_species = []
        if verbose_run_logs:
            for species in species_items:
                index_files = _resolve_species_index_files(
                    species=species,
                    index_entries=index_entries,
                    index_dir_path=index_dir_path,
                    verbose_run_logs=verbose_run_logs,
                )
                is_ambiguous = _classify_index_files(
                    species=species,
                    index_files=index_files,
                    index_available=index_available,
                    index_unavailable=index_unavailable,
                    verbose_run_logs=verbose_run_logs,
                )
                if is_ambiguous:
                    ambiguous_species.append(species)
        else:
            max_workers = _resolve_sanity_workers(args=args, num_tasks=len(species_items), verbose_run_logs=verbose_run_logs)
            if max_workers > 1:
                print(
                    'Checking index inputs for {:,} species with {:,} parallel worker(s).'.format(
                        len(species_items),
                        max_workers,
                    ),
                    flush=True,
                )
            results_by_species, failures = _run_sanity_tasks(
                task_items=species_items,
                task_fn=lambda species: _resolve_species_index_files(
                    species=species,
                    index_entries=index_entries,
                    index_dir_path=index_dir_path,
                    verbose_run_logs=False,
                ),
                max_workers=max_workers,
            )
            for species, _exc in failures:
                results_by_species[species] = []
            for species in species_items:
                is_ambiguous = _classify_index_files(
                    species=species,
                    index_files=results_by_species.get(species, []),
                    index_available=index_available,
                    index_unavailable=index_unavailable,
                    verbose_run_logs=False,
                )
                if is_ambiguous:
                    ambiguous_species.append(species)
            if ambiguous_species:
                print(
                    'Multiple possible index files detected for {:,} species. '
                    'Please keep only one index file per species.'.format(len(ambiguous_species))
                )

        if index_unavailable:
            _write_unavailable_items(
                output_dir=output_dir,
                filename="species_without_index.txt",
                item_ids=index_unavailable,
                label="species without index",
            )
        else:
            print("index found for all species in --metadata ({})".format(args.metadata))
    else:
        print("Could not find index directory ", index_dir_path, " . Did you provide the correct Path?")
        index_unavailable = list(uni_species)
        _write_unavailable_items(
            output_dir=output_dir,
            filename="species_without_index.txt",
            item_ids=index_unavailable,
            label="species without index",
        )

    return index_available, index_unavailable


def check_quant_output(args, sra_ids, output_dir):
    print("checking for quant outputs: ")
    sra_ids, num_dropped = _normalize_sra_ids(sra_ids)
    if num_dropped > 0:
        print('Ignored {:,} metadata run entr{} with missing Run ID.'.format(
            num_dropped,
            'y' if num_dropped == 1 else 'ies',
        ))
    if len(sra_ids) == 0:
        raise ValueError('No valid Run IDs were found while checking quant outputs.')
    quant_path = getattr(args, 'quant_dir', None)
    if quant_path is None:
        quant_path = os.path.join(args.out_dir, "quant")
    data_available = []
    data_unavailable = []
    verbose_run_logs, quant_run_files_map, non_dir_runs = _prepare_run_output_scan(
        args=args,
        root_path=quant_path,
        sra_ids=sra_ids,
        found_msg="amalgkit quant output folder detected. Checking presence of output files.",
        missing_msg="Could not find quant output folder {}. Have you run quant yet?".format(quant_path),
    )
    if verbose_run_logs is not None:
        if verbose_run_logs:
            for sra_id in sra_ids:
                print("\n")
                print("Looking for {}".format(sra_id))
                is_available = _check_single_quant_run(
                    sra_id=sra_id,
                    quant_path=quant_path,
                    quant_run_files_map=quant_run_files_map,
                    non_dir_runs=non_dir_runs,
                    verbose_run_logs=verbose_run_logs,
                )
                if is_available:
                    data_available.append(sra_id)
                else:
                    data_unavailable.append(sra_id)
        else:
            max_workers = _resolve_sanity_workers(args=args, num_tasks=len(sra_ids), verbose_run_logs=verbose_run_logs)
            if max_workers > 1:
                print(
                    'Checking quant outputs for {:,} runs with {:,} parallel worker(s).'.format(
                        len(sra_ids),
                        max_workers,
                    ),
                    flush=True,
                )
            results_by_sra, failures = _run_sanity_tasks(
                task_items=sra_ids,
                task_fn=lambda sra_id: _check_single_quant_run(
                    sra_id=sra_id,
                    quant_path=quant_path,
                    quant_run_files_map=quant_run_files_map,
                    non_dir_runs=non_dir_runs,
                    verbose_run_logs=False,
                ),
                max_workers=max_workers,
            )
            for sra_id, _exc in failures:
                results_by_sra[sra_id] = False
            for sra_id in sra_ids:
                is_available = bool(results_by_sra.get(sra_id, False))
                if is_available:
                    data_available.append(sra_id)
                else:
                    data_unavailable.append(sra_id)
    else:
        data_unavailable = list(sra_ids)

    if data_unavailable:
        _write_unavailable_items(
            output_dir=output_dir,
            filename="SRA_IDs_without_quant.txt",
            item_ids=data_unavailable,
            label="SRA IDs without quant output",
        )
    else:
        print("Quant outputs found for all SRA IDs in --metadata ({})".format(args.metadata))

    return data_available, data_unavailable


SANITY_CHECK_NAMES = ['getfastq', 'index', 'quant', 'merge', 'busco', 'finalize']
BUSCO_REQUIRED_COLUMNS = [
    'busco_id',
    'status',
    'sequence',
    'score',
    'length',
    'orthodb_url',
    'description',
]


def _resolve_sanity_check_path(args, subdir_name, attr_name):
    custom_path = getattr(args, attr_name, None)
    if custom_path:
        return os.path.realpath(custom_path)
    return os.path.join(os.path.realpath(args.out_dir), subdir_name)


def _resolve_metadata_path(args):
    if getattr(args, 'metadata', 'inferred') == 'inferred':
        return os.path.realpath(os.path.join(args.out_dir, 'metadata', 'metadata.tsv'))
    return os.path.realpath(args.metadata)


def _parse_csv_option(raw_value):
    if raw_value in [None, '']:
        return []
    values = []
    seen = set()
    for item in re.split(r'[\s,]+', str(raw_value).strip()):
        normalized = item.strip()
        if normalized == '':
            continue
        if normalized in seen:
            continue
        seen.add(normalized)
        values.append(normalized)
    return values


def _filter_metadata_for_sanity(args, metadata):
    run_filters = _parse_csv_option(getattr(args, 'run', None))
    species_filters = _parse_csv_option(getattr(args, 'species', None))
    filtered_df = metadata.df.copy()
    if species_filters:
        if 'scientific_name' not in filtered_df.columns:
            raise ValueError('Column "scientific_name" is required when --species is specified.')
        species_series = filtered_df.loc[:, 'scientific_name'].fillna('').astype(str).str.strip()
        filtered_df = filtered_df.loc[species_series.isin(species_filters), :].copy()
    if run_filters:
        if 'run' not in filtered_df.columns:
            raise ValueError('Column "run" is required when --run is specified.')
        run_series = filtered_df.loc[:, 'run'].fillna('').astype(str).str.strip()
        filtered_df = filtered_df.loc[run_series.isin(run_filters), :].copy()
    if filtered_df.shape[0] == 0:
        raise ValueError('No metadata entries remained after applying --run/--species filters.')
    metadata.df = filtered_df.reset_index(drop=True)
    return run_filters, species_filters


def _build_issue(check_name, severity, issue_type, target_type, target_id='', path='', message='', suggested_action=''):
    return {
        'check': str(check_name),
        'severity': str(severity),
        'issue_type': str(issue_type),
        'target_type': str(target_type),
        'target_id': '' if target_id in [None, ''] else str(target_id),
        'path': '' if path in [None, ''] else os.path.realpath(path),
        'message': '' if message in [None, ''] else str(message),
        'suggested_action': '' if suggested_action in [None, ''] else str(suggested_action),
    }


def _write_manifest(output_dir, filename, target_ids, label):
    normalized_targets = []
    seen = set()
    for target_id in target_ids:
        normalized = str(target_id).strip()
        if normalized == '':
            continue
        if normalized in seen:
            continue
        seen.add(normalized)
        normalized_targets.append(normalized)
    return _write_unavailable_items(output_dir, filename, normalized_targets, label)


def _collect_error_targets(issues, target_type):
    targets = []
    seen = set()
    for issue in issues:
        if str(issue.get('severity', '')).lower() != 'error':
            continue
        if str(issue.get('target_type', '')).lower() != str(target_type).lower():
            continue
        target_id = str(issue.get('target_id', '')).strip()
        if target_id == '':
            continue
        if target_id in seen:
            continue
        seen.add(target_id)
        targets.append(target_id)
    return targets


def _safe_get_mtime(path):
    try:
        return os.path.getmtime(path)
    except (FileNotFoundError, OSError):
        return None


def _is_stale_output(metadata_path, output_path):
    metadata_mtime = _safe_get_mtime(metadata_path)
    output_mtime = _safe_get_mtime(output_path)
    if (metadata_mtime is None) or (output_mtime is None):
        return False
    return output_mtime < metadata_mtime


def _list_root_entries(root_path):
    if (not os.path.exists(root_path)) or (not os.path.isdir(root_path)):
        return []
    with os.scandir(root_path) as entries:
        return sorted([entry.name for entry in entries if not entry.name.startswith('.')])


def _read_tsv_head(path, nrows=5, comment=None):
    return pandas.read_csv(
        path,
        sep='\t',
        header=0,
        comment=comment,
        nrows=nrows,
        low_memory=False,
    )


def _validate_nonempty_table(path, required_columns, context, comment=None, require_data_rows=True, require_non_target_columns=False):
    try:
        df = _read_tsv_head(path, comment=comment)
    except Exception as exc:
        return 'Failed to read {}: {}'.format(context, exc)
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        return 'Missing required column(s) in {}: {}'.format(context, ', '.join(missing_columns))
    if require_non_target_columns and (len(df.columns) <= len(required_columns)):
        return '{} did not include any data columns beyond {}.'.format(context, ', '.join(required_columns))
    if require_data_rows and (df.shape[0] == 0):
        return '{} did not contain any data rows.'.format(context)
    if 'target_id' in df.columns:
        target_ids = df.loc[:, 'target_id'].fillna('').astype(str).str.strip()
        if require_data_rows and target_ids.eq('').all():
            return '{} did not contain valid target_id values.'.format(context)
    return ''


def _validate_quant_run_info_json(path):
    try:
        with open(path, 'r', encoding='utf-8') as handle:
            payload = json.load(handle)
    except Exception as exc:
        return 'Failed to read quant run info JSON: {}'.format(exc)
    if 'p_pseudoaligned' not in payload:
        return 'quant run info JSON is missing "p_pseudoaligned".'
    try:
        value = float(payload.get('p_pseudoaligned'))
    except (TypeError, ValueError):
        return 'quant run info JSON has an invalid "p_pseudoaligned" value.'
    if (value < 0.0) or (value > 100.0):
        return 'quant run info JSON has out-of-range "p_pseudoaligned": {}'.format(value)
    return ''


def _validate_fastq_file(path):
    if not os.path.isfile(path):
        return 'FASTQ file not found: {}'.format(path)
    if os.path.getsize(path) <= 0:
        return 'FASTQ file is empty: {}'.format(path)
    newline_count = 0
    total_bytes = 0
    try:
        with open_fastq_binary(path) as handle:
            while True:
                chunk = handle.read(16 * 1024 * 1024)
                if chunk == b'':
                    break
                total_bytes += len(chunk)
                newline_count += chunk.count(b'\n')
    except Exception as exc:
        return 'FASTQ stream validation failed for {}: {}'.format(path, exc)
    if total_bytes == 0:
        return 'FASTQ file is empty after decompression: {}'.format(path)
    if newline_count < 4:
        return 'FASTQ file does not contain a complete record: {}'.format(path)
    if (newline_count % 4) != 0:
        return 'FASTQ line count is not divisible by 4: {}'.format(path)
    return ''


def _validate_busco_table(path):
    try:
        df = pandas.read_table(
            path,
            sep='\t',
            header=None,
            comment='#',
            names=BUSCO_REQUIRED_COLUMNS,
            dtype=str,
            low_memory=False,
            nrows=10,
        )
    except Exception as exc:
        return 'Failed to read BUSCO table: {}'.format(exc)
    if df.shape[0] == 0:
        return 'BUSCO table did not contain any data rows.'
    busco_ids = df.loc[:, 'busco_id'].fillna('').astype(str).str.lower().str.replace(r'[^a-z0-9]', '', regex=True)
    if busco_ids.eq('buscoid').all():
        return 'BUSCO table did not contain any data rows.'
    return ''


def _normalize_issue_key(issue):
    target_id = str(issue.get('target_id', '')).strip()
    if target_id != '':
        return target_id
    issue_type = str(issue.get('issue_type', 'issue')).strip()
    path = str(issue.get('path', '')).strip()
    if path != '':
        return '{}:{}'.format(issue_type, os.path.basename(path))
    return issue_type


def _build_sanity_summary_row(check_name, target_type, checked_items, issues, scanned_path, report_path='', rerun_manifest_path=''):
    checked_ids = [str(item).strip() for item in list(checked_items) if str(item).strip() != '']
    checked_id_set = set(checked_ids)
    error_checked_targets = {
        str(issue.get('target_id', '')).strip()
        for issue in issues
        if (str(issue.get('severity', '')).lower() == 'error')
        and (str(issue.get('target_id', '')).strip() in checked_id_set)
    }
    warning_targets = {
        _normalize_issue_key(issue)
        for issue in issues
        if str(issue.get('severity', '')).lower() == 'warning'
    }
    error_targets = {
        _normalize_issue_key(issue)
        for issue in issues
        if str(issue.get('severity', '')).lower() == 'error'
    }
    checked_count = len(checked_ids)
    unavailable_count = len(error_checked_targets)
    available_count = max(0, checked_count - unavailable_count)
    status = 'ok'
    if error_targets:
        status = 'error'
    elif warning_targets:
        status = 'warning'
    example_targets = sorted(list(error_targets | warning_targets))[0:5]
    return {
        'check': check_name,
        'target_type': target_type,
        'checked_count': int(checked_count),
        'available_count': int(available_count),
        'unavailable_count': int(unavailable_count),
        'warning_count': int(len(warning_targets)),
        'issue_count': int(len(issues)),
        'status': status,
        'scanned_path': os.path.realpath(scanned_path),
        'missing_report_path': '' if report_path == '' else os.path.realpath(report_path),
        'rerun_manifest_path': '' if rerun_manifest_path == '' else os.path.realpath(rerun_manifest_path),
        'example_targets': ', '.join(example_targets),
    }


def _write_sanity_summary(output_dir, summary_rows):
    summary_path = os.path.join(output_dir, 'sanity_summary.tsv')
    summary_df = pandas.DataFrame(summary_rows, columns=[
        'check',
        'target_type',
        'checked_count',
        'available_count',
        'unavailable_count',
        'warning_count',
        'issue_count',
        'status',
        'scanned_path',
        'missing_report_path',
        'rerun_manifest_path',
        'example_targets',
    ])
    summary_df.to_csv(summary_path, sep='\t', index=False)
    return summary_path


def _write_sanity_issues(output_dir, issues):
    issues_path = os.path.join(output_dir, 'sanity_issues.tsv')
    issues_df = pandas.DataFrame(issues, columns=[
        'check',
        'severity',
        'issue_type',
        'target_type',
        'target_id',
        'path',
        'message',
        'suggested_action',
    ])
    issues_df.to_csv(issues_path, sep='\t', index=False)
    return issues_path


def _load_previous_sanity_report(output_dir):
    report_path = os.path.join(output_dir, 'sanity_report.json')
    if not os.path.isfile(report_path):
        return None
    try:
        with open(report_path, 'r', encoding='utf-8') as handle:
            return json.load(handle)
    except Exception:
        return None


def _write_sanity_report_json(output_dir, payload):
    report_path = os.path.join(output_dir, 'sanity_report.json')
    with open(report_path, 'w', encoding='utf-8') as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
    return report_path


def _write_sanity_comparison(output_dir, previous_payload, summary_rows):
    if not isinstance(previous_payload, dict):
        return ''
    previous_summary = previous_payload.get('summary', [])
    previous_by_check = {
        str(row.get('check', '')).strip(): row
        for row in previous_summary
        if str(row.get('check', '')).strip() != ''
    }
    current_by_check = {
        str(row.get('check', '')).strip(): row
        for row in summary_rows
        if str(row.get('check', '')).strip() != ''
    }
    comparison_rows = []
    for check_name in SANITY_CHECK_NAMES:
        if (check_name not in previous_by_check) and (check_name not in current_by_check):
            continue
        prev_row = previous_by_check.get(check_name, {})
        curr_row = current_by_check.get(check_name, {})
        prev_unavailable = int(prev_row.get('unavailable_count', 0) or 0)
        curr_unavailable = int(curr_row.get('unavailable_count', 0) or 0)
        prev_warnings = int(prev_row.get('warning_count', 0) or 0)
        curr_warnings = int(curr_row.get('warning_count', 0) or 0)
        comparison_rows.append({
            'check': check_name,
            'previous_status': str(prev_row.get('status', 'absent')),
            'current_status': str(curr_row.get('status', 'absent')),
            'previous_unavailable_count': prev_unavailable,
            'current_unavailable_count': curr_unavailable,
            'delta_unavailable_count': int(curr_unavailable - prev_unavailable),
            'previous_warning_count': prev_warnings,
            'current_warning_count': curr_warnings,
            'delta_warning_count': int(curr_warnings - prev_warnings),
        })
    if not comparison_rows:
        return ''
    comparison_path = os.path.join(output_dir, 'sanity_comparison.tsv')
    pandas.DataFrame(comparison_rows).to_csv(comparison_path, sep='\t', index=False)
    return comparison_path


def _print_sanity_summary(summary_rows, summary_path):
    print("\nSanity summary:")
    if not summary_rows:
        print('No sanity checks were selected.')
    for row in summary_rows:
        print(
            '{}: {:,}/{:,} {} available (status={}; errors={:,}; warnings={:,}; issues={:,})'.format(
                row['check'],
                int(row['available_count']),
                int(row['checked_count']),
                row['target_type'],
                row['status'],
                int(row['unavailable_count']),
                int(row.get('warning_count', 0)),
                int(row.get('issue_count', 0)),
            )
        )
        print('scanned_path={}'.format(row['scanned_path']))
        if row['missing_report_path'] != '':
            print('missing_report_path={}'.format(row['missing_report_path']))
        if row.get('rerun_manifest_path', '') != '':
            print('rerun_manifest_path={}'.format(row['rerun_manifest_path']))
        if row.get('example_targets', '') != '':
            print('example_targets={}'.format(row['example_targets']))
    print('summary_path={}'.format(summary_path))


def _resolve_requested_sanity_checks(args):
    selected = []
    requested_checks = [value.lower() for value in _parse_csv_option(getattr(args, 'check', None))]
    if 'all' in requested_checks:
        return list(SANITY_CHECK_NAMES)
    for check_name in requested_checks:
        if check_name not in SANITY_CHECK_NAMES:
            raise ValueError(
                'Unknown sanity check "{}". Accepted values: {}.'.format(
                    check_name,
                    ', '.join(SANITY_CHECK_NAMES + ['all']),
                )
            )
        selected.append(check_name)
    if bool(getattr(args, 'all', False)):
        selected.extend(SANITY_CHECK_NAMES)
    for check_name in SANITY_CHECK_NAMES:
        if bool(getattr(args, check_name, False)):
            selected.append(check_name)
    if selected:
        return list(dict.fromkeys(selected))
    print('No specific sanity targets were selected. Running all checks (--all).')
    return list(SANITY_CHECK_NAMES)


def _resolve_strict_level(args):
    strict_level = str(getattr(args, 'strict_level', 'none')).strip().lower()
    if bool(getattr(args, 'strict', False)) and (strict_level == 'none'):
        return 'error'
    return strict_level


def _raise_if_sanity_issues(args, issues):
    strict_level = _resolve_strict_level(args)
    if strict_level == 'none':
        return
    if strict_level == 'error':
        should_fail = any(str(issue.get('severity', '')).lower() == 'error' for issue in issues)
    else:
        should_fail = len(issues) > 0
    if not should_fail:
        return
    severity_counts = {'error': 0, 'warning': 0}
    for issue in issues:
        severity = str(issue.get('severity', '')).lower()
        if severity in severity_counts:
            severity_counts[severity] += 1
    raise AmalgkitExit(
        'Sanity checks reported {:,} error(s) and {:,} warning(s) at strict_level={}.'.format(
            severity_counts['error'],
            severity_counts['warning'],
            strict_level,
        ),
        exit_code=1,
    )


def _validate_orphan_entries(check_name, root_path, expected_names, target_type):
    issues = []
    if (not os.path.exists(root_path)) or (not os.path.isdir(root_path)):
        return issues
    expected_set = set(expected_names)
    for entry_name in _list_root_entries(root_path):
        if entry_name in expected_set:
            continue
        issues.append(_build_issue(
            check_name=check_name,
            severity='warning',
            issue_type='orphan_output',
            target_type=target_type,
            target_id=entry_name,
            path=os.path.join(root_path, entry_name),
            message='Output exists under {} but is not expected from the filtered metadata.'.format(root_path),
            suggested_action='review_or_cleanup',
        ))
    return issues


def _validate_getfastq_content(metadata, sra_id, root_path, metadata_path, run_files_map):
    issues = []
    run_dir = os.path.join(root_path, sra_id)
    run_files = run_files_map.get(sra_id, set())
    try:
        sra_stat = get_sra_stat(sra_id, metadata)
    except AssertionError as exc:
        issues.append(_build_issue(
            check_name='getfastq',
            severity='error',
            issue_type='metadata_inconsistency',
            target_type='run',
            target_id=sra_id,
            path=run_dir,
            message='Metadata inconsistency while resolving getfastq files: {}'.format(exc),
            suggested_action='review_metadata',
        ))
        return issues
    try:
        ext = get_newest_intermediate_file_extension(sra_stat, run_dir, files=run_files)
    except Exception as exc:
        issues.append(_build_issue(
            check_name='getfastq',
            severity='error',
            issue_type='invalid_content',
            target_type='run',
            target_id=sra_id,
            path=run_dir,
            message='Failed to resolve getfastq intermediate extension: {}'.format(exc),
            suggested_action='rerun_getfastq',
        ))
        return issues
    if ext == '.safely_removed':
        safely_removed_files = [
            os.path.join(run_dir, name)
            for name in run_files
            if name.endswith('.safely_removed')
        ]
        if len(safely_removed_files) == 0:
            issues.append(_build_issue(
                check_name='getfastq',
                severity='error',
                issue_type='invalid_content',
                target_type='run',
                target_id=sra_id,
                path=run_dir,
                message='Run was marked as safely removed but no .safely_removed sentinel was found.',
                suggested_action='rerun_getfastq',
            ))
        return issues
    if str(sra_stat.get('layout', '')).lower() == 'paired':
        expected_paths = [
            os.path.join(run_dir, sra_id + '_1' + ext),
            os.path.join(run_dir, sra_id + '_2' + ext),
        ]
    else:
        expected_paths = [os.path.join(run_dir, sra_id + ext)]
    for expected_path in expected_paths:
        if not os.path.isfile(expected_path):
            issues.append(_build_issue(
                check_name='getfastq',
                severity='error',
                issue_type='missing_output',
                target_type='run',
                target_id=sra_id,
                path=expected_path,
                message='Expected FASTQ component was not found.',
                suggested_action='rerun_getfastq',
            ))
            continue
        validation_error = _validate_fastq_file(expected_path)
        if validation_error != '':
            issues.append(_build_issue(
                check_name='getfastq',
                severity='error',
                issue_type='invalid_content',
                target_type='run',
                target_id=sra_id,
                path=expected_path,
                message=validation_error,
                suggested_action='rerun_getfastq',
            ))
        if _is_stale_output(metadata_path, expected_path):
            issues.append(_build_issue(
                check_name='getfastq',
                severity='warning',
                issue_type='stale_output',
                target_type='run',
                target_id=sra_id,
                path=expected_path,
                message='FASTQ output is older than the metadata file.',
                suggested_action='review_or_rerun',
            ))
    return issues


def run_sanity_check_getfastq(args, metadata, uni_species, sra_ids, output_dir, metadata_path):
    _ = uni_species
    root_path = _resolve_sanity_check_path(args, 'getfastq', 'getfastq_dir')
    checked_runs = _normalize_sra_ids(sra_ids)[0]
    available, unavailable = check_getfastq_outputs(args, sra_ids, metadata, output_dir)
    issues = []
    for sra_id in unavailable:
        issues.append(_build_issue(
            check_name='getfastq',
            severity='error',
            issue_type='missing_output',
            target_type='run',
            target_id=sra_id,
            path=os.path.join(root_path, sra_id),
            message='getfastq output was not detected for this run.',
            suggested_action='rerun_getfastq',
        ))
    issues.extend(_validate_orphan_entries('getfastq', root_path, checked_runs, 'run'))
    if os.path.isdir(root_path):
        run_files_map, _non_dir_runs = _scan_target_run_dirs(root_path, set(available))
        for sra_id in available:
            issues.extend(_validate_getfastq_content(metadata, sra_id, root_path, metadata_path, run_files_map))
    rerun_manifest_path = _write_manifest(
        output_dir,
        'rerun_getfastq_run_ids.txt',
        _collect_error_targets(issues, 'run'),
        'getfastq rerun run IDs',
    )
    report_path = os.path.join(output_dir, 'SRA_IDs_without_fastq.txt') if unavailable else ''
    row = _build_sanity_summary_row('getfastq', 'runs', checked_runs, issues, root_path, report_path, rerun_manifest_path)
    return row, issues


def _validate_quant_content(sra_id, root_path, metadata_path):
    issues = []
    abundance_path = os.path.join(root_path, sra_id, sra_id + '_abundance.tsv')
    run_info_path = os.path.join(root_path, sra_id, sra_id + '_run_info.json')
    abundance_error = _validate_nonempty_table(
        abundance_path,
        required_columns=['target_id', 'length', 'eff_length', 'est_counts', 'tpm'],
        context='quant abundance table {}'.format(abundance_path),
        require_data_rows=True,
    )
    if abundance_error != '':
        issues.append(_build_issue(
            check_name='quant',
            severity='error',
            issue_type='invalid_content',
            target_type='run',
            target_id=sra_id,
            path=abundance_path,
            message=abundance_error,
            suggested_action='rerun_quant',
        ))
    run_info_error = _validate_quant_run_info_json(run_info_path)
    if run_info_error != '':
        issues.append(_build_issue(
            check_name='quant',
            severity='error',
            issue_type='invalid_content',
            target_type='run',
            target_id=sra_id,
            path=run_info_path,
            message=run_info_error,
            suggested_action='rerun_quant',
        ))
    for path in [abundance_path, run_info_path]:
        if os.path.exists(path) and _is_stale_output(metadata_path, path):
            issues.append(_build_issue(
                check_name='quant',
                severity='warning',
                issue_type='stale_output',
                target_type='run',
                target_id=sra_id,
                path=path,
                message='quant output is older than the metadata file.',
                suggested_action='review_or_rerun',
            ))
    return issues


def run_sanity_check_quant(args, metadata, uni_species, sra_ids, output_dir, metadata_path):
    _ = (metadata, uni_species)
    root_path = _resolve_sanity_check_path(args, 'quant', 'quant_dir')
    checked_runs = _normalize_sra_ids(sra_ids)[0]
    available, unavailable = check_quant_output(args, sra_ids, output_dir)
    issues = []
    for sra_id in unavailable:
        issues.append(_build_issue(
            check_name='quant',
            severity='error',
            issue_type='missing_output',
            target_type='run',
            target_id=sra_id,
            path=os.path.join(root_path, sra_id),
            message='quant output was not detected for this run.',
            suggested_action='rerun_quant',
        ))
    issues.extend(_validate_orphan_entries('quant', root_path, checked_runs, 'run'))
    if os.path.isdir(root_path):
        for sra_id in available:
            issues.extend(_validate_quant_content(sra_id, root_path, metadata_path))
    rerun_manifest_path = _write_manifest(
        output_dir,
        'rerun_quant_run_ids.txt',
        _collect_error_targets(issues, 'run'),
        'quant rerun run IDs',
    )
    report_path = os.path.join(output_dir, 'SRA_IDs_without_quant.txt') if unavailable else ''
    row = _build_sanity_summary_row('quant', 'runs', checked_runs, issues, root_path, report_path, rerun_manifest_path)
    return row, issues


def _matches_any_species_index(entry_name, species_values):
    for species in species_values:
        candidates = _get_species_prefix_candidates(species)
        fallback_prefix = _get_species_fallback_prefix(species)
        if fallback_prefix is not None:
            candidates.extend(_get_species_prefix_candidates(fallback_prefix))
        candidates = list(dict.fromkeys(candidates))
        for candidate in candidates:
            if find_species_prefixed_entries([entry_name], candidate):
                return True
    return False


def run_sanity_check_index(args, metadata, uni_species, sra_ids, output_dir, metadata_path):
    _ = (metadata, sra_ids)
    root_path = _resolve_sanity_check_path(args, 'index', 'index_dir')
    checked_species = [str(species).strip() for species in list(uni_species)]
    _available, unavailable = check_quant_index(args, uni_species, output_dir)
    issues = []
    if os.path.isdir(root_path):
        index_entries = sorted(os.listdir(root_path))
        for species in checked_species:
            index_files = _resolve_species_index_files(
                species=species,
                index_entries=index_entries,
                index_dir_path=root_path,
                verbose_run_logs=False,
            )
            if len(index_files) == 0:
                issues.append(_build_issue(
                    check_name='index',
                    severity='error',
                    issue_type='missing_output',
                    target_type='species',
                    target_id=species,
                    path=os.path.join(root_path, _normalize_species_prefix(species)),
                    message='No index file matched this species.',
                    suggested_action='rerun_index',
                ))
                continue
            if len(index_files) > 1:
                issues.append(_build_issue(
                    check_name='index',
                    severity='error',
                    issue_type='ambiguous_output',
                    target_type='species',
                    target_id=species,
                    path=index_files[0],
                    message='Multiple index files matched this species: {}'.format(', '.join(index_files)),
                    suggested_action='review_or_cleanup',
                ))
                continue
            index_path = index_files[0]
            if os.path.getsize(index_path) <= 0:
                issues.append(_build_issue(
                    check_name='index',
                    severity='error',
                    issue_type='invalid_content',
                    target_type='species',
                    target_id=species,
                    path=index_path,
                    message='Index file is empty.',
                    suggested_action='rerun_index',
                ))
            if _is_stale_output(metadata_path, index_path):
                issues.append(_build_issue(
                    check_name='index',
                    severity='warning',
                    issue_type='stale_output',
                    target_type='species',
                    target_id=species,
                    path=index_path,
                    message='Index file is older than the metadata file.',
                    suggested_action='review_or_rerun',
                ))
        for entry_name in _list_root_entries(root_path):
            entry_path = os.path.join(root_path, entry_name)
            if os.path.isdir(entry_path):
                continue
            if _matches_any_species_index(entry_name, checked_species):
                continue
            issues.append(_build_issue(
                check_name='index',
                severity='warning',
                issue_type='orphan_output',
                target_type='file',
                target_id=entry_name,
                path=entry_path,
                message='Index file does not match any species from the filtered metadata.',
                suggested_action='review_or_cleanup',
            ))
    else:
        for species in checked_species:
            issues.append(_build_issue(
                check_name='index',
                severity='error',
                issue_type='missing_output',
                target_type='species',
                target_id=species,
                path=root_path,
                message='Index directory was not found.',
                suggested_action='rerun_index',
            ))
    rerun_manifest_path = _write_manifest(
        output_dir,
        'rerun_index_species.txt',
        _collect_error_targets(issues, 'species'),
        'index rerun species',
    )
    report_path = os.path.join(output_dir, 'species_without_index.txt') if unavailable else ''
    row = _build_sanity_summary_row('index', 'species', checked_species, issues, root_path, report_path, rerun_manifest_path)
    return row, issues


def _validate_merge_species_tables(species, species_dir, metadata_path):
    issues = []
    species_tag = _normalize_species_prefix(species)
    for suffix in ['eff_length.tsv', 'est_counts.tsv', 'tpm.tsv']:
        table_path = os.path.join(species_dir, '{}_{}'.format(species_tag, suffix))
        if not os.path.isfile(table_path):
            issues.append(_build_issue(
                check_name='merge',
                severity='error',
                issue_type='missing_output',
                target_type='species',
                target_id=species,
                path=table_path,
                message='Required merged table is missing.',
                suggested_action='rerun_merge',
            ))
            continue
        table_error = _validate_nonempty_table(
            table_path,
            required_columns=['target_id'],
            context='merge table {}'.format(table_path),
            require_data_rows=True,
            require_non_target_columns=True,
        )
        if table_error != '':
            issues.append(_build_issue(
                check_name='merge',
                severity='error',
                issue_type='invalid_content',
                target_type='species',
                target_id=species,
                path=table_path,
                message=table_error,
                suggested_action='rerun_merge',
            ))
        if _is_stale_output(metadata_path, table_path):
            issues.append(_build_issue(
                check_name='merge',
                severity='warning',
                issue_type='stale_output',
                target_type='species',
                target_id=species,
                path=table_path,
                message='merge table is older than the metadata file.',
                suggested_action='review_or_rerun',
            ))
    return issues


def run_sanity_check_merge(args, metadata, uni_species, sra_ids, output_dir, metadata_path):
    _ = (metadata, sra_ids)
    root_path = _resolve_sanity_check_path(args, 'merge', 'merge_dir')
    checked_species = [str(species).strip() for species in list(uni_species)]
    species_tokens = {_normalize_species_prefix(species): species for species in checked_species}
    issues = []
    if (not os.path.exists(root_path)) or (not os.path.isdir(root_path)):
        for species in checked_species:
            issues.append(_build_issue(
                check_name='merge',
                severity='error',
                issue_type='missing_output',
                target_type='species',
                target_id=species,
                path=root_path,
                message='merge directory was not found.',
                suggested_action='rerun_merge',
            ))
    else:
        metadata_out_path = os.path.join(root_path, 'metadata.tsv')
        if not os.path.isfile(metadata_out_path):
            metadata_error = 'merge metadata.tsv was not found.'
        else:
            metadata_error = _validate_nonempty_table(
                metadata_out_path,
                required_columns=['run'],
                context='merge metadata {}'.format(metadata_out_path),
                require_data_rows=True,
            )
        if metadata_error != '':
            issues.append(_build_issue(
                check_name='merge',
                severity='error',
                issue_type='invalid_content' if os.path.isfile(metadata_out_path) else 'missing_output',
                target_type='global',
                target_id='merge',
                path=metadata_out_path,
                message=metadata_error,
                suggested_action='rerun_merge',
            ))
        elif _is_stale_output(metadata_path, metadata_out_path):
            issues.append(_build_issue(
                check_name='merge',
                severity='warning',
                issue_type='stale_output',
                target_type='global',
                target_id='merge',
                path=metadata_out_path,
                message='merge metadata.tsv is older than the input metadata file.',
                suggested_action='review_or_rerun',
            ))
        summary_pdf_path = os.path.join(root_path, 'merge_summary.pdf')
        if not os.path.isfile(summary_pdf_path):
            issues.append(_build_issue(
                check_name='merge',
                severity='warning',
                issue_type='missing_output',
                target_type='global',
                target_id='merge_summary',
                path=summary_pdf_path,
                message='merge_summary.pdf was not found.',
                suggested_action='review_merge_plots',
            ))
        for token, species in species_tokens.items():
            species_dir = os.path.join(root_path, token)
            if not os.path.isdir(species_dir):
                issues.append(_build_issue(
                    check_name='merge',
                    severity='error',
                    issue_type='missing_output',
                    target_type='species',
                    target_id=species,
                    path=species_dir,
                    message='merge species directory was not found.',
                    suggested_action='rerun_merge',
                ))
                continue
            issues.extend(_validate_merge_species_tables(species, species_dir, metadata_path))
        expected_entries = list(species_tokens.keys()) + ['metadata.tsv', 'merge_summary.pdf']
        issues.extend(_validate_orphan_entries('merge', root_path, expected_entries, 'species'))
    rerun_manifest_path = _write_manifest(
        output_dir,
        'rerun_merge_species.txt',
        _collect_error_targets(issues, 'species'),
        'merge rerun species',
    )
    row = _build_sanity_summary_row('merge', 'species', checked_species, issues, root_path, '', rerun_manifest_path)
    return row, issues


def run_sanity_check_busco(args, metadata, uni_species, sra_ids, output_dir, metadata_path):
    _ = (metadata, sra_ids)
    root_path = _resolve_sanity_check_path(args, 'busco', 'busco_dir')
    checked_species = [str(species).strip() for species in list(uni_species)]
    species_tokens = {_normalize_species_prefix(species): species for species in checked_species}
    issues = []
    if (not os.path.exists(root_path)) or (not os.path.isdir(root_path)):
        for species in checked_species:
            issues.append(_build_issue(
                check_name='busco',
                severity='error',
                issue_type='missing_output',
                target_type='species',
                target_id=species,
                path=root_path,
                message='busco directory was not found.',
                suggested_action='rerun_busco',
            ))
    else:
        completeness_pdf = os.path.join(root_path, 'busco_completeness.pdf')
        if not os.path.isfile(completeness_pdf):
            issues.append(_build_issue(
                check_name='busco',
                severity='warning',
                issue_type='missing_output',
                target_type='global',
                target_id='busco_completeness',
                path=completeness_pdf,
                message='busco_completeness.pdf was not found.',
                suggested_action='review_busco_plots',
            ))
        for token, species in species_tokens.items():
            table_path = os.path.join(root_path, '{}_busco.tsv'.format(token))
            if not os.path.isfile(table_path):
                issues.append(_build_issue(
                    check_name='busco',
                    severity='error',
                    issue_type='missing_output',
                    target_type='species',
                    target_id=species,
                    path=table_path,
                    message='Normalized BUSCO table was not found.',
                    suggested_action='rerun_busco',
                ))
            else:
                table_error = _validate_busco_table(table_path)
                if table_error != '':
                    issues.append(_build_issue(
                        check_name='busco',
                        severity='error',
                        issue_type='invalid_content',
                        target_type='species',
                        target_id=species,
                        path=table_path,
                        message=table_error,
                        suggested_action='rerun_busco',
                    ))
                if _is_stale_output(metadata_path, table_path):
                    issues.append(_build_issue(
                        check_name='busco',
                        severity='warning',
                        issue_type='stale_output',
                        target_type='species',
                        target_id=species,
                        path=table_path,
                        message='BUSCO table is older than the metadata file.',
                        suggested_action='review_or_rerun',
                    ))
            species_dir = os.path.join(root_path, token)
            if not os.path.isdir(species_dir):
                issues.append(_build_issue(
                    check_name='busco',
                    severity='warning',
                    issue_type='missing_output',
                    target_type='species',
                    target_id=species,
                    path=species_dir,
                    message='Raw BUSCO species directory was not found.',
                    suggested_action='review_busco_output',
                ))
        expected_entries = list(species_tokens.keys()) + ['busco_completeness.pdf'] + ['{}_busco.tsv'.format(token) for token in species_tokens]
        issues.extend(_validate_orphan_entries('busco', root_path, expected_entries, 'species'))
    rerun_manifest_path = _write_manifest(
        output_dir,
        'rerun_busco_species.txt',
        _collect_error_targets(issues, 'species'),
        'busco rerun species',
    )
    row = _build_sanity_summary_row('busco', 'species', checked_species, issues, root_path, '', rerun_manifest_path)
    return row, issues


def _validate_finalize_species_outputs(species, species_dir, metadata_path):
    issues = []
    species_tag = _normalize_species_prefix(species)
    expected_files = {
        '{}_metadata.tsv'.format(species_tag): ['run'],
        '{}_expression.tsv'.format(species_tag): ['target_id'],
        '{}_batch_effect_summary.tsv'.format(species_tag): ['scientific_name'],
        '{}_curation_final_summary.tsv'.format(species_tag): ['scientific_name'],
    }
    for filename, required_columns in expected_files.items():
        file_path = os.path.join(species_dir, filename)
        if not os.path.isfile(file_path):
            issues.append(_build_issue(
                check_name='finalize',
                severity='error',
                issue_type='missing_output',
                target_type='species',
                target_id=species,
                path=file_path,
                message='Required finalize output was not found.',
                suggested_action='rerun_finalize',
            ))
            continue
        file_error = _validate_nonempty_table(
            file_path,
            required_columns=required_columns,
            context='finalize output {}'.format(file_path),
            require_data_rows=True,
            require_non_target_columns=filename.endswith('_expression.tsv'),
        )
        if file_error != '':
            issues.append(_build_issue(
                check_name='finalize',
                severity='error',
                issue_type='invalid_content',
                target_type='species',
                target_id=species,
                path=file_path,
                message=file_error,
                suggested_action='rerun_finalize',
            ))
        if _is_stale_output(metadata_path, file_path):
            issues.append(_build_issue(
                check_name='finalize',
                severity='warning',
                issue_type='stale_output',
                target_type='species',
                target_id=species,
                path=file_path,
                message='finalize output is older than the metadata file.',
                suggested_action='review_or_rerun',
            ))
    return issues


def run_sanity_check_finalize(args, metadata, uni_species, sra_ids, output_dir, metadata_path):
    _ = (metadata, sra_ids)
    root_path = _resolve_sanity_check_path(args, 'finalize', 'finalize_dir')
    checked_species = [str(species).strip() for species in list(uni_species)]
    species_tokens = {_normalize_species_prefix(species): species for species in checked_species}
    issues = []
    if (not os.path.exists(root_path)) or (not os.path.isdir(root_path)):
        for species in checked_species:
            issues.append(_build_issue(
                check_name='finalize',
                severity='error',
                issue_type='missing_output',
                target_type='species',
                target_id=species,
                path=root_path,
                message='finalize directory was not found.',
                suggested_action='rerun_finalize',
            ))
    else:
        metadata_out_path = os.path.join(root_path, 'metadata.tsv')
        if not os.path.isfile(metadata_out_path):
            metadata_error = 'finalize metadata.tsv was not found.'
        else:
            metadata_error = _validate_nonempty_table(
                metadata_out_path,
                required_columns=['run'],
                context='finalize metadata {}'.format(metadata_out_path),
                require_data_rows=True,
            )
        if metadata_error != '':
            issues.append(_build_issue(
                check_name='finalize',
                severity='error',
                issue_type='invalid_content' if os.path.isfile(metadata_out_path) else 'missing_output',
                target_type='global',
                target_id='finalize',
                path=metadata_out_path,
                message=metadata_error,
                suggested_action='rerun_finalize',
            ))
        elif _is_stale_output(metadata_path, metadata_out_path):
            issues.append(_build_issue(
                check_name='finalize',
                severity='warning',
                issue_type='stale_output',
                target_type='global',
                target_id='finalize',
                path=metadata_out_path,
                message='finalize metadata.tsv is older than the input metadata file.',
                suggested_action='review_or_rerun',
            ))
        exclusion_pdf = os.path.join(root_path, 'finalize_exclusion.pdf')
        if not os.path.isfile(exclusion_pdf):
            issues.append(_build_issue(
                check_name='finalize',
                severity='warning',
                issue_type='missing_output',
                target_type='global',
                target_id='finalize_exclusion',
                path=exclusion_pdf,
                message='finalize_exclusion.pdf was not found.',
                suggested_action='review_finalize_plots',
            ))
        for token, species in species_tokens.items():
            species_dir = os.path.join(root_path, token)
            if not os.path.isdir(species_dir):
                issues.append(_build_issue(
                    check_name='finalize',
                    severity='error',
                    issue_type='missing_output',
                    target_type='species',
                    target_id=species,
                    path=species_dir,
                    message='finalize species directory was not found.',
                    suggested_action='rerun_finalize',
                ))
                continue
            issues.extend(_validate_finalize_species_outputs(species, species_dir, metadata_path))
        expected_entries = list(species_tokens.keys()) + ['metadata.tsv', 'finalize_exclusion.pdf']
        issues.extend(_validate_orphan_entries('finalize', root_path, expected_entries, 'species'))
    rerun_manifest_path = _write_manifest(
        output_dir,
        'rerun_finalize_species.txt',
        _collect_error_targets(issues, 'species'),
        'finalize rerun species',
    )
    row = _build_sanity_summary_row('finalize', 'species', checked_species, issues, root_path, '', rerun_manifest_path)
    return row, issues


def sanity_main(args):
    out_dir = os.path.realpath(args.out_dir)
    if os.path.exists(out_dir) and (not os.path.isdir(out_dir)):
        raise NotADirectoryError('Output path exists but is not a directory: {}'.format(out_dir))
    metadata = load_metadata(args)
    output_dir = os.path.join(out_dir, 'sanity')
    if os.path.exists(output_dir) and (not os.path.isdir(output_dir)):
        raise NotADirectoryError('Sanity output path exists but is not a directory: {}'.format(output_dir))
    previous_payload = _load_previous_sanity_report(output_dir) if os.path.isdir(output_dir) else None
    os.makedirs(output_dir, exist_ok=True)
    metadata_path = _resolve_metadata_path(args)
    run_filters, species_filters = _filter_metadata_for_sanity(args, metadata)
    uni_species, sra_ids = parse_metadata(args, metadata)
    requested_checks = _resolve_requested_sanity_checks(args)
    runner_by_check = {
        'getfastq': run_sanity_check_getfastq,
        'index': run_sanity_check_index,
        'quant': run_sanity_check_quant,
        'merge': run_sanity_check_merge,
        'busco': run_sanity_check_busco,
        'finalize': run_sanity_check_finalize,
    }
    summary_rows = []
    issues = []
    for check_name in requested_checks:
        row, check_issues = runner_by_check[check_name](
            args=args,
            metadata=metadata,
            uni_species=uni_species,
            sra_ids=sra_ids,
            output_dir=output_dir,
            metadata_path=metadata_path,
        )
        summary_rows.append(row)
        issues.extend(check_issues)
    summary_path = _write_sanity_summary(output_dir, summary_rows)
    issues_path = _write_sanity_issues(output_dir, issues)
    comparison_path = _write_sanity_comparison(output_dir, previous_payload, summary_rows)
    report_payload = {
        'generated_at': datetime.datetime.now().isoformat(),
        'metadata_path': metadata_path,
        'out_dir': out_dir,
        'run_filters': run_filters,
        'species_filters': species_filters,
        'requested_checks': requested_checks,
        'summary': summary_rows,
        'issues': issues,
        'summary_path': summary_path,
        'issues_path': issues_path,
        'comparison_path': comparison_path,
    }
    report_path = _write_sanity_report_json(output_dir, report_payload)
    _print_sanity_summary(summary_rows, summary_path)
    print('issues_path={}'.format(issues_path))
    if comparison_path != '':
        print('comparison_path={}'.format(comparison_path))
    print('report_path={}'.format(report_path))
    _raise_if_sanity_issues(args, issues)
