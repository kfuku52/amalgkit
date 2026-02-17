import numpy as np
from amalgkit.util import *


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
                        run_files_map[run_id] = {run_entry.name for run_entry in run_entries}
                except (FileNotFoundError, NotADirectoryError):
                    pass
    except FileNotFoundError:
        pass
    return run_files_map, non_dir_runs

def _normalize_species_prefix(species):
    return species.replace(" ", "_").replace(".", "")

def _get_species_fallback_prefix(species):
    parts = species.split(" ")
    if len(parts) <= 2:
        return None
    return _normalize_species_prefix(parts[0] + " " + parts[1])

def _should_log_per_run(args, num_runs):
    if bool(getattr(args, 'quiet', False)):
        return False
    verbose_runs = int(getattr(args, 'verbose_runs', 20))
    if verbose_runs < 0:
        return True
    return num_runs <= verbose_runs

def _print_run_log_mode(args, num_runs):
    if _should_log_per_run(args, num_runs):
        return
    verbose_runs = int(getattr(args, 'verbose_runs', 20))
    print('Per-run logs suppressed for {:,} runs (--verbose_runs={}, --quiet={}).'.format(
        num_runs,
        verbose_runs,
        bool(getattr(args, 'quiet', False)),
    ))

def _write_unavailable_items(output_dir, filename, item_ids, label):
    if not item_ids:
        return
    outpath = os.path.join(output_dir, filename)
    print("writing {} to: {}".format(label, outpath))
    with open(outpath, "w") as f:
        for item_id in item_ids:
            f.write(item_id + "\n")

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
    verbose_run_logs = _should_log_per_run(args, len(target_runs))
    _print_run_log_mode(args, len(target_runs))
    run_files_map, non_dir_runs = _scan_target_run_dirs(root_path, target_runs)
    return verbose_run_logs, run_files_map, non_dir_runs


def _print_getfastq_missing_output_message(args, sra_id):
    print("Could not find getfastq output for: ", sra_id, "\n")
    print(
        "Suggested command for rerun: getfastq -e email@adress.com --id ",
        sra_id,
        " -w ",
        args.out_dir,
        "--redo yes --gcp yes --aws yes --ncbi yes",
    )


def _check_single_getfastq_run(args, sra_id, metadata, getfastq_path, run_files_map, non_dir_runs, verbose_run_logs):
    sra_path = os.path.join(getfastq_path, sra_id)
    run_files = run_files_map.get(sra_id)
    if (run_files is None) or (sra_id in non_dir_runs):
        if verbose_run_logs:
            _print_getfastq_missing_output_message(args, sra_id)
        return False

    sra_stat = get_sra_stat(sra_id, metadata)
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
            for f in run_files
            if f.startswith(sra_id) and f.endswith(ext)
        ])
        print("Found:", files)
    return True


def parse_metadata(args, metadata):
    print("Checking essential entries from metadata file.")
    species = metadata.df.loc[:, 'scientific_name']
    uni_species = np.unique(species)
    if len(uni_species):
        print(len(uni_species), " species detected:")
        print(uni_species)
    else:
        txt = "{} species detected. Please check if --metadata ({}) has a 'scientific_name' column."
        raise ValueError(txt.format(len(uni_species), args.metadata))

    sra_ids = metadata.df.loc[:, 'run']
    uni_sra_ids = np.unique(sra_ids)
    if len(sra_ids):
        print(len(uni_sra_ids), " SRA runs detected:")
        print(uni_sra_ids)
        # check for duplicate runs
        if len(sra_ids) > len(uni_sra_ids):
            dupes = list_duplicates(sra_ids)
            raise ValueError(
                "Duplicate SRA IDs detected, where IDs should be unique. Please check these entries: {}".format(dupes)
            )
    else:
        txt = "{} SRA runs detected. Please check if --metadata ({}) has a 'run' column."
        raise ValueError(txt.format(len(uni_sra_ids), args.metadata))
    return uni_species, sra_ids


def check_getfastq_outputs(args, sra_ids, metadata, output_dir):
    print("checking for getfastq outputs: ")
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
        for sra_id in sra_ids:
            if verbose_run_logs:
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
    matched = find_prefixed_entries(index_entries, prefix)
    return [os.path.join(index_dir_path, entry) for entry in matched]

def _resolve_species_index_files(species, index_entries, index_dir_path):
    sci_name = _normalize_species_prefix(species)
    index_path = os.path.join(index_dir_path, sci_name + "*")
    print("\n")
    print("Looking for index file {} for species {}".format(index_path, species))
    index_files = _find_index_files(index_entries, index_dir_path, sci_name)
    if index_files:
        print("Found ", index_files, "!")
        return index_files

    print("could not find anything in", index_path)
    # Deprecate subspecies or variants and look again
    # I.e. if Gorilla_gorilla_gorilla.idx was not found, we look for Gorilla_gorilla.idx instead.
    fallback_prefix = _get_species_fallback_prefix(species)
    if fallback_prefix is None:
        print("Could not find any index files for ", species)
        return []
    print("Ignoring subspecies.")
    index_path = os.path.join(index_dir_path, fallback_prefix + "*")
    print("Looking for {}".format(index_path))
    index_files = _find_index_files(index_entries, index_dir_path, fallback_prefix)
    if index_files:
        print("Found ", index_files, "!")
        return index_files
    print("Could not find any index files for ", species)
    return []


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
            return index_available, index_unavailable
        index_entries = sorted(os.listdir(index_dir_path))
        for species in uni_species:
            index_files = _resolve_species_index_files(
                species=species,
                index_entries=index_entries,
                index_dir_path=index_dir_path,
            )
            if index_files:
                index_available.append(species)
            else:
                index_unavailable.append(species)
            if len(index_files) > 1:
                print("Multiple possible index files detected for ", species, ": ", index_files,
                      ". You may have to resolve ambiguity")

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

    return index_available, index_unavailable


def check_quant_output(args, sra_ids, output_dir):
    print("checking for quant outputs: ")
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
        for sra_id in sra_ids:
            if verbose_run_logs:
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


def sanity_main(args):
    metadata = load_metadata(args)
    output_dir = os.path.join(args.out_dir, 'sanity')
    os.makedirs(output_dir, exist_ok=True)
    uni_species, sra_ids = parse_metadata(args, metadata)
    if args.getfastq or args.all:
        check_getfastq_outputs(args, sra_ids, metadata, output_dir)
    if args.index or args.all:
        check_quant_index(args, uni_species, output_dir)
    if args.quant or args.all:
        check_quant_output(args, sra_ids, output_dir)
