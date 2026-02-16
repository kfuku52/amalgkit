from bisect import bisect_left

import numpy as np
from amalgkit.util import *


def list_duplicates(seq):
    seen = set()
    seen_add = seen.add
    # adds all elements it doesn't know yet to seen and all other to seen_twice
    seen_twice = set(x for x in seq if x in seen or seen_add(x))
    # turn the set into a list (as requested)
    return list(seen_twice)


def _find_prefixed_entries(sorted_entries, prefix):
    left = bisect_left(sorted_entries, prefix)
    right = bisect_left(sorted_entries, prefix + '\uffff')
    matched = []
    for i in range(left, right):
        entry = sorted_entries[i]
        if entry.startswith(prefix):
            matched.append(entry)
    return matched


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
    if os.path.exists(getfastq_path):
        print("amalgkit getfastq output folder detected. Checking presence of output files.")
        target_runs = set(sra_ids)
        verbose_run_logs = _should_log_per_run(args, len(target_runs))
        _print_run_log_mode(args, len(target_runs))
        run_files_map, non_dir_runs = _scan_target_run_dirs(getfastq_path, target_runs)
        for sra_id in sra_ids:
            if verbose_run_logs:
                print("\n")
                print("Looking for {}".format(sra_id))
            sra_path = os.path.join(getfastq_path, sra_id)
            run_files = run_files_map.get(sra_id)
            if (run_files is None) or (sra_id in non_dir_runs):
                if verbose_run_logs:
                    print("Could not find getfastq output for: ", sra_id, "\n")
                    print("Suggested command for rerun: getfastq -e email@adress.com --id ", sra_id, " -w ", args.out_dir, "--redo yes --gcp yes --aws yes --ncbi yes")
                data_unavailable.append(sra_id)
                continue
            sra_stat = get_sra_stat(sra_id, metadata)
            try:
                ext = get_newest_intermediate_file_extension(sra_stat, sra_path, files=run_files)
            except FileNotFoundError:
                if verbose_run_logs:
                    print("could not find any fastq files for ", sra_id,
                          "Please make sure amalgkit getfastq ran properly")
                data_unavailable.append(sra_id)
                continue
            if ext == 'no_extension_found':
                if verbose_run_logs:
                    print("could not find any fastq files for ", sra_id,
                          "Please make sure amalgkit getfastq ran properly")
                data_unavailable.append(sra_id)
                continue

            files = sorted([
                os.path.join(sra_path, f)
                for f in run_files
                if f.startswith(sra_id) and f.endswith(ext)
            ])
            if verbose_run_logs and (ext != '.safely_removed'):
                print("Found:", files)
            data_available.append(sra_id)

    else:
        print("Could not find getfastq output folder ", getfastq_path, ". Have you run getfastq yet?")
        data_unavailable = metadata.df['run'].tolist()

    if data_unavailable:
        print("writing SRA IDs without getfastq output to: ", os.path.join(output_dir, "SRA_IDs_without_fastq.txt"))
        with open(os.path.join(output_dir, "SRA_IDs_without_fastq.txt"), "w") as f:
            for sra_id in data_unavailable:
                f.write(sra_id + "\n")
    else:
        txt = "The getfastq output files for all SRA IDs in --metadata ({}) were found."
        print(txt.format(args.metadata))

    return data_available, data_unavailable


def check_quant_index(args, uni_species, output_dir):
    if args.index_dir:
        index_dir_path = args.index_dir
    else:
        index_dir_path = os.path.join(args.out_dir, "index")
    index_unavailable = []
    index_available = []
    if os.path.exists(index_dir_path):
        index_entries = sorted(os.listdir(index_dir_path))

        def find_index_files(prefix):
            matched = _find_prefixed_entries(index_entries, prefix)
            return [os.path.join(index_dir_path, entry) for entry in matched]

        for species in uni_species:
            sci_name = species.replace(" ", "_")
            sci_name = sci_name.replace(".", "")
            index_path = os.path.join(index_dir_path, sci_name + "*")
            print("\n")
            print("Looking for index file {} for species {}".format(index_path, species))
            index_files = find_index_files(sci_name)
            if not index_files:
                print("could not find anything in", index_path)
                sci_name = species.split(" ")
                # Deprecate subspecies or variants and look again
                # I.e. if Gorilla_gorilla_gorilla.idx was not found, we look for Gorilla_gorilla.idx instead.
                if len(sci_name) > 2:
                    sci_name = sci_name[0] + "_" + sci_name[1]
                    print("Ignoring subspecies.")
                    index_path = os.path.join(index_dir_path, sci_name + "*")
                    print("Looking for {}".format(index_path))
                    index_files = find_index_files(sci_name)
                    if index_files:
                        print("Found ", index_files, "!")
                        index_available.append(species)
                    else:
                        print("Could not find any index files for ", species)
                        index_unavailable.append(species)
                else:
                    print("Could not find any index files for ", species)
                    index_unavailable.append(species)

            else:
                print("Found ", index_files, "!")
                index_available.append(species)

            if len(index_files) > 1:
                print("Multiple possible index files detected for ", species, ": ", index_files,
                      ". You may have to resolve ambiguity")

        if index_unavailable:
            print("writing species without index to: ", os.path.join(output_dir, "species_without_index.txt"))
            with open(os.path.join(output_dir, "species_without_index.txt"), "w") as f:
                for species in index_unavailable:
                    f.write(species + "\n")
        else:
            print("index found for all species in --metadata ({})".format(args.metadata))
    else:
        print("Could not find index directory ", index_dir_path, " . Did you provide the correct Path?")

    return index_available, index_unavailable


def check_quant_output(args, sra_ids, output_dir):
    print("checking for quant outputs: ")
    quant_path = os.path.join(args.out_dir, "quant")
    data_available = []
    data_unavailable = []

    if os.path.exists(quant_path):
        print("amalgkit quant output folder detected. Checking presence of output files.")
        target_runs = set(sra_ids)
        verbose_run_logs = _should_log_per_run(args, len(target_runs))
        _print_run_log_mode(args, len(target_runs))
        quant_run_files_map, non_dir_runs = _scan_target_run_dirs(quant_path, target_runs)
        for sra_id in sra_ids:
            if verbose_run_logs:
                print("\n")
                print("Looking for {}".format(sra_id))
            sra_path = os.path.join(quant_path, sra_id)
            quant_run_files = quant_run_files_map.get(sra_id)
            if (quant_run_files is None) or (sra_id in non_dir_runs):
                if verbose_run_logs:
                    print("Could not find output folder ", sra_path, " for ", sra_id)
                data_unavailable.append(sra_id)
                continue

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
                data_available.append(sra_id)
                continue

            if verbose_run_logs and (not has_abundance):
                print(os.path.join(sra_path, abundance_file), " is missing! Please check if quant ran correctly")
            if verbose_run_logs and (not has_run_info):
                print(os.path.join(sra_path, run_info_file), " is missing! Please check if quant ran correctly")
            data_unavailable.append(sra_id)
    else:
        print("Could not find quant output folder ", quant_path, ". Have you run quant yet?")

    if data_unavailable:
        print("writing SRA IDs without quant output to: ", os.path.join(output_dir, "SRA_IDs_without_quant.txt"))
        with open(os.path.join(output_dir, "SRA_IDs_without_quant.txt"), "w") as f:
            for sra_id in data_unavailable:
                f.write(sra_id + "\n")
    else:
        print("Quant outputs found for all SRA IDs in --metadata ({})".format(args.metadata))

    return data_available, data_unavailable


def sanity_main(args):
    metadata = load_metadata(args)
    output_dir = os.path.join(args.out_dir, 'sanity')
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    uni_species, sra_ids = parse_metadata(args, metadata)
    if args.getfastq or args.all:
        check_getfastq_outputs(args, sra_ids, metadata, output_dir)
    if args.index or args.all:
        check_quant_index(args, uni_species, output_dir)
    if args.quant or args.all:
        check_quant_output(args, sra_ids, output_dir)
