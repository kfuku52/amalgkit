import pandas
import ete4

import os
from amalgkit.download_utils import (
    acquire_exclusive_lock as _acquire_exclusive_lock,
    get_ete_ncbitaxa as _get_ete_ncbitaxa,
    resolve_download_dir as _resolve_download_dir,
    resolve_ete_data_dir as _resolve_ete_data_dir,
)
from amalgkit.metadata_utils import (
    Metadata as _Metadata,
    detect_layout_from_file as _detect_layout_from_file,
    get_mapping_rate as _get_mapping_rate,
    get_metadata_row_index_by_run as _get_metadata_row_index_by_run,
    get_newest_intermediate_file_extension as _get_newest_intermediate_file_extension,
    get_sra_stat as _get_sra_stat,
    is_there_unpaired_file as _is_there_unpaired_file,
    load_metadata as _load_metadata,
    parse_bool_flags as _parse_bool_flags,
    strtobool as _strtobool,
    write_updated_metadata as _write_updated_metadata,
)
from amalgkit.orthology_utils import (
    append_unique_busco_ids as _append_unique_busco_ids,
    check_ortholog_parameter_compatibility as _check_ortholog_parameter_compatibility,
    generate_multisp_busco_table as _generate_multisp_busco_table,
    orthogroup2genecount as _orthogroup2genecount,
    parse_busco_species_name as _parse_busco_species_name,
    parse_busco_species_table as _parse_busco_species_table,
    read_busco_species_table as _read_busco_species_table,
    update_busco_meta as _update_busco_meta,
)
from amalgkit.parallel_utils import (
    is_auto_parallel_option as _is_auto_parallel_option,
    resolve_cpu_budget as _resolve_cpu_budget,
    resolve_detected_cpu_count as _resolve_detected_cpu_count,
    resolve_thread_worker_allocation as _resolve_thread_worker_allocation,
    resolve_total_core_budget as _resolve_total_core_budget,
    resolve_worker_allocation as _resolve_worker_allocation,
    run_tasks_with_optional_threads as _run_tasks_with_optional_threads,
    validate_positive_int_option as _validate_positive_int_option,
)
from amalgkit.prefix_utils import (
    find_prefixed_entries as _find_prefixed_entries,
    find_run_prefixed_entries as _find_run_prefixed_entries,
    find_species_prefixed_entries as _find_species_prefixed_entries,
)
from amalgkit.runtime_utils import (
    check_config_dir as _check_config_dir,
    cleanup_tmp_amalgkit_files as _cleanup_tmp_amalgkit_files,
    get_getfastq_run_dir as _get_getfastq_run_dir,
)

Metadata = _Metadata
get_metadata_row_index_by_run = _get_metadata_row_index_by_run
get_sra_stat = _get_sra_stat
load_metadata = _load_metadata
parse_bool_flags = _parse_bool_flags
strtobool = _strtobool
is_auto_parallel_option = _is_auto_parallel_option
resolve_cpu_budget = _resolve_cpu_budget
resolve_detected_cpu_count = _resolve_detected_cpu_count
resolve_thread_worker_allocation = _resolve_thread_worker_allocation
resolve_total_core_budget = _resolve_total_core_budget
resolve_worker_allocation = _resolve_worker_allocation
run_tasks_with_optional_threads = _run_tasks_with_optional_threads
validate_positive_int_option = _validate_positive_int_option
find_prefixed_entries = _find_prefixed_entries
find_run_prefixed_entries = _find_run_prefixed_entries
find_species_prefixed_entries = _find_species_prefixed_entries
check_config_dir = _check_config_dir
cleanup_tmp_amalgkit_files = _cleanup_tmp_amalgkit_files
get_getfastq_run_dir = _get_getfastq_run_dir

def resolve_download_dir(args):
    return _resolve_download_dir(args)

def resolve_ete_data_dir(args):
    return _resolve_ete_data_dir(args, resolve_download_dir_fn=resolve_download_dir)


acquire_exclusive_lock = _acquire_exclusive_lock


def get_ete_ncbitaxa(args=None):
    return _get_ete_ncbitaxa(
        args=args,
        acquire_exclusive_lock_fn=acquire_exclusive_lock,
        ncbitaxa_cls=ete4.NCBITaxa,
        resolve_ete_data_dir_fn=resolve_ete_data_dir,
    )

def read_config_file(file_name, dir_path):
    config_path = os.path.join(dir_path, file_name)
    if os.path.exists(config_path) and (not os.path.isfile(config_path)):
        raise IsADirectoryError('Config path exists but is not a file: {}'.format(config_path))
    try:
        df = pandas.read_csv(config_path,
                             parse_dates=False, quotechar='"', sep='\t',
                             header=None, index_col=None, skip_blank_lines=True, comment='#')
    except (FileNotFoundError, pandas.errors.EmptyDataError):
        df = pandas.DataFrame([])
    if df.shape[1]==1:
        df = df.iloc[:,0]
    return df

get_newest_intermediate_file_extension = _get_newest_intermediate_file_extension
is_there_unpaired_file = _is_there_unpaired_file
detect_layout_from_file = _detect_layout_from_file


def get_mapping_rate(metadata, quant_dir, max_workers='auto'):
    return _get_mapping_rate(
        metadata=metadata,
        quant_dir=quant_dir,
        max_workers=max_workers,
        run_tasks_with_optional_threads_fn=run_tasks_with_optional_threads,
        is_auto_parallel_option_fn=is_auto_parallel_option,
        validate_positive_int_option_fn=validate_positive_int_option,
    )


def write_updated_metadata(metadata, outpath, args, max_workers='auto'):
    return _write_updated_metadata(
        metadata=metadata,
        outpath=outpath,
        args=args,
        max_workers=max_workers,
        get_mapping_rate_fn=get_mapping_rate,
    )

orthogroup2genecount = _orthogroup2genecount
check_ortholog_parameter_compatibility = _check_ortholog_parameter_compatibility
parse_busco_species_name = _parse_busco_species_name
read_busco_species_table = _read_busco_species_table
parse_busco_species_table = _parse_busco_species_table
append_unique_busco_ids = _append_unique_busco_ids
update_busco_meta = _update_busco_meta


def generate_multisp_busco_table(dir_busco, outfile):
    return _generate_multisp_busco_table(
        dir_busco=dir_busco,
        outfile=outfile,
        run_tasks_with_optional_threads=run_tasks_with_optional_threads,
    )
