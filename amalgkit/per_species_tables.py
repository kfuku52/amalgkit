import datetime
import os
import re
import shutil
import subprocess
import sys
from types import SimpleNamespace

import pandas

from amalgkit.arg_utils import clone_namespace
from amalgkit.command_context import PerSpeciesTableContext
from amalgkit.metadata_utils import load_metadata
from amalgkit.parallel_utils import resolve_worker_allocation, run_tasks_with_optional_threads
from amalgkit.r_config import temporary_r_config
from amalgkit.runtime_utils import check_rscript
from amalgkit.subprocess_utils import run_logged_command


def validate_per_species_metadata_columns(metadata, required_columns, context):
    missing = [col for col in required_columns if col not in metadata.df.columns]
    if len(missing) > 0:
        raise ValueError(
            'Missing required metadata column(s) for {}: {}'.format(
                context,
                ', '.join(missing),
            )
        )


def get_sample_group(args, metadata):
    def normalize_sample_groups(values):
        groups = []
        for value in values:
            if value is None:
                continue
            if pandas.isna(value):
                continue
            normalized = str(value).strip()
            if normalized == '':
                continue
            groups.append(normalized)
        return list(dict.fromkeys(groups))

    def parse_sample_group_argument(sample_group_arg):
        return normalize_sample_groups(re.split(r'[,\|]+', str(sample_group_arg)))

    if args.sample_group is None:
        if 'sample_group' not in metadata.df.columns:
            txt = 'The "sample_group" column was not found in --metadata ({}). '
            txt += 'Please add this column or provide --sample_group.'
            raise ValueError(txt.format(args.metadata))
        sample_group = normalize_sample_groups(metadata.df.loc[:, 'sample_group'].tolist())
    else:
        sample_group = parse_sample_group_argument(args.sample_group)
    if (len(sample_group)==0):
        txt = 'The "sample_group" column in --metadata ({}) is not filled. Please manually edit the file. '
        txt += 'Per-species table generation recognizes samples with the same string in this column to belong to the same group.'
        raise ValueError(txt.format(args.metadata))
    print('Tissues to be included: {}'.format(', '.join(sample_group)))
    sample_group = '|'.join(sample_group)
    return sample_group

def get_completion_flag_path(per_species_dir, sp):
    return os.path.join(per_species_dir, sp, 'per_species_completion_flag.txt')

def write_completion_flag(per_species_dir, sp):
    completion_flag_path = get_completion_flag_path(per_species_dir, sp)
    with open(completion_flag_path, 'w') as f:
        f.write('amalgkit per-species table generation completed at {}\n'.format(datetime.datetime.now()))

def register_exit_status(per_species_dir, sp, exit_status, failed_species):
    if exit_status == 0:
        write_completion_flag(per_species_dir=per_species_dir, sp=sp)
    else:
        failed_species.add(sp)
        sys.stderr.write('Per-species table generation failed for species: {}\n'.format(sp))

def run_per_species_r_script(args, metadata, sp, input_dir):
    dist_method = args.dist_method
    mr_cut = args.mapping_rate
    correlation_threshold = args.correlation_threshold
    intermediate = args.plot_intermediate
    sample_group = get_sample_group(args, metadata)
    outlier_method = getattr(args, 'outlier_method', 'legacy')
    margin_threshold = float(getattr(args, 'margin_threshold', 0.0))
    robust_z_threshold = float(getattr(args, 'robust_z_threshold', -2.5))
    disable_auto_outlier_filter = bool(getattr(args, 'disable_auto_outlier_filter', False))
    dir_amalgkit_script = os.path.dirname(os.path.realpath(__file__))
    requested_script = str(getattr(args, 'r_script_name', 'prepare_tables.r'))
    if os.path.isabs(requested_script):
        r_script_path = requested_script
    else:
        r_script_path = os.path.join(dir_amalgkit_script, requested_script)
    r_script_name = os.path.basename(r_script_path)
    if not os.path.isfile(r_script_path):
        raise FileNotFoundError('R script not found: {}'.format(r_script_path))
    r_util_path = os.path.join(dir_amalgkit_script, 'util.r')
    metadata_arg = getattr(args, 'metadata', 'inferred')
    if metadata_arg != 'inferred':
        metadata_input_path = os.path.realpath(metadata_arg)
    else:
        metadata_input_path = os.path.join(input_dir, 'metadata.tsv')
    if not os.path.exists(metadata_input_path):
        raise FileNotFoundError('Metadata file not found for per-species R script: {}'.format(metadata_input_path))
    if not os.path.isfile(metadata_input_path):
        raise IsADirectoryError('Metadata path exists but is not a file: {}'.format(metadata_input_path))
    input_dir_abs = os.path.abspath(input_dir)
    len_file = os.path.join(input_dir_abs, sp, sp + '_eff_length.tsv')
    count_file_candidates = [
        os.path.join(input_dir_abs, sp, sp + '_cstmm_counts.tsv'),
        os.path.join(input_dir_abs, sp, sp + '_est_counts.tsv'),
    ]
    existing_count_files = [path for path in count_file_candidates if os.path.isfile(path)]
    if len(existing_count_files) >= 1:
        count_file = existing_count_files[0]
    else:
        count_file = count_file_candidates[0]
    has_count_file = os.path.isfile(count_file)
    has_len_file = os.path.isfile(len_file)
    if has_count_file and has_len_file:
        print("Both counts and effective length files found: {}".format(sp), flush=True)
    else:
        if not has_count_file:
            for expected_count_file in count_file_candidates:
                if os.path.exists(expected_count_file) and (not os.path.isfile(expected_count_file)):
                    sys.stderr.write('Count path exists but is not a file: {}\n'.format(expected_count_file))
                elif not os.path.exists(expected_count_file):
                    sys.stderr.write('Expected but undetected PATH of the count file: {}\n'.format(expected_count_file))
        if not has_len_file:
            if os.path.exists(len_file) and (not os.path.isfile(len_file)):
                sys.stderr.write('Effective length path exists but is not a file: {}\n'.format(len_file))
            else:
                sys.stderr.write('Expected but undetected PATH of the effective length file: {}\n'.format(len_file))
        sys.stderr.write('Skipping {}\n'.format(sp))
        print('Skipping {}'.format(sp), flush=True)
        return 1
    print("Starting Rscript ({}) to generate per-species {} tables.".format(r_script_name, args.norm), flush=True)
    config_map = {
        'est_counts_path': count_file,
        'metadata_path': metadata_input_path,
        'out_dir': os.path.realpath(args.out_dir),
        'eff_length_path': len_file,
        'dist_method': dist_method,
        'mapping_rate_cutoff': mr_cut,
        'min_dif': 0,
        'plot_intermediate': int(intermediate),
        'selected_sample_groups': sample_group,
        'sample_group_colors': args.sample_group_color,
        'transform_method': str(args.norm),
        'one_outlier_per_iteration': int(args.one_outlier_per_iter),
        'correlation_threshold': correlation_threshold,
        'batch_effect_alg': str(getattr(args, 'batch_effect_alg', 'no')),
        'clip_negative': int(getattr(args, 'clip_negative', True)),
        'maintain_zero': int(getattr(args, 'maintain_zero', True)),
        'r_util_path': os.path.realpath(r_util_path),
        'skip_curation_flag': int(getattr(args, 'skip_curation', False)),
        'outlier_method': str(outlier_method),
        'robust_margin_threshold': str(margin_threshold),
        'robust_z_threshold': str(robust_z_threshold),
        'disable_auto_outlier_filter_flag': int(disable_auto_outlier_filter),
    }
    if r_script_name == 'wsfilter.r':
        config_map['batch_effect_alg'] = 'no'
        config_map['clip_negative'] = 1
        config_map['maintain_zero'] = 1
        config_map['skip_curation_flag'] = 0
        config_map['outlier_method'] = 'robust_margin'
        config_map['disable_auto_outlier_filter_flag'] = 0
    elif r_script_name == 'finalize.r':
        ruvseq_control_genes = str(getattr(args, 'ruvseq_control_genes', 'auto'))
        ruvseq_k = str(getattr(args, 'ruvseq_k', 'auto'))
        ruvseq_k_max = int(getattr(args, 'ruvseq_k_max', 5))
        ruvseq_control_top_n = int(getattr(args, 'ruvseq_control_top_n', 1000))
        ruvseq_min_controls = int(getattr(args, 'ruvseq_min_controls', 100))
        random_seed = str(getattr(args, 'seed', 'auto'))
        sva_nsv = str(getattr(args, 'sva_nsv', 'auto'))
        sva_B = str(getattr(args, 'sva_B', 'auto'))
        sva_B_auto_max = int(getattr(args, 'sva_B_auto_max', 100))
        config_map['mapping_rate_cutoff'] = 0
        config_map['plot_intermediate'] = 0
        config_map['one_outlier_per_iteration'] = 0
        config_map['correlation_threshold'] = 0.3
        config_map['skip_curation_flag'] = 0
        config_map['outlier_method'] = 'legacy'
        config_map['robust_margin_threshold'] = 0
        config_map['robust_z_threshold'] = -2.5
        config_map['disable_auto_outlier_filter_flag'] = 1
        config_map['ruvseq_control_mode'] = ruvseq_control_genes
        config_map['ruvseq_k_setting'] = ruvseq_k
        config_map['ruvseq_k_max'] = ruvseq_k_max
        config_map['ruvseq_control_top_n'] = ruvseq_control_top_n
        config_map['ruvseq_min_controls'] = ruvseq_min_controls
        config_map['random_seed_setting'] = random_seed
        config_map['sva_nsv_setting'] = sva_nsv
        config_map['sva_B_setting'] = sva_B
        config_map['sva_B_auto_max'] = sva_B_auto_max
    elif r_script_name == 'prepare_tables.r':
        pass
    else:
        raise ValueError('Unsupported R script for per-species table pipeline: {}'.format(r_script_name))
    with temporary_r_config(config_map, prefix='amalgkit_per_species_r_') as config_path:
        r_result, _stdout_txt, _stderr_txt = run_logged_command(
            command=['Rscript', r_script_path, config_path],
            runner=subprocess.run,
            command_prefix='Rscript command',
            print_output=False,
            not_found_label='Rscript',
        )
    return r_result.returncode


def resolve_per_species_execution_context(args, context=None):
    if isinstance(context, PerSpeciesTableContext):
        return context.resolve()
    return resolve_per_species_input(args)

def resolve_per_species_input(args):
    if args.input_dir != 'inferred':
        input_dir = os.path.realpath(args.input_dir)
        print('Input_directory: {}'.format(input_dir))
        if not os.path.exists(input_dir):
            raise FileNotFoundError('Input directory not found: {}'.format(input_dir))
        if not os.path.isdir(input_dir):
            raise NotADirectoryError('Input path exists but is not a directory: {}'.format(input_dir))
        if args.metadata == 'inferred':
            input_metadata_path = os.path.join(input_dir, 'metadata.tsv')
            if not os.path.exists(input_metadata_path):
                raise FileNotFoundError(
                    'metadata.tsv not found in --input_dir: {}. '
                    'Provide --metadata explicitly or place metadata.tsv in the input directory.'.format(
                        input_metadata_path
                    )
                )
            proxy_args = SimpleNamespace(**vars(args))
            proxy_args.metadata = input_metadata_path
            metadata = load_metadata(proxy_args, batch_scope='species')
        else:
            metadata = load_metadata(args, dir_subcommand=os.path.basename(input_dir), batch_scope='species')
        return metadata, input_dir

    dir_merge = os.path.realpath(os.path.join(args.out_dir, 'merge'))
    dir_cstmm = os.path.realpath(os.path.join(args.out_dir, 'cstmm'))
    if os.path.exists(dir_cstmm):
        if not os.path.isdir(dir_cstmm):
            raise NotADirectoryError('cstmm input path exists but is not a directory: {}'.format(dir_cstmm))
        print('Subdirectory for amalgkit cstmm will be used as input: {}'.format(dir_cstmm))
        metadata = load_metadata(args, dir_subcommand='cstmm', batch_scope='species')
        return metadata, dir_cstmm
    if os.path.exists(dir_merge) and (not os.path.isdir(dir_merge)):
        raise NotADirectoryError('merge input path exists but is not a directory: {}'.format(dir_merge))
    print('Subdirectory for amalgkit merge will be used as input: {}'.format(dir_merge))
    metadata = load_metadata(args, dir_subcommand='merge', batch_scope='species')
    return metadata, dir_merge


def list_selected_species(metadata):
    validate_per_species_metadata_columns(
        metadata=metadata,
        required_columns=['scientific_name', 'exclusion'],
        context='per-species table generation',
    )
    exclusion_series = (
        metadata.df['exclusion']
        .fillna('')
        .astype(str)
        .str.strip()
        .str.lower()
    )
    is_selected = (exclusion_series == 'no')
    selected_species = (
        metadata.df
        .loc[is_selected, 'scientific_name']
        .fillna('')
        .astype(str)
        .str.strip()
    )
    selected_species = selected_species.loc[selected_species != ''].drop_duplicates().values
    return [species.replace(' ', '_') for species in selected_species]


def has_cstmm_counts_input(input_dir, selected_species):
    input_dir_real = os.path.realpath(input_dir)
    if os.path.basename(input_dir_real) == 'cstmm':
        return True
    for species in selected_species:
        cstmm_count_path = os.path.join(input_dir_real, species, species + '_cstmm_counts.tsv')
        if os.path.isfile(cstmm_count_path):
            return True
    return False


def collect_pending_species_for_tables(args, per_species_dir, selected_species):
    pending_species = []
    for species in selected_species:
        flag_path = get_completion_flag_path(per_species_dir=per_species_dir, sp=species)
        if not os.path.exists(flag_path):
            pending_species.append(species)
            continue
        if args.redo:
            print('Output file detected. Will be overwritten: {}'.format(species), flush=True)
            species_dir = os.path.join(per_species_dir, species)
            if os.path.lexists(species_dir):
                if os.path.islink(species_dir):
                    os.remove(species_dir)
                elif os.path.isdir(species_dir):
                    shutil.rmtree(species_dir)
                else:
                    raise NotADirectoryError(
                        'Per-species output path exists but is not a directory: {}'.format(species_dir)
                    )
            pending_species.append(species)
            continue
        print('Skipping. Output file detected: {}'.format(species), flush=True)
    return pending_species


def run_per_species_jobs(args, metadata, input_dir, per_species_dir, pending_species, failed_species, species_jobs):
    if (species_jobs == 1) or (len(pending_species) <= 1):
        for species in pending_species:
            print('Starting: {}'.format(species), flush=True)
            exit_status = run_per_species_r_script(args, metadata, species, input_dir)
            register_exit_status(
                per_species_dir=per_species_dir,
                sp=species,
                exit_status=exit_status,
                failed_species=failed_species,
            )
        return

    max_workers = min(species_jobs, len(pending_species))
    print(
        'Running per-species table generation for {:,} species with {:,} parallel jobs.'.format(
            len(pending_species),
            max_workers,
        ),
        flush=True,
    )
    exit_status_by_species, failures = run_tasks_with_optional_threads(
        task_items=pending_species,
        task_fn=lambda species: run_per_species_r_script(args, metadata, species, input_dir),
        max_workers=max_workers,
    )
    for species, exc in failures:
        failed_species.add(species)
        sys.stderr.write('Per-species table generation failed for species: {} ({})\n'.format(species, exc))
    for species, exit_status in exit_status_by_species.items():
        register_exit_status(
            per_species_dir=per_species_dir,
            sp=species,
            exit_status=exit_status,
            failed_species=failed_species,
        )


def generate_per_species_tables(args, context=None):
    species_jobs, _ = resolve_worker_allocation(
        requested_workers=getattr(args, 'internal_jobs', 'auto'),
        requested_threads=getattr(args, 'threads', 'auto'),
        internal_cpu_budget=getattr(args, 'internal_cpu_budget', 'auto'),
        worker_option_name='internal_jobs',
        context='per-species table generation:',
        disable_workers=(getattr(args, 'batch', None) is not None),
    )
    runtime_args = clone_namespace(args, internal_jobs=species_jobs)
    check_rscript()
    out_dir = os.path.realpath(runtime_args.out_dir)
    if os.path.exists(out_dir) and (not os.path.isdir(out_dir)):
        raise NotADirectoryError('Output path exists but is not a directory: {}'.format(out_dir))
    per_species_dir = os.path.join(out_dir, 'per_species')
    if os.path.exists(per_species_dir) and (not os.path.isdir(per_species_dir)):
        raise NotADirectoryError('Per-species path exists but is not a directory: {}'.format(per_species_dir))
    runtime_args = clone_namespace(runtime_args, out_dir=out_dir)
    metadata, input_dir = resolve_per_species_execution_context(runtime_args, context=context)
    selected_species = list_selected_species(metadata)
    if ('tpm' in str(runtime_args.norm).lower()) and has_cstmm_counts_input(input_dir, selected_species):
        txt = ("TPM and TMM are incompatible. "
               "If input data are CSTMM-normalized, "
               "please switch --norm to any of the 'fpkm' normalization methods instead.")
        raise ValueError(txt)
    os.makedirs(per_species_dir, exist_ok=True)
    failed_species = set()
    print(
        'Number of species in the selected metadata table ("exclusion"=="no"): {}'.format(len(selected_species)),
        flush=True,
    )
    pending_species = collect_pending_species_for_tables(runtime_args, per_species_dir, selected_species)
    run_per_species_jobs(
        args=runtime_args,
        metadata=metadata,
        input_dir=input_dir,
        per_species_dir=per_species_dir,
        pending_species=pending_species,
        failed_species=failed_species,
        species_jobs=species_jobs,
    )
    if len(failed_species) > 0:
        txt = 'Per-species table generation failed for {}/{} species: {}\n'
        raise RuntimeError(txt.format(len(failed_species), len(selected_species), ', '.join(sorted(failed_species))).strip())
