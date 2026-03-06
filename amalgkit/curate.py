import datetime
import os
import re
import shutil
import subprocess
import sys
import warnings
from types import SimpleNamespace

from amalgkit.r_config import temporary_r_config
from amalgkit.util import *


def validate_curate_metadata_columns(metadata, required_columns, context):
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
        txt += '`amalgkit curate` recognizes samples with the same string in this column to belong the same group.'
        raise ValueError(txt.format(args.metadata))
    print('Tissues to be included: {}'.format(', '.join(sample_group)))
    sample_group = '|'.join(sample_group)
    return sample_group

def get_curate_completion_flag_path(curate_dir, sp):
    return os.path.join(curate_dir, sp, 'curate_completion_flag.txt')

def write_curate_completion_flag(curate_dir, sp):
    file_curate_completion_flag = get_curate_completion_flag_path(curate_dir, sp)
    with open(file_curate_completion_flag, 'w') as f:
        f.write('amalgkit curate completed at {}\n'.format(datetime.datetime.now()))

def register_curate_exit_status(curate_dir, sp, exit_status, failed_species):
    if exit_status == 0:
        write_curate_completion_flag(curate_dir=curate_dir, sp=sp)
    else:
        failed_species.add(sp)
        sys.stderr.write('amalgkit curate failed for species: {}\n'.format(sp))

def run_curate_r_script(args, metadata, sp, input_dir):
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
    requested_script = str(getattr(args, 'r_script_name', 'curate.r'))
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
        path_curate_input_metadata = os.path.realpath(metadata_arg)
    else:
        path_curate_input_metadata = os.path.join(input_dir, 'metadata.tsv')
    if not os.path.exists(path_curate_input_metadata):
        raise FileNotFoundError('Metadata file not found for curate R script: {}'.format(path_curate_input_metadata))
    if not os.path.isfile(path_curate_input_metadata):
        raise IsADirectoryError('Metadata path exists but is not a file: {}'.format(path_curate_input_metadata))
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
    print("Starting Rscript ({}) to obtain curated {} values.".format(r_script_name, args.norm), flush=True)
    config_map = {
        'est_counts_path': count_file,
        'metadata_path': path_curate_input_metadata,
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
    elif r_script_name == 'curate.r':
        pass
    else:
        raise ValueError('Unsupported R script for curate pipeline: {}'.format(r_script_name))
    with temporary_r_config(config_map, prefix='amalgkit_curate_r_') as config_path:
        curate_r_result = subprocess.run(['Rscript', r_script_path, config_path], check=False)
    return curate_r_result.returncode

def resolve_curate_input(args):
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
    validate_curate_metadata_columns(
        metadata=metadata,
        required_columns=['scientific_name', 'exclusion'],
        context='curate',
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


def collect_pending_species_for_curate(args, curate_dir, selected_species):
    pending_species = []
    for species in selected_species:
        flag_path = get_curate_completion_flag_path(curate_dir=curate_dir, sp=species)
        if not os.path.exists(flag_path):
            pending_species.append(species)
            continue
        if args.redo:
            print('Output file detected. Will be overwritten: {}'.format(species), flush=True)
            species_dir = os.path.join(curate_dir, species)
            if os.path.lexists(species_dir):
                if os.path.islink(species_dir):
                    os.remove(species_dir)
                elif os.path.isdir(species_dir):
                    shutil.rmtree(species_dir)
                else:
                    raise NotADirectoryError(
                        'Curate output path exists but is not a directory: {}'.format(species_dir)
                    )
            pending_species.append(species)
            continue
        print('Skipping. Output file detected: {}'.format(species), flush=True)
    return pending_species


def run_curate_species_jobs(args, metadata, input_dir, curate_dir, pending_species, failed_species, species_jobs):
    if (species_jobs == 1) or (len(pending_species) <= 1):
        for species in pending_species:
            print('Starting: {}'.format(species), flush=True)
            exit_status = run_curate_r_script(args, metadata, species, input_dir)
            register_curate_exit_status(
                curate_dir=curate_dir,
                sp=species,
                exit_status=exit_status,
                failed_species=failed_species,
            )
        return

    max_workers = min(species_jobs, len(pending_species))
    print('Running curate for {:,} species with {:,} parallel jobs.'.format(len(pending_species), max_workers), flush=True)
    exit_status_by_species, failures = run_tasks_with_optional_threads(
        task_items=pending_species,
        task_fn=lambda species: run_curate_r_script(args, metadata, species, input_dir),
        max_workers=max_workers,
    )
    for species, exc in failures:
        failed_species.add(species)
        sys.stderr.write('amalgkit curate failed for species: {} ({})\n'.format(species, exc))
    for species, exit_status in exit_status_by_species.items():
        register_curate_exit_status(
            curate_dir=curate_dir,
            sp=species,
            exit_status=exit_status,
            failed_species=failed_species,
        )


def curate_main(args):
    species_jobs, _ = resolve_worker_allocation(
        requested_workers=getattr(args, 'internal_jobs', 'auto'),
        requested_threads=getattr(args, 'threads', 'auto'),
        internal_cpu_budget=getattr(args, 'internal_cpu_budget', 'auto'),
        worker_option_name='internal_jobs',
        context='curate:',
        disable_workers=(getattr(args, 'batch', None) is not None),
    )
    args.internal_jobs = species_jobs
    check_rscript()
    out_dir = os.path.realpath(args.out_dir)
    if os.path.exists(out_dir) and (not os.path.isdir(out_dir)):
        raise NotADirectoryError('Output path exists but is not a directory: {}'.format(out_dir))
    curate_dir = os.path.join(out_dir, 'curate')
    if os.path.exists(curate_dir) and (not os.path.isdir(curate_dir)):
        raise NotADirectoryError('Curate path exists but is not a directory: {}'.format(curate_dir))
    metadata = getattr(args, '_resolved_metadata', None)
    input_dir = getattr(args, '_resolved_input_dir', None)
    if (metadata is None) or (input_dir is None):
        metadata, input_dir = resolve_curate_input(args)
    selected_species = list_selected_species(metadata)
    if ('tpm' in str(args.norm).lower()) and has_cstmm_counts_input(input_dir, selected_species):
        txt = ("TPM and TMM are incompatible. "
               "If input data are CSTMM-normalized, "
               "please switch --norm to any of the 'fpkm' normalization methods instead.")
        raise ValueError(txt)
    os.makedirs(curate_dir, exist_ok=True)
    failed_species = set()
    print(
        'Number of species in the selected metadata table ("exclusion"=="no"): {}'.format(len(selected_species)),
        flush=True,
    )
    pending_species = collect_pending_species_for_curate(args, curate_dir, selected_species)
    run_curate_species_jobs(
        args=args,
        metadata=metadata,
        input_dir=input_dir,
        curate_dir=curate_dir,
        pending_species=pending_species,
        failed_species=failed_species,
        species_jobs=species_jobs,
    )
    if len(failed_species) > 0:
        txt = 'amalgkit curate failed for {}/{} species: {}\n'
        raise RuntimeError(txt.format(len(failed_species), len(selected_species), ', '.join(sorted(failed_species))).strip())
