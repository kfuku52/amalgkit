import shutil
import re

from amalgkit.util import *

def resolve_select_config_dir(args):
    if args.config_dir == 'inferred':
        return os.path.join(args.out_dir, 'config')
    return args.config_dir


def filter_metadata_by_sample_group(metadata, sample_group_arg):
    if sample_group_arg is None:
        return metadata
    txt = '{}: Extracting pre-selected sample_group entries: {}'
    print(txt.format(datetime.datetime.now(), sample_group_arg), flush=True)
    if 'sample_group' not in metadata.df.columns:
        raise ValueError(
            'Column "sample_group" is required in metadata when --sample_group is specified.'
        )
    selected_sample_groups = [token.strip() for token in re.split(r'[,\|]+', str(sample_group_arg))]
    selected_sample_groups = [token for token in selected_sample_groups if token != '']
    if len(selected_sample_groups) == 0:
        raise ValueError('No sample_group was selected. Provide non-empty --sample_group value.')
    sample_groups = metadata.df['sample_group'].fillna('').astype(str).str.strip()
    metadata.df = metadata.df.loc[sample_groups.isin(selected_sample_groups), :].reset_index(drop=True)
    return metadata


def apply_select_filters(metadata, args, dir_config):
    metadata.nspot_cutoff(args.min_nspots)
    metadata.mark_missing_rank(args.mark_missing_rank)
    metadata.mark_redundant_biosample(args.mark_redundant_biosamples)
    metadata.remove_specialchars()
    metadata.group_attributes(dir_config)
    metadata.mark_exclude_keywords(dir_config)
    metadata.mark_treatment_terms(dir_config)
    metadata.label_sampled_data(args.max_sample)
    metadata.reorder(omit_misc=True)
    return metadata


def write_select_outputs(path_metadata_original, path_metadata_table, metadata_dir, metadata, metadata_original_df=None):
    if os.path.exists(metadata_dir) and (not os.path.isdir(metadata_dir)):
        raise NotADirectoryError('Output metadata path exists but is not a directory: {}'.format(metadata_dir))
    for output_path in [path_metadata_original, path_metadata_table]:
        if os.path.exists(output_path) and (not os.path.isfile(output_path)):
            raise NotADirectoryError('Output path exists but is not a file: {}'.format(output_path))
    os.makedirs(metadata_dir, exist_ok=True)
    if os.path.exists(path_metadata_original):
        print('Original metadata copy exists and will not be renewed: {}'.format(path_metadata_original), flush=True)
    else:
        print('Original metadata copy does not exist. Creating: {}'.format(path_metadata_original), flush=True)
        if os.path.exists(path_metadata_table):
            shutil.copyfile(path_metadata_table, path_metadata_original)
        else:
            if metadata_original_df is None:
                metadata_original_df = metadata.df
            metadata_original_df.to_csv(path_metadata_original, sep='\t', index=False)
    print('Updating metadata table at: {}'.format(path_metadata_table), flush=True)
    metadata.df.to_csv(path_metadata_table, sep='\t', index=False)
    pivot_specs = [
        ('pivot_qualified.tsv', False),
        ('pivot_selected.tsv', True),
    ]
    for filename, sampled_only in pivot_specs:
        pivot_df = metadata.pivot(n_sp_cutoff=0, qualified_only=True, sampled_only=sampled_only)
        pivot_df.to_csv(os.path.join(metadata_dir, filename), sep='\t')

def select_main(args):
    out_dir = os.path.realpath(args.out_dir)
    if os.path.exists(out_dir) and (not os.path.isdir(out_dir)):
        raise NotADirectoryError('Output path exists but is not a directory: {}'.format(out_dir))
    args.out_dir = out_dir
    metadata_dir = os.path.join(out_dir, 'metadata')
    dir_config = resolve_select_config_dir(args)
    check_config_dir(dir_path=dir_config, mode='select')
    path_metadata_table = os.path.join(metadata_dir, 'metadata.tsv')
    path_metadata_original = os.path.join(metadata_dir, 'metadata_original.tsv')

    metadata = load_metadata(args)
    metadata_original_df = metadata.df.copy(deep=True)
    metadata = filter_metadata_by_sample_group(metadata, args.sample_group)
    metadata = apply_select_filters(metadata, args, dir_config)
    write_select_outputs(
        path_metadata_original=path_metadata_original,
        path_metadata_table=path_metadata_table,
        metadata_dir=metadata_dir,
        metadata=metadata,
        metadata_original_df=metadata_original_df,
    )
