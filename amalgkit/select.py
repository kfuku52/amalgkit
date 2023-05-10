import shutil

from amalgkit.util import *

def write_select_outputs(path_metadata_original, path_metadata_table, metadata_dir, metadata):
    if os.path.exists(path_metadata_original):
        print('Original metadata copy exists and will not be renewed: {}'.format(path_metadata_original), flush=True)
    else:
        print('Original metadata copy does not exist. Creating: {}'.format(path_metadata_original), flush=True)
        shutil.copyfile(path_metadata_table, path_metadata_original)
    print('Updating metadata table at: {}'.format(path_metadata_table), flush=True)
    metadata.df.to_csv(path_metadata_table, sep='\t', index=False)
    sra_qualified_pivot = metadata.pivot(n_sp_cutoff=0, qualified_only=True, sampled_only=False)
    sra_qualified_pivot.to_csv(os.path.join(metadata_dir, 'pivot_qualified.tsv'), sep='\t')
    sra_selected_pivot = metadata.pivot(n_sp_cutoff=0, qualified_only=True, sampled_only=True)
    sra_selected_pivot.to_csv(os.path.join(metadata_dir, 'pivot_selected.tsv'), sep='\t')

def select_main(args):
    metadata_dir = os.path.join(args.out_dir, 'metadata')
    if args.config_dir=='inferred':
        dir_config = os.path.join(args.out_dir, 'config')
    else:
        dir_config = args.config_dir
    check_config_dir(dir_path=dir_config, mode='select')
    path_metadata_table = os.path.join(metadata_dir, 'metadata.tsv')
    path_metadata_original = os.path.join(metadata_dir, 'metadata_original.tsv')

    metadata = load_metadata(args)
    metadata.nspot_cutoff(args.min_nspots)
    if args.mark_redundant_biosamples:
        metadata.mark_redundant_biosample()
    metadata.remove_specialchars()
    metadata.group_attributes()
    metadata.mark_exclude_keywords()
    metadata.mark_treatment_terms()
    metadata.label_sampled_data(args.max_sample)
    metadata.reorder(omit_misc=True)
    write_select_outputs(path_metadata_original, path_metadata_table, metadata_dir, metadata)
