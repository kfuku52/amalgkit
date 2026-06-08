"""
amalgkit dataset: Extract bundled datasets, initialize empty workspaces,
and export packaged select rule sets.

This module provides functionality to extract pre-bundled datasets for
testing and demonstration purposes, initialize an empty workspace
scaffold, and export packaged ``select_rules.tsv`` files for use with
``amalgkit select``.
"""

import os
import shutil

from amalgkit.exceptions import AmalgkitExit

try:
    import importlib.resources as ir
except ImportError:
    import importlib_resources as ir

SELECT_RULES_FILENAME = 'select_rules.tsv'
SPECIES_TSV_FILENAME = 'species.tsv'
ORGAN_TERMS_TSV_FILENAME = 'organ_terms.tsv'
WORKSPACE_README_FILENAME = 'WORKSPACE_README.md'
DEFAULT_INIT_RULE_SET = 'base'
WORKSPACE_DIRS = {
    'fasta': 'fasta',
    'private_fastq': 'private_fastq',
    'downloads': 'downloads',
    'download_locks': os.path.join('downloads', 'locks'),
    'metadata': 'metadata',
    'metadata_specieswise': 'metadata_specieswise',
}
ROOT_LEVEL_DATASET_FILE_TYPES = {'rules'}
INIT_TEMPLATE_TEXT = {
    'species_tsv': 'scientific_name\tspecies_token\n',
    'organ_terms_tsv': 'sample_group\ttitle_terms\n',
}
WORKSPACE_README_TEXT = """# Amalgkit Workspace Scaffold

This workspace was initialized by `amalgkit dataset --name init`.

## Files
- `species.tsv`: input table for `amalgkit metadata --species_tsv`. Required column: `scientific_name`. Optional column: `species_token`.
- `organ_terms.tsv`: optional table for `amalgkit metadata --species_tsv --mode title_union` or `title_split`. Required columns: `sample_group`, `title_terms`. Separate multiple title terms with `;`.
- `select_rules.tsv`: default rule set loaded by `amalgkit select --out_dir ...`.

## Directories
- `fasta/`: reference FASTA files.
- `private_fastq/`: local FASTQ files for `amalgkit integrate`.
- `downloads/`: shared cache and download directory.
- `metadata/`: merged metadata outputs.
- `metadata_specieswise/`: per-species metadata outputs for batch workflows.

## Quick Start
1. Fill `species.tsv` if you want species-wise metadata retrieval.
2. Optionally fill `organ_terms.tsv` for title-based batch queries.
3. Review `select_rules.tsv`.
4. Run `amalgkit metadata --out_dir ./ --species_tsv ./species.tsv`.
5. Run `amalgkit select --out_dir ./`.
"""

# Available datasets and their descriptions
DATASETS = {
    'init': {
        'description': 'Initialize an empty workspace scaffold with template TSVs and a default select_rules.tsv',
        'species': [],
        'files': {},
    },
    'yeast': {
        'description': 'Two yeast species (S. cerevisiae + S. pombe) with small BUSCO-focused test FASTAs and matching BUSCO tables (~3 MB)',
        'species': ['Saccharomyces_cerevisiae', 'Schizosaccharomyces_pombe'],
        'files': {
            'fasta': [
                'Saccharomyces_cerevisiae.fa.gz',
                'Schizosaccharomyces_pombe.fa.gz',
            ],
            'busco': [
                'Saccharomyces_cerevisiae_busco.tsv',
                'Schizosaccharomyces_pombe_busco.tsv',
            ],
            'rules': [
                'select_rules.tsv',
            ],
        },
    },
}


def get_dataset_dir():
    """Return the path to the bundled datasets directory."""
    return os.path.join(os.path.dirname(__file__), 'datasets')


def get_rule_set_root():
    """Return the importlib.resources root for bundled select rule sets."""
    return ir.files('amalgkit.select_rule_sets')


def list_available_rule_sets():
    """Return available bundled select rule set names."""
    root = get_rule_set_root()
    rule_sets = []
    for entry in root.iterdir():
        if not entry.is_dir():
            continue
        if entry.name.startswith('_'):
            continue
        if entry.joinpath(SELECT_RULES_FILENAME).is_file():
            rule_sets.append(entry.name)
    return sorted(rule_sets)


def validate_rule_set_name(name):
    if name in list_available_rule_sets():
        return
    available = ', '.join(list_available_rule_sets())
    raise ValueError(f'Unknown rule set "{name}". Available: {available}')


def list_dataset_assets():
    """List available datasets and bundled select rule sets."""
    print('Available datasets:')
    for name, info in DATASETS.items():
        print(f'  {name}: {info["description"]}')
        if len(info.get('species', [])) > 0:
            print(f'    Species: {", ".join(info["species"])}')
    print('Available rule sets:')
    for name in list_available_rule_sets():
        print(f'  {name}')


def validate_dataset_name(name):
    if name in DATASETS:
        return
    available = ', '.join(DATASETS.keys())
    raise ValueError(f'Unknown dataset "{name}". Available: {available}')


def resolve_dataset_source_dir(name):
    dataset_src = os.path.join(get_dataset_dir(), name)
    if os.path.isdir(dataset_src):
        return dataset_src
    raise FileNotFoundError(f'Dataset directory not found: {dataset_src}')


def build_extracted_dirs(out_dir, extra_dir_names=None):
    out_dir = os.path.realpath(out_dir)
    if os.path.exists(out_dir) and (not os.path.isdir(out_dir)):
        raise NotADirectoryError('Output path exists but is not a directory: {}'.format(out_dir))
    extracted_dirs = {
        key: os.path.join(out_dir, relative_path)
        for key, relative_path in WORKSPACE_DIRS.items()
    }
    for dir_name in list(extra_dir_names or []):
        if dir_name in extracted_dirs:
            continue
        if dir_name in ROOT_LEVEL_DATASET_FILE_TYPES:
            continue
        extracted_dirs[dir_name] = os.path.join(out_dir, dir_name)
    for path_dir in extracted_dirs.values():
        os.makedirs(path_dir, exist_ok=True)
    return extracted_dirs


def validate_dataset_source_files(dataset, dataset_src):
    missing_sources = []
    for file_type, files in dataset['files'].items():
        for filename in files:
            src = os.path.join(dataset_src, filename)
            if not os.path.isfile(src):
                missing_sources.append(src)
    if missing_sources:
        raise FileNotFoundError(
            'Dataset source file(s) not found: {}'.format(', '.join(missing_sources))
        )


def build_select_rules_output_path(out_dir):
    out_dir = os.path.realpath(out_dir)
    if os.path.exists(out_dir) and (not os.path.isdir(out_dir)):
        raise NotADirectoryError('Output path exists but is not a directory: {}'.format(out_dir))
    os.makedirs(out_dir, exist_ok=True)
    path_rules = os.path.join(out_dir, SELECT_RULES_FILENAME)
    if os.path.exists(path_rules):
        if not os.path.isfile(path_rules):
            raise IsADirectoryError(
                'Select rules output path exists but is not a file: {}'.format(path_rules)
            )
    return path_rules


def preflight_rule_set_export(rule_set_name, out_dir, overwrite=False):
    validate_rule_set_name(rule_set_name)
    out_dir = os.path.realpath(out_dir)
    if os.path.exists(out_dir) and (not os.path.isdir(out_dir)):
        raise NotADirectoryError('Output path exists but is not a directory: {}'.format(out_dir))
    path_rules = os.path.join(out_dir, SELECT_RULES_FILENAME)
    if not os.path.lexists(path_rules):
        return path_rules
    if not os.path.isfile(path_rules):
        raise IsADirectoryError(
            'Select rules output path exists but is not a file: {}'.format(path_rules)
        )
    if overwrite:
        return path_rules
    print('Output select rules already exists: {}'.format(path_rules))
    raise AmalgkitExit('--overwrite is set to "no". Exiting.', exit_code=0, use_stderr=False)


def build_workspace_file_output_path(out_dir, filename, label):
    out_dir = os.path.realpath(out_dir)
    if os.path.exists(out_dir) and (not os.path.isdir(out_dir)):
        raise NotADirectoryError('Output path exists but is not a directory: {}'.format(out_dir))
    os.makedirs(out_dir, exist_ok=True)
    outpath = os.path.join(out_dir, filename)
    if os.path.exists(outpath) and (not os.path.isfile(outpath)):
        raise IsADirectoryError(
            '{} output path exists but is not a file: {}'.format(label, outpath)
        )
    return outpath


def write_workspace_template(out_dir, filename, text, label, overwrite=False):
    outpath = build_workspace_file_output_path(out_dir, filename, label)
    if os.path.exists(outpath) and (not overwrite):
        print('  Skipping (exists): {}'.format(outpath))
        return outpath
    if os.path.exists(outpath) and overwrite:
        print('  Overwriting: {}'.format(outpath))
    else:
        print('  Writing: {}'.format(outpath))
    with open(outpath, 'w', encoding='utf-8') as handle:
        handle.write(text)
    return outpath


def copy_dataset_files(dataset, dataset_src, extracted_dirs, out_dir, overwrite=False, skip_file_types=None):
    out_dir = os.path.realpath(out_dir)
    skip_file_types = set(skip_file_types or [])
    for file_type, files in dataset['files'].items():
        if file_type in skip_file_types:
            continue
        dest_dir = extracted_dirs.get(file_type)
        if dest_dir is None:
            if file_type not in ROOT_LEVEL_DATASET_FILE_TYPES:
                continue
        for filename in files:
            src = os.path.join(dataset_src, filename)
            if file_type in ROOT_LEVEL_DATASET_FILE_TYPES:
                dst = os.path.join(out_dir, filename)
                dest_label = out_dir
            else:
                dst = os.path.join(dest_dir, filename)
                dest_label = dest_dir
            if os.path.exists(dst) and (not os.path.isfile(dst)):
                raise NotADirectoryError(
                    'Destination path exists but is not a file: {}'.format(dst)
                )
            if os.path.exists(dst) and not overwrite:
                print(f'  Skipping (exists): {dst}')
                continue
            shutil.copy2(src, dst)
            print(f'  Copied: {filename} -> {dest_label}/')


def _copy_rule_set_to_output(rule_set_name, path_rules):
    validate_rule_set_name(rule_set_name)
    rule_root = get_rule_set_root().joinpath(rule_set_name)
    rules_file = rule_root.joinpath(SELECT_RULES_FILENAME)
    if not rules_file.is_file():
        raise FileNotFoundError(
            'Bundled select_rules.tsv not found for rule set "{}".'.format(rule_set_name)
        )
    print('Copying from {} to {}'.format(rules_file, path_rules))
    with open(path_rules, mode='wb') as f:
        f.write(rules_file.read_bytes())
    return path_rules


def export_rule_set(rule_set_name, out_dir, overwrite=False):
    path_rules = build_select_rules_output_path(out_dir)
    if os.path.exists(path_rules):
        print('Output select rules already exists: {}'.format(path_rules))
        if overwrite:
            print('--overwrite is set to "yes". The select rules file will be overwritten: {}'.format(path_rules))
        else:
            raise AmalgkitExit('--overwrite is set to "no". Exiting.', exit_code=0, use_stderr=False)
    return _copy_rule_set_to_output(rule_set_name, path_rules)


def ensure_rule_set(rule_set_name, out_dir, overwrite=False):
    path_rules = build_select_rules_output_path(out_dir)
    if os.path.exists(path_rules) and (not overwrite):
        print('  Skipping (exists): {}'.format(path_rules))
        return path_rules
    if os.path.exists(path_rules) and overwrite:
        print('--overwrite is set to "yes". The select rules file will be overwritten: {}'.format(path_rules))
    return _copy_rule_set_to_output(rule_set_name, path_rules)


def extract_dataset(name, out_dir, overwrite=False, skip_file_types=None):
    """
    Extract a bundled dataset to the specified output directory.

    Parameters
    ----------
    name : str
        Name of the dataset (e.g., 'yeast')
    out_dir : str
        Output directory path
    overwrite : bool
        Whether to overwrite existing files

    Returns
    -------
    dict
        Paths to extracted directories and root-level files.
    """
    validate_dataset_name(name)
    if name == 'init':
        raise ValueError('extract_dataset("init") is not supported. Use initialize_workspace().')
    dataset = DATASETS[name]
    dataset_src = resolve_dataset_source_dir(name)
    validate_dataset_source_files(dataset, dataset_src)
    extracted_dirs = build_extracted_dirs(out_dir, extra_dir_names=dataset['files'].keys())
    copy_dataset_files(
        dataset,
        dataset_src,
        extracted_dirs,
        out_dir=out_dir,
        overwrite=overwrite,
        skip_file_types=skip_file_types,
    )
    extracted_paths = dict(extracted_dirs)
    if ('rules' in dataset['files']) and ('rules' not in set(skip_file_types or [])):
        extracted_paths['select_rules_tsv'] = os.path.join(
            os.path.realpath(out_dir),
            SELECT_RULES_FILENAME,
        )
    return extracted_paths


def initialize_workspace(out_dir, overwrite=False, default_rule_set_name=DEFAULT_INIT_RULE_SET):
    workspace_dirs = build_extracted_dirs(out_dir)
    workspace_readme_path = write_workspace_template(
        out_dir=out_dir,
        filename=WORKSPACE_README_FILENAME,
        text=WORKSPACE_README_TEXT,
        label=WORKSPACE_README_FILENAME,
        overwrite=overwrite,
    )
    species_tsv_path = write_workspace_template(
        out_dir=out_dir,
        filename=SPECIES_TSV_FILENAME,
        text=INIT_TEMPLATE_TEXT['species_tsv'],
        label='species.tsv',
        overwrite=overwrite,
    )
    organ_terms_tsv_path = write_workspace_template(
        out_dir=out_dir,
        filename=ORGAN_TERMS_TSV_FILENAME,
        text=INIT_TEMPLATE_TEXT['organ_terms_tsv'],
        label='organ_terms.tsv',
        overwrite=overwrite,
    )
    select_rules_path = ''
    if default_rule_set_name not in [None, '']:
        select_rules_path = ensure_rule_set(
            rule_set_name=default_rule_set_name,
            out_dir=out_dir,
            overwrite=overwrite,
        )
    return {
        'workspace_dirs': workspace_dirs,
        'workspace_readme': workspace_readme_path,
        'species_tsv': species_tsv_path,
        'organ_terms_tsv': organ_terms_tsv_path,
        'select_rules_tsv': select_rules_path,
    }


def dataset_main(args):
    """Main entry point for the dataset command."""

    if args.list:
        list_dataset_assets()
        return

    if args.name is None and args.rule_set is None:
        raise ValueError('Please specify --name, --rule_set, or use --list to see available assets.')

    if args.name is not None:
        validate_dataset_name(args.name)
        if args.name != 'init':
            dataset_src = resolve_dataset_source_dir(args.name)
            validate_dataset_source_files(DATASETS[args.name], dataset_src)
    if args.rule_set is not None:
        preflight_rule_set_export(args.rule_set, args.out_dir, overwrite=args.overwrite)

    dataset_paths = None
    init_paths = None
    rules_path = None

    if args.name is not None:
        print(f'Extracting dataset "{args.name}" to: {args.out_dir}')
        if args.name == 'init':
            default_rule_set_name = None if args.rule_set is not None else DEFAULT_INIT_RULE_SET
            init_paths = initialize_workspace(
                out_dir=args.out_dir,
                overwrite=args.overwrite,
                default_rule_set_name=default_rule_set_name,
            )
        else:
            dataset_paths = extract_dataset(
                name=args.name,
                out_dir=args.out_dir,
                overwrite=args.overwrite,
                skip_file_types={'rules'} if args.rule_set is not None else None,
            )

    if args.rule_set is not None:
        print(f'Exporting rule set "{args.rule_set}" to: {args.out_dir}')
        rules_path = export_rule_set(
            rule_set_name=args.rule_set,
            out_dir=args.out_dir,
            overwrite=args.overwrite,
        )

    print('')
    print('Completed successfully!')
    if dataset_paths is not None:
        print(f'  FASTA files:      {dataset_paths["fasta"]}')
        print(f'  BUSCO files:      {dataset_paths["busco"]}')
        if 'select_rules_tsv' in dataset_paths:
            print(f'  Dataset rules:    {dataset_paths["select_rules_tsv"]}')
        if args.name == 'yeast':
            print('  Note: the yeast dataset uses small BUSCO-focused test FASTAs, so BUSCO completeness is intentionally lower than a full gene set.')
    if init_paths is not None:
        print(f'  FASTA dir:        {init_paths["workspace_dirs"]["fasta"]}')
        print(f'  private_fastq:    {init_paths["workspace_dirs"]["private_fastq"]}')
        print(f'  downloads:        {init_paths["workspace_dirs"]["downloads"]}')
        print(f'  metadata dir:     {init_paths["workspace_dirs"]["metadata"]}')
        print(f'  Workspace README: {init_paths["workspace_readme"]}')
        print(f'  species TSV:      {init_paths["species_tsv"]}')
        print(f'  organ terms TSV:  {init_paths["organ_terms_tsv"]}')
        if init_paths['select_rules_tsv'] != '':
            print(f'  Select rules:     {init_paths["select_rules_tsv"]}')
    if rules_path is not None:
        print(f'  Exported rules:   {rules_path}')
    print('')
    print('Example usage:')
    if init_paths is not None:
        print(f'  amalgkit metadata --out_dir {args.out_dir} --species_tsv {init_paths["species_tsv"]}')
        print(f'  amalgkit select --out_dir {args.out_dir}')
    elif dataset_paths is not None:
        print(f'  amalgkit dataset --out_dir {args.out_dir} --rule_set base --overwrite yes')
    else:
        print(f'  amalgkit select --out_dir {args.out_dir} --select_rules_tsv {rules_path}')
