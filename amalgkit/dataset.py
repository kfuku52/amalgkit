"""
amalgkit dataset: Extract bundled datasets and bundled select rule sets.

This module provides functionality to extract pre-bundled datasets for
testing and demonstration purposes, and to export packaged
``select_rules.tsv`` files for use with ``amalgkit select``.
"""

import os
import shutil

from amalgkit.exceptions import AmalgkitExit

try:
    import importlib.resources as ir
except ImportError:
    import importlib_resources as ir

SELECT_RULES_FILENAME = 'select_rules.tsv'

# Available datasets and their descriptions
DATASETS = {
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


def build_extracted_dirs(out_dir):
    out_dir = os.path.realpath(out_dir)
    if os.path.exists(out_dir) and (not os.path.isdir(out_dir)):
        raise NotADirectoryError('Output path exists but is not a directory: {}'.format(out_dir))
    extracted_dirs = {
        'fasta': os.path.join(out_dir, 'fasta'),
        'busco': os.path.join(out_dir, 'busco'),
        'rules': os.path.join(out_dir, 'rules'),
    }
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


def copy_dataset_files(dataset, dataset_src, extracted_dirs, overwrite=False):
    for file_type, files in dataset['files'].items():
        dest_dir = extracted_dirs.get(file_type)
        if dest_dir is None:
            continue
        for filename in files:
            src = os.path.join(dataset_src, filename)
            dst = os.path.join(dest_dir, filename)
            if os.path.exists(dst) and (not os.path.isfile(dst)):
                raise NotADirectoryError(
                    'Destination path exists but is not a file: {}'.format(dst)
                )
            if os.path.exists(dst) and not overwrite:
                print(f'  Skipping (exists): {dst}')
                continue
            shutil.copy2(src, dst)
            print(f'  Copied: {filename} -> {dest_dir}/')


def export_rule_set(rule_set_name, out_dir, overwrite=False):
    validate_rule_set_name(rule_set_name)
    path_rules = build_select_rules_output_path(out_dir)
    if os.path.exists(path_rules):
        print('Output select rules already exists: {}'.format(path_rules))
        if overwrite:
            print('--overwrite is set to "yes". The select rules file will be overwritten: {}'.format(path_rules))
        else:
            raise AmalgkitExit('--overwrite is set to "no". Exiting.', exit_code=0, use_stderr=False)
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


def extract_dataset(name, out_dir, overwrite=False):
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
        Paths to extracted directories: {'fasta': ..., 'busco': ..., 'rules': ...}
    """
    validate_dataset_name(name)
    dataset = DATASETS[name]
    dataset_src = resolve_dataset_source_dir(name)
    validate_dataset_source_files(dataset, dataset_src)
    extracted_dirs = build_extracted_dirs(out_dir)
    copy_dataset_files(dataset, dataset_src, extracted_dirs, overwrite=overwrite)
    return extracted_dirs


def dataset_main(args):
    """Main entry point for the dataset command."""

    if args.list:
        list_dataset_assets()
        return

    if args.name is None and args.rule_set is None:
        raise ValueError('Please specify --name, --rule_set, or use --list to see available assets.')

    dataset_paths = None
    rules_path = None

    if args.name is not None:
        print(f'Extracting dataset "{args.name}" to: {args.out_dir}')
        dataset_paths = extract_dataset(
            name=args.name,
            out_dir=args.out_dir,
            overwrite=args.overwrite,
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
        print(f'  Dataset rules:    {dataset_paths["rules"]}')
        if args.name == 'yeast':
            print('  Note: the yeast dataset uses small BUSCO-focused test FASTAs, so BUSCO completeness is intentionally lower than a full gene set.')
    if rules_path is not None:
        print(f'  Exported rules:   {rules_path}')
    print('')
    print('Example usage:')
    if dataset_paths is not None:
        print(f'  amalgkit dataset --out_dir {args.out_dir} --rule_set base --overwrite yes')
    else:
        print(f'  amalgkit select --out_dir {args.out_dir} --select_rules_tsv {rules_path}')
