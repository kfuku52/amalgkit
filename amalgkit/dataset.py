"""
amalgkit dataset: Extract bundled test datasets.

This module provides functionality to extract pre-bundled datasets
for testing and demonstration purposes.
"""

import os
import shutil
import sys

# Available datasets and their descriptions
DATASETS = {
    'yeast': {
        'description': 'Two yeast species (S. cerevisiae + S. pombe) with BUSCO genes only (~3 MB)',
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
            'config': [
                'control_term.config',
                'exclude_keyword.config',
                'group_attribute.config',
            ],
        },
    },
}


def get_dataset_dir():
    """Return the path to the bundled datasets directory."""
    return os.path.join(os.path.dirname(__file__), 'datasets')


def list_datasets():
    """List available datasets with descriptions."""
    print('Available datasets:')
    for name, info in DATASETS.items():
        print(f'  {name}: {info["description"]}')
        print(f'    Species: {", ".join(info["species"])}')


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
        Paths to extracted directories: {'fasta': ..., 'busco': ..., 'config': ...}
    """
    if name not in DATASETS:
        available = ', '.join(DATASETS.keys())
        sys.stderr.write(f'Error: Unknown dataset "{name}". Available: {available}\n')
        sys.exit(1)

    dataset = DATASETS[name]
    dataset_src = os.path.join(get_dataset_dir(), name)

    if not os.path.isdir(dataset_src):
        sys.stderr.write(f'Error: Dataset directory not found: {dataset_src}\n')
        sys.exit(1)

    # Create output directories
    out_dir = os.path.realpath(out_dir)
    fasta_dir = os.path.join(out_dir, 'fasta')
    busco_dir = os.path.join(out_dir, 'busco')
    config_dir = os.path.join(out_dir, 'config')

    for d in [fasta_dir, busco_dir, config_dir]:
        os.makedirs(d, exist_ok=True)

    # Copy files
    copied = {'fasta': [], 'busco': [], 'config': []}

    for file_type, files in dataset['files'].items():
        if file_type == 'fasta':
            dest_dir = fasta_dir
        elif file_type == 'busco':
            dest_dir = busco_dir
        elif file_type == 'config':
            dest_dir = config_dir
        else:
            continue

        for filename in files:
            src = os.path.join(dataset_src, filename)
            dst = os.path.join(dest_dir, filename)

            if os.path.exists(dst) and not overwrite:
                print(f'  Skipping (exists): {dst}')
                continue

            if not os.path.exists(src):
                sys.stderr.write(f'Warning: Source file not found: {src}\n')
                continue

            shutil.copy2(src, dst)
            copied[file_type].append(dst)
            print(f'  Copied: {filename} -> {dest_dir}/')

    return {
        'fasta': fasta_dir,
        'busco': busco_dir,
        'config': config_dir,
    }


def dataset_main(args):
    """Main entry point for the dataset command."""

    if args.list:
        list_datasets()
        return

    if args.name is None:
        sys.stderr.write('Error: Please specify --name or use --list to see available datasets.\n')
        sys.exit(1)

    print(f'Extracting dataset "{args.name}" to: {args.out_dir}')

    paths = extract_dataset(
        name=args.name,
        out_dir=args.out_dir,
        overwrite=args.overwrite,
    )

    print('')
    print('Dataset extracted successfully!')
    print(f'  FASTA files:  {paths["fasta"]}')
    print(f'  BUSCO files:  {paths["busco"]}')
    print(f'  Config files: {paths["config"]}')
    print('')
    print('Example usage with this dataset:')
    print(f'  amalgkit config --out_dir {args.out_dir} --config base --overwrite yes')
    print(f'  cp {paths["config"]}/*.config {args.out_dir}/config_base/')
