import os
import re
import warnings

import numpy
import pandas

from amalgkit.parallel_utils import run_tasks_with_optional_threads as default_run_tasks_with_optional_threads


def orthogroup2genecount(file_orthogroup, file_genecount, spp):
    df = pandas.read_csv(file_orthogroup, sep='\t', header=0, low_memory=False)
    if 'busco_id' not in df.columns:
        raise ValueError('Column "busco_id" is required in orthogroup table: {}'.format(file_orthogroup))
    missing_species = [sp for sp in spp if sp not in df.columns]
    if len(missing_species) > 0:
        raise ValueError(
            'Species column(s) not found in orthogroup table ({}): {}'.format(
                file_orthogroup,
                ', '.join(missing_species),
            )
        )
    orthogroup_df = pandas.DataFrame({'orthogroup_id': df['busco_id'].to_numpy()})
    is_spp = df.columns.isin(spp)
    df = df.loc[:, is_spp].fillna('').replace('-', '').astype(str)
    values = df.to_numpy(dtype=str)
    non_empty = (values != '')
    comma_counts = numpy.char.count(values, ',').astype(int)
    gc_values = non_empty.astype(int) + (comma_counts * non_empty)
    gc = pandas.DataFrame(gc_values, index=df.index, columns=df.columns)
    gc = pandas.concat([orthogroup_df, gc], axis=1)
    col_order = ['orthogroup_id'] + [col for col in gc.columns if col != 'orthogroup_id']
    gc = gc[col_order]
    gc.to_csv(file_genecount, index=False, sep='\t')


def check_ortholog_parameter_compatibility(args):
    orthogroup_table = getattr(args, 'orthogroup_table', None)
    dir_busco = getattr(args, 'dir_busco', None)
    if isinstance(orthogroup_table, str):
        orthogroup_table = orthogroup_table.strip()
        if orthogroup_table == '':
            orthogroup_table = None
    if isinstance(dir_busco, str):
        dir_busco = dir_busco.strip()
        if dir_busco == '':
            dir_busco = None
    if (orthogroup_table is None) and (dir_busco is None):
        raise ValueError('One of --orthogroup_table and --dir_busco should be specified.')
    if (orthogroup_table is not None) and (dir_busco is not None):
        raise ValueError('Only one of --orthogroup_table and --dir_busco should be specified.')
    return orthogroup_table, dir_busco


BUSCO_TABLE_COLUMNS = ['busco_id', 'status', 'sequence', 'score', 'length', 'orthodb_url', 'description']
BUSCO_TABLE_USE_COLUMNS = ['busco_id', 'sequence', 'orthodb_url', 'description']
BUSCO_SPECIES_SUFFIX_PATTERN = re.compile(r'\.tsv(?:\.gz)?$', re.IGNORECASE)
BUSCO_SPECIES_CLEANUP_PATTERN = re.compile(r'(_busco|_full_table.*)$', re.IGNORECASE)
BUSCO_SPECIES_MATCH_PATTERN = re.compile(r'^([^_]+_[^_]+)')


def parse_busco_species_name(species_infile):
    species_colname = BUSCO_SPECIES_SUFFIX_PATTERN.sub('', species_infile)
    species_colname = BUSCO_SPECIES_CLEANUP_PATTERN.sub('', species_colname)
    matched = BUSCO_SPECIES_MATCH_PATTERN.match(species_colname)
    if matched is not None:
        species_colname = matched.group(1)
    return species_colname


def read_busco_species_table(path_to_table):
    tmp_table = pandas.read_table(
        path_to_table,
        sep='\t',
        header=None,
        comment='#',
        names=BUSCO_TABLE_COLUMNS,
        usecols=BUSCO_TABLE_USE_COLUMNS,
    )
    busco_id_key = (
        tmp_table.loc[:, 'busco_id']
        .fillna('')
        .astype(str)
        .str.lower()
        .str.replace(r'[^a-z0-9]', '', regex=True)
    )
    tmp_table = tmp_table.loc[busco_id_key != 'buscoid', :].copy()
    tmp_table.loc[:, 'sequence'] = tmp_table.loc[:, 'sequence'].str.replace(r':[-\.0-9]*$', '', regex=True)
    for col in ['sequence', 'orthodb_url', 'description']:
        tmp_table[col] = tmp_table[col].fillna('').astype(str)
        tmp_table.loc[(tmp_table[col] == ''), col] = '-'
    return tmp_table


def parse_busco_species_table(dir_busco, species_infile):
    path_to_table = os.path.join(dir_busco, species_infile)
    if not os.path.isfile(path_to_table):
        if os.path.exists(path_to_table):
            warnings.warn('full_table.tsv path exists but is not a file. Skipping: {}'.format(species_infile))
        else:
            warnings.warn('full_table.tsv does not exist. Skipping: {}'.format(species_infile))
        return None
    tmp_table = read_busco_species_table(path_to_table)
    meta_rows = tmp_table.loc[:, ['busco_id', 'orthodb_url', 'description']].drop_duplicates(
        subset=['busco_id'],
        keep='first',
        inplace=False,
    )
    grouped = tmp_table.loc[:, ['busco_id', 'sequence']].groupby('busco_id', sort=False)['sequence'].agg(','.join)
    species_colname = parse_busco_species_name(species_infile)
    return species_colname, meta_rows, grouped


def append_unique_busco_ids(busco_id_seen, busco_id_order, busco_ids):
    for busco_id in busco_ids:
        if busco_id in busco_id_seen:
            continue
        busco_id_seen.add(busco_id)
        busco_id_order.append(busco_id)


def update_busco_meta(busco_meta, busco_ids, urls, descriptions):
    for busco_id, orthodb_url, description in zip(busco_ids, urls, descriptions):
        if busco_id not in busco_meta:
            busco_meta[busco_id] = {
                'orthodb_url': orthodb_url,
                'description': description,
            }
            continue
        if (busco_meta[busco_id]['orthodb_url'] == '-') and (orthodb_url != '-'):
            busco_meta[busco_id]['orthodb_url'] = orthodb_url
        if (busco_meta[busco_id]['description'] == '-') and (description != '-'):
            busco_meta[busco_id]['description'] = description


def generate_multisp_busco_table(dir_busco, outfile, run_tasks_with_optional_threads=None):
    if run_tasks_with_optional_threads is None:
        run_tasks_with_optional_threads = default_run_tasks_with_optional_threads
    dir_busco = os.path.realpath(dir_busco)
    if not os.path.exists(dir_busco):
        raise FileNotFoundError('BUSCO directory not found: {}'.format(dir_busco))
    if not os.path.isdir(dir_busco):
        raise NotADirectoryError('BUSCO path exists but is not a directory: {}'.format(dir_busco))
    print('Generating multi-species BUSCO table.', flush=True)
    species_infiles = [
        f for f in os.listdir(path=dir_busco)
        if BUSCO_SPECIES_SUFFIX_PATTERN.search(f)
        and os.path.isfile(os.path.join(dir_busco, f))
    ]
    species_infiles = sorted(species_infiles)
    print('BUSCO full tables for {} species were detected at: {}'.format(len(species_infiles), dir_busco), flush=True)
    if len(species_infiles) == 0:
        raise FileNotFoundError('No BUSCO full table file (.tsv) was detected in: {}'.format(dir_busco))

    busco_id_seen = set()
    busco_id_order = []
    busco_meta = dict()
    species_series = dict()
    species_order = []

    max_workers = min(8, len(species_infiles))
    parsed_by_file, parse_failures = run_tasks_with_optional_threads(
        task_items=species_infiles,
        task_fn=lambda species_infile: parse_busco_species_table(dir_busco, species_infile),
        max_workers=max_workers,
    )
    for species_infile, exc in parse_failures:
        warnings.warn('Failed to parse BUSCO table {}: {}'.format(species_infile, exc))
        parsed_by_file[species_infile] = None
    parsed_results = [parsed_by_file.get(species_infile) for species_infile in species_infiles]

    species_source_file = dict()
    for species_infile, parsed in zip(species_infiles, parsed_results):
        if parsed is None:
            continue
        species_colname, meta_rows, grouped = parsed
        existing_infile = species_source_file.get(species_colname)
        if (existing_infile is not None) and (existing_infile != species_infile):
            raise ValueError(
                'Duplicate species label was detected across BUSCO tables: {} ({} vs {}). '
                'Rename input BUSCO table files to avoid species label collision.'.format(
                    species_colname,
                    existing_infile,
                    species_infile,
                )
            )
        species_source_file[species_colname] = species_infile
        busco_ids = meta_rows['busco_id'].to_numpy()
        urls = meta_rows['orthodb_url'].to_numpy()
        descriptions = meta_rows['description'].to_numpy()
        append_unique_busco_ids(busco_id_seen, busco_id_order, busco_ids)
        update_busco_meta(busco_meta, busco_ids, urls, descriptions)
        species_order.append(species_colname)
        species_series[species_colname] = grouped

    if len(species_series) == 0:
        raise ValueError('Failed to parse any BUSCO table under: {}'.format(dir_busco))

    merged_table = pandas.DataFrame({'busco_id': busco_id_order})
    merged_table['orthodb_url'] = merged_table['busco_id'].map(
        lambda bid: busco_meta.get(bid, {}).get('orthodb_url', '-')
    )
    merged_table['description'] = merged_table['busco_id'].map(
        lambda bid: busco_meta.get(bid, {}).get('description', '-')
    )
    for species_colname in species_order:
        merged_table[species_colname] = merged_table['busco_id'].map(species_series[species_colname])
    merged_table.to_csv(outfile, sep='\t', index=None, doublequote=False)
