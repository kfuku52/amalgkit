import glob
import os
import re
import shlex
import shutil
import subprocess
import sys

import pandas

from amalgkit.util import load_metadata


REQUIRED_COLUMNS = [
    'busco_id',
    'status',
    'sequence',
    'score',
    'length',
    'orthodb_url',
    'description',
]


def normalize_busco_columns(df):
    def keyify(name):
        return re.sub(r'[^a-z0-9]', '', str(name).lower())

    lookup = {keyify(col): col for col in df.columns}
    required_keys = {
        'buscoid': 'busco_id',
        'status': 'status',
        'sequence': 'sequence',
        'score': 'score',
        'length': 'length',
        'orthodburl': 'orthodb_url',
        'description': 'description',
    }
    missing = [target for key, target in required_keys.items() if key not in lookup]
    if missing:
        raise ValueError('Missing required BUSCO columns: {}'.format(', '.join(missing)))
    renamed = {lookup[key]: target for key, target in required_keys.items()}
    df = df.rename(columns=renamed)
    return df


def normalize_busco_table(src_path, dest_path):
    header = None
    with open(src_path, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                continue
            candidate = line.lstrip('#').strip()
            if candidate.lower().startswith('busco'):
                header = candidate.split('\t')
    if header:
        df = pandas.read_table(src_path, sep='\t', header=None, comment='#', names=header, dtype=str, low_memory=False)
    else:
        df = pandas.read_table(src_path, sep='\t', header=0, comment='#', dtype=str, low_memory=False)
    if df.shape[0] == 0:
        raise ValueError('BUSCO table is empty: {}'.format(src_path))
    if not header and df.shape[1] == len(REQUIRED_COLUMNS):
        df.columns = REQUIRED_COLUMNS
    df = normalize_busco_columns(df)
    df = df.loc[:, REQUIRED_COLUMNS]
    with open(dest_path, 'w') as f:
        f.write('# Busco id\tStatus\tSequence\tScore\tLength\tOrthoDB url\tDescription\n')
        df.to_csv(f, sep='\t', index=False, header=False)


def find_full_table(output_dir):
    patterns = [
        '**/full_table*.tsv',
        '**/full_table*.tsv.gz',
        '**/*full_table*.tsv',
        '**/*full_table*.tsv.gz',
    ]
    matches = []
    for pattern in patterns:
        matches.extend(glob.glob(os.path.join(output_dir, pattern), recursive=True))
    matches = sorted(set(matches))
    if not matches:
        raise FileNotFoundError('No full_table.tsv found under: {}'.format(output_dir))
    if len(matches) > 1:
        raise ValueError('Multiple BUSCO full_table files detected: {}'.format(', '.join(matches)))
    return matches[0]


def resolve_fasta_dir(args):
    if args.fasta_dir == 'inferred':
        return os.path.join(args.out_dir, 'fasta')
    return os.path.realpath(args.fasta_dir)


def resolve_species_fasta(sci_name, fasta_dir):
    sci_name = sci_name.replace(' ', '_')
    fasta_files = []
    for ext in ['*.fa', '*.fasta', '*.fa.gz', '*.fasta.gz']:
        fasta_files.extend(glob.glob(os.path.join(fasta_dir, sci_name + ext)))
    if len(fasta_files) > 1:
        raise ValueError('Found multiple reference fasta files for {}: {}'.format(sci_name, ', '.join(fasta_files)))
    if len(fasta_files) == 0:
        raise FileNotFoundError('Could not find reference fasta file for {} in: {}'.format(sci_name, fasta_dir))
    return fasta_files[0]


def select_tool(args):
    tool = args.tool
    if tool == 'auto':
        if shutil.which(args.compleasm_exe):
            return 'compleasm'
        if shutil.which(args.busco_exe):
            return 'busco'
        raise FileNotFoundError('Neither compleasm nor busco was found on PATH.')
    if tool == 'compleasm':
        if not shutil.which(args.compleasm_exe):
            raise FileNotFoundError('compleasm executable not found: {}'.format(args.compleasm_exe))
        return 'compleasm'
    if tool == 'busco':
        if not shutil.which(args.busco_exe):
            raise FileNotFoundError('busco executable not found: {}'.format(args.busco_exe))
        return 'busco'
    raise ValueError('Unknown tool: {}'.format(tool))


def ensure_clean_dir(path, redo):
    if os.path.exists(path) and redo:
        shutil.rmtree(path)
    os.makedirs(path, exist_ok=True)


def run_command(cmd):
    print('Command: {}'.format(' '.join(cmd)), flush=True)
    out = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    print(out.stdout.decode('utf8'))
    print(out.stderr.decode('utf8'))
    if out.returncode != 0:
        raise RuntimeError('Command failed with exit code {}: {}'.format(out.returncode, ' '.join(cmd)))


def run_busco(fasta_path, sci_name, output_root, args, extra_args):
    out_name = sci_name.replace(' ', '_')
    cmd = [
        args.busco_exe,
        '-i', fasta_path,
        '-o', out_name,
        '-l', args.lineage,
        '-m', 'transcriptome',
        '--out_path', output_root,
        '--cpu', str(args.threads),
    ]
    if args.redo:
        cmd.append('--force')
    cmd.extend(extra_args)
    run_command(cmd)
    return os.path.join(output_root, out_name)


def run_compleasm(fasta_path, sci_name, output_root, args, extra_args):
    out_dir = os.path.join(output_root, sci_name.replace(' ', '_'))
    ensure_clean_dir(out_dir, args.redo)
    cmd = [
        args.compleasm_exe,
        'run',
        '-i', fasta_path,
        '-o', out_dir,
        '-l', args.lineage,
        '-t', str(args.threads),
        '--mode', 'transcriptome',
    ]
    cmd.extend(extra_args)
    run_command(cmd)
    return out_dir


def collect_species(args, metadata):
    if args.fasta is not None:
        if args.species is None:
            raise ValueError('--species is required when --fasta is provided.')
        return [args.species], {args.species: os.path.realpath(args.fasta)}
    fasta_dir = resolve_fasta_dir(args)
    if not os.path.exists(fasta_dir):
        raise FileNotFoundError('FASTA directory not found: {}'.format(fasta_dir))
    species = metadata.df.loc[:, 'scientific_name'].dropna().unique().tolist()
    fasta_map = {}
    for sp in species:
        fasta_map[sp] = resolve_species_fasta(sp, fasta_dir)
    return species, fasta_map


def busco_main(args):
    if not args.lineage:
        raise ValueError('--lineage is required.')
    tool = select_tool(args)
    extra_args = shlex.split(args.tool_args) if args.tool_args else []
    if args.fasta is not None:
        metadata = None
    else:
        metadata = load_metadata(args)
    species, fasta_map = collect_species(args, metadata)
    busco_dir = os.path.join(os.path.realpath(args.out_dir), 'busco')
    os.makedirs(busco_dir, exist_ok=True)
    for sp in species:
        print('Processing species: {}'.format(sp), flush=True)
        fasta_path = fasta_map[sp]
        if tool == 'busco':
            output_root = busco_dir
            ensure_clean_dir(os.path.join(output_root, sp.replace(' ', '_')), args.redo)
            tool_out_dir = run_busco(fasta_path, sp, output_root, args, extra_args)
        else:
            tool_out_dir = run_compleasm(fasta_path, sp, busco_dir, args, extra_args)
        full_table = find_full_table(tool_out_dir)
        out_table = os.path.join(busco_dir, sp.replace(' ', '_') + '_busco.tsv')
        normalize_busco_table(full_table, out_table)
        print('BUSCO table written: {}'.format(out_table), flush=True)
