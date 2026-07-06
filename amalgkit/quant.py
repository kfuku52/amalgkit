import json
import os
import re
import shlex
import subprocess
import sys
import tempfile

import numpy
import pandas

from amalgkit.arg_utils import clone_namespace
from amalgkit.command_context import PrefetchedDirEntries, QuantRuntimeContext
from amalgkit.download_utils import acquire_exclusive_lock
from amalgkit.metadata_utils import (
    get_metadata_row_index_by_run,
    get_newest_intermediate_file_extension,
    get_sra_stat,
    is_private_fastq_scientific_name_placeholder,
    load_metadata,
)
from amalgkit.parallel_utils import (
    is_auto_parallel_option,
    resolve_detected_cpu_count,
    resolve_thread_worker_allocation,
    run_tasks_with_optional_threads,
)
from amalgkit.prefix_utils import find_run_prefixed_entries, find_species_prefixed_entries
from amalgkit.subprocess_utils import probe_dependency_command, run_logged_command

INDEX_BUILD_LOCK_POLL_SECONDS = 5
INDEX_BUILD_LOCK_TIMEOUT_SECONDS = 3600

INDEX_FASTA_SUFFIXES = ('.fa', '.fasta', '.fa.gz', '.fasta.gz')
KALLISTO_INDEX_SUFFIX = '.idx'
OARFISH_INDEX_SUFFIX = '.mmi'
FASTA_REFERENCE_STEM_SUFFIXES = ('_for_kallisto_index',)
QUANT_BACKENDS = ('auto', 'kallisto', 'oarfish')
OARFISH_SEQ_TECHS = ('auto', 'ont-cdna', 'ont-drna', 'pac-bio', 'pac-bio-hifi')
SHORT_READ_SPOT_LENGTH_THRESHOLD = 1000

PLATFORM_REGEX_BY_FAMILY = {
    'ont': (
        r'oxford[\s_/-]*nanopore',
        r'\bnanopore\b',
        r'\bpromethion\b',
        r'\bgridion\b',
        r'\bminion\b',
        r'\bont\b',
    ),
    'pacbio': (
        r'\bpacbio\b',
        r'pacific[\s_/-]*biosciences',
        r'\bsmrt\b',
        r'\bsequel\b',
        r'\brevio\b',
        r'\brs[\s_-]*ii\b',
    ),
    'illumina': (
        r'\billumina\b',
        r'\bnovaseq\b',
        r'\bhiseq\b',
        r'\bnextseq\b',
        r'\bmiseq\b',
        r'\biseq\b',
        r'genome[\s_-]*analyzer',
    ),
}

OARFISH_PACBIO_HIFI_PATTERNS = (
    r'\bhifi\b',
    r'\bccs\b',
    r'\biso[\s_-]*seq\b',
    r'\bisoseq\b',
)
OARFISH_ONT_DIRECT_RNA_PATTERNS = (
    r'\bdirect[\s_-]*rna\b',
    r'\bd[\s_-]*rna\b',
    r'\bdrna\b',
)


def purge_existing_quant_outputs(sra_id, output_dir):
    stale_names = [
        sra_id + '_abundance.tsv',
        sra_id + '_run_info.json',
        sra_id + '_abundance.h5',
        'abundance.tsv',
        'run_info.json',
        'abundance.h5',
        sra_id + '.quant',
        sra_id + '.meta_info.json',
        sra_id + '.ambig_info.tsv',
        sra_id + '.infreps.pq',
        sra_id + '.prob',
        sra_id + '.prob.lz4',
    ]
    for stale_name in stale_names:
        stale_path = os.path.join(output_dir, stale_name)
        if not os.path.exists(stale_path):
            continue
        if not os.path.isfile(stale_path):
            raise IsADirectoryError('Quant output path exists but is not a file: {}'.format(stale_path))
        os.remove(stale_path)


def quant_output_exists(sra_id, output_dir):
    abundance_path = os.path.join(output_dir, sra_id + '_abundance.tsv')
    run_info_path = os.path.join(output_dir, sra_id + '_run_info.json')
    has_abundance = os.path.isfile(abundance_path)
    has_run_info = os.path.isfile(run_info_path)
    if has_abundance and has_run_info:
        print('Output files detected: {}, {}'.format(abundance_path, run_info_path))
        return True
    if not has_abundance:
        print('Output file was not detected: {}'.format(abundance_path))
    if not has_run_info:
        print('Output file was not detected: {}'.format(run_info_path))
    return False


def _normalize_metadata_text(value):
    if pandas.isna(value):
        return ''
    text = str(value).strip()
    if text.lower() in {'', 'nan', 'none'}:
        return ''
    return re.sub(r'\s+', ' ', text)


def _normalize_search_text(text):
    text = _normalize_metadata_text(text).lower()
    return re.sub(r'[^0-9a-z]+', ' ', text).strip()


def _metadata_row_index(metadata, sra_id):
    return get_metadata_row_index_by_run(metadata, sra_id)


def _metadata_value(metadata, sra_id, column_name, default=''):
    if column_name not in metadata.df.columns:
        return default
    return metadata.df.at[_metadata_row_index(metadata, sra_id), column_name]


def _metadata_numeric_value(metadata, sra_id, column_name):
    if column_name not in metadata.df.columns:
        return numpy.nan
    return pandas.to_numeric(metadata.df.at[_metadata_row_index(metadata, sra_id), column_name], errors='coerce')


def _collect_metadata_text_fragments(metadata, sra_id):
    fragments = []
    for column_name in [
        'platform',
        'instrument',
        'protocol',
        'lib_strategy',
        'lib_selection',
        'lib_source',
        'lib_name',
        'exp_title',
        'design',
        'sample_title',
        'sample_description',
    ]:
        text = _normalize_metadata_text(_metadata_value(metadata, sra_id, column_name, default=''))
        if text != '':
            fragments.append(text)
    return fragments


def _describe_metadata_backend_context(metadata, sra_id):
    values = []
    for column_name in ['platform', 'instrument', 'lib_layout', 'spot_length', 'total_spots', 'total_bases']:
        if column_name in metadata.df.columns:
            values.append('{}={}'.format(column_name, _metadata_value(metadata, sra_id, column_name, default='')))
    if not values:
        return 'no platform-related metadata columns were available'
    return ', '.join(values)


def _estimate_spot_length_from_metadata(metadata, sra_id):
    spot_length = _metadata_numeric_value(metadata, sra_id, 'spot_length')
    if numpy.isfinite(spot_length) and (float(spot_length) > 0):
        return float(spot_length)
    total_spots = _metadata_numeric_value(metadata, sra_id, 'total_spots')
    total_bases = _metadata_numeric_value(metadata, sra_id, 'total_bases')
    if (
        numpy.isfinite(total_spots)
        and numpy.isfinite(total_bases)
        and (float(total_spots) > 0)
        and (float(total_bases) > 0)
    ):
        return float(total_bases) / float(total_spots)
    return numpy.nan


def _matches_search_patterns(search_text, patterns):
    if search_text == '':
        return False
    for pattern in patterns:
        if re.search(pattern, search_text):
            return True
    return False


def infer_platform_family(metadata, sra_id):
    search_text = _normalize_search_text(' '.join(_collect_metadata_text_fragments(metadata, sra_id)))
    for family in ['ont', 'pacbio', 'illumina']:
        if _matches_search_patterns(search_text, PLATFORM_REGEX_BY_FAMILY[family]):
            return family
    layout = _normalize_metadata_text(_metadata_value(metadata, sra_id, 'lib_layout', default='')).lower()
    if layout == 'paired':
        return 'illumina'
    spot_length = _estimate_spot_length_from_metadata(metadata, sra_id)
    if numpy.isfinite(spot_length):
        if float(spot_length) <= SHORT_READ_SPOT_LENGTH_THRESHOLD:
            return 'illumina'
        return 'long-read-unknown'
    return 'unknown'


def resolve_quant_backend(args, metadata, sra_id):
    requested_backend = str(getattr(args, 'quant_backend', 'kallisto')).strip().lower()
    if requested_backend not in QUANT_BACKENDS:
        raise ValueError('Unsupported quant backend: {}'.format(requested_backend))
    if requested_backend in {'kallisto', 'oarfish'}:
        return requested_backend
    platform_family = infer_platform_family(metadata, sra_id)
    if platform_family == 'illumina':
        return 'kallisto'
    if platform_family in {'ont', 'pacbio', 'long-read-unknown'}:
        return 'oarfish'
    raise ValueError(
        'Could not infer quant backend for run {} from metadata ({}). '
        'Populate metadata platform/instrument fields or set --quant_backend explicitly.'.format(
            sra_id,
            _describe_metadata_backend_context(metadata, sra_id),
        )
    )


def resolve_oarfish_seq_tech(args, metadata, sra_id):
    requested_seq_tech = str(getattr(args, 'oarfish_seq_tech', 'auto')).strip().lower()
    if requested_seq_tech not in OARFISH_SEQ_TECHS:
        raise ValueError('Unsupported oarfish sequencing technology preset: {}'.format(requested_seq_tech))
    if requested_seq_tech != 'auto':
        return requested_seq_tech
    platform_family = infer_platform_family(metadata, sra_id)
    search_text = _normalize_search_text(' '.join(_collect_metadata_text_fragments(metadata, sra_id)))
    if platform_family == 'ont':
        if _matches_search_patterns(search_text, OARFISH_ONT_DIRECT_RNA_PATTERNS):
            return 'ont-drna'
        return 'ont-cdna'
    if platform_family == 'pacbio':
        if _matches_search_patterns(search_text, OARFISH_PACBIO_HIFI_PATTERNS):
            return 'pac-bio-hifi'
        return 'pac-bio'
    raise ValueError(
        'Could not infer oarfish sequencing technology for run {} from metadata ({}). '
        'Populate metadata platform/instrument fields or set --oarfish_seq_tech explicitly.'.format(
            sra_id,
            _describe_metadata_backend_context(metadata, sra_id),
        )
    )


def resolve_quant_backends_for_tasks(args, metadata, tasks):
    backend_by_run = {}
    oarfish_seq_tech_by_run = {}
    for sra_id, _sci_name in tasks:
        backend = resolve_quant_backend(args, metadata, sra_id)
        backend_by_run[sra_id] = backend
        if backend == 'oarfish':
            oarfish_seq_tech_by_run[sra_id] = resolve_oarfish_seq_tech(args, metadata, sra_id)
    return backend_by_run, oarfish_seq_tech_by_run


def resolve_nominal_length_for_kallisto(metadata, sra_id, sra_stat):
    nominal_length = sra_stat.get('nominal_length', numpy.nan)
    if numpy.isnan(pandas.to_numeric(nominal_length, errors='coerce')):
        try:
            idx = get_metadata_row_index_by_run(metadata, sra_id)
            nominal_length = metadata.df.at[idx, 'nominal_length']
        except AssertionError:
            nominal_length = numpy.nan
    nominal_length = pandas.to_numeric(nominal_length, errors='coerce')
    if numpy.isnan(nominal_length) or nominal_length <= 0:
        print("Could not find nominal length in metadata. Assuming fragment length.")
        nominal_length = 200
    elif nominal_length < 200:
        print('Nominal length in metadata is unusually small ({}). Setting it to 200.'.format(nominal_length))
        nominal_length = 200
    print("Fragment length set to: {}".format(nominal_length))
    fragment_sd = nominal_length / 10
    print("Fragment length standard deviation set to: {}".format(fragment_sd))
    return nominal_length, fragment_sd


def parse_quant_option_args(option_string, option_name):
    if not option_string:
        return []
    try:
        return shlex.split(option_string)
    except ValueError as e:
        raise ValueError('Invalid {} string: {}'.format(option_name, option_string)) from e


def build_kallisto_quant_command(args, in_files, lib_layout, output_dir, index, nominal_length=None, fragment_sd=None):
    extra_option_args = parse_quant_option_args(getattr(args, 'kallisto_options', None), '--kallisto_options')
    if lib_layout == 'single':
        if len(in_files) != 1:
            txt = "Library layout: {} and expected 1 input file. " \
                  "Received {} input file[s]. Please check your inputs and metadata."
            raise ValueError(txt.format(lib_layout, len(in_files)))
        return [
            'kallisto', 'quant', '--threads', str(args.threads), '--index', index, '-o', output_dir,
            '--single', '-l', str(nominal_length), '-s', str(fragment_sd),
        ] + extra_option_args + [in_files[0]]
    if lib_layout == 'paired':
        if len(in_files) != 2:
            txt = "Library layout: {} and expected 2 input files. " \
                  "Received {} input file[s]. Please check your inputs and metadata."
            raise ValueError(txt.format(lib_layout, len(in_files)))
        return [
            'kallisto', 'quant', '--threads', str(args.threads), '-i', index, '-o', output_dir,
        ] + extra_option_args + [in_files[0], in_files[1]]
    raise ValueError("Unsupported library layout: {}. Expected 'single' or 'paired'.".format(lib_layout))


def build_oarfish_quant_command(args, in_files, output_prefix, index, seq_tech):
    extra_option_args = parse_quant_option_args(getattr(args, 'oarfish_options', None), '--oarfish_options')
    if len(in_files) != 1:
        raise ValueError(
            'oarfish requires exactly one input read file per run. Received {} input file(s).'.format(len(in_files))
        )
    return [
        'oarfish', '-j', str(args.threads),
        '--reads', in_files[0],
        '--index', index,
        '--seq-tech', seq_tech,
    ] + extra_option_args + ['-o', output_prefix]


def rename_kallisto_outputs(output_dir, sra_id):
    output_specs = [
        ('run_info.json', sra_id + '_run_info.json', True),
        ('abundance.tsv', sra_id + '_abundance.tsv', True),
        ('abundance.h5', sra_id + '_abundance.h5', False),
    ]
    for src_name, dst_name, is_required in output_specs:
        src_path = os.path.join(output_dir, src_name)
        dst_path = os.path.join(output_dir, dst_name)
        if not os.path.exists(src_path):
            if is_required:
                raise FileNotFoundError('kallisto output file was not generated: {}'.format(src_path))
            continue
        if not os.path.isfile(src_path):
            raise IsADirectoryError('kallisto output path exists but is not a file: {}'.format(src_path))
        if os.path.exists(dst_path) and (not os.path.isfile(dst_path)):
            raise IsADirectoryError('kallisto renamed output path exists but is not a file: {}'.format(dst_path))
        os.replace(src_path, dst_path)


def _compute_compatibility_tpm(counts, lengths):
    counts = numpy.asarray(counts, dtype=float)
    lengths = numpy.asarray(lengths, dtype=float)
    with numpy.errstate(divide='ignore', invalid='ignore'):
        rate = counts / lengths
    rate[~numpy.isfinite(rate)] = 0.0
    denom = float(numpy.nansum(rate))
    if denom <= 0:
        return numpy.zeros(rate.shape[0], dtype=float)
    return rate * 1e6 / denom


def _normalize_oarfish_quant_columns(quant_df):
    rename_map = {}
    for column_name in quant_df.columns:
        normalized = re.sub(r'[^0-9a-z]+', '', str(column_name).strip().lower())
        if normalized in {'tname', 'targetid', 'targetname', 'name'}:
            rename_map[column_name] = 'target_id'
        elif normalized in {'len', 'length'}:
            rename_map[column_name] = 'length'
        elif normalized in {'numreads', 'numread', 'estcounts', 'estcount'}:
            rename_map[column_name] = 'est_counts'
        elif normalized == 'tpm':
            rename_map[column_name] = 'tpm'
    return quant_df.rename(columns=rename_map)


def adapt_oarfish_outputs(output_dir, sra_id, sra_stat, output_prefix, seq_tech):
    quant_path = output_prefix + '.quant'
    meta_info_path = output_prefix + '.meta_info.json'
    if not os.path.exists(quant_path):
        raise FileNotFoundError('oarfish output file was not generated: {}'.format(quant_path))
    if not os.path.isfile(quant_path):
        raise IsADirectoryError('oarfish output path exists but is not a file: {}'.format(quant_path))
    if not os.path.exists(meta_info_path):
        raise FileNotFoundError('oarfish output file was not generated: {}'.format(meta_info_path))
    if not os.path.isfile(meta_info_path):
        raise IsADirectoryError('oarfish output path exists but is not a file: {}'.format(meta_info_path))

    quant_df = pandas.read_csv(quant_path, sep='\t', header=0)
    quant_df = _normalize_oarfish_quant_columns(quant_df)
    required_columns = ['target_id', 'length', 'est_counts']
    missing_columns = [col for col in required_columns if col not in quant_df.columns]
    if missing_columns:
        raise ValueError(
            'oarfish quant output is missing required column(s): {}. Found columns: {}'.format(
                ', '.join(missing_columns),
                ', '.join([str(col) for col in quant_df.columns.tolist()]),
            )
        )

    target_ids = quant_df['target_id'].astype(str)
    lengths = pandas.to_numeric(quant_df['length'], errors='coerce').fillna(0.0)
    est_counts = pandas.to_numeric(quant_df['est_counts'], errors='coerce').fillna(0.0)
    if 'tpm' in quant_df.columns:
        tpm = pandas.to_numeric(quant_df['tpm'], errors='coerce').fillna(0.0)
    else:
        tpm = pandas.Series(_compute_compatibility_tpm(est_counts.to_numpy(), lengths.to_numpy()))

    abundance_df = pandas.DataFrame({
        'target_id': target_ids,
        'length': lengths,
        'eff_length': lengths,
        'est_counts': est_counts,
        'tpm': tpm,
    })
    abundance_path = os.path.join(output_dir, sra_id + '_abundance.tsv')
    abundance_df.to_csv(abundance_path, sep='\t', index=False)

    with open(meta_info_path) as meta_handle:
        meta_info = json.load(meta_handle)
    if not isinstance(meta_info, dict):
        meta_info = {'oarfish_meta_info': meta_info}
    mapped_reads = float(est_counts.sum())
    total_reads = float(sra_stat['total_spot'])
    p_pseudoaligned = 0.0
    if total_reads > 0:
        p_pseudoaligned = mapped_reads / total_reads * 100.0
    p_pseudoaligned = min(max(float(p_pseudoaligned), 0.0), 100.0)
    run_info = dict(meta_info)
    run_info['quant_backend'] = 'oarfish'
    run_info['oarfish_seq_tech'] = seq_tech
    run_info['num_processed'] = total_reads
    run_info['num_pseudoaligned'] = mapped_reads
    run_info['p_pseudoaligned'] = p_pseudoaligned
    run_info_path = os.path.join(output_dir, sra_id + '_run_info.json')
    with open(run_info_path, 'w') as run_info_handle:
        json.dump(run_info, run_info_handle, indent=2, sort_keys=True)


def call_kallisto(args, in_files, metadata, sra_stat, output_dir, index):
    sra_id = sra_stat['sra_id']
    lib_layout = sra_stat['layout']
    if lib_layout == 'single':
        print('Single end reads detected. Proceeding in single mode')
        nominal_length, fragment_sd = resolve_nominal_length_for_kallisto(metadata, sra_id, sra_stat)
        kallisto_cmd = build_kallisto_quant_command(
            args=args,
            in_files=in_files,
            lib_layout=lib_layout,
            output_dir=output_dir,
            index=index,
            nominal_length=nominal_length,
            fragment_sd=fragment_sd,
        )
    else:
        if lib_layout == 'paired':
            print('Paired-end reads detected. Running in paired read mode.')
        kallisto_cmd = build_kallisto_quant_command(
            args=args,
            in_files=in_files,
            lib_layout=lib_layout,
            output_dir=output_dir,
            index=index,
        )

    kallisto_out, stdout_txt, stderr_txt = run_logged_command(
        command=kallisto_cmd,
        runner=subprocess.run,
        print_command=True,
        print_output=True,
        stdout_label='kallisto quant stdout:',
        stderr_label='kallisto quant stderr:',
    )
    if kallisto_out.returncode != 0:
        sys.stderr.write('kallisto did not finish safely.\n')
        if 'Zero reads pseudoaligned' in stderr_txt:
            sys.stderr.write(
                'No reads are mapped to the reference. This sample will be excluded by downstream filtering/final export.'
            )
        raise RuntimeError(
            'kallisto quant failed with exit code {} for {}.'.format(kallisto_out.returncode, sra_id)
        )

    rename_kallisto_outputs(output_dir=output_dir, sra_id=sra_id)
    return kallisto_out


def call_oarfish(args, in_files, metadata, sra_stat, output_dir, index, seq_tech):
    _ = metadata
    sra_id = sra_stat['sra_id']
    if sra_stat['layout'] != 'single':
        raise ValueError(
            'oarfish backend currently supports single-end long-read runs only. '
            'Run {} was labeled as {}.'.format(sra_id, sra_stat['layout'])
        )
    output_prefix = os.path.join(output_dir, sra_id)
    oarfish_cmd = build_oarfish_quant_command(
        args=args,
        in_files=in_files,
        output_prefix=output_prefix,
        index=index,
        seq_tech=seq_tech,
    )
    oarfish_out, _stdout_txt, _stderr_txt = run_logged_command(
        command=oarfish_cmd,
        runner=subprocess.run,
        print_command=True,
        print_output=True,
        stdout_label='oarfish quant stdout:',
        stderr_label='oarfish quant stderr:',
    )
    if oarfish_out.returncode != 0:
        raise RuntimeError(
            'oarfish quant failed with exit code {} for {}.'.format(oarfish_out.returncode, sra_id)
        )
    adapt_oarfish_outputs(
        output_dir=output_dir,
        sra_id=sra_id,
        sra_stat=sra_stat,
        output_prefix=output_prefix,
        seq_tech=seq_tech,
    )
    return oarfish_out


def check_kallisto_dependency():
    probe_dependency_command(
        command=['kallisto', 'version'],
        label='kallisto',
        runner=subprocess.run,
    )


def check_oarfish_dependency():
    probe_dependency_command(
        command=['oarfish', '--version'],
        label='oarfish',
        runner=subprocess.run,
    )


def check_quant_dependencies(backend_by_run):
    backends = set([backend for backend in backend_by_run.values() if backend != ''])
    if 'kallisto' in backends:
        check_kallisto_dependency()
    if 'oarfish' in backends:
        check_oarfish_dependency()


def list_getfastq_run_files(output_dir):
    return list_dir_entries(output_dir)


def list_dir_entries(path_dir):
    try:
        with os.scandir(path_dir) as entries:
            return {
                entry.name
                for entry in entries
                if entry.is_file()
            }
    except FileNotFoundError:
        return set()


def _normalize_index_suffixes(suffixes):
    if suffixes is None:
        return None
    return tuple([str(suffix).lower() for suffix in suffixes])


def _strip_filename_suffix(filename, suffixes):
    filename = str(filename)
    filename_lower = filename.lower()
    for suffix in sorted(suffixes, key=len, reverse=True):
        suffix_lower = str(suffix).lower()
        if filename_lower.endswith(suffix_lower):
            return filename[:-len(suffix_lower)]
    return filename


def _normalize_species_stem_from_entry(entry_name, suffixes, trailing_suffixes=()):
    stem = _strip_filename_suffix(os.path.basename(str(entry_name)), suffixes)
    stem_lower = stem.lower()
    for trailing_suffix in trailing_suffixes:
        trailing_suffix_lower = str(trailing_suffix).lower()
        if stem_lower.endswith(trailing_suffix_lower):
            stem = stem[:-len(trailing_suffix_lower)]
            break
    return _normalize_species_identifier(stem)


def _find_species_prefixed_files(
    path_dir,
    species_prefix,
    entries=None,
    suffixes=None,
    exact_stem_suffixes=(),
):
    if entries is None:
        entries = list_dir_entries(path_dir)
    normalized_suffixes = _normalize_index_suffixes(suffixes)
    matched_paths = [
        os.path.join(path_dir, entry)
        for entry in find_species_prefixed_entries(entries, species_prefix)
        if (
            (normalized_suffixes is None)
            or str(entry).lower().endswith(normalized_suffixes)
        )
        and os.path.isfile(os.path.join(path_dir, entry))
    ]
    normalized_prefix = _normalize_species_identifier(species_prefix)
    exact_paths = [
        path
        for path in matched_paths
        if _normalize_species_stem_from_entry(
            path,
            suffixes=suffixes or (),
            trailing_suffixes=exact_stem_suffixes,
        ) == normalized_prefix
    ]
    return exact_paths, matched_paths


def _is_nonempty_regular_file(path):
    try:
        return os.path.isfile(path) and os.path.getsize(path) > 0
    except OSError:
        return False


def _remove_empty_regular_file(path):
    try:
        if os.path.isfile(path) and os.path.getsize(path) == 0:
            os.remove(path)
            print('Removed empty index file before rebuilding: {}'.format(path), flush=True)
    except OSError as exc:
        raise RuntimeError('Failed to remove empty index file {}: {}'.format(path, exc)) from exc


def find_species_index_files(index_dir, sci_name, entries=None, suffixes=(KALLISTO_INDEX_SUFFIX,)):
    if entries is None:
        entries = list_dir_entries(index_dir)
    normalized_suffixes = _normalize_index_suffixes(suffixes)
    normalized_prefix = _normalize_species_identifier(sci_name)
    return [
        os.path.join(index_dir, entry)
        for entry in entries
        if (
            (normalized_suffixes is None)
            or str(entry).lower().endswith(normalized_suffixes)
        )
        and (
            _normalize_species_stem_from_entry(
                entry,
                suffixes=suffixes or (),
            ) == normalized_prefix
        )
        and _is_nonempty_regular_file(os.path.join(index_dir, entry))
    ]


def find_species_fasta_files(path_fasta_dir, sci_name, entries=None, alias_names=None):
    _, matched_files = _find_species_fasta_files(
        path_fasta_dir,
        sci_name,
        entries=entries,
        alias_names=alias_names,
    )
    return matched_files


def _find_species_fasta_files(path_fasta_dir, sci_name, entries=None, alias_names=None):
    _ = alias_names
    if os.path.exists(path_fasta_dir) and (not os.path.isdir(path_fasta_dir)):
        raise NotADirectoryError('Fasta path exists but is not a directory: {}'.format(path_fasta_dir))
    if entries is None:
        entries = list_dir_entries(path_fasta_dir)
    species_prefixes = _build_species_identifier_candidates_from_values(sci_name)
    for species_prefix in species_prefixes:
        exact_matches, _matched_paths = _find_species_prefixed_files(
            path_fasta_dir,
            species_prefix,
            entries=entries,
            suffixes=INDEX_FASTA_SUFFIXES,
            exact_stem_suffixes=FASTA_REFERENCE_STEM_SUFFIXES,
        )
        if exact_matches:
            return species_prefix, exact_matches
    return None, []


def _normalize_fasta_candidate_key(fasta_file):
    filename = os.path.basename(fasta_file)
    filename_lower = filename.lower()
    for suffix in sorted(INDEX_FASTA_SUFFIXES, key=len, reverse=True):
        if filename_lower.endswith(suffix):
            return filename_lower[:-len(suffix)]
    return filename_lower


def _resolve_duplicate_fasta_candidates(fasta_files):
    preferred_by_key = {}
    for fasta_file in fasta_files:
        key = _normalize_fasta_candidate_key(fasta_file)
        existing = preferred_by_key.get(key)
        if existing is None:
            preferred_by_key[key] = fasta_file
            continue
        existing_is_gz = str(existing).lower().endswith('.gz')
        current_is_gz = str(fasta_file).lower().endswith('.gz')
        if existing_is_gz and (not current_is_gz):
            preferred_by_key[key] = fasta_file
    return sorted(preferred_by_key.values())


def check_layout_mismatch(sra_stat, output_dir, files=None):
    if files is None:
        files = list_getfastq_run_files(output_dir)
    if sra_stat['layout'] == 'paired':
        fastq_files = [
            f for f in find_run_prefixed_entries(files, sra_stat['sra_id'])
            if '.fastq' in f
        ]
        if len(fastq_files) == 1:
            sys.stderr.write('Single-end fastq was detected even though layout = {}\n'.format(sra_stat['layout']))
            sys.stderr.write('This sample will be treated as single-end sequencing.\n')
            sra_stat['layout'] = 'single'
    return sra_stat


def resolve_input_fastq_files(sra_stat, output_dir_getfastq, ext, files=None):
    if files is None:
        files = list_getfastq_run_files(output_dir_getfastq)
    sra_id = sra_stat['sra_id']
    if sra_stat['layout'] == 'paired':
        pair1 = sra_id + '_1' + ext
        pair2 = sra_id + '_2' + ext
        if (pair1 in files) and (pair2 in files):
            return [
                os.path.join(output_dir_getfastq, pair1),
                os.path.join(output_dir_getfastq, pair2),
            ]
    elif sra_stat['layout'] == 'single':
        single = sra_id + ext
        if single in files:
            return [os.path.join(output_dir_getfastq, single)]

    matched = sorted([
        f for f in find_run_prefixed_entries(files, sra_id)
        if f.endswith(ext)
    ])
    return [os.path.join(output_dir_getfastq, f) for f in matched]


def run_quant(args, metadata, sra_id, index, runtime_context=None, backend=None, oarfish_seq_tech=None):
    output_dir = os.path.join(args.out_dir, 'quant', sra_id)
    if os.path.exists(output_dir) and (not os.path.isdir(output_dir)):
        raise NotADirectoryError('Quant run output path exists but is not a directory: {}'.format(output_dir))
    os.makedirs(output_dir, exist_ok=True)
    if args.redo:
        purge_existing_quant_outputs(sra_id=sra_id, output_dir=output_dir)
    is_quant_output_available = quant_output_exists(sra_id, output_dir)
    if is_quant_output_available:
        if args.redo:
            print('The output will be overwritten. Set "--redo no" to not overwrite results.')
        else:
            print('Continued. The output will not be overwritten. If you want to overwrite the results, set "--redo yes".')
            return
    output_dir_getfastq = os.path.join(args.out_dir, 'getfastq', sra_id)
    if os.path.exists(output_dir_getfastq) and (not os.path.isdir(output_dir_getfastq)):
        raise NotADirectoryError(
            'getfastq run path exists but is not a directory: {}'.format(output_dir_getfastq)
        )
    sra_stat = get_sra_stat(sra_id, metadata, num_bp_per_sra=None)
    if backend is None:
        if isinstance(runtime_context, QuantRuntimeContext) and (sra_id in runtime_context.quant_backend_by_run):
            backend = runtime_context.quant_backend_by_run[sra_id]
        else:
            backend = resolve_quant_backend(args, metadata, sra_id)
    if (backend == 'oarfish') and (oarfish_seq_tech is None):
        if isinstance(runtime_context, QuantRuntimeContext) and (sra_id in runtime_context.oarfish_seq_tech_by_run):
            oarfish_seq_tech = runtime_context.oarfish_seq_tech_by_run[sra_id]
        else:
            oarfish_seq_tech = resolve_oarfish_seq_tech(args, metadata, sra_id)
    run_files = None
    if isinstance(runtime_context, QuantRuntimeContext):
        run_files = runtime_context.run_files_by_run.get(sra_id, None)
    if run_files is None:
        run_files = list_getfastq_run_files(output_dir_getfastq)
    sra_stat = check_layout_mismatch(sra_stat, output_dir_getfastq, files=run_files)
    sra_stat['getfastq_sra_dir'] = output_dir_getfastq
    ext = get_newest_intermediate_file_extension(sra_stat, work_dir=output_dir_getfastq, files=run_files)
    if ext == '.safely_removed':
        print('These files have been safe-deleted. If you wish to re-obtain the .fastq file(s), run: getfastq --id ', sra_id, ' -w ', args.out_dir)
        print('Skipping.')
        return
    if ext == 'no_extension_found':
        sys.stderr.write('getfastq output not found in: {}, layout = {}\n'.format(sra_stat['getfastq_sra_dir'], sra_stat['layout']))
        txt = 'Exiting. If you wish to obtain the .fastq file(s), run: getfastq --id {}\n'
        sys.stderr.write(txt.format(sra_stat['sra_id']))
        raise FileNotFoundError('getfastq output not found for {}.'.format(sra_stat['sra_id']))
    in_files = resolve_input_fastq_files(sra_stat, output_dir_getfastq, ext, files=run_files)
    if len(in_files) == 0:
        run_files = list_getfastq_run_files(output_dir_getfastq)
        in_files = resolve_input_fastq_files(sra_stat, output_dir_getfastq, ext, files=run_files)
    if not in_files:
        raise FileNotFoundError('{}: Fastq file not found. Check {}'.format(sra_id, output_dir_getfastq))
    print('Input fastq detected:', ', '.join(in_files))
    print('Quant backend selected for {}: {}'.format(sra_id, backend))
    if backend == 'kallisto':
        call_kallisto(args, in_files, metadata, sra_stat, output_dir, index)
    elif backend == 'oarfish':
        print('oarfish seq-tech selected for {}: {}'.format(sra_id, oarfish_seq_tech))
        call_oarfish(args, in_files, metadata, sra_stat, output_dir, index, oarfish_seq_tech)
    else:
        raise ValueError('Unsupported quant backend: {}'.format(backend))
    if (args.clean_fastq and quant_output_exists(sra_id, output_dir)):
        print('Safe-deleting getfastq files.', flush=True)
        for in_file in in_files:
            print('Output file detected. Safely removing fastq:', in_file)
            os.remove(in_file)
            with open(in_file + '.safely_removed', 'w') as handle:
                handle.write('This fastq file was safely removed after `amalgkit quant`.')
    else:
        print('Skipping the deletion of getfastq files.', flush=True)


def _resolve_index_lock_options(args):
    lock_poll_seconds = int(getattr(args, 'index_lock_poll', INDEX_BUILD_LOCK_POLL_SECONDS))
    lock_timeout_seconds = int(getattr(args, 'index_lock_timeout', INDEX_BUILD_LOCK_TIMEOUT_SECONDS))
    if lock_poll_seconds <= 0:
        raise ValueError('--index_lock_poll must be > 0 (seconds).')
    if lock_timeout_seconds <= 0:
        raise ValueError('--index_lock_timeout must be > 0 (seconds).')
    return lock_poll_seconds, lock_timeout_seconds


def _resolve_index_dir(args):
    if args.index_dir is not None:
        index_dir = args.index_dir
    else:
        index_dir = os.path.join(args.out_dir, 'index')
    if (not os.path.exists(index_dir)) and args.build_index:
        os.makedirs(index_dir, exist_ok=True)
    if not os.path.exists(index_dir):
        raise FileNotFoundError('Could not find index folder at: {}'.format(index_dir))
    if not os.path.isdir(index_dir):
        raise NotADirectoryError('Index path exists but is not a directory: {}'.format(index_dir))
    return index_dir


def _normalize_species_identifier(text):
    normalized = str(text).strip()
    if normalized == '':
        return ''
    normalized = re.sub(r'\s+', '_', normalized)
    normalized = re.sub(r'_+', '_', normalized)
    return normalized


def _collect_species_identifier_values(sci_name, alias_names=None):
    _ = alias_names
    values = []

    def add(value):
        text = _normalize_metadata_text(value)
        if (text != '') and (text not in values):
            values.append(text)

    add(sci_name)
    return values


def _build_species_identifier_candidates(sci_name):
    normalized = _normalize_species_identifier(sci_name)
    if normalized == '':
        return []
    return [normalized]


def _build_species_identifier_candidates_from_values(sci_name, alias_names=None):
    _ = alias_names
    return _build_species_identifier_candidates(sci_name)


def _resolve_task_species_identifier_values(sra_id, sci_name, runtime_context=None):
    _ = (sra_id, runtime_context)
    return _collect_species_identifier_values(sci_name)


def _resolve_quant_species_identifier_values(metadata, tasks):
    _ = metadata
    species_identifier_values_by_run = {}
    for sra_id, sci_name in tasks:
        species_identifier_values_by_run[sra_id] = _collect_species_identifier_values(sci_name)
    return species_identifier_values_by_run


def _resolve_index_suffix(backend):
    if backend == 'kallisto':
        return KALLISTO_INDEX_SUFFIX
    if backend == 'oarfish':
        return OARFISH_INDEX_SUFFIX
    raise ValueError('Unsupported quant backend: {}'.format(backend))


def _build_index_stem(sci_name, backend='kallisto', oarfish_seq_tech=None):
    normalized_sci_name = _normalize_species_identifier(sci_name)
    if backend == 'kallisto':
        return normalized_sci_name
    if oarfish_seq_tech in {None, ''}:
        raise ValueError('oarfish index stem requires sequencing technology.')
    return '{}.{}'.format(normalized_sci_name, oarfish_seq_tech)


def _build_index_cache_key(backend, sci_name, oarfish_seq_tech=None, alias_names=None):
    _ = alias_names
    normalized_sci_name = _normalize_species_identifier(sci_name)
    if backend == 'kallisto':
        return normalized_sci_name
    return (backend, normalized_sci_name, str(oarfish_seq_tech))


def _find_single_index_file(index_dir, sci_name, entries=None, backend='kallisto', oarfish_seq_tech=None, alias_names=None):
    _ = alias_names
    backend_label = 'Kallisto' if backend == 'kallisto' else 'oarfish'
    index_suffix = _resolve_index_suffix(backend)
    species_prefixes = _build_species_identifier_candidates_from_values(sci_name)
    for species_prefix in species_prefixes:
        prefix = _build_index_stem(species_prefix, backend=backend, oarfish_seq_tech=oarfish_seq_tech)
        index_files, _matched_paths = _find_species_prefixed_files(
            index_dir,
            prefix,
            entries=entries,
            suffixes=(index_suffix,),
        )
        index_files = [path for path in index_files if _is_nonempty_regular_file(path)]
        if len(index_files) > 1:
            raise ValueError(
                'Found multiple {} index files for species. Please make sure there is only one index file for this species.'.format(
                    backend_label
                )
            )
        if len(index_files) == 1:
            index_file = index_files[0]
            print('{} index file found: {}'.format(backend_label, index_file), flush=True)
            return index_file
    return None


def _resolve_prefetched_index_entries(index_dir, runtime_context=None):
    if isinstance(runtime_context, QuantRuntimeContext):
        return runtime_context.prefetched_index.resolve_entries(index_dir)
    return None


def _find_single_fasta_match(args, sci_name, runtime_context=None, alias_names=None):
    if args.fasta_dir == 'inferred':
        path_fasta_dir = os.path.join(args.out_dir, 'fasta')
    else:
        path_fasta_dir = args.fasta_dir
    prefetched_entries = None
    if isinstance(runtime_context, QuantRuntimeContext):
        prefetched_entries = runtime_context.prefetched_fasta.resolve_entries(path_fasta_dir)
    try:
        matched_prefix, fasta_files = _find_species_fasta_files(
            path_fasta_dir=path_fasta_dir,
            sci_name=sci_name,
            entries=prefetched_entries,
            alias_names=alias_names,
        )
    except TypeError as exc:
        if ('unexpected keyword argument' not in str(exc)) or ('alias_names' not in str(exc)):
            raise
        matched_prefix, fasta_files = _find_species_fasta_files(
            path_fasta_dir=path_fasta_dir,
            sci_name=sci_name,
            entries=prefetched_entries,
        )
    fasta_files = _resolve_duplicate_fasta_candidates(fasta_files)
    if len(fasta_files) > 1:
        txt = 'Found multiple reference fasta files for this species: {}\n'
        txt += 'Please make sure there is only one index file for this species.\n{}'
        raise ValueError(txt.format(sci_name, ', '.join(fasta_files)))
    if len(fasta_files) == 0:
        txt = 'Could not find reference fasta file for this species: {}\n'.format(sci_name)
        txt += 'If the reference fasta file is correctly placed, the column "scientific_name" of the --metadata file may need to be edited.'
        raise FileNotFoundError(txt)
    return matched_prefix, fasta_files[0]


def _resolve_single_fasta_file(args, sci_name, runtime_context=None, alias_names=None):
    _matched_prefix, fasta_file = _find_single_fasta_match(
        args,
        sci_name,
        runtime_context=runtime_context,
        alias_names=alias_names,
    )
    return fasta_file


def _build_kallisto_index(index_path, fasta_file, sci_name):
    index_path = os.path.realpath(index_path)
    fasta_file = os.path.realpath(fasta_file)
    index_dir = os.path.dirname(index_path)
    print('Reference fasta file found: {}'.format(fasta_file), flush=True)
    print('Building index: {}'.format(index_path), flush=True)
    kallisto_build_cmd = ['kallisto', 'index', '-i', index_path, fasta_file]

    with tempfile.TemporaryDirectory(prefix='amalgkit_kallisto_index_', dir=index_dir) as build_cwd:
        print('Using isolated kallisto index work directory: {}'.format(build_cwd), flush=True)

        def run_in_build_cwd(command, stdout, stderr):
            try:
                return subprocess.run(command, stdout=stdout, stderr=stderr, cwd=build_cwd)
            except TypeError as exc:
                if 'cwd' not in str(exc):
                    raise
                return subprocess.run(command, stdout=stdout, stderr=stderr)

        index_out, _stdout_txt, _stderr_txt = run_logged_command(
            command=kallisto_build_cmd,
            runner=run_in_build_cwd,
            print_command=True,
            print_output=True,
            stdout_label='kallisto index stdout:',
            stderr_label='kallisto index stderr:',
        )
    if index_out.returncode != 0:
        raise RuntimeError('kallisto index failed for {}.'.format(sci_name))
    if not _is_nonempty_regular_file(index_path):
        raise RuntimeError('Index file was not generated: {}'.format(index_path))


def _build_oarfish_index(args, index_path, fasta_file, sci_name, seq_tech):
    print('Reference fasta file found: {}'.format(fasta_file), flush=True)
    print('Building index: {}'.format(index_path), flush=True)
    oarfish_build_cmd = [
        'oarfish', '-j', str(args.threads),
        '--annotated', fasta_file,
        '--seq-tech', seq_tech,
        '--only-index',
        '--index-out', index_path,
    ]
    index_out, _stdout_txt, _stderr_txt = run_logged_command(
        command=oarfish_build_cmd,
        runner=subprocess.run,
        print_command=True,
        print_output=True,
        stdout_label='oarfish index stdout:',
        stderr_label='oarfish index stderr:',
    )
    if index_out.returncode != 0:
        raise RuntimeError('oarfish index failed for {}.'.format(sci_name))
    if not _is_nonempty_regular_file(index_path):
        raise RuntimeError('Index file was not generated: {}'.format(index_path))


def get_index(args, sci_name, runtime_context=None, backend=None, oarfish_seq_tech=None, alias_names=None):
    _ = alias_names
    backend = str(backend if backend is not None else getattr(args, 'quant_backend', 'kallisto')).strip().lower()
    if backend == 'auto':
        backend = 'kallisto'
    lock_poll_seconds, lock_timeout_seconds = _resolve_index_lock_options(args)
    index_dir = _resolve_index_dir(args)
    build_sci_name = sci_name
    index_stem = _build_index_stem(build_sci_name, backend=backend, oarfish_seq_tech=oarfish_seq_tech)
    index_suffix = _resolve_index_suffix(backend)
    index_path = os.path.join(index_dir, index_stem + index_suffix)
    lock_path = os.path.join(index_dir, '.{}.lock'.format(os.path.basename(index_path)))
    prefetched_index_lookup_entries = _resolve_prefetched_index_entries(index_dir, runtime_context=runtime_context)
    use_prefetched_index = prefetched_index_lookup_entries is not None

    if use_prefetched_index:
        index = _find_single_index_file(
            index_dir=index_dir,
            sci_name=sci_name,
            entries=prefetched_index_lookup_entries,
            backend=backend,
            oarfish_seq_tech=oarfish_seq_tech,
        )
    else:
        index = _find_single_index_file(
            index_dir=index_dir,
            sci_name=sci_name,
            backend=backend,
            oarfish_seq_tech=oarfish_seq_tech,
        )
    if index is not None:
        return index

    if not args.build_index:
        sys.stderr.write('No index file was found in: {}\n'.format(index_dir))
        sys.stderr.write('Try --fasta_dir PATH and --build_index yes\n')
        raise FileNotFoundError('Could not find index file.')

    with acquire_exclusive_lock(
        lock_path=lock_path,
        lock_label='Index lock',
        poll_seconds=lock_poll_seconds,
        timeout_seconds=lock_timeout_seconds,
    ):
        index = _find_single_index_file(
            index_dir=index_dir,
            sci_name=sci_name,
            entries=prefetched_index_lookup_entries,
            backend=backend,
            oarfish_seq_tech=oarfish_seq_tech,
        )
        if index is not None:
            print('Detected completed index after lock acquisition: {}'.format(index), flush=True)
            return index

        print('--build_index set. Building index for {}'.format(sci_name))
        _remove_empty_regular_file(index_path)
        matched_prefix, fasta_file = _find_single_fasta_match(
            args,
            sci_name,
            runtime_context=runtime_context,
        )
        if backend == 'kallisto':
            _build_kallisto_index(index_path=index_path, fasta_file=fasta_file, sci_name=build_sci_name)
        elif backend == 'oarfish':
            _build_oarfish_index(
                args=args,
                index_path=index_path,
                fasta_file=fasta_file,
                sci_name=build_sci_name,
                seq_tech=oarfish_seq_tech,
            )
        else:
            raise ValueError('Unsupported quant backend: {}'.format(backend))
        print('Index file found: {}'.format(index_path), flush=True)
        return index_path


def _get_index_for_backend(args, sci_name, runtime_context=None, backend='kallisto', oarfish_seq_tech=None, alias_names=None):
    _ = alias_names
    if (backend == 'kallisto') and (oarfish_seq_tech in {None, ''}):
        return get_index(args, sci_name, runtime_context=runtime_context)
    return get_index(
        args,
        sci_name,
        runtime_context=runtime_context,
        backend=backend,
        oarfish_seq_tech=oarfish_seq_tech,
    )


def run_quant_for_sra(args, metadata, sra_id, sci_name, runtime_context=None):
    print('')
    print('Species: {}'.format(sci_name))
    print('SRA Run ID: {}'.format(sra_id))
    normalized_sci_name = _normalize_species_identifier(sci_name)
    backend = None
    oarfish_seq_tech = None
    index_cache = None
    if isinstance(runtime_context, QuantRuntimeContext):
        backend = runtime_context.quant_backend_by_run.get(sra_id, None)
        oarfish_seq_tech = runtime_context.oarfish_seq_tech_by_run.get(sra_id, None)
        index_cache = runtime_context.resolved_index_cache
    if backend is None:
        backend = resolve_quant_backend(args, metadata, sra_id)
    if (backend == 'oarfish') and (oarfish_seq_tech is None):
        oarfish_seq_tech = resolve_oarfish_seq_tech(args, metadata, sra_id)
    index_cache_key = _build_index_cache_key(
        backend,
        normalized_sci_name,
        oarfish_seq_tech=oarfish_seq_tech,
    )
    if isinstance(index_cache, dict) and (index_cache_key in index_cache):
        index = index_cache[index_cache_key]
        print('Using pre-resolved index: {}'.format(index))
    else:
        print('Looking for index folder in {}'.format(args.out_dir))
        index = _get_index_for_backend(
            args,
            normalized_sci_name,
            runtime_context=runtime_context,
            backend=backend,
            oarfish_seq_tech=oarfish_seq_tech,
        )
    run_quant(
        args,
        metadata,
        sra_id,
        index,
        runtime_context=runtime_context,
        backend=backend,
        oarfish_seq_tech=oarfish_seq_tech,
    )


def _resolve_task_backend_info(args, tasks, runtime_context=None):
    task_info = []
    for sra_id, sci_name in tasks:
        backend = 'kallisto'
        oarfish_seq_tech = None
        if isinstance(runtime_context, QuantRuntimeContext):
            backend = runtime_context.quant_backend_by_run.get(sra_id, backend)
            oarfish_seq_tech = runtime_context.oarfish_seq_tech_by_run.get(sra_id, None)
        elif str(getattr(args, 'quant_backend', 'kallisto')).strip().lower() == 'oarfish':
            backend = 'oarfish'
            requested_seq_tech = str(getattr(args, 'oarfish_seq_tech', 'auto')).strip().lower()
            if requested_seq_tech == 'auto':
                raise ValueError(
                    'pre_resolve_species_indices requires metadata-backed oarfish seq-tech resolution '
                    'or an explicit --oarfish_seq_tech override.'
                )
            oarfish_seq_tech = requested_seq_tech
        task_info.append((
            sra_id,
            _normalize_species_identifier(sci_name),
            backend,
            oarfish_seq_tech,
        ))
    return task_info


def pre_resolve_species_indices(args, tasks, runtime_context=None):
    if runtime_context is None:
        runtime_context = QuantRuntimeContext()
    target_info = _resolve_task_backend_info(args, tasks, runtime_context=runtime_context)
    unique_targets = sorted(set([
        (backend, sci_name, '' if oarfish_seq_tech is None else oarfish_seq_tech)
        for _sra_id, sci_name, backend, oarfish_seq_tech in target_info
    ]))
    if len(unique_targets) == 0:
        return {}
    prefetched_fasta_entries = None
    prefetched_fasta_dir = None
    if getattr(args, 'build_index', False):
        if args.fasta_dir == 'inferred':
            path_fasta_dir = os.path.join(args.out_dir, 'fasta')
        else:
            path_fasta_dir = args.fasta_dir
        if os.path.exists(path_fasta_dir) and (not os.path.isdir(path_fasta_dir)):
            raise NotADirectoryError('Fasta path exists but is not a directory: {}'.format(path_fasta_dir))
        prefetched_fasta_entries = list_dir_entries(path_fasta_dir)
        prefetched_fasta_dir = path_fasta_dir
    runtime_context.prefetched_fasta = PrefetchedDirEntries.from_entries(
        entries=prefetched_fasta_entries,
        path_dir=prefetched_fasta_dir,
    )
    if args.index_dir is not None:
        index_dir = args.index_dir
    else:
        index_dir = os.path.join(args.out_dir, 'index')
    if os.path.exists(index_dir) and (not os.path.isdir(index_dir)):
        raise NotADirectoryError('Index path exists but is not a directory: {}'.format(index_dir))
    if (not os.path.exists(index_dir)) and args.build_index:
        os.makedirs(index_dir, exist_ok=True)
    prefetched_index_entries = list_dir_entries(index_dir)
    runtime_context.prefetched_index = PrefetchedDirEntries.from_entries(
        entries=prefetched_index_entries,
        path_dir=index_dir,
    )
    print('Resolving quant indices for {:,} target(s).'.format(len(unique_targets)), flush=True)
    requested_jobs = getattr(args, 'internal_jobs', 'auto')
    if is_auto_parallel_option(requested_jobs):
        requested_jobs = resolve_detected_cpu_count()
    try:
        max_workers = min(max(1, int(requested_jobs)), len(unique_targets))
    except (TypeError, ValueError):
        max_workers = 1

    def resolve_one_target(target):
        backend, sci_name, oarfish_seq_tech = target
        seq_tech_arg = None if oarfish_seq_tech == '' else oarfish_seq_tech
        return _get_index_for_backend(
            args,
            sci_name,
            runtime_context=runtime_context,
            backend=backend,
            oarfish_seq_tech=seq_tech_arg,
        )

    if max_workers <= 1:
        resolved = dict()
        for backend, sci_name, oarfish_seq_tech in unique_targets:
            print('Pre-resolving index for {} ({})'.format(sci_name, backend), flush=True)
            cache_key = _build_index_cache_key(
                backend,
                sci_name,
                oarfish_seq_tech=None if oarfish_seq_tech == '' else oarfish_seq_tech,
            )
            resolved[cache_key] = resolve_one_target((backend, sci_name, oarfish_seq_tech))
        return resolved

    print('Pre-resolving indices with {:,} parallel jobs.'.format(max_workers), flush=True)
    resolved_by_target, failures = run_tasks_with_optional_threads(
        task_items=unique_targets,
        task_fn=resolve_one_target,
        max_workers=max_workers,
    )
    if failures:
        details = '; '.join([
            '{}:{}: {}'.format(target[0], target[1], err)
            for target, err in failures
        ])
        raise RuntimeError(
            'Failed to pre-resolve index for {}/{} target(s). {}'.format(
                len(failures),
                len(unique_targets),
                details,
            )
        )
    resolved = dict()
    for backend, sci_name, oarfish_seq_tech in unique_targets:
        cache_key = _build_index_cache_key(
            backend,
            sci_name,
            oarfish_seq_tech=None if oarfish_seq_tech == '' else oarfish_seq_tech,
        )
        resolved[cache_key] = resolved_by_target[(backend, sci_name, oarfish_seq_tech)]
    return resolved


def prefetch_getfastq_run_files(args, tasks):
    run_ids = sorted(set([sra_id for sra_id, _ in tasks]))
    if len(run_ids) == 0:
        return {}
    getfastq_root = os.path.join(args.out_dir, 'getfastq')
    prefetched = {}
    for run_id in run_ids:
        run_dir = os.path.join(getfastq_root, run_id)
        try:
            with os.scandir(run_dir) as run_entries:
                prefetched[run_id] = {
                    run_entry.name
                    for run_entry in run_entries
                    if run_entry.is_file()
                }
        except (FileNotFoundError, NotADirectoryError):
            continue
    return prefetched


def prepare_quant_runtime_context(args, tasks, metadata=None, backend_by_run=None, oarfish_seq_tech_by_run=None):
    runtime_context = QuantRuntimeContext()
    runtime_context.run_files_by_run = prefetch_getfastq_run_files(args, tasks)
    if backend_by_run is None or oarfish_seq_tech_by_run is None:
        if metadata is None:
            backend_by_run = {sra_id: 'kallisto' for sra_id, _sci_name in tasks}
            oarfish_seq_tech_by_run = {}
        else:
            backend_by_run, oarfish_seq_tech_by_run = resolve_quant_backends_for_tasks(args, metadata, tasks)
    runtime_context.quant_backend_by_run = dict(backend_by_run)
    runtime_context.oarfish_seq_tech_by_run = dict(oarfish_seq_tech_by_run)
    runtime_context.species_identifier_values_by_run = _resolve_quant_species_identifier_values(metadata, tasks)
    runtime_context.resolved_index_cache = pre_resolve_species_indices(
        args,
        tasks,
        runtime_context=runtime_context,
    )
    return runtime_context


def build_quant_tasks(metadata):
    required_columns = ['run', 'scientific_name']
    missing_columns = [col for col in required_columns if col not in metadata.df.columns]
    if len(missing_columns) > 0:
        raise ValueError('Missing required metadata column(s) for quant: {}'.format(', '.join(missing_columns)))
    runs = metadata.df['run'].fillna('').astype(str).str.strip()
    species = metadata.df['scientific_name'].fillna('').astype(str).str.strip()
    metadata.df['run'] = runs
    metadata.df['scientific_name'] = species
    missing_species_runs = []
    placeholder_species_runs = []
    missing_run_count = 0
    duplicate_runs = []
    seen_runs = set()
    tasks = []
    for run_id, sci_name in zip(runs.tolist(), species.tolist()):
        if run_id == '':
            missing_run_count += 1
            continue
        if sci_name == '':
            missing_species_runs.append(run_id)
            continue
        if is_private_fastq_scientific_name_placeholder(sci_name):
            placeholder_species_runs.append(run_id)
            continue
        if run_id in seen_runs:
            duplicate_runs.append(run_id)
            continue
        seen_runs.add(run_id)
        tasks.append((run_id, sci_name))
    if missing_run_count > 0:
        raise ValueError('Missing run ID in metadata for {:,} row(s).'.format(missing_run_count))
    if len(missing_species_runs) > 0:
        raise ValueError('Missing scientific_name in metadata for run(s): {}'.format(', '.join(missing_species_runs)))
    if len(placeholder_species_runs) > 0:
        raise ValueError(
            'Placeholder scientific_name from amalgkit integrate was found for run(s): {}. '
            'Edit the "scientific_name" column before running quant.'.format(
                ', '.join(placeholder_species_runs)
            )
        )
    if len(duplicate_runs) > 0:
        duplicate_runs = list(dict.fromkeys(duplicate_runs))
        raise ValueError('Duplicate run ID in metadata for run(s): {}'.format(', '.join(duplicate_runs)))
    if len(tasks) == 0:
        raise ValueError('No valid run/scientific_name entries were found in metadata.')
    return tasks


def quant_main(args):
    threads, jobs, _ = resolve_thread_worker_allocation(
        requested_threads=getattr(args, 'threads', 'auto'),
        requested_workers=getattr(args, 'internal_jobs', 'auto'),
        internal_cpu_budget=getattr(args, 'internal_cpu_budget', 'auto'),
        worker_option_name='internal_jobs',
        context='quant:',
        disable_workers=(getattr(args, 'batch', None) is not None),
    )
    out_dir = os.path.realpath(args.out_dir)
    if os.path.exists(out_dir) and (not os.path.isdir(out_dir)):
        raise NotADirectoryError('Output path exists but is not a directory: {}'.format(out_dir))
    quant_dir = os.path.join(out_dir, 'quant')
    if os.path.exists(quant_dir) and (not os.path.isdir(quant_dir)):
        raise NotADirectoryError('Quant path exists but is not a directory: {}'.format(quant_dir))
    runtime_args = clone_namespace(args, threads=threads, internal_jobs=jobs, out_dir=out_dir)
    metadata = load_metadata(runtime_args)
    tasks = build_quant_tasks(metadata)
    backend_by_run, oarfish_seq_tech_by_run = resolve_quant_backends_for_tasks(runtime_args, metadata, tasks)
    check_quant_dependencies(backend_by_run)
    runtime_context = prepare_quant_runtime_context(
        runtime_args,
        tasks,
        metadata=metadata,
        backend_by_run=backend_by_run,
        oarfish_seq_tech_by_run=oarfish_seq_tech_by_run,
    )
    if (jobs == 1) or (len(tasks) <= 1):
        for sra_id, sci_name in tasks:
            run_quant_for_sra(runtime_args, metadata, sra_id, sci_name, runtime_context=runtime_context)
        return

    max_workers = min(jobs, len(tasks))
    print('Running quant for {:,} SRA runs with {:,} parallel jobs.'.format(len(tasks), max_workers), flush=True)
    _, failures = run_tasks_with_optional_threads(
        task_items=tasks,
        task_fn=lambda task: run_quant_for_sra(
            runtime_args,
            metadata,
            task[0],
            task[1],
            runtime_context=runtime_context,
        ),
        max_workers=max_workers,
    )
    if failures:
        details = '; '.join(['{}: {}'.format(task[0], err) for task, err in failures])
        raise RuntimeError('quant failed for {}/{} SRA runs. {}'.format(len(failures), len(tasks), details))
