import datetime
import json
import os
import re
import sys
import time  # noqa: F401
import xml.etree.ElementTree as ET  # noqa: F401

import pandas

from Bio import Entrez

from amalgkit.__init__ import __version__
from amalgkit.arg_utils import clone_namespace
from amalgkit.exceptions import AmalgkitExit
from amalgkit.metadata_utils import Metadata
from amalgkit.output_utils import atomic_output_path, sanitize_dataframe_for_tsv
from amalgkit.parallel_utils import (
    is_auto_parallel_option,
    resolve_worker_allocation,
    run_tasks_with_optional_threads,
)
from amalgkit.sra import (
    esearch_sra_with_retry as _esearch_sra_with_retry,
    fetch_sra_xml as _fetch_sra_xml,
    fetch_sra_xml_chunk as _fetch_sra_xml_chunk,
    inspect_accession_search_mismatches as _inspect_accession_search_mismatches,
    iter_sra_xml_chunks as _iter_sra_xml_chunks,
    merge_xml_chunk as _merge_xml_chunk,
    raise_if_xml_has_error as _raise_if_xml_has_error,
    search_sra_record_ids as _search_sra_record_ids,
)

BASE_QUERY_TEMPLATE = (
    '("platform illumina"[Properties]) AND '
    '("type rnaseq"[Filter]) AND '
    '("sra biosample"[Filter]) AND '
    '("{species}"[Organism])'
)
DEFAULT_ORGAN_TERMS = {
    'flower': ['flower'],
    'leaf': ['leaf'],
    'root': ['root'],
}
METADATA_AUTO_MAX_SPECIES_JOBS = 4


def merge_xml_chunk(root, chunk):
    return _merge_xml_chunk(root, chunk)


def esearch_sra_with_retry(search_term):
    return _esearch_sra_with_retry(search_term)


def fetch_sra_xml_chunk(record_ids, start, end, retmax, max_retry=10):
    return _fetch_sra_xml_chunk(record_ids, start, end, retmax, max_retry=max_retry, verbose=True)


def raise_if_xml_has_error(root):
    return _raise_if_xml_has_error(root)


def search_sra_record_ids(search_term):
    return _search_sra_record_ids(search_term, verbose=True)


def inspect_accession_search_mismatches(search_term):
    return _inspect_accession_search_mismatches(search_term)


def iter_sra_xml_chunks(record_ids, retmax=1000):
    for chunk in _iter_sra_xml_chunks(
        record_ids=record_ids,
        retmax=retmax,
        verbose=True,
        timestamp_logs=True,
        progress_label='Retrieving SRA XML',
    ):
        raise_if_xml_has_error(chunk)
        yield chunk


def fetch_sra_xml(search_term, retmax=1000):
    return _fetch_sra_xml(
        search_term=search_term,
        retmax=retmax,
        verbose=True,
        timestamp_logs=True,
        progress_label='Retrieving SRA XML',
    )


def _normalize_search_string(search_term_raw):
    if search_term_raw is None:
        return ''
    return str(search_term_raw).strip()


def _realpath_or_empty(path):
    if path in [None, '']:
        return ''
    return os.path.realpath(path)


def _ensure_directory(path, label):
    real_path = os.path.realpath(path)
    if os.path.exists(real_path) and (not os.path.isdir(real_path)):
        raise NotADirectoryError('{} exists but is not a directory: {}'.format(label, real_path))
    os.makedirs(real_path, exist_ok=True)
    return real_path


def _read_json_if_exists(path):
    if (path is None) or (not os.path.exists(path)):
        return None
    with open(path, 'r', encoding='utf-8') as handle:
        return json.load(handle)


def _write_json_atomic(payload, outpath):
    with atomic_output_path(outpath=outpath, suffix='.json') as tmp_path:
        with open(tmp_path, 'w', encoding='utf-8') as handle:
            json.dump(payload, handle, indent=2, sort_keys=True, ensure_ascii=True)
            handle.write('\n')


def _write_tsv_atomic(df, outpath):
    sanitized = sanitize_dataframe_for_tsv(df)
    with atomic_output_path(outpath=outpath, suffix='.tsv') as tmp_path:
        sanitized.to_csv(tmp_path, sep='\t', index=False)


def _validate_metadata_dataframe(df, context='metadata'):
    if not isinstance(df, pandas.DataFrame):
        raise TypeError('Expected pandas.DataFrame for {}.'.format(context))
    required_columns = ['run', 'scientific_name']
    missing_columns = [col for col in required_columns if col not in df.columns]
    if len(missing_columns) > 0:
        raise ValueError(
            'Missing required metadata column(s) for {}: {}'.format(
                context,
                ', '.join(missing_columns),
            )
        )
    validated = df.copy()
    validated['run'] = validated['run'].fillna('').astype(str).str.strip()
    validated['scientific_name'] = validated['scientific_name'].fillna('').astype(str).str.strip()
    missing_run_count = int((validated['run'] == '').sum())
    if missing_run_count > 0:
        raise ValueError('Missing run ID in metadata for {} row(s) [{}].'.format(missing_run_count, context))
    missing_species_runs = validated.loc[
        (validated['run'] != '') & (validated['scientific_name'] == ''),
        'run',
    ].tolist()
    if len(missing_species_runs) > 0:
        raise ValueError(
            'Missing scientific_name in metadata for run(s) [{}]: {}'.format(
                context,
                ', '.join(missing_species_runs),
            )
        )
    nonempty_runs = validated.loc[validated['run'] != '', 'run']
    duplicate_runs = nonempty_runs[nonempty_runs.duplicated()].unique().tolist()
    if len(duplicate_runs) > 0:
        raise ValueError(
            'Duplicate run ID in metadata for run(s) [{}]: {}'.format(
                context,
                ', '.join(duplicate_runs),
            )
        )
    return validated


def _drop_rows_missing_run(df, context='metadata'):
    if not isinstance(df, pandas.DataFrame):
        raise TypeError('Expected pandas.DataFrame for {}.'.format(context))
    if 'run' not in df.columns:
        return df.copy(), {
            'missing_run_drop_count': 0,
            'missing_run_drop_examples': [],
        }
    cleaned = df.copy()
    cleaned['run'] = cleaned['run'].fillna('').astype(str).str.strip()
    missing_mask = (cleaned['run'] == '')
    missing_count = int(missing_mask.sum())
    examples = []
    if missing_count > 0:
        example_columns = [
            col for col in ['scientific_name', 'biosample', 'experiment', 'sample_title']
            if col in cleaned.columns
        ]
        for _, row in cleaned.loc[missing_mask, example_columns].head(10).iterrows():
            example = {}
            for col in example_columns:
                value = row[col]
                if pandas.isna(value):
                    value = ''
                example[col] = str(value).strip()
            examples.append(example)
        sys.stderr.write(
            'Warning: Dropping {} metadata row(s) without run ID [{}].\n'.format(
                missing_count,
                context,
            )
        )
        cleaned = cleaned.loc[~missing_mask, :].reset_index(drop=True)
    return cleaned, {
        'missing_run_drop_count': missing_count,
        'missing_run_drop_examples': examples,
    }


def _write_metadata_tsv_with_validation(df, outpath, context, source_metadata=None):
    validated = _validate_metadata_dataframe(df=df, context=context)
    metadata = Metadata.from_DataFrame(validated)
    if source_metadata is not None:
        for attr_name in [
            'missing_run_drop_count',
            'missing_run_drop_examples',
            'sample_attribute_collision_count',
            'sample_attribute_collision_examples',
        ]:
            if hasattr(source_metadata, attr_name):
                setattr(metadata, attr_name, getattr(source_metadata, attr_name))
    sanitized = sanitize_dataframe_for_tsv(metadata.df)
    with atomic_output_path(outpath=outpath, suffix='.tsv') as tmp_path:
        sanitized.to_csv(tmp_path, sep='\t', index=False)
        reloaded = pandas.read_csv(
            tmp_path,
            sep='\t',
            dtype=str,
            keep_default_na=False,
            low_memory=False,
        )
        if list(reloaded.columns) != list(sanitized.columns):
            raise ValueError(
                'Read-back metadata columns changed after TSV write [{}].'.format(context)
            )
        if reloaded.shape[0] != sanitized.shape[0]:
            raise ValueError(
                'Read-back metadata row count changed after TSV write [{}].'.format(context)
            )
        _validate_metadata_dataframe(df=reloaded, context='{} read-back'.format(context))
    return metadata


def _normalized_text_series(df, column_name):
    if column_name not in df.columns:
        return pandas.Series([''] * len(df), index=df.index, dtype='object')
    return df[column_name].fillna('').astype(str).str.strip()


def _build_metadata_summary(df, top_n=10):
    rows = []

    def add_metric(name, value):
        rows.append(
            {
                'section': 'metric',
                'name': name,
                'value': str(value),
                'count': int(value) if isinstance(value, (int, bool)) else '',
                'rank': '',
            }
        )

    row_count = int(df.shape[0])
    runs = _normalized_text_series(df, 'run')
    species = _normalized_text_series(df, 'scientific_name')
    bioprojects = _normalized_text_series(df, 'bioproject')
    biosamples = _normalized_text_series(df, 'biosample')
    add_metric('row_count', row_count)
    add_metric('unique_run_count', int((runs != '').sum()))
    add_metric('unique_scientific_name_count', int(species[species != ''].nunique()))
    add_metric('unique_bioproject_count', int(bioprojects[bioprojects != ''].nunique()))
    add_metric('unique_biosample_count', int(biosamples[biosamples != ''].nunique()))
    add_metric('missing_run_count', int((runs == '').sum()))
    add_metric('missing_scientific_name_count', int(((runs != '') & (species == '')).sum()))
    for column_name in ['tissue', 'sample_group', 'source_name', 'sample_title']:
        values = _normalized_text_series(df, column_name)
        add_metric('nonempty_{}_count'.format(column_name), int((values != '').sum()))
        counts = values[values != ''].value_counts().head(top_n)
        for rank, (value, count) in enumerate(counts.items(), start=1):
            rows.append(
                {
                    'section': 'top_value',
                    'name': column_name,
                    'value': value,
                    'count': int(count),
                    'rank': rank,
                }
            )
    return pandas.DataFrame.from_records(
        rows,
        columns=['section', 'name', 'value', 'count', 'rank'],
    )


def _build_query_info(
    args,
    metadata,
    metadata_path,
    summary_path,
    query_info_path,
    search_string=None,
    record_ids=None,
    species_name=None,
    query_label=None,
    mode='single',
    used_cached_metadata=False,
    previous_info=None,
    kind='query',
    extra=None,
):
    payload = {}
    if isinstance(previous_info, dict):
        payload.update(previous_info)
    missing_run_drop_count = payload.get('missing_run_drop_count')
    missing_run_drop_examples = payload.get('missing_run_drop_examples')
    metadata_missing_run_drop_count = getattr(metadata, 'missing_run_drop_count', None)
    metadata_missing_run_drop_examples = getattr(metadata, 'missing_run_drop_examples', None)
    if (not used_cached_metadata) or (missing_run_drop_count is None):
        missing_run_drop_count = metadata_missing_run_drop_count
    elif metadata_missing_run_drop_count not in [None, 0]:
        missing_run_drop_count = metadata_missing_run_drop_count
    if (not used_cached_metadata) or (missing_run_drop_examples is None):
        missing_run_drop_examples = metadata_missing_run_drop_examples
    elif metadata_missing_run_drop_examples not in [None, []]:
        missing_run_drop_examples = metadata_missing_run_drop_examples
    collision_count = payload.get('sample_attribute_collision_count')
    collision_examples = payload.get('sample_attribute_collision_examples')
    metadata_collision_count = getattr(metadata, 'sample_attribute_collision_count', None)
    metadata_collision_examples = getattr(metadata, 'sample_attribute_collision_examples', None)
    if (not used_cached_metadata) or (collision_count is None):
        collision_count = metadata_collision_count
    elif metadata_collision_count not in [None, 0]:
        collision_count = metadata_collision_count
    if (not used_cached_metadata) or (collision_examples is None):
        collision_examples = metadata_collision_examples
    elif metadata_collision_examples not in [None, []]:
        collision_examples = metadata_collision_examples
    record_id_count = payload.get('record_id_count')
    if record_ids is not None:
        record_id_count = len(record_ids)
    payload.update(
        {
            'amalgkit_version': __version__,
            'created_at': datetime.datetime.now().isoformat(timespec='seconds'),
            'kind': kind,
            'mode': mode,
            'species': species_name,
            'query_label': query_label,
            'search_string': search_string,
            'record_id_count': record_id_count,
            'metadata_row_count': int(metadata.df.shape[0]),
            'resolve_names': bool(getattr(args, 'resolve_names', False)),
            'redo': bool(getattr(args, 'redo', False)),
            'used_cached_metadata': bool(used_cached_metadata),
            'metadata_path': _realpath_or_empty(metadata_path),
            'summary_path': _realpath_or_empty(summary_path),
            'query_info_path': _realpath_or_empty(query_info_path),
            'missing_run_drop_count': 0 if missing_run_drop_count in [None, ''] else int(missing_run_drop_count),
            'missing_run_drop_examples': [] if missing_run_drop_examples in [None, ''] else missing_run_drop_examples,
            'sample_attribute_collision_count': collision_count,
            'sample_attribute_collision_examples': collision_examples,
        }
    )
    if extra is not None:
        payload.update(extra)
    return payload


def _metadata_output_paths(out_dir):
    out_dir = _ensure_directory(out_dir, 'Output path')
    metadata_dir = _ensure_directory(os.path.join(out_dir, 'metadata'), 'Metadata path')
    metadata_path = os.path.join(metadata_dir, 'metadata.tsv')
    if os.path.exists(metadata_path) and (not os.path.isfile(metadata_path)):
        raise NotADirectoryError('Output path exists but is not a file: {}'.format(metadata_path))
    return {
        'out_dir': out_dir,
        'metadata_dir': metadata_dir,
        'metadata_path': metadata_path,
        'summary_path': os.path.join(metadata_dir, 'metadata.summary.tsv'),
        'query_info_path': os.path.join(metadata_dir, 'metadata.query_info.json'),
    }


def _read_metadata_from_path(metadata_path, context):
    df = pandas.read_csv(
        metadata_path,
        sep='\t',
        dtype=str,
        keep_default_na=False,
        low_memory=False,
    )
    validated = _validate_metadata_dataframe(df=df, context=context)
    return Metadata.from_DataFrame(validated)


def _prepare_single_metadata(metadata, args):
    if 'tissue' not in metadata.df.columns:
        metadata.df['tissue'] = ''
    metadata.df['tissue'] = metadata.df['tissue'].fillna('').astype(str)
    metadata.df.loc[:, 'sample_group'] = metadata.df.loc[:, 'tissue'].replace('nan', '').str.lower()
    if 'taxid' not in metadata.df.columns:
        metadata.df['taxid'] = pandas.Series([pandas.NA] * metadata.df.shape[0], dtype='Int64')
    else:
        metadata.df['taxid'] = pandas.to_numeric(metadata.df['taxid'], errors='coerce').astype('Int64')
    metadata.add_standard_rank_taxids(args=args)
    if getattr(args, 'resolve_names', False):
        metadata.resolve_scientific_names(args=args)
    metadata.reorder(omit_misc=False)
    metadata.df, missing_run_stats = _drop_rows_missing_run(
        df=metadata.df,
        context='metadata generation',
    )
    metadata.missing_run_drop_count = missing_run_stats['missing_run_drop_count']
    metadata.missing_run_drop_examples = missing_run_stats['missing_run_drop_examples']
    metadata.df = _validate_metadata_dataframe(df=metadata.df, context='metadata generation')
    return metadata


def _build_zero_record_guidance(search_term):
    diagnostics = inspect_accession_search_mismatches(search_term)
    if len(diagnostics) == 0:
        return ''
    lines = [
        'Diagnostic: the full --search_string returned zero SRA records, '
        'but accession-specific lookups found the following existing record(s):\n'
    ]
    for diagnostic in diagnostics:
        fields = []
        for key in ['scientific_name', 'platform', 'library_strategy', 'library_source', 'library_selection']:
            value = diagnostic.get(key, '')
            if value not in [None, '']:
                fields.append('{}={}'.format(key, value))
        line = '  {} exists in SRA'.format(diagnostic.get('accession', 'accession'))
        if len(fields) > 0:
            line += ' ({})'.format('; '.join(fields))
        line += '.'
        lines.append(line + '\n')
    lines.append(
        'Diagnostic: this usually means the accession does not satisfy additional query clauses such as '
        '[Platform] or [Strategy]. Consider relaxing those filters if this run is intended.\n'
    )
    return ''.join(lines)


def _run_single_query(
    args,
    out_dir=None,
    search_string=None,
    species_name=None,
    query_label=None,
    mode='single',
    allow_cached=False,
):
    paths = _metadata_output_paths(out_dir if out_dir is not None else args.out_dir)
    metadata_path = paths['metadata_path']
    if os.path.exists(metadata_path) and (not getattr(args, 'redo', False)) and (not allow_cached):
        raise AmalgkitExit(
            'Output file already exists (set --redo yes to overwrite): {}'.format(metadata_path),
            exit_code=0,
            use_stderr=False,
        )

    search_term = _normalize_search_string(
        search_string if search_string is not None else getattr(args, 'search_string', None)
    )
    if search_term == '':
        raise ValueError('--search_string is required.')

    used_cached_metadata = bool(os.path.exists(metadata_path) and (not getattr(args, 'redo', False)) and allow_cached)
    previous_info = _read_json_if_exists(paths['query_info_path'])

    if used_cached_metadata:
        metadata = _read_metadata_from_path(
            metadata_path=metadata_path,
            context='cached metadata {}'.format(species_name if species_name else search_term),
        )
        record_ids = None
    else:
        print('Entrez search term:', search_term)
        record_ids = search_sra_record_ids(search_term)
        metadata = Metadata.from_xml_roots(iter_sra_xml_chunks(record_ids=record_ids))
        metadata = _prepare_single_metadata(metadata=metadata, args=args)
        metadata = _write_metadata_tsv_with_validation(
            df=metadata.df,
            outpath=metadata_path,
            context='metadata {}'.format(species_name if species_name else search_term),
            source_metadata=metadata,
        )

    summary_df = _build_metadata_summary(metadata.df)
    _write_tsv_atomic(summary_df, paths['summary_path'])
    query_info = _build_query_info(
        args=args,
        metadata=metadata,
        metadata_path=metadata_path,
        summary_path=paths['summary_path'],
        query_info_path=paths['query_info_path'],
        search_string=search_term,
        record_ids=record_ids,
        species_name=species_name,
        query_label=query_label,
        mode=mode,
        used_cached_metadata=used_cached_metadata,
        previous_info=previous_info,
    )
    _write_json_atomic(query_info, paths['query_info_path'])
    return {
        'metadata': metadata,
        'paths': paths,
        'query_info': query_info,
        'search_string': search_term,
        'query_label': query_label,
    }


def _sanitize_species_token(species_name):
    token = re.sub(r'[^0-9A-Za-z._-]+', '_', species_name.strip())
    token = re.sub(r'_+', '_', token).strip('_.')
    if token == '':
        raise ValueError('Could not derive a filesystem-safe token from species name: {}'.format(species_name))
    return token


def _load_species_list(species_tsv, limit=None):
    species_df = pandas.read_csv(
        species_tsv,
        sep='\t',
        dtype=str,
        keep_default_na=False,
        low_memory=False,
    )
    if 'scientific_name' not in species_df.columns:
        raise ValueError('species_tsv must contain a scientific_name column.')
    species = []
    seen = set()
    for raw_name in species_df['scientific_name'].tolist():
        normalized = str(raw_name).strip()
        if (normalized == '') or (normalized in seen):
            continue
        seen.add(normalized)
        species.append(normalized)
        if (limit is not None) and (len(species) >= limit):
            break
    if len(species) == 0:
        raise ValueError('No scientific_name entries were found in species_tsv.')
    return species


def _split_title_terms(raw_text):
    return [term.strip() for term in re.split(r'[;,]', str(raw_text)) if term.strip() != '']


def _load_organ_terms(organ_terms_tsv):
    if organ_terms_tsv in [None, '']:
        return {group: list(terms) for group, terms in DEFAULT_ORGAN_TERMS.items()}
    organ_terms_df = pandas.read_csv(
        organ_terms_tsv,
        sep='\t',
        dtype=str,
        keep_default_na=False,
        low_memory=False,
    )
    required_columns = ['sample_group', 'title_terms']
    missing_columns = [col for col in required_columns if col not in organ_terms_df.columns]
    if len(missing_columns) > 0:
        raise ValueError(
            'organ_terms_tsv is missing required column(s): {}'.format(', '.join(missing_columns))
        )
    organ_terms = {}
    for _, row in organ_terms_df.iterrows():
        sample_group = str(row['sample_group']).strip()
        title_terms = _split_title_terms(row['title_terms'])
        if (sample_group == '') or (len(title_terms) == 0):
            continue
        organ_terms[sample_group] = title_terms
    missing_groups = [group for group in ['flower', 'leaf', 'root'] if group not in organ_terms]
    if len(missing_groups) > 0:
        raise ValueError(
            'organ_terms_tsv is missing sample_group(s): {}'.format(', '.join(missing_groups))
        )
    return organ_terms


def _build_title_union(terms):
    return ' OR '.join(['"{}"[Title]'.format(term.replace('"', '')) for term in terms])


def _build_species_queries(species_name, mode, title_terms, organ_terms):
    cleaned_species = species_name.replace('"', '')
    base_query = BASE_QUERY_TEMPLATE.format(species=cleaned_species)
    if mode == 'base':
        return [('base', base_query)]
    if organ_terms is not None:
        if mode == 'title_union':
            merged_terms = []
            for group in ['flower', 'leaf', 'root']:
                merged_terms.extend(organ_terms[group])
            return [('title_union', '({}) AND ({})'.format(base_query, _build_title_union(merged_terms)))]
        if mode == 'title_split':
            return [
                (group, '({}) AND ({})'.format(base_query, _build_title_union(organ_terms[group])))
                for group in ['flower', 'leaf', 'root']
            ]
    parsed_terms = _split_title_terms(title_terms)
    if len(parsed_terms) == 0:
        raise ValueError('--title_terms did not contain any usable terms.')
    if mode == 'title_union':
        return [('title_union', '({}) AND ({})'.format(base_query, _build_title_union(parsed_terms)))]
    if mode == 'title_split':
        return [
            (term, '({}) AND ({})'.format(base_query, _build_title_union([term])))
            for term in parsed_terms
        ]
    raise ValueError('Unknown metadata mode: {}'.format(mode))


def _resolve_metadata_species_jobs(args, task_count):
    species_jobs, _ = resolve_worker_allocation(
        requested_workers=getattr(args, 'internal_jobs', 'auto'),
        requested_threads=getattr(args, 'threads', 'auto'),
        internal_cpu_budget=getattr(args, 'internal_cpu_budget', 'auto'),
        worker_option_name='internal_jobs',
        context='metadata:',
    )
    if is_auto_parallel_option(getattr(args, 'internal_jobs', 'auto')):
        capped_jobs = min(species_jobs, METADATA_AUTO_MAX_SPECIES_JOBS)
        if capped_jobs < species_jobs:
            print(
                'metadata: capping auto --internal_jobs from {} to {} to avoid overloading Entrez.'.format(
                    species_jobs,
                    capped_jobs,
                ),
                flush=True,
            )
        species_jobs = capped_jobs
    if species_jobs > int(task_count):
        species_jobs = int(task_count)
    if species_jobs <= 0:
        species_jobs = 1
    return species_jobs


def _merge_metadata_results(query_results):
    if len(query_results) == 0:
        return Metadata()
    frames = [result['metadata'].df for result in query_results]
    merged_df = pandas.concat(frames, ignore_index=True, sort=False)
    if 'run' in merged_df.columns:
        merged_df['run'] = merged_df['run'].fillna('').astype(str).str.strip()
        merged_df = merged_df.drop_duplicates(subset=['run'], keep='first')
    merged_df = merged_df.reset_index(drop=True)
    return Metadata.from_DataFrame(merged_df)


def _write_species_merged_metadata(args, species_name, species_dir, species_token, query_results):
    merged_metadata = _merge_metadata_results(query_results)
    merged_path = os.path.join(species_dir, '{}.metadata.tsv'.format(species_token))
    merged_query_info_path = os.path.join(species_dir, '{}.query_info.json'.format(species_token))
    merged_summary_path = os.path.join(species_dir, '{}.summary.tsv'.format(species_token))
    if os.path.exists(merged_path) and (not os.path.isfile(merged_path)):
        raise NotADirectoryError('Output path exists but is not a file: {}'.format(merged_path))
    merged_metadata = _write_metadata_tsv_with_validation(
        df=merged_metadata.df,
        outpath=merged_path,
        context='merged metadata {}'.format(species_name),
        source_metadata=merged_metadata,
    )
    merged_summary = _build_metadata_summary(merged_metadata.df)
    _write_tsv_atomic(merged_summary, merged_summary_path)
    source_query_infos = [result['query_info'] for result in query_results]
    source_record_id_total = sum(
        [
            int(info['record_id_count'])
            for info in source_query_infos
            if info.get('record_id_count') is not None
        ]
    )
    source_missing_run_drop_total = sum(
        [
            int(info.get('missing_run_drop_count', 0))
            for info in source_query_infos
        ]
    )
    source_missing_run_drop_examples = []
    for info in source_query_infos:
        for example in info.get('missing_run_drop_examples', []):
            if len(source_missing_run_drop_examples) >= 10:
                break
            source_missing_run_drop_examples.append(example)
        if len(source_missing_run_drop_examples) >= 10:
            break
    merged_query_info = _build_query_info(
        args=args,
        metadata=merged_metadata,
        metadata_path=merged_path,
        summary_path=merged_summary_path,
        query_info_path=merged_query_info_path,
        search_string=None,
        record_ids=None,
        species_name=species_name,
        query_label=None,
        mode=getattr(args, 'mode', 'base'),
        used_cached_metadata=any(info.get('used_cached_metadata', False) for info in source_query_infos),
        previous_info=_read_json_if_exists(merged_query_info_path),
        kind='merged',
        extra={
            'source_query_labels': [result['query_label'] for result in query_results],
            'source_query_info_paths': [result['paths']['query_info_path'] for result in query_results],
            'source_metadata_paths': [result['paths']['metadata_path'] for result in query_results],
            'source_record_id_count_total': source_record_id_total,
            'missing_run_drop_count': source_missing_run_drop_total,
            'missing_run_drop_examples': source_missing_run_drop_examples,
        },
    )
    _write_json_atomic(merged_query_info, merged_query_info_path)
    print('  merged_rows={}'.format(merged_metadata.df.shape[0]), flush=True)


def _run_species_batch_task(
    task,
    args,
    species_count,
    mode,
    organ_terms,
    metadata_specieswise_dir,
):
    index, species_name = task
    print('[{}/{}] {}'.format(index, species_count, species_name), flush=True)
    species_token = _sanitize_species_token(species_name)
    species_dir = _ensure_directory(
        os.path.join(metadata_specieswise_dir, species_token),
        'species metadata path',
    )
    query_results = []
    for query_label, search_string in _build_species_queries(
        species_name=species_name,
        mode=mode,
        title_terms=getattr(args, 'title_terms', 'flower,leaf,root'),
        organ_terms=organ_terms,
    ):
        query_out_dir = os.path.join(species_dir, query_label)
        query_args = clone_namespace(
            args,
            out_dir=query_out_dir,
            search_string=search_string,
            species_tsv=None,
        )
        query_results.append(
            _run_single_query(
                args=query_args,
                out_dir=query_out_dir,
                search_string=search_string,
                species_name=species_name,
                query_label=query_label,
                mode=mode,
                allow_cached=True,
            )
        )
    if getattr(args, 'merge', True):
        _write_species_merged_metadata(
            args=args,
            species_name=species_name,
            species_dir=species_dir,
            species_token=species_token,
            query_results=query_results,
        )
    return species_name


def _run_species_batch(args):
    species_list = _load_species_list(
        species_tsv=getattr(args, 'species_tsv', None),
        limit=getattr(args, 'species_limit', None),
    )
    mode = getattr(args, 'mode', 'base')
    organ_terms = None
    if mode != 'base':
        organ_terms = _load_organ_terms(getattr(args, 'organ_terms_tsv', None))
    metadata_specieswise_dir = _ensure_directory(
        os.path.join(args.out_dir, 'metadata_specieswise'),
        'metadata_specieswise path',
    )
    task_items = list(enumerate(species_list, start=1))
    species_jobs = _resolve_metadata_species_jobs(args=args, task_count=len(task_items))
    if (species_jobs == 1) or (len(task_items) <= 1):
        for task in task_items:
            _run_species_batch_task(
                task=task,
                args=args,
                species_count=len(species_list),
                mode=mode,
                organ_terms=organ_terms,
                metadata_specieswise_dir=metadata_specieswise_dir,
            )
        return
    print(
        'metadata: running {:,} species queries with {:,} parallel job(s).'.format(
            len(task_items),
            species_jobs,
        ),
        flush=True,
    )
    _, failures = run_tasks_with_optional_threads(
        task_items=task_items,
        task_fn=lambda task: _run_species_batch_task(
            task=task,
            args=args,
            species_count=len(species_list),
            mode=mode,
            organ_terms=organ_terms,
            metadata_specieswise_dir=metadata_specieswise_dir,
        ),
        max_workers=species_jobs,
    )
    if failures:
        details = '; '.join(
            ['{}: {}'.format(task[1] if isinstance(task, tuple) and len(task) >= 2 else task, err) for task, err in failures]
        )
        raise RuntimeError(
            'metadata batch failed for {}/{} species. {}'.format(
                len(failures),
                len(task_items),
                details,
            )
        )


def metadata_main(args):
    Entrez.email = getattr(args, 'entrez_email', '')
    search_term = _normalize_search_string(getattr(args, 'search_string', None))
    species_tsv = getattr(args, 'species_tsv', None)
    if species_tsv not in [None, '']:
        if search_term != '':
            raise ValueError('Use either --search_string or --species_tsv, not both.')
        _run_species_batch(args)
        return
    result = _run_single_query(args=args, out_dir=args.out_dir, search_string=search_term, mode='single')
    if result['metadata'].df.shape[0] == 0:
        query_info = result.get('query_info', {})
        record_id_count = query_info.get('record_id_count')
        try:
            record_id_count = int(record_id_count)
        except (TypeError, ValueError):
            record_id_count = None
        if (record_id_count == 0) and (search_term != ''):
            sys.stderr.write(_build_zero_record_guidance(search_term))
        txt = (
            'No entry was found/survived in the metadata processing. '
            'Please reconsider the --search_string specification.\n'
        )
        sys.stderr.write(txt)
