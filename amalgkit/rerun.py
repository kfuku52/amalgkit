import json
import fcntl
import hashlib
import os
import shlex
import shutil
import tempfile
from contextlib import contextmanager
from types import SimpleNamespace

import pandas

from amalgkit.arg_utils import clone_namespace
from amalgkit.command_context import PerSpeciesTableContext
from amalgkit.exceptions import AmalgkitExit
from amalgkit.filter_utils import (
    copy_per_species_pdfs,
    load_merged_per_species_metadata,
    merge_metadata_by_run,
    save_exclusion_plot_pdf,
    staged_output_dir,
)
from amalgkit.finalize import _build_per_species_args, _copy_species_tables
from amalgkit.getfastq import getfastq_main
from amalgkit.merge import (
    collect_valid_run_ids,
    generate_merge_plot_pdfs,
    merge_fastp_stats_into_metadata,
    merge_species_quant_tables,
    scan_quant_abundance_paths,
)
from amalgkit.metadata_utils import Metadata, load_metadata, write_updated_metadata
from amalgkit.output_utils import atomic_output_path
from amalgkit.per_species_tables import generate_per_species_tables, resolve_per_species_input
from amalgkit.quant import (
    build_quant_tasks,
    check_quant_dependencies,
    prepare_quant_runtime_context,
    quant_main,
    resolve_quant_backends_for_tasks,
)
from amalgkit.runtime_utils import (
    resolve_species_token,
    safe_join_component,
    validate_safe_path_component,
    validate_unique_species_tokens,
)
from amalgkit.busco import (
    build_busco_species_token_map,
    collect_species as collect_busco_species,
    generate_busco_species_plot,
    load_metadata_if_needed_for_busco,
    process_species_busco,
    select_tool as select_busco_tool,
)


RERUN_CHECK_NAMES = ['getfastq', 'index', 'quant', 'merge', 'busco', 'finalize']
RERUN_MANIFEST_FILENAME = 'rerun_manifest.json'
MERGE_ROOT_OUTPUT_FILES = [
    'metadata.tsv',
    'merge_mapping_rate.pdf',
    'merge_total_spots.pdf',
    'merge_total_bases.pdf',
    'merge_library_layout.pdf',
    'merge_mean_expression_boxplot.pdf',
    'merge_fastp_duplication_rate_histogram.pdf',
    'merge_fastp_insert_size_peak_histogram.pdf',
    'merge_exclusion.pdf',
]
FINALIZE_ROOT_OUTPUT_FILES = [
    'metadata.tsv',
    'finalize_exclusion.pdf',
]
MERGE_REQUIRED_SPECIES_SUFFIXES = (
    'eff_length.tsv',
    'est_counts.tsv',
    'tpm.tsv',
)
FINALIZE_REQUIRED_SPECIES_SUFFIXES = (
    'metadata.tsv',
    'expression.tsv',
    'batch_effect_summary.tsv',
    'curation_final_summary.tsv',
)


def _parse_csv_option(value):
    if value in [None, '']:
        return []
    return [
        token.strip()
        for token in str(value).split(',')
        if token.strip() != ''
    ]


def _read_json(path):
    with open(path, 'r', encoding='utf-8') as handle:
        return json.load(handle)


def _write_json(path, payload):
    parent_dir = os.path.dirname(path)
    if parent_dir != '':
        os.makedirs(parent_dir, exist_ok=True)
    if os.path.lexists(path) and os.path.islink(path):
        raise ValueError('Refusing symbolic-link JSON output path: {}'.format(path))
    if os.path.exists(path) and (not os.path.isfile(path)):
        raise IsADirectoryError('JSON output path exists but is not a file: {}'.format(path))
    with atomic_output_path(path, suffix='.json') as temporary_path:
        with open(temporary_path, 'w', encoding='utf-8') as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)


def _resolve_report_path(args):
    if getattr(args, 'report', 'inferred') != 'inferred':
        report_path = os.path.realpath(args.report)
    else:
        report_path = os.path.realpath(os.path.join(args.out_dir, 'sanity', 'sanity_report.json'))
    if not os.path.exists(report_path):
        raise FileNotFoundError('Sanity report not found: {}'.format(report_path))
    if not os.path.isfile(report_path):
        raise IsADirectoryError('Sanity report path exists but is not a file: {}'.format(report_path))
    return report_path


def _resolve_execution_out_dir(args, report_payload):
    requested_absolute = os.path.abspath(args.out_dir)
    if os.path.lexists(requested_absolute) and os.path.islink(requested_absolute):
        raise ValueError(
            'Refusing symbolic-link rerun output root: {}'.format(requested_absolute)
        )
    requested = os.path.realpath(requested_absolute)
    report_out_dir = str(report_payload.get('out_dir', '')).strip()
    if report_out_dir == '':
        return requested
    report_out_absolute = os.path.abspath(report_out_dir)
    if os.path.lexists(report_out_absolute) and os.path.islink(report_out_absolute):
        raise ValueError(
            'Refusing symbolic-link rerun report output root: {}'.format(
                report_out_absolute
            )
        )
    report_out_dir = os.path.realpath(report_out_absolute)
    if requested == os.path.realpath('./'):
        return report_out_dir
    return requested


def _resolve_metadata_path(args, report_payload, out_dir):
    requested = str(getattr(args, 'metadata', 'report')).strip()
    if requested == 'report':
        report_metadata = str(report_payload.get('metadata_path', '')).strip()
        if report_metadata == '':
            raise ValueError('Sanity report did not include metadata_path. Specify --metadata explicitly.')
        return os.path.realpath(report_metadata)
    if requested == 'inferred':
        return os.path.realpath(os.path.join(out_dir, 'metadata', 'metadata.tsv'))
    return os.path.realpath(requested)


def _resolve_manifest_path(args, out_dir):
    requested = str(getattr(args, 'manifest', 'inferred')).strip()
    if requested.lower() == 'none':
        return None
    if requested == 'inferred':
        manifest_path = os.path.abspath(
            os.path.join(out_dir, 'sanity', RERUN_MANIFEST_FILENAME)
        )
    else:
        manifest_path = os.path.abspath(requested)
    manifest_parent = os.path.dirname(manifest_path)
    if os.path.lexists(manifest_parent) and os.path.islink(manifest_parent):
        raise ValueError(
            'Refusing symbolic-link rerun manifest directory: {}'.format(
                manifest_parent
            )
        )
    if os.path.lexists(manifest_path) and os.path.islink(manifest_path):
        raise ValueError(
            'Refusing symbolic-link rerun manifest path: {}'.format(manifest_path)
        )
    return os.path.join(
        os.path.realpath(manifest_parent),
        os.path.basename(manifest_path),
    )


def _load_metadata_from_path(metadata_path, out_dir):
    proxy_args = SimpleNamespace(metadata=metadata_path, out_dir=out_dir, batch=None)
    return load_metadata(proxy_args)


def _resolve_requested_checks(args, report_payload):
    requested = [value.lower() for value in _parse_csv_option(getattr(args, 'check', None))]
    if requested:
        if 'all' in requested:
            return list(RERUN_CHECK_NAMES)
        invalid = [name for name in requested if name not in RERUN_CHECK_NAMES]
        if invalid:
            raise ValueError(
                'Unknown rerun check(s): {}. Accepted values: {}.'.format(
                    ', '.join(invalid),
                    ', '.join(RERUN_CHECK_NAMES + ['all']),
                )
            )
        return list(dict.fromkeys(requested))
    summary = report_payload.get('summary', [])
    summary_checks = []
    for row in summary:
        check_name = str(row.get('check', '')).strip().lower()
        if check_name in RERUN_CHECK_NAMES:
            summary_checks.append(check_name)
    if summary_checks:
        return list(dict.fromkeys(summary_checks))
    return list(RERUN_CHECK_NAMES)


def _collect_issue_targets(report_payload, check_name, target_type, include_warnings=False):
    targets = []
    seen = set()
    severities = {'error', 'warning'} if include_warnings else {'error'}
    for issue in report_payload.get('issues', []):
        if str(issue.get('check', '')).strip().lower() != check_name:
            continue
        if str(issue.get('severity', '')).strip().lower() not in severities:
            continue
        if str(issue.get('target_type', '')).strip().lower() != target_type:
            continue
        target_id = str(issue.get('target_id', '')).strip()
        if (target_id == '') or (target_id in seen):
            continue
        seen.add(target_id)
        targets.append(target_id)
    return targets


def _has_global_issue(report_payload, check_name, include_warnings=False):
    severities = {'error', 'warning'} if include_warnings else {'error'}
    for issue in report_payload.get('issues', []):
        if str(issue.get('check', '')).strip().lower() != check_name:
            continue
        if str(issue.get('severity', '')).strip().lower() not in severities:
            continue
        if str(issue.get('target_type', '')).strip().lower() == 'global':
            return True
    return False


def _filter_runs_by_user_filters(metadata, run_ids, allowed_runs, allowed_species):
    if len(run_ids) == 0:
        return []
    run_series = metadata.df.loc[:, 'run'].fillna('').astype(str).str.strip() if 'run' in metadata.df.columns else pandas.Series('', index=metadata.df.index)
    species_series = metadata.df.loc[:, 'scientific_name'].fillna('').astype(str).str.strip() if 'scientific_name' in metadata.df.columns else pandas.Series('', index=metadata.df.index)
    run_to_species = {}
    for run_id, species in zip(run_series.tolist(), species_series.tolist()):
        if run_id not in run_to_species:
            run_to_species[run_id] = species
    filtered = []
    for run_id in run_ids:
        if run_id not in run_to_species:
            continue
        if allowed_runs and (run_id not in allowed_runs):
            continue
        if allowed_species and (run_to_species.get(run_id, '') not in allowed_species):
            continue
        filtered.append(run_id)
    return filtered


def _filter_species_by_user_filters(species_names, allowed_species, available_species=None):
    available_set = None if available_species is None else set(available_species)
    filtered = []
    for species in species_names:
        if (available_set is not None) and (species not in available_set):
            continue
        if allowed_species and (species not in allowed_species):
            continue
        filtered.append(species)
    return filtered


def _subset_metadata(metadata, run_ids=None, species_names=None):
    df = metadata.df.copy(deep=True)
    if species_names:
        if 'scientific_name' not in df.columns:
            raise ValueError('scientific_name column is required for species-based rerun.')
        species_series = df.loc[:, 'scientific_name'].fillna('').astype(str).str.strip()
        df = df.loc[species_series.isin(species_names), :].copy()
    if run_ids:
        if 'run' not in df.columns:
            raise ValueError('run column is required for run-based rerun.')
        run_series = df.loc[:, 'run'].fillna('').astype(str).str.strip()
        df = df.loc[run_series.isin(run_ids), :].copy()
    if df.shape[0] == 0:
        raise ValueError('No metadata rows remained after applying rerun target filters.')
    return Metadata.from_DataFrame(df.reset_index(drop=True))


def _write_temp_metadata(df):
    handle = tempfile.NamedTemporaryFile(prefix='amalgkit_rerun_', suffix='.tsv', delete=False)
    handle.close()
    df.to_csv(handle.name, sep='\t', index=False)
    return handle.name


def _remove_path_if_exists(path):
    if not os.path.lexists(path):
        return
    if os.path.islink(path) or os.path.isfile(path):
        os.remove(path)
        return
    shutil.rmtree(path)


def _resolve_rerun_redo(args):
    return bool(getattr(args, 'redo', True))


def _normalize_relative_paths(relative_paths):
    normalized = []
    seen = set()
    for rel_path in relative_paths:
        rel_path = str(rel_path).strip()
        if (rel_path == '') or (rel_path in seen):
            continue
        rel_path = validate_safe_path_component(rel_path, label='rerun relative path')
        seen.add(rel_path)
        normalized.append(rel_path)
    return normalized


def _preview_path_list(paths, limit=5):
    if len(paths) <= limit:
        return ', '.join(paths)
    return '{}, ... (+{} more)'.format(', '.join(paths[:limit]), len(paths) - limit)


def _preflight_rerun_overwrite(target_root, relative_paths, redo):
    if redo:
        return
    existing_paths = []
    for rel_path in _normalize_relative_paths(relative_paths):
        target_path = safe_join_component(target_root, rel_path, label='rerun relative path')
        if os.path.lexists(target_path):
            existing_paths.append(target_path)
    if len(existing_paths) == 0:
        return
    print('Existing rerun outputs detected: {}'.format(_preview_path_list(existing_paths)), flush=True)
    raise AmalgkitExit('--redo is set to "no". Exiting.', exit_code=0, use_stderr=False)


def _normalize_species_dir_name(species, explicit_token=None):
    return resolve_species_token(
        species,
        explicit_token,
        label='species_token',
    )


def _build_rerun_species_token_map(metadata):
    species_names = metadata.df['scientific_name'].fillna('').astype(str).str.strip()
    if 'species_token' in metadata.df.columns:
        explicit_tokens = metadata.df['species_token'].fillna('').astype(str).str.strip()
    else:
        explicit_tokens = pandas.Series('', index=metadata.df.index)
    validate_unique_species_tokens(
        list(zip(species_names.tolist(), explicit_tokens.tolist())),
        context='rerun output',
    )
    token_by_species = {}
    for species, explicit_token in zip(species_names.tolist(), explicit_tokens.tolist()):
        if species == '':
            continue
        token = _normalize_species_dir_name(species, explicit_token)
        previous_token = token_by_species.get(species)
        if previous_token is not None and previous_token != token:
            raise ValueError(
                'Scientific name has conflicting species_token values: {} ({}, {})'.format(
                    species,
                    previous_token,
                    token,
                )
            )
        token_by_species[species] = token
    return token_by_species


def _create_staging_root(target_root, prefix):
    absolute_target_root = os.path.abspath(target_root)
    if os.path.lexists(absolute_target_root) and os.path.islink(absolute_target_root):
        raise ValueError(
            'Refusing to stage through symbolic-link rerun output root: {}'.format(
                absolute_target_root
            )
        )
    parent_dir = os.path.dirname(os.path.realpath(absolute_target_root))
    if parent_dir != '':
        os.makedirs(parent_dir, exist_ok=True)
        return tempfile.mkdtemp(prefix=prefix, dir=parent_dir)
    return tempfile.mkdtemp(prefix=prefix)


BACKUP_TEMP_PREFIX = 'amalgkit_rerun_backup_'
CLEANUP_TEMP_PREFIX = 'amalgkit_rerun_cleanup_'
TRANSACTION_MANIFEST_FILENAME = '.amalgkit-rerun-transaction.json'
TRANSACTION_SCHEMA_VERSION = 1
TRANSACTION_LOCK_PREFIX = '.amalgkit-rerun-transaction-'


def _fsync_directory(path):
    """Durably flush directory-entry changes for a directory."""
    if path in [None, '']:
        return
    fd = os.open(path, os.O_RDONLY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def _fsync_file(path):
    if os.path.islink(path) or not os.path.isfile(path):
        return
    fd = os.open(path, os.O_RDONLY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def _fsync_tree(path):
    """Flush staged file contents and the directory entries that contain them."""
    if os.path.islink(path):
        raise ValueError('Refusing to fsync symbolic-link rerun path: {}'.format(path))
    if os.path.isfile(path):
        _fsync_file(path)
        return
    if not os.path.isdir(path):
        raise FileNotFoundError('Rerun staged path disappeared before commit: {}'.format(path))
    for current_root, dirnames, filenames in os.walk(path, topdown=False, followlinks=False):
        for filename in filenames:
            child_path = os.path.join(current_root, filename)
            if not os.path.islink(child_path):
                _fsync_file(child_path)
        for dirname in dirnames:
            child_path = os.path.join(current_root, dirname)
            if not os.path.islink(child_path):
                _fsync_directory(child_path)
        _fsync_directory(current_root)


def _ensure_directory_durable(path):
    path = os.path.abspath(path)
    missing_paths = []
    cursor = path
    while not os.path.exists(cursor):
        missing_paths.append(cursor)
        parent = os.path.dirname(cursor)
        if parent == cursor:
            break
        cursor = parent
    if os.path.lexists(path) and (os.path.islink(path) or not os.path.isdir(path)):
        raise NotADirectoryError('Rerun directory path is not a regular directory: {}'.format(path))
    os.makedirs(path, exist_ok=True)
    for created_path in reversed(missing_paths):
        _fsync_directory(created_path)
        _fsync_directory(os.path.dirname(created_path))


def _durable_rename(source_path, destination_path):
    source_parent = os.path.dirname(source_path)
    destination_parent = os.path.dirname(destination_path)
    os.rename(source_path, destination_path)
    _fsync_directory(source_parent)
    if destination_parent != source_parent:
        _fsync_directory(destination_parent)


def _durable_remove_path(path):
    if not os.path.lexists(path):
        return
    parent_dir = os.path.dirname(path)
    _remove_path_if_exists(path)
    _fsync_directory(parent_dir)


def _transaction_lock_path(target_root):
    target_root = os.path.realpath(target_root)
    target_digest = hashlib.sha256(target_root.encode('utf-8')).hexdigest()[:16]
    return os.path.join(
        os.path.dirname(target_root),
        '{}{}.lock'.format(TRANSACTION_LOCK_PREFIX, target_digest),
    )


@contextmanager
def _target_transaction_lock(target_root):
    lock_path = _transaction_lock_path(target_root)
    _ensure_directory_durable(os.path.dirname(lock_path))
    flags = os.O_RDWR | os.O_CREAT | getattr(os, 'O_NOFOLLOW', 0)
    fd = os.open(lock_path, flags, 0o600)
    try:
        fcntl.flock(fd, fcntl.LOCK_EX)
        yield
    finally:
        fcntl.flock(fd, fcntl.LOCK_UN)
        os.close(fd)


def _transaction_manifest_path(backup_root):
    return os.path.join(backup_root, TRANSACTION_MANIFEST_FILENAME)


def _write_transaction_manifest(backup_root, payload):
    manifest_path = _transaction_manifest_path(backup_root)
    temporary_path = manifest_path + '.tmp'
    with open(temporary_path, 'w', encoding='utf-8') as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write('\n')
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary_path, manifest_path)
    _fsync_directory(backup_root)


def _read_transaction_manifest(backup_root):
    manifest_path = _transaction_manifest_path(backup_root)
    if not os.path.isfile(manifest_path) or os.path.islink(manifest_path):
        raise RuntimeError(
            'Unjournaled rerun backup was preserved for manual recovery: {}'.format(backup_root)
        )
    try:
        with open(manifest_path, 'r', encoding='utf-8') as handle:
            payload = json.load(handle)
    except (OSError, ValueError) as exc:
        raise RuntimeError(
            'Unreadable rerun transaction was preserved for manual recovery: {}'.format(
                backup_root
            )
        ) from exc
    if payload.get('schema_version') != TRANSACTION_SCHEMA_VERSION:
        raise RuntimeError(
            'Unsupported rerun transaction was preserved for manual recovery: {}'.format(
                backup_root
            )
        )
    if payload.get('state') not in {
        'prepared',
        'installing',
        'rolling_back',
        'rolled_back',
        'committed',
    }:
        raise RuntimeError('Invalid rerun transaction state in: {}'.format(backup_root))
    target_root = payload.get('target_root')
    if not isinstance(target_root, str) or target_root == '':
        raise RuntimeError('Invalid rerun transaction target in: {}'.format(backup_root))
    entries = payload.get('paths')
    if not isinstance(entries, list):
        raise RuntimeError('Invalid rerun transaction path list in: {}'.format(backup_root))
    seen_paths = set()
    for entry in entries:
        if (
            not isinstance(entry, dict)
            or not isinstance(entry.get('path'), str)
            or not isinstance(entry.get('had_original'), bool)
        ):
            raise RuntimeError('Invalid rerun transaction path entry in: {}'.format(backup_root))
        rel_path = validate_safe_path_component(entry['path'], label='rerun transaction path')
        if rel_path in seen_paths:
            raise RuntimeError('Duplicate rerun transaction path in: {}'.format(backup_root))
        seen_paths.add(rel_path)
    return payload


def _cleanup_directory_path(backup_root):
    backup_name = os.path.basename(backup_root)
    suffix = backup_name[len(BACKUP_TEMP_PREFIX):]
    return os.path.join(os.path.dirname(backup_root), CLEANUP_TEMP_PREFIX + suffix)


def _retire_transaction_dir(backup_root):
    """Move a resolved transaction out of the recovery namespace before deletion."""
    if not os.path.lexists(backup_root):
        return
    cleanup_root = _cleanup_directory_path(backup_root)
    _durable_rename(backup_root, cleanup_root)
    try:
        shutil.rmtree(cleanup_root)
    except OSError:
        return
    _fsync_directory(os.path.dirname(cleanup_root))


def _sweep_resolved_cleanup_dirs(parent_dir):
    if parent_dir in [None, ''] or not os.path.isdir(parent_dir):
        return
    for name in os.listdir(parent_dir):
        if not name.startswith(CLEANUP_TEMP_PREFIX):
            continue
        candidate = os.path.join(parent_dir, name)
        if os.path.islink(candidate):
            os.remove(candidate)
            continue
        if os.path.isdir(candidate):
            shutil.rmtree(candidate, ignore_errors=True)
    _fsync_directory(parent_dir)


def _recover_interrupted_transaction(backup_root, expected_target_root):
    payload = _read_transaction_manifest(backup_root)
    target_root = os.path.realpath(payload['target_root'])
    if target_root != os.path.realpath(expected_target_root):
        return False
    state = payload['state']
    if state in {'committed', 'rolled_back'}:
        _retire_transaction_dir(backup_root)
        return True

    entries = payload['paths']
    for entry in entries:
        if not entry['had_original']:
            continue
        rel_path = entry['path']
        target_path = safe_join_component(target_root, rel_path, label='rerun transaction path')
        backup_path = safe_join_component(backup_root, rel_path, label='rerun transaction path')
        if state == 'installing' and not os.path.lexists(backup_path):
            raise RuntimeError(
                'Incomplete rerun transaction is missing its required backup: {}'.format(
                    backup_path
                )
            )
        if state in {'prepared', 'rolling_back'} and not (
            os.path.lexists(backup_path) or os.path.lexists(target_path)
        ):
            raise RuntimeError(
                'Incomplete rerun transaction cannot restore missing output: {}'.format(
                    target_path
                )
            )

    payload['state'] = 'rolling_back'
    _write_transaction_manifest(backup_root, payload)
    for entry in entries:
        rel_path = entry['path']
        target_path = safe_join_component(target_root, rel_path, label='rerun transaction path')
        backup_path = safe_join_component(backup_root, rel_path, label='rerun transaction path')
        if entry['had_original']:
            if not os.path.lexists(backup_path):
                continue
            _durable_remove_path(target_path)
            _ensure_directory_durable(os.path.dirname(target_path))
            _durable_rename(backup_path, target_path)
        elif os.path.lexists(target_path):
            _durable_remove_path(target_path)
    _fsync_directory(target_root)
    payload['state'] = 'rolled_back'
    _write_transaction_manifest(backup_root, payload)
    _retire_transaction_dir(backup_root)
    return True


def _recover_interrupted_transactions(parent_dir, target_root):
    _sweep_resolved_cleanup_dirs(parent_dir)
    if parent_dir in [None, ''] or not os.path.isdir(parent_dir):
        return
    for name in sorted(os.listdir(parent_dir)):
        if not name.startswith(BACKUP_TEMP_PREFIX):
            continue
        backup_root = os.path.join(parent_dir, name)
        if os.path.islink(backup_root) or not os.path.isdir(backup_root):
            raise RuntimeError(
                'Invalid rerun backup was preserved for manual recovery: {}'.format(backup_root)
            )
        _recover_interrupted_transaction(backup_root, target_root)


def _recover_target_transactions(target_root):
    absolute_target_root = os.path.abspath(target_root)
    if os.path.lexists(absolute_target_root) and os.path.islink(absolute_target_root):
        raise ValueError(
            'Refusing to recover through symbolic-link rerun output root: {}'.format(
                absolute_target_root
            )
        )
    target_root = os.path.realpath(absolute_target_root)
    with _target_transaction_lock(target_root):
        _recover_interrupted_transactions(os.path.dirname(target_root), target_root)


def _commit_staged_paths(target_root, staged_root, relative_paths):
    absolute_target_root = os.path.abspath(target_root)
    absolute_staged_root = os.path.abspath(staged_root)
    if os.path.lexists(absolute_target_root) and os.path.islink(absolute_target_root):
        raise ValueError(
            'Refusing to commit through symbolic-link rerun output root: {}'.format(
                absolute_target_root
            )
        )
    if os.path.lexists(absolute_staged_root) and os.path.islink(absolute_staged_root):
        raise ValueError(
            'Refusing symbolic-link rerun staging root: {}'.format(absolute_staged_root)
        )
    target_root = os.path.realpath(absolute_target_root)
    staged_root = os.path.realpath(absolute_staged_root)
    with _target_transaction_lock(target_root):
        _commit_staged_paths_locked(target_root, staged_root, relative_paths)


def _commit_staged_paths_locked(target_root, staged_root, relative_paths):
    parent_dir = os.path.dirname(target_root)
    if parent_dir != '':
        _ensure_directory_durable(parent_dir)
    _ensure_directory_durable(target_root)
    _recover_interrupted_transactions(parent_dir, target_root)
    requested_paths = _normalize_relative_paths(relative_paths)
    replace_paths = [
        rel_path
        for rel_path in requested_paths
        if os.path.lexists(
            safe_join_component(staged_root, rel_path, label='rerun relative path')
        )
    ]
    if len(replace_paths) == 0:
        return
    for rel_path in replace_paths:
        staged_path = safe_join_component(staged_root, rel_path, label='rerun relative path')
        _fsync_tree(staged_path)
    backup_root = tempfile.mkdtemp(
        prefix=BACKUP_TEMP_PREFIX,
        dir=parent_dir if parent_dir != '' else None,
    )
    _fsync_directory(parent_dir)
    transaction = {
        'schema_version': TRANSACTION_SCHEMA_VERSION,
        'state': 'prepared',
        'target_root': target_root,
        'paths': [
            {
                'path': rel_path,
                'had_original': os.path.lexists(
                    safe_join_component(target_root, rel_path, label='rerun relative path')
                ),
            }
            for rel_path in replace_paths
        ],
    }
    try:
        _write_transaction_manifest(backup_root, transaction)
        for entry in transaction['paths']:
            if not entry['had_original']:
                continue
            rel_path = entry['path']
            target_path = safe_join_component(target_root, rel_path, label='rerun relative path')
            backup_path = safe_join_component(backup_root, rel_path, label='rerun relative path')
            _ensure_directory_durable(os.path.dirname(backup_path))
            _durable_rename(target_path, backup_path)
        transaction['state'] = 'installing'
        _write_transaction_manifest(backup_root, transaction)
        for rel_path in replace_paths:
            staged_path = safe_join_component(staged_root, rel_path, label='rerun relative path')
            target_path = safe_join_component(target_root, rel_path, label='rerun relative path')
            _ensure_directory_durable(os.path.dirname(target_path))
            _durable_rename(staged_path, target_path)
        _fsync_directory(target_root)
        _fsync_directory(backup_root)
        _fsync_directory(parent_dir)
        transaction['state'] = 'committed'
        _write_transaction_manifest(backup_root, transaction)
    except Exception:
        if os.path.isfile(_transaction_manifest_path(backup_root)):
            _recover_interrupted_transaction(backup_root, target_root)
        else:
            shutil.rmtree(backup_root, ignore_errors=True)
            _fsync_directory(parent_dir)
        raise
    _retire_transaction_dir(backup_root)


def _validate_staged_species_outputs(staged_root, species_tokens, required_suffixes, label):
    for token in species_tokens:
        species_dir = safe_join_component(
            staged_root,
            token,
            label='{} species token'.format(label),
        )
        if os.path.islink(species_dir) or not os.path.isdir(species_dir):
            raise RuntimeError(
                '{} rerun did not generate the required species directory: {}'.format(
                    label,
                    species_dir,
                )
            )
        for suffix in required_suffixes:
            required_path = os.path.join(
                species_dir,
                '{}_{}'.format(token, suffix),
            )
            if (
                os.path.islink(required_path)
                or not os.path.isfile(required_path)
                or os.path.getsize(required_path) <= 0
            ):
                raise RuntimeError(
                    '{} rerun did not generate required output: {}'.format(
                        label,
                        required_path,
                    )
                )


def _list_unique_species(metadata):
    if 'scientific_name' not in metadata.df.columns:
        return []
    species = (
        metadata.df.loc[:, 'scientific_name']
        .fillna('')
        .astype(str)
        .str.strip()
    )
    return [value for value in dict.fromkeys(species.tolist()) if value != '']


def _list_selected_species(metadata):
    if 'scientific_name' not in metadata.df.columns:
        return []
    species_series = metadata.df.loc[:, 'scientific_name'].fillna('').astype(str).str.strip()
    if 'exclusion' not in metadata.df.columns:
        return [value for value in dict.fromkeys(species_series.tolist()) if value != '']
    exclusion_series = metadata.df.loc[:, 'exclusion'].fillna('').astype(str).str.strip().str.lower()
    selected = species_series.loc[exclusion_series.eq('no')]
    return [value for value in dict.fromkeys(selected.tolist()) if value != '']


def _copy_existing_merge_species_dirs(
    source_merge_dir,
    stage_payload_dir,
    species_names,
    skip_species=None,
    species_token_by_name=None,
):
    if not os.path.isdir(source_merge_dir):
        return
    species_token_by_name = species_token_by_name or {}

    def resolve_token(species):
        token = species_token_by_name.get(species)
        if token is None:
            token = _normalize_species_dir_name(species)
        return token

    skip_tokens = {
        resolve_token(species)
        for species in list(skip_species or [])
    }
    species_tokens = _normalize_relative_paths([
        resolve_token(species)
        for species in list(species_names or [])
    ])
    for token in species_tokens:
        if token in skip_tokens:
            continue
        source_dir = os.path.join(source_merge_dir, token)
        if not os.path.lexists(source_dir):
            continue
        if not os.path.isdir(source_dir):
            raise NotADirectoryError('Merge species path exists but is not a directory: {}'.format(source_dir))
        destination_dir = os.path.join(stage_payload_dir, token)
        if os.path.lexists(destination_dir):
            continue
        shutil.copytree(source_dir, destination_dir)


def _load_existing_metadata_table(path):
    if not os.path.isfile(path):
        return None
    return pandas.read_csv(path, sep='\t', low_memory=False)


def _print_targets(label, targets):
    if len(targets) == 0:
        print('{}: none'.format(label), flush=True)
        return
    print('{} ({}): {}'.format(label, len(targets), ', '.join(targets)), flush=True)


def _resolve_global_species_targets(check_name, metadata, allowed_species):
    if check_name == 'merge':
        return _filter_species_by_user_filters(_list_unique_species(metadata), allowed_species)
    if check_name == 'finalize':
        return _filter_species_by_user_filters(_list_selected_species(metadata), allowed_species)
    return []


def _build_rerun_plan(report_payload, metadata, requested_checks, allowed_runs, allowed_species, include_warnings):
    plan = []
    for check_name in requested_checks:
        item = {
            'check': check_name,
            'target_runs': [],
            'target_species': [],
            'has_global_issue': False,
            'rebuild_summary_only': False,
            'selection_basis': '',
            'will_execute': False,
            'skip_reason': '',
        }
        if check_name in ['getfastq', 'quant']:
            run_targets = _collect_issue_targets(
                report_payload=report_payload,
                check_name=check_name,
                target_type='run',
                include_warnings=include_warnings,
            )
            run_targets = _filter_runs_by_user_filters(metadata, run_targets, allowed_runs, allowed_species)
            item['target_runs'] = run_targets
            item['selection_basis'] = 'run_issues'
            item['will_execute'] = len(run_targets) > 0
            if not item['will_execute']:
                item['skip_reason'] = 'no matching targets'
            plan.append(item)
            continue

        species_targets = _collect_issue_targets(
            report_payload=report_payload,
            check_name=check_name,
            target_type='species',
            include_warnings=include_warnings,
        )
        species_targets = _filter_species_by_user_filters(
            species_targets,
            allowed_species,
            available_species=_list_unique_species(metadata),
        )
        has_global_issue = _has_global_issue(
            report_payload=report_payload,
            check_name=check_name,
            include_warnings=include_warnings,
        )
        item['has_global_issue'] = has_global_issue
        item['target_species'] = species_targets
        if check_name == 'busco' and has_global_issue and (len(species_targets) == 0):
            item['rebuild_summary_only'] = True
            item['selection_basis'] = 'global_issue_summary_rebuild'
            item['will_execute'] = True
        elif has_global_issue and check_name in {'merge', 'finalize'}:
            item['target_species'] = _resolve_global_species_targets(check_name, metadata, allowed_species)
            item['selection_basis'] = 'global_issue_all_species'
            item['will_execute'] = len(item['target_species']) > 0
            if not item['will_execute']:
                item['skip_reason'] = 'no matching targets'
        elif has_global_issue and (len(species_targets) == 0):
            item['target_species'] = _resolve_global_species_targets(check_name, metadata, allowed_species)
            item['selection_basis'] = 'global_issue_all_species'
            item['will_execute'] = len(item['target_species']) > 0
            if not item['will_execute']:
                item['skip_reason'] = 'no matching targets'
        else:
            item['selection_basis'] = 'species_issues'
            item['will_execute'] = len(item['target_species']) > 0
            if not item['will_execute']:
                item['skip_reason'] = 'no matching targets'
        plan.append(item)
    return plan


RERUN_MANIFEST_SCHEMA_VERSION = 1


def _now_iso_for_runtime_log():
    """Wall-clock run time, used only for the separate non-reproducible runtime log."""
    import time
    return time.strftime('%Y-%m-%dT%H:%M:%S', time.localtime())


def _write_rerun_runtime_log(out_dir, started_at_iso):
    """Record the non-reproducible wall-clock run time outside the compared manifest."""
    manifest_dir = os.path.join(out_dir, 'sanity')
    os.makedirs(manifest_dir, exist_ok=True)
    runtime_log_path = os.path.join(manifest_dir, 'rerun_runtime.log')
    with open(runtime_log_path, 'w', encoding='utf-8') as handle:
        handle.write('amalgkit rerun started at {}\n'.format(started_at_iso))
    return runtime_log_path


def _build_rerun_manifest(report_path, out_dir, metadata_path, requested_checks, allowed_runs, allowed_species, include_warnings, dry_run, plan):
    return {
        'schema_version': RERUN_MANIFEST_SCHEMA_VERSION,
        'report_path': report_path,
        'out_dir': out_dir,
        'metadata_path': metadata_path,
        'requested_checks': requested_checks,
        'filters': {
            'run': sorted(allowed_runs),
            'species': sorted(allowed_species),
            'include_warnings': include_warnings,
        },
        'dry_run': dry_run,
        'checks': plan,
    }


def rerun_getfastq_check(args, metadata, target_runs, dry_run=False):
    print('rerun getfastq', flush=True)
    _print_targets('target_runs', target_runs)
    if len(target_runs) == 0:
        return
    if dry_run:
        return
    subset_metadata = _subset_metadata(metadata, run_ids=target_runs)
    temp_metadata = _write_temp_metadata(subset_metadata.df)
    try:
        redo = _resolve_rerun_redo(args)
        runtime_args = clone_namespace(
            args,
            metadata=temp_metadata,
            redo=redo,
            batch=None,
            id=None,
            id_list=None,
        )
        getfastq_main(runtime_args)
    finally:
        if os.path.exists(temp_metadata):
            os.remove(temp_metadata)


def rerun_quant_check(args, metadata, target_runs, dry_run=False):
    print('rerun quant', flush=True)
    _print_targets('target_runs', target_runs)
    if len(target_runs) == 0:
        return
    if dry_run:
        return
    subset_metadata = _subset_metadata(metadata, run_ids=target_runs)
    temp_metadata = _write_temp_metadata(subset_metadata.df)
    try:
        redo = _resolve_rerun_redo(args)
        runtime_args = clone_namespace(
            args,
            metadata=temp_metadata,
            redo=redo,
            batch=None,
        )
        quant_main(runtime_args)
    finally:
        if os.path.exists(temp_metadata):
            os.remove(temp_metadata)


def rerun_index_check(args, metadata, target_species, dry_run=False):
    print('rerun index', flush=True)
    _print_targets('target_species', target_species)
    if len(target_species) == 0:
        return
    if dry_run:
        return
    subset_metadata = _subset_metadata(metadata, species_names=target_species)
    temp_metadata = _write_temp_metadata(subset_metadata.df)
    try:
        redo = _resolve_rerun_redo(args)
        runtime_args = clone_namespace(
            args,
            metadata=temp_metadata,
            build_index=True,
            clean_fastq=False,
            redo=redo,
            batch=None,
        )
        subset_loaded = load_metadata(runtime_args)
        tasks = build_quant_tasks(subset_loaded)
        backend_by_run, oarfish_seq_tech_by_run = resolve_quant_backends_for_tasks(runtime_args, subset_loaded, tasks)
        check_quant_dependencies(backend_by_run)
        prepare_quant_runtime_context(
            runtime_args,
            tasks,
            metadata=subset_loaded,
            backend_by_run=backend_by_run,
            oarfish_seq_tech_by_run=oarfish_seq_tech_by_run,
        )
    finally:
        if os.path.exists(temp_metadata):
            os.remove(temp_metadata)


def rerun_merge_check(args, metadata, target_species, force_all_species=False, dry_run=False):
    print('rerun merge', flush=True)
    if force_all_species and (len(target_species) == 0):
        target_species = _list_unique_species(metadata)
    _print_targets('target_species', target_species)
    if len(target_species) == 0:
        return
    if dry_run:
        return
    merge_dir = os.path.join(args.out_dir, 'merge')
    _recover_target_transactions(merge_dir)
    quant_dir = os.path.join(args.out_dir, 'quant')
    species_token_by_name = _build_rerun_species_token_map(metadata)
    target_species_dirs = [species_token_by_name[species] for species in target_species]
    redo = _resolve_rerun_redo(args)
    _preflight_rerun_overwrite(
        target_root=merge_dir,
        relative_paths=target_species_dirs + MERGE_ROOT_OUTPUT_FILES,
        redo=redo,
    )
    stage_root = _create_staging_root(merge_dir, prefix='amalgkit_rerun_merge_')
    stage_payload_dir = os.path.join(stage_root, 'payload')
    os.makedirs(stage_payload_dir, exist_ok=True)
    run_abundance_paths = scan_quant_abundance_paths(
        quant_dir=quant_dir,
        target_runs=set(collect_valid_run_ids(metadata.df.loc[:, 'run'].values)) if 'run' in metadata.df.columns else None,
    )
    try:
        for species in target_species:
            merge_species_quant_tables(
                sp=species,
                metadata=metadata,
                quant_dir=quant_dir,
                merge_dir=stage_payload_dir,
                run_abundance_paths=run_abundance_paths,
            )
        _copy_existing_merge_species_dirs(
            source_merge_dir=merge_dir,
            stage_payload_dir=stage_payload_dir,
            species_names=_list_unique_species(metadata),
            skip_species=target_species,
            species_token_by_name=species_token_by_name,
        )
        refreshed_metadata = merge_fastp_stats_into_metadata(
            Metadata.from_DataFrame(metadata.df.copy(deep=True)),
            args.out_dir,
            max_workers=getattr(args, 'threads', 'auto'),
        )
        merge_metadata_path = os.path.join(stage_payload_dir, 'metadata.tsv')
        write_updated_metadata(
            refreshed_metadata,
            merge_metadata_path,
            args,
            max_workers=getattr(args, 'threads', 'auto'),
        )
        generate_merge_plot_pdfs(
            merge_dir=stage_payload_dir,
            path_metadata_merge=merge_metadata_path,
        )
        _validate_staged_species_outputs(
            staged_root=stage_payload_dir,
            species_tokens=target_species_dirs,
            required_suffixes=MERGE_REQUIRED_SPECIES_SUFFIXES,
            label='merge',
        )
        _commit_staged_paths(
            target_root=merge_dir,
            staged_root=stage_payload_dir,
            relative_paths=target_species_dirs + MERGE_ROOT_OUTPUT_FILES,
        )
    finally:
        shutil.rmtree(stage_root, ignore_errors=True)


def rerun_busco_check(args, metadata, target_species, rebuild_summary_only=False, dry_run=False):
    print('rerun busco', flush=True)
    _print_targets('target_species', target_species)
    if (len(target_species) == 0) and (not rebuild_summary_only):
        return
    if dry_run:
        return
    busco_dir = os.path.join(args.out_dir, 'busco')
    all_species = _list_unique_species(metadata)
    species_token_by_name = _build_rerun_species_token_map(metadata)
    target_relative_paths = []
    for species in target_species:
        token = species_token_by_name[species]
        target_relative_paths.extend([token, token + '_busco.tsv'])
    target_relative_paths.append('busco_completeness.pdf')
    redo = _resolve_rerun_redo(args)
    _preflight_rerun_overwrite(
        target_root=busco_dir,
        relative_paths=target_relative_paths,
        redo=redo,
    )
    runtime_args = clone_namespace(args, redo=redo)
    with staged_output_dir(busco_dir, redo=True, prefix='amalgkit_rerun_busco_') as stage_busco_dir:
        if os.path.isdir(busco_dir):
            shutil.copytree(busco_dir, stage_busco_dir, dirs_exist_ok=True)
        if len(target_species) > 0:
            tool = select_busco_tool(runtime_args)
            full_metadata = load_metadata_if_needed_for_busco(runtime_args)
            subset_metadata = None
            if full_metadata is not None:
                subset_metadata = _subset_metadata(full_metadata, species_names=target_species)
            species, fasta_map = collect_busco_species(runtime_args, subset_metadata)
            busco_token_by_name = build_busco_species_token_map(
                subset_metadata,
                species,
            )
            extra_args = shlex.split(runtime_args.tool_args) if getattr(runtime_args, 'tool_args', None) else []
            stage_args = clone_namespace(runtime_args, redo=True)
            for species_name in species:
                process_species_busco(
                    sp=species_name,
                    fasta_path=fasta_map[species_name],
                    busco_dir=stage_busco_dir,
                    tool=tool,
                    args=stage_args,
                    extra_args=extra_args,
                    species_token=busco_token_by_name[species_name],
                )
        summary_path = os.path.join(stage_busco_dir, 'busco_completeness.pdf')
        generate_busco_species_plot(
            busco_dir=stage_busco_dir,
            species_order=all_species,
            out_path=summary_path,
            species_name_by_token={
                token: species_name
                for species_name, token in species_token_by_name.items()
            },
        )
        if (
            not os.path.isfile(summary_path)
            or os.path.getsize(summary_path) <= 0
        ):
            raise RuntimeError(
                'BUSCO rerun did not generate busco_completeness.pdf.'
            )


def rerun_finalize_check(args, metadata, target_species, force_all_species=False, dry_run=False):
    print('rerun finalize', flush=True)
    if force_all_species and (len(target_species) == 0):
        target_species = _list_selected_species(metadata)
    _print_targets('target_species', target_species)
    if len(target_species) == 0:
        return
    if dry_run:
        return
    finalize_dir = os.path.join(args.out_dir, 'finalize')
    _recover_target_transactions(finalize_dir)
    species_token_by_name = _build_rerun_species_token_map(metadata)
    target_species_dirs = [species_token_by_name[species] for species in target_species]
    redo = _resolve_rerun_redo(args)
    _preflight_rerun_overwrite(
        target_root=finalize_dir,
        relative_paths=target_species_dirs + FINALIZE_ROOT_OUTPUT_FILES,
        redo=redo,
    )
    resolve_args = clone_namespace(args)
    _resolved_metadata, input_dir = resolve_per_species_input(resolve_args)
    filtered_metadata = metadata if len(target_species) == 0 else _subset_metadata(metadata, species_names=target_species)
    stage_root = _create_staging_root(finalize_dir, prefix='amalgkit_rerun_finalize_')
    tmp_out_dir = os.path.join(stage_root, 'generated')
    stage_payload_dir = os.path.join(stage_root, 'payload')
    os.makedirs(tmp_out_dir, exist_ok=True)
    os.makedirs(stage_payload_dir, exist_ok=True)
    try:
        per_species_args = _build_per_species_args(args=resolve_args, input_dir=input_dir, tmp_out_dir=tmp_out_dir)
        per_species_args = clone_namespace(per_species_args, redo=redo)
        generate_per_species_tables(
            per_species_args,
            context=PerSpeciesTableContext(metadata=filtered_metadata, input_dir=input_dir),
        )
        per_species_dir = os.path.join(tmp_out_dir, 'per_species')
        _copy_species_tables(
            per_species_dir=per_species_dir,
            finalize_dir=stage_payload_dir,
            batch_effect_alg=per_species_args.batch_effect_alg,
        )
        copy_per_species_pdfs(
            per_species_dir=per_species_dir,
            dst_dir=stage_payload_dir,
            species_subset=target_species_dirs,
        )
        updated_metadata = load_merged_per_species_metadata(per_species_dir=per_species_dir)
        existing_finalize_metadata = _load_existing_metadata_table(os.path.join(finalize_dir, 'metadata.tsv'))
        source_df = metadata.df if existing_finalize_metadata is None else existing_finalize_metadata
        merged_metadata = merge_metadata_by_run(source_df, updated_metadata)
        merged_metadata.to_csv(os.path.join(stage_payload_dir, 'metadata.tsv'), sep='\t', index=False)
        save_exclusion_plot_pdf(
            df_metadata=merged_metadata,
            out_pdf_path=os.path.join(stage_payload_dir, 'finalize_exclusion.pdf'),
            y_label='Sample count',
            font_size=8,
        )
        _validate_staged_species_outputs(
            staged_root=stage_payload_dir,
            species_tokens=target_species_dirs,
            required_suffixes=FINALIZE_REQUIRED_SPECIES_SUFFIXES,
            label='finalize',
        )
        _commit_staged_paths(
            target_root=finalize_dir,
            staged_root=stage_payload_dir,
            relative_paths=target_species_dirs + FINALIZE_ROOT_OUTPUT_FILES,
        )
    finally:
        shutil.rmtree(stage_root, ignore_errors=True)


def rerun_main(args):
    report_path = _resolve_report_path(args)
    report_payload = _read_json(report_path)
    out_dir = _resolve_execution_out_dir(args, report_payload)
    metadata_path = _resolve_metadata_path(args, report_payload, out_dir)
    runtime_args = clone_namespace(
        args,
        out_dir=out_dir,
        metadata=metadata_path,
    )
    metadata = _load_metadata_from_path(metadata_path, out_dir)
    _build_rerun_species_token_map(metadata)
    allowed_runs = set(_parse_csv_option(getattr(args, 'run', None)))
    allowed_species = set(_parse_csv_option(getattr(args, 'species', None)))
    include_warnings = bool(getattr(args, 'include_warnings', False))
    dry_run = bool(getattr(args, 'dry_run', False))
    requested_checks = _resolve_requested_checks(args, report_payload)
    manifest_path = _resolve_manifest_path(args, out_dir)
    plan = _build_rerun_plan(
        report_payload=report_payload,
        metadata=metadata,
        requested_checks=requested_checks,
        allowed_runs=allowed_runs,
        allowed_species=allowed_species,
        include_warnings=include_warnings,
    )

    print('Using sanity report: {}'.format(report_path), flush=True)
    print('Using out_dir: {}'.format(out_dir), flush=True)
    print('Using metadata: {}'.format(metadata_path), flush=True)
    print('Requested rerun checks: {}'.format(', '.join(requested_checks)), flush=True)
    if manifest_path is not None:
        manifest_payload = _build_rerun_manifest(
            report_path=report_path,
            out_dir=out_dir,
            metadata_path=metadata_path,
            requested_checks=requested_checks,
            allowed_runs=allowed_runs,
            allowed_species=allowed_species,
            include_warnings=include_warnings,
            dry_run=dry_run,
            plan=plan,
        )
        runtime_log_path = _write_rerun_runtime_log(out_dir, _now_iso_for_runtime_log())
        manifest_payload['runtime_log_path'] = runtime_log_path
        _write_json(manifest_path, manifest_payload)
        print('Wrote rerun manifest: {}'.format(manifest_path), flush=True)

    for item in plan:
        check_name = item['check']
        if not item['will_execute']:
            print('rerun {}: {}'.format(check_name, item['skip_reason']), flush=True)
            continue
        if check_name == 'getfastq':
            rerun_getfastq_check(runtime_args, metadata, item['target_runs'], dry_run=dry_run)
        elif check_name == 'quant':
            rerun_quant_check(runtime_args, metadata, item['target_runs'], dry_run=dry_run)
        elif check_name == 'index':
            rerun_index_check(runtime_args, metadata, item['target_species'], dry_run=dry_run)
        elif check_name == 'merge':
            rerun_merge_check(runtime_args, metadata, item['target_species'], dry_run=dry_run)
        elif check_name == 'busco':
            rerun_busco_check(
                runtime_args,
                metadata,
                item['target_species'],
                rebuild_summary_only=item['rebuild_summary_only'],
                dry_run=dry_run,
            )
        elif check_name == 'finalize':
            rerun_finalize_check(runtime_args, metadata, item['target_species'], dry_run=dry_run)
        else:
            raise AmalgkitExit('Unsupported rerun check: {}'.format(check_name), exit_code=1)
