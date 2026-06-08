import numpy
import pandas


def _load_pycombat_seq():
    try:
        from inmoose.pycombat import pycombat_seq
    except ImportError as exc:
        raise ImportError(
            'Python Combat-seq backend requires the "inmoose" package.'
        ) from exc
    return pycombat_seq


def _align_metadata_to_counts(counts_df, metadata_df):
    if 'run' not in metadata_df.columns:
        raise ValueError('Missing required metadata column: run')
    run_indexed = metadata_df.copy()
    run_indexed['run'] = run_indexed['run'].astype(str)
    missing_runs = [run_id for run_id in counts_df.columns if run_id not in set(run_indexed['run'])]
    if missing_runs:
        raise ValueError('Metadata is missing rows for runs: {}'.format(', '.join(missing_runs)))
    aligned = run_indexed.drop_duplicates(subset=['run'], keep='first').set_index('run')
    return aligned.loc[list(counts_df.columns), :].reset_index()


def _coerce_corrected_matrix(corrected, index, columns):
    if isinstance(corrected, pandas.DataFrame):
        corrected_df = corrected.copy()
        corrected_df.index = index
        corrected_df.columns = columns
        return corrected_df
    values = numpy.asarray(corrected)
    if values.shape != (len(index), len(columns)):
        raise ValueError(
            'Unexpected Combat-seq output shape: expected {}x{}, got {}x{}.'.format(
                len(index),
                len(columns),
                values.shape[0] if values.ndim >= 1 else 0,
                values.shape[1] if values.ndim >= 2 else 0,
            )
        )
    return pandas.DataFrame(values, index=index, columns=columns)


def run_combatseq_backend(
    counts_df,
    metadata_df,
    batch_column='bioproject',
    sample_group_column='sample_group',
):
    if counts_df.shape[1] == 0:
        raise ValueError('counts_df must contain at least one sample column.')
    aligned_metadata = _align_metadata_to_counts(counts_df=counts_df, metadata_df=metadata_df)
    if batch_column not in aligned_metadata.columns:
        raise ValueError('Missing required metadata column: {}'.format(batch_column))

    batch_labels = aligned_metadata.loc[:, batch_column].fillna('').astype(str).str.strip()
    if (batch_labels == '').any():
        raise ValueError('Batch column contains empty values: {}'.format(batch_column))
    batch_sizes = batch_labels.value_counts()

    corrected_run_ids = []
    uncorrected_run_ids = []
    for run_id, batch_label in zip(counts_df.columns, batch_labels):
        if int(batch_sizes[batch_label]) > 1:
            corrected_run_ids.append(str(run_id))
        else:
            uncorrected_run_ids.append(str(run_id))

    if len(corrected_run_ids) == 0:
        summary = {
            'backend': 'combatseq',
            'method': 'all_singleton',
            'skip_reason': 'combatseq_all_singleton',
            'stable': None,
            'corrected_run_ids': [],
            'uncorrected_run_ids': uncorrected_run_ids,
            'batch_column': batch_column,
            'group_model_used': False,
            'group_fallback_used': False,
        }
        return counts_df.copy(), summary

    pycombat_seq = _load_pycombat_seq()
    combat_counts = counts_df.loc[:, corrected_run_ids]
    combat_metadata = aligned_metadata.set_index('run').loc[corrected_run_ids, :]
    combat_batches = combat_metadata.loc[:, batch_column].astype(str).tolist()

    group_model_used = False
    group_fallback_used = False
    method = 'no_group'
    group_error_message = ''

    def run_without_group():
        corrected = pycombat_seq(
            counts=combat_counts,
            batch=combat_batches,
        )
        return _coerce_corrected_matrix(
            corrected=corrected,
            index=combat_counts.index,
            columns=combat_counts.columns,
        )

    if sample_group_column in combat_metadata.columns:
        sample_groups = combat_metadata.loc[:, sample_group_column].fillna('').astype(str).str.strip()
        if ((sample_groups != '').all()) and (sample_groups.nunique() > 1):
            covariates = pandas.DataFrame(
                {sample_group_column: sample_groups.tolist()},
                index=combat_counts.columns,
            )
            try:
                corrected = pycombat_seq(
                    counts=combat_counts,
                    batch=combat_batches,
                    covar_mod=covariates,
                )
                group_model_used = True
                method = 'group'
                corrected_df = _coerce_corrected_matrix(
                    corrected=corrected,
                    index=combat_counts.index,
                    columns=combat_counts.columns,
                )
            except Exception as exc:  # pragma: no cover - exercised via integration tests
                group_fallback_used = True
                group_error_message = str(exc)
                corrected_df = run_without_group()
        else:
            corrected_df = run_without_group()
    else:
        corrected_df = run_without_group()

    corrected_full = counts_df.copy()
    corrected_full.loc[:, corrected_run_ids] = corrected_df.loc[:, corrected_run_ids]
    summary = {
        'backend': 'combatseq',
        'method': method,
        'skip_reason': 'combatseq_singleton_kept' if len(uncorrected_run_ids) > 0 else '',
        'stable': None,
        'corrected_run_ids': corrected_run_ids,
        'uncorrected_run_ids': uncorrected_run_ids,
        'batch_column': batch_column,
        'group_model_used': group_model_used,
        'group_fallback_used': group_fallback_used,
        'group_error_message': group_error_message,
    }
    return corrected_full, summary
