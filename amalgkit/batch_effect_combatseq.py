import numpy
import pandas
import warnings

from amalgkit.batch_effect_common import align_metadata_to_counts


def _load_pycombat_seq():
    try:
        from inmoose.pycombat import pycombat_seq
    except ImportError as exc:
        raise ImportError(
            'Python Combat-seq backend requires the "inmoose" package. '
            'Install it with: pip install "amalgkit[combatseq]"'
        ) from exc
    return pycombat_seq


def _align_metadata_to_counts(counts_df, metadata_df):
    return align_metadata_to_counts(counts_df=counts_df, metadata_df=metadata_df)


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

    combat_counts = counts_df.loc[:, corrected_run_ids]
    combat_metadata = aligned_metadata.set_index('run').loc[corrected_run_ids, :]
    combat_batches = combat_metadata.loc[:, batch_column].astype(str).tolist()
    if len(set(combat_batches)) < 2:
        summary = {
            'backend': 'combatseq',
            'method': 'insufficient_batches',
            'skip_reason': 'combatseq_insufficient_batches',
            'stable': None,
            'corrected_run_ids': [],
            'uncorrected_run_ids': [str(run_id) for run_id in counts_df.columns],
            'batch_column': batch_column,
            'group_model_used': False,
            'group_fallback_used': False,
        }
        return counts_df.copy(), summary

    pycombat_seq = _load_pycombat_seq()

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
            except (FloatingPointError, ValueError, numpy.linalg.LinAlgError) as exc:
                group_fallback_used = True
                group_error_message = str(exc)
                warnings.warn(
                    'Combat-seq could not fit the sample-group covariate model and '
                    'fell back to batch-only correction: {}'.format(exc),
                    stacklevel=2,
                )
                corrected_df = run_without_group()
        else:
            corrected_df = run_without_group()
    else:
        corrected_df = run_without_group()

    corrected_full = counts_df.copy()
    for run_id in corrected_run_ids:
        values = corrected_df[run_id]
        original_dtype = counts_df[run_id].dtype
        if pandas.api.types.is_integer_dtype(original_dtype):
            bounds = numpy.iinfo(getattr(original_dtype, 'numpy_dtype', original_dtype))
            array = values.to_numpy(dtype=float)
            if (
                numpy.isfinite(array).all()
                and (array == numpy.floor(array)).all()
                and (array >= bounds.min).all()
                and (array < int(bounds.max) + 1).all()
            ):
                # Retain the integer count contract only for lossless casts.
                values = values.astype(original_dtype)
        # Replacing the column supports pandas 3 without truncating fractional
        # results or changing the untouched singleton-batch columns.
        corrected_full[run_id] = values
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
