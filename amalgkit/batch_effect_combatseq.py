import numpy
import pandas

from amalgkit.batch_effect_common import align_metadata_to_counts
from amalgkit.batch_effect_contract import BatchModelError, batch_backend


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
        if not corrected.index.equals(index) or not corrected.columns.equals(columns):
            raise RuntimeError('ComBat-seq returned reordered or changed IDs.')
        return corrected.copy()
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


@batch_backend('combatseq', factors=False)
def run_combatseq_backend(
    counts_df,
    metadata_df,
    batch_column='bioproject',
    sample_group_column='sample_group',
    protected_design=None,
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

    diagnostics = {'group_model_used': False, 'group_fallback_used': False}
    if (batch_sizes <= 1).any():
        reason = 'combatseq_all_singleton' if (batch_sizes <= 1).all() else 'combatseq_singleton_batch'
        raise BatchModelError(reason, 'ComBat-seq requires at least two samples in every batch.', diagnostics)
    if len(batch_sizes) < 2:
        raise BatchModelError('combatseq_insufficient_batches', 'ComBat-seq requires at least two batches.', diagnostics)
    if protected_design.diagnostics['batch_design_confounded']:
        raise BatchModelError('combatseq_confounded_design', 'Biological covariates and batch are not identifiable.', diagnostics)
    if protected_design.diagnostics['batch_design_rank'] >= counts_df.shape[1]:
        raise BatchModelError('combatseq_no_residual_df', 'No residual degrees of freedom for ComBat-seq.', diagnostics)
    if counts_df.shape[0] == 0:
        return counts_df.copy(), dict(backend='combatseq', method='empty', skip_reason='no_expressed_genes',
                                     corrected_run_ids=[], uncorrected_run_ids=list(counts_df.columns))
    pycombat_seq = _load_pycombat_seq()
    covariates = protected_design.matrix
    group_model_used = protected_design.diagnostics['group_protected'] and metadata_df[sample_group_column].nunique() > 1
    method = 'group' if group_model_used else 'no_group'
    call = {'counts': counts_df, 'batch': batch_labels.tolist()}
    if covariates.shape[1] > 1:
        # InMoose accepts a prebuilt patsy design and bypasses formula parsing.
        # This also preserves continuous covariates and unusual category names.
        # Keep the intercept: InMoose uses this matrix for within-batch
        # dispersion fits before removing its intercept from the joint design.
        from patsy import DesignInfo, DesignMatrix
        call['covar_mod'] = DesignMatrix(covariates.to_numpy(dtype=float),
                                         DesignInfo(list(covariates.columns)))
    try:
        corrected = pycombat_seq(**call)
    except (FloatingPointError, ValueError, numpy.linalg.LinAlgError) as exc:
        raise BatchModelError('combatseq_fit_failed', str(exc), dict(diagnostics, group_error_message=str(exc))) from exc
    corrected_df = _coerce_corrected_matrix(corrected, counts_df.index, counts_df.columns)
    if not numpy.isfinite(corrected_df.to_numpy(dtype=float)).all():
        raise BatchModelError('combatseq_nonfinite_fit', 'ComBat-seq returned nonfinite values.', diagnostics)
    corrected_run_ids = list(counts_df.columns)
    uncorrected_run_ids = []
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
        'group_fallback_used': False,
        'group_error_message': '',
    }
    return corrected_full, summary
