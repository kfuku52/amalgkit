"""Shared input, biological-design and recoverable-failure contract.

Only BatchModelError is recoverable. Invalid inputs and unexpected exceptions
must remain errors, even when the requested failure policy is ``skip``.
"""
from dataclasses import dataclass
from functools import wraps
import inspect
from importlib.metadata import PackageNotFoundError, version
import sys

import numpy
import pandas

from amalgkit.batch_effect_common import BatchEffectResult, align_metadata_to_counts


class BatchModelError(RuntimeError):
    def __init__(self, reason, message, diagnostics=None):
        super().__init__(message)
        self.reason = reason
        self.diagnostics = {} if diagnostics is None else diagnostics


@dataclass
class BiologicalDesign:
    matrix: pandas.DataFrame
    diagnostics: dict


def orthogonal_basis(matrix, reference_scale=None):
    values = numpy.asarray(matrix, dtype=float)
    left, singular, _ = numpy.linalg.svd(values, full_matrices=False)
    scale = max(1.0, float(singular[0]) if singular.size else 0.0)
    if reference_scale is not None:
        scale = max(scale, float(reference_scale))
    tolerance = numpy.finfo(float).eps * max(values.shape, default=1) * scale
    return left[:, singular > tolerance]


def removal_basis(factors, design):
    factors = numpy.asarray(factors, dtype=float)
    basis = orthogonal_basis(design)
    residual = factors - basis @ (basis.T @ factors)
    return orthogonal_basis(residual, reference_scale=numpy.linalg.norm(factors))


def matrix_change(before, after, operation, scale):
    """Describe an actual postprocessing step, on its stated numeric scale."""
    before, after = numpy.asarray(before, dtype=float), numpy.asarray(after, dtype=float)
    changed = before != after
    # Equal log(0) cells are unchanged, not an undefined (-inf)-(-inf) delta.
    delta = numpy.zeros_like(before)
    numpy.subtract(after, before, out=delta, where=changed)
    numpy.abs(delta, out=delta)
    return {
        'operation': operation, 'scale': scale,
        'changed_cells': int(changed.sum()),
        'changed_genes': int(changed.any(axis=1).sum()),
        'changed_runs': int(changed.any(axis=0).sum()),
        'absolute_change_sum': float(delta.sum()),
        'maximum_absolute_change': float(delta.max()) if delta.size else 0.0,
    }


def build_biological_design(metadata, sample_group_column='sample_group',
                            categorical_covariates=(), continuous_covariates=(),
                            protect_group=True, batch_column='bioproject'):
    categorical = list(categorical_covariates or ())
    continuous = list(continuous_covariates or ())
    if protect_group:
        categorical.insert(0, sample_group_column)
    if len(set(categorical + continuous)) != len(categorical + continuous):
        raise ValueError('Protected covariates must be unique and have one explicit type.')
    matrix = pandas.DataFrame({'Intercept': numpy.ones(len(metadata))}, index=metadata['run'])
    encoding = {}
    for column in categorical + continuous:
        if column not in metadata:
            raise ValueError('Missing required metadata column: {}'.format(column))
        source = metadata[column]
        if source.isna().any() or source.astype(str).str.strip().eq('').any():
            raise ValueError('Protected covariate contains missing values: {}'.format(column))
        if column in categorical:
            labels = source.astype(str).str.strip()
            levels = sorted(labels.unique())
            encoding[column] = {'type': 'categorical', 'levels': levels, 'reference': levels[0]}
            for level in levels[1:]:
                name = '{}[{}]'.format(column, level)
                if name in matrix:
                    raise ValueError('Ambiguous protected-design column: {}'.format(name))
                matrix[name] = labels.eq(level).to_numpy(dtype=float)
        else:
            if column in matrix:
                raise ValueError('Ambiguous protected-design column: {}'.format(column))
            values = pandas.to_numeric(source, errors='raise').to_numpy(dtype=float)
            if not numpy.isfinite(values).all():
                raise ValueError('Nonfinite continuous covariate: {}'.format(column))
            center, scale = float(values.mean()), float(values.std())
            encoding[column] = {'type': 'continuous', 'center': center, 'scale': scale}
            matrix[column] = (values - center) / scale if scale > 0 else numpy.zeros(len(values))
    singular = numpy.linalg.svd(matrix.to_numpy(), compute_uv=False)
    rank = orthogonal_basis(matrix).shape[1]
    diagnostics = {
        'design_columns': list(matrix.columns), 'design_rank': rank,
        'design_residual_df': len(metadata) - rank,
        'design_singular_values': singular.tolist(), 'design_encoding': encoding,
        'design_run_ids': metadata['run'].tolist(), 'group_protected': bool(protect_group),
        'design_matrix': matrix.to_numpy(dtype=float).tolist(),
        'batch_column': batch_column,
    }
    if batch_column in metadata:
        batch = metadata[batch_column].fillna('').astype(str).str.strip()
        diagnostics['batch_labels_complete'] = bool(batch.ne('').all())
        if batch.ne('').all():
            batch_design = pandas.get_dummies(batch, drop_first=True, dtype=float).to_numpy()
            combined = numpy.column_stack([matrix.to_numpy(), batch_design])
            diagnostics['batch_design_columns'] = combined.shape[1]
            diagnostics['batch_design_rank'] = orthogonal_basis(combined).shape[1]
            diagnostics['batch_design_confounded'] = diagnostics['batch_design_rank'] < combined.shape[1]
            if sample_group_column in metadata:
                table = pandas.crosstab(metadata[sample_group_column], batch)
                diagnostics['group_batch_counts'] = table.to_dict()
    return BiologicalDesign(matrix, diagnostics)


def batch_backend(name, factors=True):
    """Add the same policy to direct, runner and finalize backend calls."""
    def decorate(function):
        signature = inspect.signature(function)

        @wraps(function)
        def run(*args, failure_policy='skip', categorical_covariates=(),
                continuous_covariates=(), protect_group=True, **kwargs):
            if failure_policy not in {'skip', 'error'}:
                raise ValueError('batch failure policy must be skip or error.')
            bound = signature.bind(*args, **kwargs)
            bound.apply_defaults()
            counts = bound.arguments['counts_df']
            if counts.shape[1] == 0:
                raise ValueError('Batch correction requires at least one sample.')
            if counts.index.has_duplicates or counts.index.isna().any():
                raise ValueError('Expression gene IDs must be unique and nonmissing.')
            values = counts.to_numpy(dtype=float)
            if not numpy.isfinite(values).all():
                raise ValueError('Batch correction input must contain only finite values.')
            scale = bound.arguments.get('input_scale', 'counts')
            if scale == 'counts' and numpy.any(values < 0):
                raise ValueError('Count-scale input must not contain negative values.')
            if name == 'combatseq' and numpy.any(values != numpy.floor(values)):
                raise ValueError('ComBat-seq requires integer raw counts; do not pass normalized or fractional counts.')
            metadata = align_metadata_to_counts(counts, bound.arguments['metadata_df'])
            design = build_biological_design(
                metadata, bound.arguments.get('sample_group_column', 'sample_group'),
                categorical_covariates, continuous_covariates, protect_group,
                bound.arguments.get('batch_column', 'bioproject'),
            )
            bound.arguments['protected_design'] = design
            bound.arguments['metadata_df'] = metadata
            common = {
                'schema_version': 2, 'batch_failure_policy': failure_policy,
                'input_scale': scale, 'design': design.diagnostics,
                'design_rank': design.diagnostics['design_rank'],
                'design_residual_df': design.diagnostics['design_residual_df'],
                'batch_design_confounded': design.diagnostics.get('batch_design_confounded'),
                'requested_parameters': {key: value for key, value in bound.arguments.items()
                                         if key not in {'counts_df', 'metadata_df', 'protected_design', 'control_gene_ids'}},
                'package_versions': {},
            }
            for package in ('amalgkit', 'numpy', 'scipy', 'pandas', 'statsmodels', 'inmoose'):
                try:
                    common['package_versions'][package] = version(package)
                except PackageNotFoundError:
                    continue
            try:
                if design.diagnostics['design_rank'] < design.matrix.shape[1]:
                    raise BatchModelError('protected_design_rank_deficient', 'Protected design is rank deficient.')
                if counts.shape[1] < 2:
                    raise BatchModelError('single_sample', 'Batch correction requires more than one sample.')
                try:
                    result = function(*bound.args, **bound.kwargs)
                except (FloatingPointError, numpy.linalg.LinAlgError) as exc:
                    raise BatchModelError(name + '_numerical_failure', str(exc)) from exc
            except BatchModelError as exc:
                summary = BatchEffectResult(backend=name, method='skipped').to_jsonable()
                summary.update(exc.diagnostics)
                summary.update(common)
                summary.update({
                    'backend': name, 'method': 'skipped', 'status': 'skipped',
                    'skip_reason': exc.reason, 'error_message': str(exc), 'stable': False,
                    'corrected_run_ids': [], 'uncorrected_run_ids': list(counts.columns),
                })
                if failure_policy == 'error':
                    exc.diagnostics = dict(summary, status='error', method='error')
                    raise
                print('Batch correction skipped ({}): {}: {}'.format(name, exc.reason, exc), file=sys.stderr)
                if factors:
                    return counts.copy(), pandas.DataFrame(index=counts.columns), summary
                return counts.copy(), summary
            corrected, summary = result[0], dict(result[-1])
            if not corrected.index.equals(counts.index) or not corrected.columns.equals(counts.columns):
                raise RuntimeError('Batch backend changed expression IDs or their order.')
            if not numpy.isfinite(corrected.to_numpy(dtype=float)).all():
                raise RuntimeError('Batch backend returned nonfinite expression values.')
            summary.update(common)
            summary.setdefault('status', 'corrected' if summary.get('corrected_run_ids') else 'not_needed')
            if factors:
                factor_df = result[1]
                basis = removal_basis(factor_df.to_numpy(dtype=float), design.matrix.to_numpy())
                summary['removal_basis'] = basis.tolist()
                summary['factor_columns'] = list(factor_df.columns)
                summary['factor_values'] = factor_df.to_numpy(dtype=float).tolist()
            return (*result[:-1], summary)
        return run
    return decorate


def backend_options(args):
    return {
        'failure_policy': getattr(args, 'batch_failure_policy', 'skip'),
        'categorical_covariates': getattr(args, 'batch_categorical_covariates', ()) or (),
        'continuous_covariates': getattr(args, 'batch_continuous_covariates', ()) or (),
    }


def read_control_gene_ids(path):
    if path is None:
        return None
    with open(path, encoding='utf-8') as handle:
        ids = [line.strip() for line in handle if line.strip()]
    if not ids or len(set(ids)) != len(ids) or any('\t' in value for value in ids):
        raise ValueError('Control file must contain one unique gene ID per line, without a header.')
    return ids
