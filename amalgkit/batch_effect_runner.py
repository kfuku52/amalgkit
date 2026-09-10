import argparse
import sys

import numpy

from amalgkit.batch_effect_common import BatchEffectResult
from amalgkit.batch_effect_contract import backend_options, read_control_gene_ids
from amalgkit.batch_effect_combatseq import run_combatseq_backend
from amalgkit.batch_effect_io import (
    read_backend_summary_dcf,
    read_backend_summary_json,
    read_expression_matrix_tsv,
    read_metadata_tsv,
    write_backend_summary_dcf,
    write_backend_summary_json,
    write_expression_matrix_tsv,
)
from amalgkit.batch_effect_latent_glm import run_latent_glm_backend
from amalgkit.batch_effect_ruvseq import run_ruvseq_backend
from amalgkit.batch_effect_sva import run_sva_backend


SUPPORTED_BACKENDS = ('sva', 'combatseq', 'ruvseq', 'latent_glm', 'latent_loglinear')


def build_parser():
    parser = argparse.ArgumentParser(
        description='Internal AMALGKIT helper for Python batch-effect backends.'
    )
    parser.add_argument('--backend', choices=SUPPORTED_BACKENDS, required=True)
    parser.add_argument('--batch_failure_policy', choices=['skip', 'error'], default='skip')
    parser.add_argument('--batch_categorical_covariates', nargs='+', default=[])
    parser.add_argument('--batch_continuous_covariates', nargs='+', default=[])
    parser.add_argument('--combatseq_group_model', choices=['protect', 'batch_only'], default='protect')
    parser.add_argument('--sva_irw_iterations', type=int, default=5)
    parser.add_argument('--sva_estimation_method', choices=['be', 'leek', 'be_then_leek'], default='be')
    parser.add_argument('--ruvseq_k_selection', choices=['manual', 'legacy'], default='manual')
    parser.add_argument('--ruvseq_control_file', default=None)
    parser.add_argument('--latent_k_selection', choices=['manual', 'legacy'], default='manual')
    parser.add_argument('--counts_tsv', required=True)
    parser.add_argument('--metadata_tsv', required=True)
    parser.add_argument('--options_json', required=False, default=None)
    parser.add_argument('--out_summary_json', required=False, default=None)
    parser.add_argument('--out_summary_dcf', required=False, default=None)
    parser.add_argument('--out_counts_tsv', required=False, default=None)
    parser.add_argument('--out_sv_tsv', required=False, default=None)
    parser.add_argument('--sva_nsv', required=False, default='auto')
    parser.add_argument('--sva_B', '--sva_nsv_permutations', dest='sva_B', required=False, default='auto')
    parser.add_argument('--sva_B_auto_max', required=False, default='100')
    parser.add_argument(
        '--sva_input_scale',
        choices=('counts', 'transformed'),
        required=False,
        default='counts',
    )
    parser.add_argument('--batch_column', required=False, default='bioproject')
    parser.add_argument('--sample_group_column', required=False, default='sample_group')
    parser.add_argument('--ruvseq_control_mode', required=False, default='auto')
    parser.add_argument('--ruvseq_k', required=False, default='auto')
    parser.add_argument('--ruvseq_k_max', required=False, default='5')
    parser.add_argument('--ruvseq_control_top_n', required=False, default='1000')
    parser.add_argument('--ruvseq_min_controls', required=False, default='100')
    parser.add_argument('--latent_weighting', '--latent_family', dest='latent_family', required=False, default='dispersion')
    parser.add_argument('--latent_k', required=False, default='auto')
    parser.add_argument('--latent_k_max', required=False, default='5')
    parser.add_argument('--latent_max_iter', required=False, default='200')
    parser.add_argument('--latent_tol', required=False, default='1e-5')
    parser.add_argument('--random_seed', required=False, default='0')
    return parser


def _write_failure_summary(args, counts, metadata, exc, skip_reason):
    summary = BatchEffectResult(
        backend=args.backend,
        method='error',
        skip_reason=getattr(exc, 'reason', skip_reason),
        extra={
            **getattr(exc, 'diagnostics', {}),
            'counts_shape': [int(counts.shape[0]), int(counts.shape[1])],
            'metadata_rows': int(metadata.shape[0]),
            'error_message': str(exc),
            'status': 'error',
            'batch_failure_policy': args.batch_failure_policy,
        },
    )
    if args.out_summary_json is not None:
        write_backend_summary_json(summary, args.out_summary_json)
    if args.out_summary_dcf is not None:
        write_backend_summary_dcf(summary, args.out_summary_dcf)
    sys.stderr.write('ERROR: {}\n'.format(exc))
    return 1


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    counts = read_expression_matrix_tsv(args.counts_tsv)
    metadata = read_metadata_tsv(args.metadata_tsv)
    _ = args.options_json
    if args.backend == 'sva':
        try:
            corrected_df, sv_df, summary = run_sva_backend(
                counts_df=counts,
                metadata_df=metadata,
                nsv_setting=args.sva_nsv,
                B_setting=args.sva_B,
                B_auto_max=args.sva_B_auto_max,
                sample_group_column=args.sample_group_column,
                random_seed=args.random_seed,
                input_scale=args.sva_input_scale,
                irw_iterations=args.sva_irw_iterations,
                estimation_method=args.sva_estimation_method,
                **backend_options(args),
            )
        except (ValueError, RuntimeError, numpy.linalg.LinAlgError) as exc:
            return _write_failure_summary(
                args=args,
                counts=counts,
                metadata=metadata,
                exc=exc,
                skip_reason='sva_fit_failed',
            )
        if args.out_counts_tsv is not None:
            write_expression_matrix_tsv(corrected_df, args.out_counts_tsv)
        if (args.out_sv_tsv is not None) and (sv_df is not None):
            write_expression_matrix_tsv(sv_df, args.out_sv_tsv)
        if args.out_summary_json is not None:
            write_backend_summary_json(summary, args.out_summary_json)
        if args.out_summary_dcf is not None:
            write_backend_summary_dcf(summary, args.out_summary_dcf)
        return 0
    if args.backend == 'combatseq':
        try:
            corrected_df, summary = run_combatseq_backend(
                counts_df=counts,
                metadata_df=metadata,
                batch_column=args.batch_column,
                sample_group_column=args.sample_group_column,
                protect_group=args.combatseq_group_model == 'protect',
                **backend_options(args),
            )
        except (ImportError, ValueError, RuntimeError) as exc:
            return _write_failure_summary(
                args=args,
                counts=counts,
                metadata=metadata,
                exc=exc,
                skip_reason='backend_dependency_missing' if isinstance(exc, ImportError) else 'combatseq_fit_failed',
            )
        if args.out_counts_tsv is not None:
            write_expression_matrix_tsv(corrected_df, args.out_counts_tsv)
        if args.out_summary_json is not None:
            write_backend_summary_json(summary, args.out_summary_json)
        if args.out_summary_dcf is not None:
            write_backend_summary_dcf(summary, args.out_summary_dcf)
        return 0
    if args.backend == 'ruvseq':
        try:
            corrected_df, w_df, summary = run_ruvseq_backend(
                counts_df=counts,
                metadata_df=metadata,
                control_mode=args.ruvseq_control_mode,
                k_setting=args.ruvseq_k,
                k_max=args.ruvseq_k_max,
                top_n=args.ruvseq_control_top_n,
                min_controls=args.ruvseq_min_controls,
                k_selection=args.ruvseq_k_selection,
                control_gene_ids=read_control_gene_ids(args.ruvseq_control_file),
                **backend_options(args),
                batch_column=args.batch_column,
                sample_group_column=args.sample_group_column,
            )
        except (ValueError, RuntimeError, numpy.linalg.LinAlgError) as exc:
            return _write_failure_summary(
                args=args,
                counts=counts,
                metadata=metadata,
                exc=exc,
                skip_reason='ruvseq_fit_failed',
            )
        if args.out_counts_tsv is not None:
            write_expression_matrix_tsv(corrected_df, args.out_counts_tsv)
        if (args.out_sv_tsv is not None) and (w_df is not None):
            write_expression_matrix_tsv(w_df, args.out_sv_tsv)
        if args.out_summary_json is not None:
            write_backend_summary_json(summary, args.out_summary_json)
        if args.out_summary_dcf is not None:
            write_backend_summary_dcf(summary, args.out_summary_dcf)
        return 0
    if args.backend in {'latent_glm', 'latent_loglinear'}:
        try:
            corrected_df, latent_df, summary = run_latent_glm_backend(
                counts_df=counts,
                metadata_df=metadata,
                family=args.latent_family,
                k_setting=args.latent_k,
                k_max=args.latent_k_max,
                sample_group_column=args.sample_group_column,
                max_iter=args.latent_max_iter,
                tol=float(args.latent_tol),
                k_selection=args.latent_k_selection,
                **backend_options(args),
            )
        except (ValueError, RuntimeError, numpy.linalg.LinAlgError) as exc:
            return _write_failure_summary(
                args=args,
                counts=counts,
                metadata=metadata,
                exc=exc,
                skip_reason='latent_glm_fit_failed',
            )
        if args.out_counts_tsv is not None:
            write_expression_matrix_tsv(corrected_df, args.out_counts_tsv)
        if (args.out_sv_tsv is not None) and (latent_df is not None):
            write_expression_matrix_tsv(latent_df, args.out_sv_tsv)
        if args.out_summary_json is not None:
            write_backend_summary_json(summary, args.out_summary_json)
        if args.out_summary_dcf is not None:
            write_backend_summary_dcf(summary, args.out_summary_dcf)
        return 0
    summary = BatchEffectResult(
        backend=args.backend,
        method='unimplemented',
        skip_reason='backend_not_implemented',
        extra={
            'counts_shape': [int(counts.shape[0]), int(counts.shape[1])],
            'metadata_rows': int(metadata.shape[0]),
        },
    )
    if args.out_summary_json is not None:
        write_backend_summary_json(summary, args.out_summary_json)
    if args.out_summary_dcf is not None:
        write_backend_summary_dcf(summary, args.out_summary_dcf)
    sys.stderr.write(
        'ERROR: Python batch backend "{}" is not implemented yet.\n'.format(args.backend)
    )
    return 1


__all__ = [
    'SUPPORTED_BACKENDS',
    'build_parser',
    'main',
    'read_backend_summary_json',
    'read_backend_summary_dcf',
]


if __name__ == '__main__':
    raise SystemExit(main())
