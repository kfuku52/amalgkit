import pytest

from amalgkit.batch_effect_common import BatchEffectResult
from amalgkit.batch_effect_io import (
    read_backend_summary_dcf,
    write_backend_summary_dcf,
)


def test_backend_summary_dcf_round_trip_restores_known_field_types(tmp_path):
    summary = BatchEffectResult(
        backend='sva',
        method='irw',
        stable=True,
        corrected_run_ids=['RUN1', 'RUN2'],
        uncorrected_run_ids=[],
        resolved_sva_nsv=2,
        resolved_sva_B=10,
        resolved_latent_k=None,
        latent_family=None,
        ruv_residual_method='least_squares',
        ruv_pvalue_method='one_way_anova',
        ruv_fallback_used=True,
        ruv_fallback_reason='statsmodels_unavailable',
        ruv_nb_fallback_genes=0,
        ruv_anova_failure_genes=1,
        latent_iterations=12,
        latent_objective=1.25,
        negative_values_before_clip=3,
        negative_values_after_clip=0,
        extra={
            'sva_stable': False,
            'trace_B': [5, 10],
            'trace_nsv': [1, 2],
            'trace_method': ['be', 'leek'],
            'counts_shape': [100, 2],
            'metadata_rows': 2,
            'group_model_used': True,
            'group_fallback_used': False,
            'ruv_baseline_score': 0.75,
            'latent_converged': True,
            'custom_text': '007',
        },
    )
    path = tmp_path / 'summary.dcf'

    write_backend_summary_dcf(summary, path)
    observed = read_backend_summary_dcf(path)

    assert observed == summary.to_jsonable()


def test_backend_summary_dcf_keeps_unknown_fields_as_strings(tmp_path):
    path = tmp_path / 'summary.dcf'
    write_backend_summary_dcf(
        {
            'custom_int': 2,
            'custom_bool': True,
            'custom_list': ['A', 'B'],
            'custom_empty': '',
        },
        path,
    )

    observed = read_backend_summary_dcf(path)

    assert observed == {
        'custom_int': '2',
        'custom_bool': 'TRUE',
        'custom_list': 'A|B',
        'custom_empty': '',
    }


@pytest.mark.parametrize(
    ('line', 'message'),
    [
        ('stable: maybe\n', 'Invalid boolean value'),
        ('resolved_sva_nsv: two\n', 'Invalid int value'),
        ('trace_B: 5|ten\n', 'Invalid integer-list value'),
        ('latent_objective: many\n', 'Invalid float value'),
    ],
)
def test_backend_summary_dcf_rejects_invalid_known_field_values(tmp_path, line, message):
    path = tmp_path / 'summary.dcf'
    path.write_text(line, encoding='utf-8')

    with pytest.raises(ValueError, match=message):
        read_backend_summary_dcf(path)
