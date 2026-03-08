from amalgkit.batch_effect_sva import (
    SVAEstimate,
    resolve_sva_B_value,
    resolve_sva_parameters,
)


def test_resolve_sva_B_value_honors_valid_manual_setting():
    assert resolve_sva_B_value(B_setting='25', sample_count=8, auto_max=100) == 25


def test_resolve_sva_B_value_falls_back_to_auto_for_invalid_manual_setting():
    assert resolve_sva_B_value(B_setting='invalid', sample_count=10, auto_max=100) == 12


def test_resolve_sva_B_value_uses_small_sample_default_cap():
    assert resolve_sva_B_value(B_setting='auto', sample_count=1, auto_max=100) == 20


def test_resolve_sva_parameters_manual_nsv_clamps_to_maximum():
    resolved = resolve_sva_parameters(
        num_samples=4,
        design_columns=2,
        nsv_setting='99',
        B_setting='auto',
    )
    assert resolved.nsv == 1
    assert resolved.method == 'manual_nsv'
    assert resolved.stable is True


def test_resolve_sva_parameters_returns_max_nsv_zero_without_estimator():
    resolved = resolve_sva_parameters(
        num_samples=3,
        design_columns=2,
        nsv_setting='auto',
        B_setting='auto',
    )
    assert resolved.nsv == 0
    assert resolved.method == 'max_nsv_zero'
    assert resolved.stable is True


def test_resolve_sva_parameters_uses_manual_B_estimate_when_requested():
    resolved = resolve_sva_parameters(
        num_samples=8,
        design_columns=2,
        nsv_setting='auto',
        B_setting='15',
        estimate_nsv_at_B=lambda B, max_nsv: SVAEstimate(nsv=max_nsv + 5, method='be'),
    )
    assert resolved.B == 15
    assert resolved.nsv == 5
    assert resolved.method == 'manual_B_be'
    assert resolved.trace_B == [15]
    assert resolved.trace_nsv == [5]


def test_resolve_sva_parameters_auto_B_stops_when_two_consecutive_estimates_match():
    estimates = {
        15: SVAEstimate(nsv=3, method='be'),
        23: SVAEstimate(nsv=2, method='be'),
        30: SVAEstimate(nsv=2, method='be'),
    }
    resolved = resolve_sva_parameters(
        num_samples=8,
        design_columns=2,
        nsv_setting='auto',
        B_setting='auto',
        B_auto_max=30,
        estimate_nsv_at_B=lambda B, max_nsv: estimates[B],
    )
    assert resolved.B == 30
    assert resolved.nsv == 2
    assert resolved.method == 'auto_be'
    assert resolved.stable is True
    assert resolved.trace_B == [15, 23, 30]
    assert resolved.trace_nsv == [3, 2, 2]


def test_resolve_sva_parameters_auto_B_defaults_to_zero_when_estimates_fail():
    resolved = resolve_sva_parameters(
        num_samples=8,
        design_columns=2,
        nsv_setting='auto',
        B_setting='auto',
        B_auto_max=30,
        estimate_nsv_at_B=lambda B, max_nsv: SVAEstimate(nsv=None, method='failed'),
    )
    assert resolved.nsv == 0
    assert resolved.method == 'auto_failed'
    assert resolved.stable is False
