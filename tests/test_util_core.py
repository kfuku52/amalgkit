from contextlib import nullcontext
import threading
import time

import numpy
import pandas
import pytest

from amalgkit.parallel_utils import raise_task_failures

from amalgkit.util import (
    strtobool,
    parse_bool_flags,
    run_tasks_with_optional_threads,
    validate_positive_int_option,
    resolve_cpu_budget,
    resolve_thread_worker_allocation,
    resolve_worker_allocation,
    find_prefixed_entries,
    find_species_prefixed_entries,
    find_run_prefixed_entries,
)


# ---------------------------------------------------------------------------
# strtobool
# ---------------------------------------------------------------------------


class TestStrtobool:
    @pytest.mark.parametrize("val", ["y", "yes", "Yes", "YES", "t", "true", "True", "on", "1"])
    def test_true_values(self, val):
        assert strtobool(val) is True

    @pytest.mark.parametrize("val", ["n", "no", "No", "NO", "f", "false", "False", "off", "0"])
    def test_false_values(self, val):
        assert strtobool(val) is False

    def test_invalid_value(self):
        with pytest.raises(ValueError):
            strtobool("maybe")


class TestParseBoolFlags:
    def test_fills_missing_values_with_default(self):
        result = parse_bool_flags((v for v in ['yes', None, '']), column_name='is_sampled', default='no')
        assert result.tolist() == [True, False, False]

    def test_raises_for_invalid_values(self):
        with pytest.raises(ValueError, match='invalid boolean flag'):
            parse_bool_flags(['yes', 'maybe'], column_name='is_sampled', default='no')

    @pytest.mark.parametrize('as_generator', [False, True])
    def test_preserves_boolean_tokens_mixed_values_and_order(self, as_generator):
        true_tokens = ['y', 'YES', ' t ', 'True', 'ON', '1', True, 1]
        false_tokens = ['n', 'NO', ' f ', 'False', 'OFF', '0', False, 0]
        values = [value for pair in zip(true_tokens, false_tokens) for value in pair] * 3
        result = parse_bool_flags((v for v in values) if as_generator else values)
        assert isinstance(result, numpy.ndarray)
        assert result.dtype == numpy.dtype(bool)
        assert result.shape == (len(values),)
        assert result.tolist() == [True, False] * len(true_tokens) * 3

    @pytest.mark.parametrize('default, expected', [('no', False), ('yes', True), (1, True)])
    def test_missing_values_with_string_or_numeric_default(self, default, expected):
        result = parse_bool_flags([None, numpy.nan, pandas.NA], default=default)
        assert result.tolist() == [expected] * 3

    @pytest.mark.parametrize('default', [1, True, False])
    def test_empty_string_preserves_nonstring_default_error(self, default):
        with pytest.raises(AttributeError, match='lower'):
            parse_bool_flags([''], default=default)

    def test_invalid_values_report_original_sorted_unique_text(self):
        values = [' MAYBE ', 'maybe', 'NA', '0001', ' MAYBE ', None, '']
        with pytest.raises(ValueError) as error:
            parse_bool_flags((v for v in values), column_name='is_sampled', default='no')
        assert str(error.value) == (
            'Column "is_sampled" contains invalid boolean flag(s):  MAYBE , 0001, NA, maybe'
        )

    def test_empty_iterable_is_one_dimensional_boolean_array(self):
        result = parse_bool_flags(iter(()), default='no')
        assert result.dtype == numpy.dtype(bool)
        assert result.shape == (0,)

    @pytest.mark.parametrize('default', [True, False, numpy.bool_(True), numpy.bool_(False)])
    def test_boolean_default_fills_float_missing_values_without_warning(self, default):
        result = parse_bool_flags([numpy.nan], default=default)
        assert result.tolist() == [bool(default)]

    @pytest.mark.parametrize('default', [1, numpy.int64(1)])
    @pytest.mark.parametrize('future_downcasting', [False, True])
    def test_numeric_default_preserves_active_object_inference(self, default, future_downcasting):
        pandas_major, pandas_minor = (int(v) for v in pandas.__version__.split('.')[:2])
        has_downcast_option = pandas_major == 2 and pandas_minor >= 2
        context = (
            pandas.option_context('future.no_silent_downcasting', future_downcasting)
            if has_downcast_option else nullcontext()
        )
        legacy_inference = pandas_major == 2 and not (has_downcast_option and future_downcasting)
        expected_invalid = '1.0, <NA>' if legacy_inference else '1.0'
        with context:
            with pytest.raises(ValueError) as error:
                parse_bool_flags([1.0, pandas.NA], column_name='is_sampled', default=default)
            assert str(error.value) == (
                'Column "is_sampled" contains invalid boolean flag(s): ' + expected_invalid
            )
            if has_downcast_option:
                assert pandas.get_option('future.no_silent_downcasting') is future_downcasting

class TestRunTasksWithOptionalThreads:
    def test_empty_tasks(self):
        results, failures = run_tasks_with_optional_threads([], lambda x: x, max_workers=4)
        assert results == {}
        assert failures == []

    def test_serial_collects_successes_and_failures(self):
        def worker(x):
            if x == 2:
                raise RuntimeError('boom')
            return x * 10

        results, failures = run_tasks_with_optional_threads([1, 2, 3], worker, max_workers=1)

        assert results == {1: 10, 3: 30}
        assert len(failures) == 1
        assert failures[0][0] == 2
        assert isinstance(failures[0][1], RuntimeError)

    def test_parallel_collects_successes_and_failures(self):
        def worker(x):
            if x == 3:
                raise ValueError('bad')
            return x + 1

        results, failures = run_tasks_with_optional_threads([1, 2, 3, 4], worker, max_workers=3)

        assert results == {1: 2, 2: 3, 4: 5}
        assert len(failures) == 1
        assert failures[0][0] == 3
        assert isinstance(failures[0][1], ValueError)

    def test_parallel_converts_system_exit_to_failure(self):
        def worker(x):
            if x == 2:
                raise SystemExit(2)
            return x

        results, failures = run_tasks_with_optional_threads([1, 2], worker, max_workers=2)

        assert results == {1: 1}
        assert len(failures) == 1
        assert failures[0][0] == 2
        assert isinstance(failures[0][1], RuntimeError)
        assert 'Task requested exit with code 2.' in str(failures[0][1])

    def test_accepts_string_max_workers(self):
        results, failures = run_tasks_with_optional_threads([1, 2], lambda x: x * 2, max_workers='2')
        assert results == {1: 2, 2: 4}
        assert failures == []

    def test_auto_max_workers_falls_back_to_serial(self):
        results, failures = run_tasks_with_optional_threads([1, 2], lambda x: x + 1, max_workers='auto')
        assert results == {1: 2, 2: 3}
        assert failures == []

    def test_invalid_max_workers_raises_clear_error(self):
        with pytest.raises(ValueError, match='max_workers must be an integer'):
            run_tasks_with_optional_threads([1], lambda x: x, max_workers='two')

    def test_stop_on_failure_does_not_start_more_tasks_and_keeps_inflight_results(self):
        started = []
        started_lock = threading.Lock()
        both_started = threading.Barrier(2)

        def worker(x):
            with started_lock:
                started.append(x)
            both_started.wait(timeout=1)
            if x == 1:
                raise RuntimeError('stop')
            time.sleep(0.02)
            return x

        results, failures = run_tasks_with_optional_threads(
            [1, 2, 3, 4],
            worker,
            max_workers=2,
            stop_scheduling_on_failure=True,
        )

        assert len(failures) == 1
        assert failures[0][0] == 1
        assert set(started).issubset({1, 2})
        assert results == {2: 2}

    def test_fail_fast_remains_a_backwards_compatible_alias(self):
        results, failures = run_tasks_with_optional_threads(
            [1, 2, 3],
            lambda value: (_ for _ in ()).throw(RuntimeError('stop')) if value == 1 else value,
            max_workers=1,
            fail_fast=True,
        )

        assert results == {}
        assert [task for task, _exc in failures] == [1]

    def test_parallel_failure_order_follows_input_order(self):
        def worker(value):
            if value == 1:
                time.sleep(0.02)
            raise ValueError('bad {}'.format(value))

        _results, failures = run_tasks_with_optional_threads(
            [1, 2],
            worker,
            max_workers=2,
        )

        assert [task for task, _exc in failures] == [1, 2]

    def test_summary_error_retains_original_exceptions_as_cause(self):
        original = ValueError('root cause')

        with pytest.raises(RuntimeError, match='workflow failed') as exc_info:
            raise_task_failures([('run-1', original)], 'workflow failed')

        assert exc_info.value.failures == (('run-1', original),)
        assert isinstance(exc_info.value.__cause__, ExceptionGroup)
        assert exc_info.value.__cause__.exceptions == (original,)
        assert "Parallel task: 'run-1'" in original.__notes__

class TestValidatePositiveIntOption:
    def test_accepts_positive(self):
        assert validate_positive_int_option(3, 'internal_jobs') == 3
        assert validate_positive_int_option('2', 'internal_jobs') == 2

    def test_rejects_nonpositive(self):
        with pytest.raises(ValueError, match='--internal_jobs must be > 0'):
            validate_positive_int_option(0, 'internal_jobs')


class TestCpuBudgetHelpers:
    def test_resolve_cpu_budget_auto_uses_os_cpu_count(self, monkeypatch):
        monkeypatch.setattr('amalgkit.util.os.cpu_count', lambda: 12)
        assert resolve_cpu_budget(0) == 12

    def test_resolve_cpu_budget_auto_fallbacks_to_one(self, monkeypatch):
        monkeypatch.setattr('amalgkit.util.os.cpu_count', lambda: None)
        assert resolve_cpu_budget(0) == 1

    def test_resolve_cpu_budget_rejects_negative(self):
        with pytest.raises(ValueError, match='--internal_cpu_budget must be >= 0'):
            resolve_cpu_budget(-1)

    def test_resolve_thread_worker_allocation_caps_workers(self):
        threads, workers, budget = resolve_thread_worker_allocation(
            requested_threads=4,
            requested_workers=4,
            internal_cpu_budget=8,
            worker_option_name='internal_jobs',
        )
        assert threads == 1
        assert workers == 4
        assert budget == 4

    def test_resolve_thread_worker_allocation_caps_threads(self):
        threads, workers, budget = resolve_thread_worker_allocation(
            requested_threads=16,
            requested_workers=2,
            internal_cpu_budget=8,
            worker_option_name='internal_jobs',
        )
        assert threads == 4
        assert workers == 2
        assert budget == 8

    def test_resolve_thread_worker_allocation_respects_total_core_budget(self):
        threads, workers, budget = resolve_thread_worker_allocation(
            requested_threads=8,
            requested_workers='auto',
            internal_cpu_budget=64,
            worker_option_name='internal_jobs',
        )
        assert threads == 1
        assert workers == 8
        assert budget == 8

    def test_resolve_worker_allocation_caps_workers(self):
        workers, budget = resolve_worker_allocation(
            requested_workers=10,
            internal_cpu_budget=3,
            worker_option_name='internal_jobs',
        )
        assert workers == 3
        assert budget == 3

class TestFindPrefixedEntries:
    def test_set_entries_returns_sorted_output(self):
        entries = {'Homo_sapiens_b.idx', 'Mus_musculus.idx', 'Homo_sapiens_a.idx'}
        out = find_prefixed_entries(entries, 'Homo_sapiens')
        assert out == ['Homo_sapiens_a.idx', 'Homo_sapiens_b.idx']

    def test_no_match(self):
        entries = ['Arabidopsis_thaliana.fa']
        out = find_prefixed_entries(entries, 'Homo_sapiens')
        assert out == []

class TestFindSpeciesPrefixedEntries:
    def test_rejects_similar_species_prefix(self):
        entries = ['Homo_sapiens.fa', 'Homo_sapiens2.fa', 'Homo_sapiens_k31.idx']
        out = find_species_prefixed_entries(entries, 'Homo_sapiens')
        assert out == ['Homo_sapiens.fa', 'Homo_sapiens_k31.idx']


class TestFindRunPrefixedEntries:
    def test_rejects_similar_run_prefix(self):
        entries = ['SRR001.fastq.gz', 'SRR0010.fastq.gz', 'SRR001-legacy.fastq.gz', 'SRR001_1.fastq.gz']
        out = find_run_prefixed_entries(entries, 'SRR001')
        assert out == ['SRR001.fastq.gz', 'SRR001_1.fastq.gz']

# ---------------------------------------------------------------------------
# Metadata class
# ---------------------------------------------------------------------------
