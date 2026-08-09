import threading
import time

import pytest

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
        result = parse_bool_flags(['yes', None, ''], column_name='is_sampled', default='no')
        assert result.tolist() == [True, False, False]

    def test_raises_for_invalid_values(self):
        with pytest.raises(ValueError, match='invalid boolean flag'):
            parse_bool_flags(['yes', 'maybe'], column_name='is_sampled', default='no')

    def test_accepts_generator_input(self):
        values = (v for v in ['yes', 'no', None])
        result = parse_bool_flags(values, column_name='is_sampled', default='no')
        assert result.tolist() == [True, False, False]

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
        with pytest.raises(ValueError, match='--internal_jobs must be > 0'):
            validate_positive_int_option(-1, 'internal_jobs')

class TestFindPrefixedEntries:
    def test_sorted_list_entries(self):
        entries = ['Homo_sapiens.fa', 'Homo_sapiens.idx', 'Mus_musculus.fa']
        out = find_prefixed_entries(entries, 'Homo_sapiens')
        assert out == ['Homo_sapiens.fa', 'Homo_sapiens.idx']

    def test_set_entries_returns_sorted_output(self):
        entries = {'Homo_sapiens_b.idx', 'Mus_musculus.idx', 'Homo_sapiens_a.idx'}
        out = find_prefixed_entries(entries, 'Homo_sapiens')
        assert out == ['Homo_sapiens_a.idx', 'Homo_sapiens_b.idx']

    def test_no_match(self):
        entries = ['Arabidopsis_thaliana.fa']
        out = find_prefixed_entries(entries, 'Homo_sapiens')
        assert out == []

    def test_unsorted_list_entries(self):
        entries = ['Aardvark.fa', 'Homo_sapiens.fa', 'Mus_musculus.fa', 'Homo_sapiens.idx']
        out = find_prefixed_entries(entries, 'Homo_sapiens')
        assert out == ['Homo_sapiens.fa', 'Homo_sapiens.idx']


class TestFindSpeciesPrefixedEntries:
    def test_rejects_similar_species_prefix(self):
        entries = ['Homo_sapiens.fa', 'Homo_sapiens2.fa', 'Homo_sapiens_k31.idx']
        out = find_species_prefixed_entries(entries, 'Homo_sapiens')
        assert out == ['Homo_sapiens.fa', 'Homo_sapiens_k31.idx']


class TestFindRunPrefixedEntries:
    def test_rejects_similar_run_prefix(self):
        entries = ['SRR001.fastq.gz', 'SRR0010.fastq.gz', 'SRR001_1.fastq.gz']
        out = find_run_prefixed_entries(entries, 'SRR001')
        assert out == ['SRR001.fastq.gz', 'SRR001_1.fastq.gz']

    def test_rejects_hyphen_suffix_variants(self):
        entries = ['SRR001.fastq.gz', 'SRR001-legacy.fastq.gz', 'SRR001_1.fastq.gz']
        out = find_run_prefixed_entries(entries, 'SRR001')
        assert out == ['SRR001.fastq.gz', 'SRR001_1.fastq.gz']


# ---------------------------------------------------------------------------
# Metadata class
# ---------------------------------------------------------------------------
