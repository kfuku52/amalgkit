import itertools
import os
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, as_completed, wait
from typing import Any

from amalgkit.logging_utils import get_logger


def _sort_failures_by_task_order(failures, task_order):
    order_by_task: dict[Any, int] = {}
    for index, task in enumerate(task_order):
        order_by_task.setdefault(task, index)
    return sorted(failures, key=lambda failure: order_by_task.get(failure[0], len(task_order)))


def _finish_task_batch(results, failures, tasks):
    ordered_failures = _sort_failures_by_task_order(failures, tasks)
    get_logger("parallel").debug(
        "task_batch_end",
        extra={
            "event": "task_batch_end",
            "task_count": len(tasks),
            "failure_count": len(ordered_failures),
        },
    )
    return results, ordered_failures


def raise_task_failures(failures, message, error_type=RuntimeError):
    """Raise a summary error while retaining every original task traceback.

    The user-facing exception remains the requested legacy type and message.
    Its cause is an ``ExceptionGroup`` containing the original exceptions, so
    debug output and programmatic callers can inspect every failed task.
    """
    failures = list(failures)
    if not failures:
        return
    original_exceptions = []
    for task, exc in failures:
        exc.add_note(f"Parallel task: {task!r}")
        original_exceptions.append(exc)
    grouped = ExceptionGroup(message, original_exceptions)
    summary_error = error_type(message)
    summary_error.failures = tuple(failures)
    raise summary_error from grouped


def run_tasks_with_optional_threads(
    task_items,
    task_fn,
    max_workers=1,
    fail_fast=False,
    stop_scheduling_on_failure=None,
):
    """Run tasks and return completed results plus failures.

    ``stop_scheduling_on_failure`` stops submitting new work after the first
    observed failure. Already-running threads cannot be terminated safely, so
    they are allowed to finish and their results or failures are retained.
    ``fail_fast`` is kept as a backwards-compatible alias for this behavior.
    """
    tasks = list(task_items)
    get_logger("parallel").debug(
        "task_batch_start",
        extra={"event": "task_batch_start", "task_count": len(tasks)},
    )
    results: dict[Any, Any] = {}
    failures: list[tuple[Any, Exception]] = []
    if stop_scheduling_on_failure is None:
        stop_after_failure = bool(fail_fast)
    else:
        stop_after_failure = bool(stop_scheduling_on_failure)
        if bool(fail_fast) and not stop_after_failure:
            raise ValueError("fail_fast=True conflicts with stop_scheduling_on_failure=False.")
    if len(tasks) == 0:
        return _finish_task_batch(results, failures, tasks)
    if max_workers is None:
        worker_limit = 1
    elif is_auto_parallel_option(max_workers):
        worker_limit = 1
    else:
        try:
            worker_limit = int(max_workers)
        except (TypeError, ValueError) as exc:
            raise ValueError('max_workers must be an integer >= 1, None, or "auto".') from exc
        if worker_limit <= 1:
            worker_limit = 1
    if (worker_limit <= 1) or (len(tasks) <= 1):
        for task in tasks:
            try:
                results[task] = task_fn(task)
            except SystemExit as exc:
                failures.append((task, RuntimeError(f"Task requested exit with code {exc.code}.")))
            except Exception as exc:
                failures.append((task, exc))
            if stop_after_failure and failures:
                break
        return _finish_task_batch(results, failures, tasks)
    worker_count = min(worker_limit, len(tasks))
    if stop_after_failure:
        task_iter = iter(tasks)
        with ThreadPoolExecutor(max_workers=worker_count) as executor:
            futures = {executor.submit(task_fn, task): task for task in itertools.islice(task_iter, worker_count)}
            stop_submitting = False
            while futures:
                completed, _pending = wait(futures, return_when=FIRST_COMPLETED)
                completed_without_failure = 0
                for future in completed:
                    task = futures.pop(future)
                    try:
                        results[task] = future.result()
                        completed_without_failure += 1
                    except SystemExit as exc:
                        failures.append((task, RuntimeError(f"Task requested exit with code {exc.code}.")))
                    except Exception as exc:
                        failures.append((task, exc))
                if failures and not stop_submitting:
                    stop_submitting = True
                    for future, task in list(futures.items()):
                        if future.cancel():
                            futures.pop(future)
                if not stop_submitting:
                    for task in itertools.islice(task_iter, completed_without_failure):
                        futures[executor.submit(task_fn, task)] = task
        return _finish_task_batch(results, failures, tasks)
    with ThreadPoolExecutor(max_workers=worker_count) as executor:
        futures = {executor.submit(task_fn, task): task for task in tasks}
        for future in as_completed(futures):
            task = futures[future]
            try:
                results[task] = future.result()
            except SystemExit as exc:
                failures.append((task, RuntimeError(f"Task requested exit with code {exc.code}.")))
            except Exception as exc:
                failures.append((task, exc))
    return _finish_task_batch(results, failures, tasks)


def validate_positive_int_option(value, option_name):
    if option_name in ("jobs", "species_jobs"):
        option_name = "internal_jobs"
    if is_auto_parallel_option(value):
        raise ValueError(f"--{option_name} must be > 0.")
    int_value = int(value)
    if int_value <= 0:
        raise ValueError(f"--{option_name} must be > 0.")
    return int_value


def is_auto_parallel_option(value):
    if value is None:
        return True
    if isinstance(value, str):
        return value.strip().lower() in ("", "auto")
    return False


def resolve_detected_cpu_count():
    detected = os.cpu_count()
    if (detected is None) or (detected <= 0):
        return 1
    return int(detected)


def resolve_cpu_budget(internal_cpu_budget="auto"):
    if is_auto_parallel_option(internal_cpu_budget):
        return resolve_detected_cpu_count()
    internal_cpu_budget = int(internal_cpu_budget)
    if internal_cpu_budget < 0:
        raise ValueError("--internal_cpu_budget must be >= 0.")
    if internal_cpu_budget == 0:
        return resolve_detected_cpu_count()
    return internal_cpu_budget


def resolve_total_core_budget(
    requested_threads="auto",
    internal_cpu_budget="auto",
    context="",
):
    if is_auto_parallel_option(requested_threads):
        requested_total = resolve_detected_cpu_count()
    else:
        requested_total = validate_positive_int_option(requested_threads, "threads")
    budget_cap = resolve_cpu_budget(internal_cpu_budget=internal_cpu_budget)
    budget = min(requested_total, budget_cap)
    if budget < requested_total:
        print(
            "{} reducing total cores from {} to {} to fit --internal_cpu_budget {}.".format(
                context if context else "CPU budget:",
                requested_total,
                budget,
                budget_cap,
            ),
            flush=True,
        )
    return budget


def resolve_thread_worker_allocation(
    requested_threads="auto",
    requested_workers="auto",
    internal_cpu_budget="auto",
    worker_option_name="internal_jobs",
    context="",
    disable_workers=False,
    task_count=None,
):
    budget = resolve_total_core_budget(
        requested_threads=requested_threads,
        internal_cpu_budget=internal_cpu_budget,
        context=context,
    )
    if disable_workers:
        if not is_auto_parallel_option(requested_workers):
            workers = validate_positive_int_option(requested_workers, worker_option_name)
            if workers != 1:
                print(
                    "{} --batch is set. Forcing --{} to 1.".format(
                        context if context else "Parallel:",
                        worker_option_name,
                    ),
                    flush=True,
                )
        effective_workers = 1
        effective_threads = budget
    else:
        if is_auto_parallel_option(requested_workers):
            workers = budget
        else:
            workers = validate_positive_int_option(requested_workers, worker_option_name)
        effective_workers = min(workers, budget)
        if effective_workers < workers:
            print(
                "{} reducing --{} from {} to {} to fit total core budget {}.".format(
                    context if context else "CPU budget:",
                    worker_option_name,
                    workers,
                    effective_workers,
                    budget,
                ),
                flush=True,
            )
        if task_count is not None:
            task_limit = max(1, int(task_count))
            if effective_workers > task_limit:
                print(
                    "{} reducing --{} from {} to {} to fit task count {}.".format(
                        context if context else "CPU budget:",
                        worker_option_name,
                        effective_workers,
                        task_limit,
                        task_limit,
                    ),
                    flush=True,
                )
                effective_workers = task_limit
        effective_threads = max(1, budget // effective_workers)
    print(
        "{} effective parallelism: {} x {} = {} core(s) max.".format(
            context if context else "CPU budget:",
            effective_workers,
            effective_threads,
            effective_workers * effective_threads,
        ),
        flush=True,
    )
    return effective_threads, effective_workers, budget


def resolve_worker_allocation(
    requested_workers="auto",
    requested_threads="auto",
    internal_cpu_budget="auto",
    worker_option_name="internal_jobs",
    context="",
    disable_workers=False,
):
    budget = resolve_total_core_budget(
        requested_threads=requested_threads,
        internal_cpu_budget=internal_cpu_budget,
        context=context,
    )
    if disable_workers:
        if not is_auto_parallel_option(requested_workers):
            workers = validate_positive_int_option(requested_workers, worker_option_name)
            if workers != 1:
                print(
                    "{} --batch is set. Forcing --{} to 1.".format(
                        context if context else "CPU budget:",
                        worker_option_name,
                    ),
                    flush=True,
                )
        effective_workers = 1
    else:
        if is_auto_parallel_option(requested_workers):
            workers = budget
        else:
            workers = validate_positive_int_option(requested_workers, worker_option_name)
        effective_workers = min(workers, budget)
        if effective_workers < workers:
            print(
                "{} reducing --{} from {} to {} to fit total core budget {}.".format(
                    context if context else "CPU budget:",
                    worker_option_name,
                    workers,
                    effective_workers,
                    budget,
                ),
                flush=True,
            )
    print(
        "{} effective parallel workers: {} (total core budget {}).".format(
            context if context else "CPU budget:",
            effective_workers,
            budget,
        ),
        flush=True,
    )
    return effective_workers, budget
