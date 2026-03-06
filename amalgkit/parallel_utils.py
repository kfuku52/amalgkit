import os
from concurrent.futures import ThreadPoolExecutor, as_completed


def run_tasks_with_optional_threads(task_items, task_fn, max_workers=1):
    tasks = list(task_items)
    results = dict()
    failures = list()
    if len(tasks) == 0:
        return results, failures
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
                failures.append((task, RuntimeError('Task requested exit with code {}.'.format(exc.code))))
            except Exception as exc:
                failures.append((task, exc))
        return results, failures
    worker_count = min(worker_limit, len(tasks))
    with ThreadPoolExecutor(max_workers=worker_count) as executor:
        futures = {
            executor.submit(task_fn, task): task
            for task in tasks
        }
        for future in as_completed(futures):
            task = futures[future]
            try:
                results[task] = future.result()
            except SystemExit as exc:
                failures.append((task, RuntimeError('Task requested exit with code {}.'.format(exc.code))))
            except Exception as exc:
                failures.append((task, exc))
    return results, failures


def validate_positive_int_option(value, option_name):
    if option_name in ('jobs', 'species_jobs'):
        option_name = 'internal_jobs'
    if is_auto_parallel_option(value):
        raise ValueError('--{} must be > 0.'.format(option_name))
    int_value = int(value)
    if int_value <= 0:
        raise ValueError('--{} must be > 0.'.format(option_name))
    return int_value


def is_auto_parallel_option(value):
    if value is None:
        return True
    if isinstance(value, str):
        return value.strip().lower() in ('', 'auto')
    return False


def resolve_detected_cpu_count():
    detected = os.cpu_count()
    if (detected is None) or (detected <= 0):
        return 1
    return int(detected)


def resolve_cpu_budget(internal_cpu_budget='auto'):
    if is_auto_parallel_option(internal_cpu_budget):
        return resolve_detected_cpu_count()
    internal_cpu_budget = int(internal_cpu_budget)
    if internal_cpu_budget < 0:
        raise ValueError('--internal_cpu_budget must be >= 0.')
    if internal_cpu_budget == 0:
        return resolve_detected_cpu_count()
    return internal_cpu_budget


def resolve_total_core_budget(
    requested_threads='auto',
    internal_cpu_budget='auto',
    context='',
):
    if is_auto_parallel_option(requested_threads):
        requested_total = resolve_detected_cpu_count()
    else:
        requested_total = validate_positive_int_option(requested_threads, 'threads')
    budget_cap = resolve_cpu_budget(internal_cpu_budget=internal_cpu_budget)
    budget = min(requested_total, budget_cap)
    if budget < requested_total:
        print(
            '{} reducing total cores from {} to {} to fit --internal_cpu_budget {}.'.format(
                context if context else 'CPU budget:',
                requested_total,
                budget,
                budget_cap,
            ),
            flush=True,
        )
    return budget


def resolve_thread_worker_allocation(
    requested_threads='auto',
    requested_workers='auto',
    internal_cpu_budget='auto',
    worker_option_name='internal_jobs',
    context='',
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
                    '{} --batch is set. Forcing --{} to 1.'.format(
                        context if context else 'Parallel:',
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
                '{} reducing --{} from {} to {} to fit total core budget {}.'.format(
                    context if context else 'CPU budget:',
                    worker_option_name,
                    workers,
                    effective_workers,
                    budget,
                ),
                flush=True,
            )
        effective_threads = max(1, budget // effective_workers)
    print(
        '{} effective parallelism: {} x {} = {} core(s) max.'.format(
            context if context else 'CPU budget:',
            effective_workers,
            effective_threads,
            effective_workers * effective_threads,
        ),
        flush=True,
    )
    return effective_threads, effective_workers, budget


def resolve_worker_allocation(
    requested_workers='auto',
    requested_threads='auto',
    internal_cpu_budget='auto',
    worker_option_name='internal_jobs',
    context='',
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
                    '{} --batch is set. Forcing --{} to 1.'.format(
                        context if context else 'CPU budget:',
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
                '{} reducing --{} from {} to {} to fit total core budget {}.'.format(
                    context if context else 'CPU budget:',
                    worker_option_name,
                    workers,
                    effective_workers,
                    budget,
                ),
                flush=True,
            )
    print(
        '{} effective parallel workers: {} (total core budget {}).'.format(
            context if context else 'CPU budget:',
            effective_workers,
            budget,
        ),
        flush=True,
    )
    return effective_workers, budget
