import subprocess


def format_command(command):
    return ' '.join([str(part) for part in command])


def decode_command_output(output):
    if output is None:
        return ''
    if isinstance(output, str):
        return output
    return output.decode('utf8', errors='replace')


def print_command_output(stdout_txt, stderr_txt, stdout_label=None, stderr_label=None):
    if (stdout_label is None) and (stderr_label is None):
        print(stdout_txt, flush=True)
        print(stderr_txt, flush=True)
        return
    if stdout_label is not None:
        print(stdout_label, flush=True)
        print(stdout_txt, flush=True)
    elif stdout_txt != '':
        print(stdout_txt, flush=True)
    if stderr_label is not None:
        print(stderr_label, flush=True)
        print(stderr_txt, flush=True)
    elif stderr_txt != '':
        print(stderr_txt, flush=True)


def run_logged_command(
    command,
    runner=subprocess.run,
    print_command=True,
    command_prefix='Command',
    print_output=False,
    stdout_label=None,
    stderr_label=None,
    not_found_label=None,
):
    command_txt = format_command(command)
    if print_command:
        print('{}: {}'.format(command_prefix, command_txt), flush=True)
    try:
        result = runner(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except FileNotFoundError as exc:
        if not_found_label is None:
            raise
        raise FileNotFoundError(
            '{} executable not found: {}'.format(not_found_label, command[0])
        ) from exc
    stdout_txt = decode_command_output(result.stdout)
    stderr_txt = decode_command_output(result.stderr)
    if print_output:
        print_command_output(
            stdout_txt=stdout_txt,
            stderr_txt=stderr_txt,
            stdout_label=stdout_label,
            stderr_label=stderr_label,
        )
    return result, stdout_txt, stderr_txt


def run_checked_command(
    command,
    runner=subprocess.run,
    print_command=True,
    command_prefix='Command',
    print_output=False,
    stdout_label=None,
    stderr_label=None,
    not_found_label=None,
    failure_message=None,
):
    result, stdout_txt, stderr_txt = run_logged_command(
        command=command,
        runner=runner,
        print_command=print_command,
        command_prefix=command_prefix,
        print_output=print_output,
        stdout_label=stdout_label,
        stderr_label=stderr_label,
        not_found_label=not_found_label,
    )
    if result.returncode != 0:
        command_txt = format_command(command)
        if callable(failure_message):
            message = failure_message(result, stdout_txt, stderr_txt, command_txt)
        elif failure_message is not None:
            message = failure_message
        else:
            message = 'Command failed with exit code {}: {}'.format(result.returncode, command_txt)
        raise RuntimeError(message)
    return result, stdout_txt, stderr_txt


def run_logged_check_call(
    command,
    runner=subprocess.check_call,
    print_command=True,
    command_prefix='Command',
    not_found_label=None,
):
    command_txt = format_command(command)
    if print_command:
        print('{}: {}'.format(command_prefix, command_txt), flush=True)
    try:
        return runner(command)
    except FileNotFoundError as exc:
        if not_found_label is None:
            raise
        raise FileNotFoundError(
            '{} executable not found: {}'.format(not_found_label, command[0])
        ) from exc


def probe_dependency_command(command, label, runner=subprocess.run):
    return run_checked_command(
        command=command,
        runner=runner,
        print_command=False,
        print_output=False,
        not_found_label=label,
        failure_message=lambda result, _stdout, _stderr, command_txt: (
            '{} dependency probe failed with exit code {}: {}'.format(
                label,
                result.returncode,
                command_txt,
            )
        ),
    )
