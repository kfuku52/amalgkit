import os
import sys

from amalgkit.__init__ import __version__
from amalgkit.cli_parser import build_parser
from amalgkit.cli_utils import (
    build_timed_command_handler,
    print_runtime_banner,
    resolve_active_command,
    should_print_runtime_banner,
)
from amalgkit.exceptions import AmalgkitExit

SCRIPT_DIR = os.path.realpath(os.path.dirname(__file__))
SCRIPT_PARENT_DIR = os.path.realpath(os.path.dirname(SCRIPT_DIR))


def sanitize_sys_path():
    cleaned = []
    for path_entry in sys.path:
        real_entry = os.path.realpath(path_entry if path_entry != '' else os.getcwd())
        if real_entry == SCRIPT_DIR:
            continue
        cleaned.append(path_entry)
    sys.path[:] = cleaned
    if SCRIPT_PARENT_DIR not in [os.path.realpath(p if p != '' else os.getcwd()) for p in sys.path]:
        sys.path.insert(0, SCRIPT_PARENT_DIR)


COMMAND_SPECS = [
    ('metadata', 'amalgkit.metadata', 'metadata_main'),
    ('select', 'amalgkit.select', 'select_main'),
    ('getfastq', 'amalgkit.getfastq', 'getfastq_main'),
    ('quant', 'amalgkit.quant', 'quant_main'),
    ('merge', 'amalgkit.merge', 'merge_main'),
    ('busco', 'amalgkit.busco', 'busco_main'),
    ('cstmm', 'amalgkit.cstmm', 'cstmm_main'),
    ('wsfilter', 'amalgkit.wsfilter', 'wsfilter_main'),
    ('csfilter', 'amalgkit.csfilter', 'csfilter_main'),
    ('finalize', 'amalgkit.finalize', 'finalize_main'),
    ('sanity', 'amalgkit.sanity', 'sanity_main'),
    ('integrate', 'amalgkit.integrate', 'integrate_main'),
    ('config', 'amalgkit.config', 'config_main'),
    ('dataset', 'amalgkit.dataset', 'dataset_main'),
]


def build_main_parser():
    command_handlers = {
        name: build_timed_command_handler(name, module_name, function_name)
        for name, module_name, function_name in COMMAND_SPECS
    }
    command_names = list(command_handlers)
    return build_parser(
        command_handlers=command_handlers,
        command_names=command_names,
        version=__version__,
        prog='amalgkit',
    )


def main(argv=None):
    sanitize_sys_path()
    parser = build_main_parser()
    if argv is None:
        argv = sys.argv
    args = parser.parse_args(argv[1:])
    active_command = resolve_active_command(argv)
    try:
        if should_print_runtime_banner(active_command):
            print_runtime_banner(argv)
        if hasattr(args, 'handler'):
            args.handler(args)
        else:
            parser.print_help()
        return 0
    except AmalgkitExit as exc:
        if exc.message != '':
            if exc.use_stderr:
                sys.stderr.write(exc.message.rstrip('\n') + '\n')
            else:
                print(exc.message)
        return exc.exit_code
    except KeyboardInterrupt:
        sys.stderr.write('Interrupted.\n')
        return 130
    except Exception as exc:
        sys.stderr.write('ERROR: {}\n'.format(exc))
        return 1
