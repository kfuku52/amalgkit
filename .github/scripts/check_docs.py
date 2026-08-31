"""Check documentation against the local CLI and source tree, without network I/O.

This checks syntax, literal defaults and local links, not biological validity or
external tools. Small executable workflow examples live in the test suite.
"""

import argparse
import contextlib
import io
import re
import shlex
import sys
from pathlib import Path
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from amalgkit.main import build_main_parser  # noqa: E402


def shell_commands(text):
    for block in re.finditer(r'^```(?:bash|sh|shell)\s*\n(.*?)^```', text, re.M | re.S):
        first_line = text.count('\n', 0, block.start(1)) + 1
        pending = ''
        start_line = first_line
        for offset, line in enumerate(block.group(1).splitlines()):
            if not pending:
                start_line = first_line + offset
            continued = line.rstrip().endswith('\\')
            pending += (line.rstrip()[:-1] if continued else line).strip() + ' '
            if continued:
                continue
            if pending.lstrip().startswith('amalgkit '):
                yield start_line, pending.strip()
            pending = ''


def local_link_target(path, url):
    parsed = urlsplit(url)
    location = unquote(parsed.path)
    wiki_prefix = '/kfuku52/amalgkit/wiki'
    if parsed.netloc == 'github.com' and (location == wiki_prefix or location.startswith(wiki_prefix + '/')):
        page = location[len(wiki_prefix):].strip('/') or 'Home'
        return ROOT / '.wiki' / (page + '.md'), unquote(parsed.fragment)
    repo_prefix = '/kfuku52/amalgkit/'
    if parsed.netloc == 'github.com' and location.startswith(repo_prefix):
        for kind in ('blob/master/', 'tree/master/'):
            if location.startswith(repo_prefix + kind):
                return ROOT / location[len(repo_prefix + kind):], unquote(parsed.fragment)
    if parsed.scheme or parsed.netloc:
        return None
    target = path.parent / location if location else path
    if path.parent == ROOT / '.wiki' and location and not target.exists():
        target = Path(str(target) + '.md')
    return target, unquote(parsed.fragment)


def heading_anchors(text):
    used = {}
    anchors = set()
    for heading in re.findall(r'^#{1,6}\s+(.+?)\s*#*$', text, re.M):
        heading = re.sub(r'\[([^]]+)\]\([^)]+\)', r'\1', heading)
        slug = re.sub(r'[^\w\-\s]', '', heading.lower()).replace(' ', '-')
        suffix = used.get(slug, 0)
        used[slug] = suffix + 1
        anchors.add(slug + ('-' + str(suffix) if suffix else ''))
    return anchors


def check_documents(paths):
    parser = build_main_parser()
    commands = next(action.choices for action in parser._actions if isinstance(action, argparse._SubParsersAction))
    errors = []
    counts = {'commands': 0, 'defaults': 0, 'links': 0}
    substitutions = {'$SLURM_CPUS_PER_TASK': '2', '$SLURM_ARRAY_TASK_ID': '1'}
    for path in paths:
        text = path.read_text(encoding='utf-8')
        historical = path.name.startswith('Release-notes-')
        for line, command in ([] if historical else shell_commands(text)):
            capture = io.StringIO()
            try:
                args = [substitutions.get(value, value) for value in shlex.split(command)[1:]]
                with contextlib.redirect_stdout(capture), contextlib.redirect_stderr(capture):
                    parser.parse_args(args)
            except SystemExit as exc:
                if exc.code:
                    errors.append(f'{path.relative_to(ROOT)}:{line}: {capture.getvalue().strip()}')
            except ValueError as exc:
                errors.append(f'{path.relative_to(ROOT)}:{line}: {exc}')
            counts['commands'] += 1
        command_name = path.stem.removeprefix('amalgkit-')
        if command_name in commands and not historical:
            options = {option: action for action in commands[command_name]._actions for option in action.option_strings}
            default_column = None
            for line_number, line in enumerate(text.splitlines(), 1):
                if not line.startswith('|'):
                    default_column = None
                    continue
                cells = [cell.strip() for cell in line.strip('|').split('|')]
                if 'Default' in cells:
                    default_column = cells.index('Default')
                    continue
                if default_column is None or len(cells) <= default_column:
                    continue
                option = re.fullmatch(r'`(--[\w_]+)`', cells[0])
                value = re.fullmatch(r'`([\w.\-]+)`', cells[default_column])
                if option and value and option[1] in options:
                    default = options[option[1]].default
                    if str(default) != value[1]:
                        errors.append(f'{path.relative_to(ROOT)}:{line_number}: {option[1]} default '
                                      f'is {default!r}, documented as {value[1]!r}')
                    counts['defaults'] += 1
        for match in re.finditer(r'\[[^\]\n]*\]\(([^\s)]+)\)', text):
            resolved = local_link_target(path, match[1])
            if resolved is None:
                continue
            target, anchor = resolved
            counts['links'] += 1
            line = text.count('\n', 0, match.start()) + 1
            if not target.exists():
                errors.append(f'{path.relative_to(ROOT)}:{line}: missing link target {match[1]}')
            elif anchor and target.suffix == '.md' and anchor not in heading_anchors(target.read_text(encoding='utf-8')):
                errors.append(f'{path.relative_to(ROOT)}:{line}: missing heading {match[1]}')
    return errors, counts


def main():
    paths = sorted(ROOT.glob('*.md')) + sorted((ROOT / '.wiki').glob('*.md'))
    paths += [ROOT / path for path in ('tests/README.md', 'tests/reference/README.md', 'benchmarks/README.md')]
    paths = [path for path in paths if path.name != 'AGENTS.md']
    errors, counts = check_documents(paths)
    for error in errors:
        print(error, file=sys.stderr)
    print(f'Documentation checked: {len(paths)} files, {counts}, {len(errors)} errors')
    return bool(errors)


if __name__ == '__main__':
    raise SystemExit(main())
