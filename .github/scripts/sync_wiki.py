"""Compare or stage the canonical .wiki pages in a reviewed Wiki checkout.

Never commits, pushes, deletes public-only pages, or uses credentials. Fetch/pull
and review the Wiki before --apply; publish with a normal non-force git push.
"""

import argparse
import ast
import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WIKI_REMOTES = {
    'https://github.com/kfuku52/amalgkit.wiki.git',
    'git@github.com:kfuku52/amalgkit.wiki.git',
    'ssh://git@github.com/kfuku52/amalgkit.wiki.git',
}


def git(directory, *arguments):
    executable = shutil.which('git')
    if executable is None:
        raise FileNotFoundError('git is required to inspect the source and Wiki checkouts.')
    return subprocess.run(  # noqa: S603 - resolved git executable, argument vector, no shell
        [executable, '-C', str(directory), *arguments], check=True, capture_output=True, text=True,
    ).stdout.strip()


def expected_pages(source=ROOT):
    version = None
    for node in ast.parse((source / 'amalgkit/__init__.py').read_text()).body:
        if isinstance(node, ast.Assign) and any(isinstance(target, ast.Name) and target.id == '__version__' for target in node.targets):
            version = ast.literal_eval(node.value)
    if not isinstance(version, str):
        raise ValueError('Cannot resolve the source project version.')
    commit = git(source, 'rev-parse', 'HEAD')
    pages = {}
    for path in sorted((source / '.wiki').glob('*.md')):
        if path.is_symlink() or not path.is_file():
            raise ValueError(f'Wiki source must be a regular file: {path}')
        pages[path.name] = path.read_bytes()
    if 'Home.md' not in pages:
        raise ValueError('Canonical Wiki must contain Home.md.')
    if '_Footer.md' in pages:
        raise ValueError('_Footer.md is generated from the published source version and commit.')
    pages['_Footer.md'] = (
        f'Documentation for **AMALGKIT {version}** · '
        f'[source `{commit[:7]}`](https://github.com/kfuku52/amalgkit/tree/{commit}/.wiki).\n\n'
        'These pages track the default branch, including patch updates. '
        'Historical release notes describe their named version.\n'
    ).encode('utf-8')
    return pages


def compare_pages(wiki_dir, pages):
    differences = []
    for name, content in pages.items():
        target = wiki_dir / name
        if target.is_symlink() or (target.exists() and not target.is_file()):
            raise ValueError(f'Wiki target must be a regular file: {target}')
        if not target.exists():
            differences.append(('missing', name))
        elif target.read_bytes() != content:
            differences.append(('changed', name))
    differences += [('public-only', path.name) for path in sorted(wiki_dir.glob('*.md')) if path.name not in pages]
    return differences


def synchronize(wiki_dir, *, source=ROOT, apply=False, expected_head=None):
    if wiki_dir.is_symlink() or not wiki_dir.is_dir():
        raise ValueError('Wiki checkout must be an existing directory, not a symbolic link.')
    pages = expected_pages(source)
    differences = compare_pages(wiki_dir, pages)
    if not apply:
        return differences
    if not expected_head or git(wiki_dir, 'rev-parse', 'HEAD') != expected_head:
        raise ValueError('--apply requires --expect-head matching the full, reviewed Wiki commit.')
    if git(wiki_dir, 'remote', 'get-url', 'origin') not in WIKI_REMOTES:
        raise ValueError('Refusing to modify a checkout whose origin is not the AMALGKIT Wiki.')
    if git(source, 'status', '--porcelain') or git(wiki_dir, 'status', '--porcelain'):
        raise ValueError('Commit source changes first; preserve/review changes in a dirty Wiki checkout before applying.')
    if any(status == 'public-only' for status, _ in differences):
        raise ValueError('Public-only Wiki pages need explicit review; no files were changed or deleted.')
    for _, name in differences:
        (wiki_dir / name).write_bytes(pages[name])
    return differences


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--wiki-dir', type=Path, required=True)
    parser.add_argument('--apply', action='store_true', help='Stage file changes only; never commit or push.')
    parser.add_argument('--expect-head', help='Full Wiki commit reviewed before applying changes.')
    args = parser.parse_args()
    try:
        differences = synchronize(args.wiki_dir, apply=args.apply, expected_head=args.expect_head)
    except (ValueError, OSError, subprocess.CalledProcessError) as exc:
        parser.exit(2, f'{exc}\n')
    for status, name in differences:
        print(f'{status}: {name}')
    if args.apply:
        print('Wiki files prepared. Review the git diff, commit and push normally, then rerun this script without --apply.')
    else:
        print('Wiki matches the source version/commit.' if not differences else 'Wiki is not synchronized.')
    return 0 if args.apply else bool(differences)


if __name__ == '__main__':
    raise SystemExit(main())
