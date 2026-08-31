"""Check that documentation validation and publication safeguards catch drift."""

import runpy
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]


def test_documentation_checker_reports_bad_cli_default_and_links(tmp_path):
    tool = runpy.run_path(str(ROOT / '.github/scripts/check_docs.py'))
    check = tool['check_documents']
    check.__globals__['ROOT'] = tmp_path
    wiki = tmp_path / '.wiki'
    wiki.mkdir()
    page = wiki / 'amalgkit-finalize.md'
    page.write_text(
        '```bash\namalgkit finalize --no_such_option yes\n```\n'
        '| Option | Default |\n| --- | --- |\n| `--seed` | `auto` |\n'
        '[Missing](./Missing)\n[Bad heading](./Target#absent)\n')
    (wiki / 'Target.md').write_text('# Present\n')
    errors, counts = check([page])
    assert len(errors) == 4
    assert any('--seed default' in error for error in errors)
    assert counts == {'commands': 1, 'defaults': 1, 'links': 2}


@pytest.fixture
def wiki_checkouts(tmp_path):
    tool = runpy.run_path(str(ROOT / '.github/scripts/sync_wiki.py'))
    git = tool['git']
    source = tmp_path / 'source'
    wiki = tmp_path / 'published'
    for directory in (source, wiki):
        directory.mkdir()
        git(directory, 'init', '-q')
        git(directory, 'config', 'user.name', 'Docs test')
        git(directory, 'config', 'user.email', 'docs-test@example.invalid')
        git(directory, 'config', 'commit.gpgsign', 'false')
    (source / '.wiki').mkdir()
    (source / '.wiki' / 'Home.md').write_text('# Current documentation\n')
    (source / 'amalgkit').mkdir()
    (source / 'amalgkit' / '__init__.py').write_text('__version__ = "1.2.3"\n')
    (wiki / 'Home.md').write_text('# Old documentation\n')
    for directory in (source, wiki):
        git(directory, 'add', '.')
        git(directory, 'commit', '-qm', 'fixture')
    git(wiki, 'remote', 'add', 'origin', 'https://github.com/kfuku52/amalgkit.wiki.git')
    return tool, source, wiki, git(wiki, 'rev-parse', 'HEAD')


def test_wiki_check_is_read_only_and_apply_stamps_source(wiki_checkouts):
    tool, source, wiki, head = wiki_checkouts
    synchronize = tool['synchronize']
    assert synchronize(wiki, source=source) == [('changed', 'Home.md'), ('missing', '_Footer.md')]
    assert (wiki / 'Home.md').read_text() == '# Old documentation\n'
    assert not (wiki / '_Footer.md').exists()
    synchronize(wiki, source=source, apply=True, expected_head=head)
    assert not synchronize(wiki, source=source)
    footer = (wiki / '_Footer.md').read_text()
    assert '1.2.3' in footer
    assert tool['git'](source, 'rev-parse', 'HEAD') in footer


@pytest.mark.parametrize('problem', ['wrong_head', 'dirty_source', 'dirty_wiki', 'public_only', 'wrong_remote', 'symlink'])
def test_wiki_apply_preserves_unreviewed_content(wiki_checkouts, problem, tmp_path):
    tool, source, wiki, head = wiki_checkouts
    if problem == 'wrong_head':
        head = '0' * 40
    elif problem == 'dirty_source':
        (source / '.wiki/Home.md').write_text('uncommitted source')
    elif problem == 'dirty_wiki':
        (wiki / 'untracked.txt').write_text('unreviewed notes')
    elif problem == 'public_only':
        (wiki / 'Extra.md').write_text('public-only page')
        tool['git'](wiki, 'add', '.')
        tool['git'](wiki, 'commit', '-qm', 'extra page')
        head = tool['git'](wiki, 'rev-parse', 'HEAD')
    elif problem == 'wrong_remote':
        tool['git'](wiki, 'remote', 'set-url', 'origin', 'https://example.invalid/not-the-wiki.git')
    else:
        target = tmp_path / 'valuable.md'
        target.write_text('must not be overwritten')
        (wiki / 'Home.md').unlink()
        (wiki / 'Home.md').symlink_to(target)
    original = (wiki / 'Home.md').read_bytes()
    with pytest.raises(ValueError):
        tool['synchronize'](wiki, source=source, apply=True, expected_head=head)
    assert (wiki / 'Home.md').read_bytes() == original
    assert not (wiki / '_Footer.md').exists()
