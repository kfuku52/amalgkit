import os

import pytest

from amalgkit.filter_utils import staged_output_dir


def test_staged_output_dir_commits_new_output(tmp_path):
    target_dir = tmp_path / 'wsfilter'

    with staged_output_dir(str(target_dir), redo=False, prefix='test_stage_') as stage_dir:
        with open(os.path.join(stage_dir, 'result.txt'), 'w', encoding='utf-8') as handle:
            handle.write('new output')

    assert target_dir.is_dir()
    assert (target_dir / 'result.txt').read_text(encoding='utf-8') == 'new output'


def test_staged_output_dir_restores_existing_output_when_commit_fails(tmp_path, monkeypatch):
    target_dir = tmp_path / 'finalize'
    target_dir.mkdir()
    (target_dir / 'result.txt').write_text('old output', encoding='utf-8')
    real_rename = os.rename
    rename_calls = {'count': 0}

    def flaky_rename(src, dst):
        rename_calls['count'] += 1
        if rename_calls['count'] == 2:
            raise OSError('simulated commit failure')
        return real_rename(src, dst)

    monkeypatch.setattr('amalgkit.filter_utils.os.rename', flaky_rename)

    with pytest.raises(OSError, match='simulated commit failure'):
        with staged_output_dir(str(target_dir), redo=True, prefix='test_stage_') as stage_dir:
            with open(os.path.join(stage_dir, 'result.txt'), 'w', encoding='utf-8') as handle:
                handle.write('new output')

    assert target_dir.is_dir()
    assert (target_dir / 'result.txt').read_text(encoding='utf-8') == 'old output'
