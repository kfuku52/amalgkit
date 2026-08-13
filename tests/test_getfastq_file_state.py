"""Direct tests for amalgkit.getfastq_file_state.

Note: the lane brief named this class ``GetfastqFileState``; the production
module exports it as ``RunFileState``.  These tests exercise the real public
name so they keep working if the class is ever renamed back.
"""

from amalgkit.getfastq_file_state import RunFileState, list_run_dir_files


def test_initializes_from_optional_iterable():
    state = RunFileState('/some/run', files=['a.fastq', 'b.fastq'])
    assert state.to_set() == {'a.fastq', 'b.fastq'}


def test_has_add_and_discard_update_inventory():
    state = RunFileState('/some/run', files=[])
    assert not state.has('c.fastq')
    state.add('c.fastq')
    assert state.has('c.fastq')
    state.discard('c.fastq')
    assert not state.has('c.fastq')
    # Discarding a missing name must not raise.
    state.discard('never-there.fastq')


def test_to_set_returns_a_copy_not_the_internal_set():
    state = RunFileState('/some/run', files=['a'])
    snapshot = state.to_set()
    snapshot.add('mutated')
    assert state.to_set() == {'a'}


def test_path_joins_work_dir_and_filename():
    state = RunFileState('/runs/sample1', files=[])
    assert state.path('x.fastq') == '/runs/sample1/x.fastq'


def test_list_run_dir_files_missing_dir_returns_empty_set(tmp_path):
    missing = tmp_path / 'does_not_exist'
    assert list_run_dir_files(str(missing)) == set()


def test_list_run_dir_files_returns_regular_files_only(tmp_path):
    (tmp_path / 'a.fastq').write_text('data')
    (tmp_path / 'b.fastq').write_text('data')
    (tmp_path / 'sub').mkdir()
    assert list_run_dir_files(str(tmp_path)) == {'a.fastq', 'b.fastq'}
