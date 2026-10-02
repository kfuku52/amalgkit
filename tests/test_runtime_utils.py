import pytest

from amalgkit.runtime_utils import (
    build_species_token_map,
    resolve_species_token,
    safe_join_existing_input_component,
    safe_join_component,
    validate_safe_path_component,
    validate_unique_species_tokens,
)


@pytest.mark.parametrize(
    'value',
    [
        '',
        '.',
        '..',
        '../escape',
        'nested/name',
        r'nested\name',
        '/absolute',
        'line\nbreak',
    ],
)
def test_validate_safe_path_component_rejects_unsafe_values(value):
    with pytest.raises(ValueError):
        validate_safe_path_component(value, label='test component')


def test_safe_join_component_stays_within_root(tmp_path):
    assert safe_join_component(str(tmp_path), 'SRR001') == str(tmp_path / 'SRR001')


def test_safe_join_component_rejects_symbolic_link_leaf(tmp_path):
    outside = tmp_path / 'outside'
    outside.mkdir()
    (tmp_path / 'SRR001').symlink_to(outside, target_is_directory=True)

    with pytest.raises(ValueError, match='symbolic-link'):
        safe_join_component(str(tmp_path), 'SRR001', label='run ID')


def test_safe_join_component_rejects_symbolic_link_root(tmp_path):
    real_root = tmp_path / 'real'
    real_root.mkdir()
    linked_root = tmp_path / 'linked'
    linked_root.symlink_to(real_root, target_is_directory=True)

    with pytest.raises(ValueError, match='symbolic-link output root'):
        safe_join_component(str(linked_root), 'SRR001', label='run ID')


def test_safe_join_existing_input_component_accepts_symbolic_link_root(tmp_path):
    real_root = tmp_path / 'real'
    run_dir = real_root / 'SRR001'
    run_dir.mkdir(parents=True)
    linked_root = tmp_path / 'linked'
    linked_root.symlink_to(real_root, target_is_directory=True)

    observed = safe_join_existing_input_component(
        str(linked_root),
        'SRR001',
        label='run ID',
    )

    assert observed == str(run_dir)


def test_safe_join_existing_input_component_rejects_symbolic_link_leaf(tmp_path):
    real_root = tmp_path / 'real'
    real_root.mkdir()
    outside = tmp_path / 'outside'
    outside.mkdir()
    (real_root / 'SRR001').symlink_to(outside, target_is_directory=True)

    with pytest.raises(ValueError, match='symbolic-link run ID'):
        safe_join_existing_input_component(
            str(real_root),
            'SRR001',
            label='run ID',
        )


def test_safe_join_existing_input_component_rejects_dangling_symbolic_link_root(tmp_path):
    missing_root = tmp_path / 'missing'
    linked_root = tmp_path / 'linked'
    linked_root.symlink_to(missing_root, target_is_directory=True)

    with pytest.raises(NotADirectoryError, match='existing directory'):
        safe_join_existing_input_component(
            str(linked_root),
            'SRR001',
            label='run ID',
        )


def test_explicit_species_token_must_be_canonical():
    with pytest.raises(ValueError, match='ASCII letters'):
        resolve_species_token('Species A', 'token (x)')


def test_build_species_token_map_rejects_conflicting_tokens_for_one_species():
    with pytest.raises(ValueError, match='conflicting species_token'):
        build_species_token_map(
            ['Homo sapiens', 'Homo sapiens'],
            ['human', 'Homo_sapiens'],
        )


@pytest.mark.parametrize(
    'pairs',
    [
        [('Alpha', None), ('alpha', None)],
        [('Species A', 'Token'), ('Species B', 'token')],
        [('Species A', '\u00e9'), ('Species B', 'e\u0301')],
    ],
)
def test_species_token_collision_detection_is_filesystem_conservative(pairs):
    with pytest.raises(ValueError):
        validate_unique_species_tokens(pairs)


@pytest.mark.parametrize('ids', [['R1', 'r1'], ['caf\u00e9', 'cafe\u0301']])
def test_metadata_consumers_reject_colliding_run_ids(tmp_path, ids):
    import pandas
    from types import SimpleNamespace
    from amalgkit.getfastq import collect_valid_run_ids
    from amalgkit.metadata_utils import Metadata, load_metadata
    from amalgkit.quant import build_quant_tasks
    from amalgkit.runtime_utils import validate_unique_run_ids
    frame = pandas.DataFrame({'run': ids, 'scientific_name': ['Species A'] * 2})
    path = tmp_path / 'metadata.tsv'
    frame.to_csv(path, sep='\t', index=False)
    for operation in [lambda: validate_unique_run_ids(ids), lambda: collect_valid_run_ids(ids),
                      lambda: load_metadata(SimpleNamespace(metadata=str(path))),
                      lambda: build_quant_tasks(Metadata.from_DataFrame(frame))]:
        with pytest.raises(ValueError, match='Colliding run IDs'):
            operation()


def test_run_directory_alias_is_rejected_before_output_reuse(tmp_path):
    from types import SimpleNamespace
    from amalgkit.runtime_utils import get_getfastq_run_dir
    root = tmp_path / 'getfastq'
    (root / 'R1').mkdir(parents=True)
    for operation in [lambda: safe_join_component(str(root), 'r1', label='run ID'),
                      lambda: safe_join_existing_input_component(str(root), 'r1', label='run ID'),
                      lambda: get_getfastq_run_dir(SimpleNamespace(out_dir=str(tmp_path)), 'r1')]:
        with pytest.raises(ValueError, match='Colliding run ID paths'):
            operation()


@pytest.mark.parametrize('ids', [['R1', 'r1'], ['caf\u00e9', 'cafe\u0301']])
def test_integrate_rejects_colliding_discovered_runs(ids):
    from amalgkit.integrate import assign_unique_run_ids
    with pytest.raises(ValueError, match='Colliding run IDs'):
        assign_unique_run_ids({((), run): {} for run in ids})
