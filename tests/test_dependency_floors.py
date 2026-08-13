import re
import tomllib
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]

PIN_PATTERN = re.compile(r"([A-Za-z0-9_.-]+)==(.+)")


def _normalize_name(value):
    return re.sub(r"[-_.]+", "-", value).lower()


def _parse_pin_file(rel_path):
    """Parse a ``name==version`` requirements file into {normalized_name: version}."""
    pins = {}
    for raw_line in (REPO_ROOT / rel_path).read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if line == "" or line.startswith("#"):
            continue
        match = PIN_PATTERN.fullmatch(line)
        assert match is not None, f"invalid pin line in {rel_path}: {line!r}"
        name, version = match.groups()
        pins[_normalize_name(name)] = version
    return pins


def _parse_declared_floors():
    project = tomllib.loads((REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    declared = {}
    for requirement in project["project"]["dependencies"]:
        match = re.fullmatch(r"([A-Za-z0-9_.-]+)>=(.+)", requirement)
        assert match is not None, requirement
        declared[_normalize_name(match.group(1))] = match.group(2)
    return declared


def test_minimum_requirements_match_runtime_dependency_floors():
    declared = _parse_declared_floors()
    pinned = _parse_pin_file("requirements/minimum.txt")

    assert pinned == declared


def test_minimum_constraints_file_is_valid_and_conflict_free():
    declared = _parse_declared_floors()
    pinned = _parse_pin_file("requirements/minimum.txt")
    constrained = _parse_pin_file("requirements/minimum-constraints.txt")

    # The constraints file exists to keep the floors installable by pinning
    # transitive dependencies; it must not be empty.
    assert constrained, "minimum-constraints.txt must not be empty"

    # A transitive override must not also be pinned as a direct floor with a
    # different version (that would make the two files contradict each other).
    assert constrained.keys().isdisjoint(pinned)

    # Overrides apply to transitive dependencies only, never to a declared
    # direct floor from pyproject.toml.
    assert constrained.keys().isdisjoint(declared)
