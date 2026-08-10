import re
import tomllib
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


def _normalize_name(value):
    return re.sub(r"[-_.]+", "-", value).lower()


def test_minimum_requirements_match_runtime_dependency_floors():
    project = tomllib.loads((REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    declared = {}
    for requirement in project["project"]["dependencies"]:
        match = re.fullmatch(r"([A-Za-z0-9_.-]+)>=(.+)", requirement)
        assert match is not None, requirement
        declared[_normalize_name(match.group(1))] = match.group(2)

    pinned = {}
    for raw_line in (REPO_ROOT / "requirements" / "minimum.txt").read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if line == "" or line.startswith("#"):
            continue
        name, version = line.split("==", maxsplit=1)
        pinned[_normalize_name(name)] = version

    assert pinned == declared
