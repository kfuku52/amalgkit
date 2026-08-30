"""Use one boundary catalog for local and CI formatting/type checks."""

import subprocess
import sys
import tomllib
from pathlib import Path

root = Path(__file__).resolve().parents[2]
with (root / 'pyproject.toml').open('rb') as handle:
    boundaries = tomllib.load(handle)['tool']['mypy']['files']
for arguments in (['ruff', 'check', '.'], ['ruff', 'format', '--check', *boundaries], ['mypy']):
    subprocess.run(  # noqa: S603 - fixed QA modules and repository-owned paths
        [sys.executable, '-m', *arguments], cwd=root, check=True,
    )
