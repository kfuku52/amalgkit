# Contributing

AMALGKIT supports Python 3.11 through 3.14. Create an isolated environment and
install the package with its test and quality extras:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e ".[test,quality]"
```

Before opening a pull request, run the same checks used by CI:

```bash
ruff check .
python -m pytest -q
python -m build
```

Add a regression test for bug fixes. Keep external bioinformatics tools out of
unit tests by injecting or mocking command runners and network clients.

Package metadata and dependency floors live in `pyproject.toml`. When changing
the supported Python or dependency range, update the CI matrix and README in the
same pull request.

## Releases

The package version is defined once in `amalgkit/__init__.py`. After the full
test workflow succeeds for a push to `master`, the release workflow creates the
matching annotated Git tag and GitHub Release if that version has not already
been released. Do not reuse or move an existing version tag.
