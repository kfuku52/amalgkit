# Contributing

AMALGKIT supports Python 3.11 through 3.14. Create an isolated environment and
install the package with its test and quality extras:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip "setuptools>=83"
python -m pip install -e ".[test,quality]"
```

Before opening a pull request, run the same checks used by CI:

```bash
ruff check .
mypy
python -m pytest -q -n 2 --cov=amalgkit --cov-branch --cov-fail-under=75
python -m build
```

For the shortest development feedback loop, skip integration tests, deliberate
wall-clock waits, real PDF rendering, and optional-dependency coverage:

```bash
python -m pytest -q -n 2 -m "not integration and not slow and not optional_dependency"
```

Run `python -m pytest -q` when debugging order or process-isolation issues in a
single worker. Tests without a marker are expected to be deterministic and
fast. Use `integration` for public workflows spanning multiple production
components, `slow` for intentional process waits or real rendering, and
`optional_dependency` for tests that require an optional project extra.

Add a regression test for bug fixes. Keep external bioinformatics tools out of
unit tests by injecting or mocking command runners and network clients.

See [ARCHITECTURE.md](ARCHITECTURE.md) before moving responsibilities between
command modules. Migrated boundary modules are formatted and typed; add a newly
migrated module to the CI format/mypy lists only after those checks pass. For a
performance-sensitive change, run:

```bash
python benchmarks/benchmark_core.py --quick --output benchmark-results.json
```

Package metadata and dependency floors live in `pyproject.toml`. When changing
the supported Python or dependency range, update the CI matrix and README in the
same pull request.

## Releases

The package version is defined once in `amalgkit/__init__.py`. After the full
test workflow succeeds for a push to `master`, the release workflow creates the
matching annotated Git tag and GitHub Release only when the patch component is
zero (for example, `0.17.0` or `1.0.0`). Eligible releases include both wheel
and source distributions. Patch-only versions remain available from the default
branch without a tag or GitHub Release. Do not manually tag patch-only versions,
or reuse or move an existing version tag.

The Bioconda recipe follows these major and minor release tags. Do not update
the Bioconda recipe for patch-only versions; patch fixes remain available from
the default branch until the next major or minor release.
