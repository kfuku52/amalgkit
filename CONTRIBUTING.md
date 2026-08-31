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
python .github/scripts/check_quality.py
python -m pytest -q -n 2 --cov=amalgkit --cov-branch --cov-fail-under=75
python -m build
```

For the shortest development feedback loop, skip integration tests, deliberate
wall-clock waits, real PDF rendering, and optional-dependency coverage:

```bash
python -m pytest -q -n 2 -m "not integration and not slow and not optional_dependency"
```

CI installs into isolated environments with cached `uv` downloads/wheels. For
the same local install path, `uv venv --python 3.14` followed by
`uv pip install -e ".[test,quality]"` is an alternative to pip. This is a tooling
choice, not a new runtime dependency or a dependency lock.

When running several test processes, set `OMP_NUM_THREADS=1`,
`OPENBLAS_NUM_THREADS=1` and `MKL_NUM_THREADS=1` to avoid nested native thread
pools. CI uses these settings. Test the optional Combat-seq extra separately
with `python -m pytest -q -m optional_dependency` after installing
`.[combatseq,test]`; missing optional packages must not be mistaken for coverage.

Run `python -m pytest -q` when debugging order or process-isolation issues in a
single worker. Tests without a marker are expected to be deterministic and
fast. Use `integration` for public workflows spanning multiple production
components, `slow` for intentional process waits or real rendering, and
`optional_dependency` for tests that require an optional project extra.

Add a regression test for bug fixes. Keep external bioinformatics tools out of
unit tests by injecting or mocking command runners and network clients.

See [ARCHITECTURE.md](ARCHITECTURE.md) before moving responsibilities between
command modules. Migrated boundary modules are formatted and typed; add a newly
migrated module to `[tool.mypy].files` only after those checks pass. For a
performance-sensitive change, run:

```bash
python benchmarks/benchmark_core.py --quick --output benchmark-results.json
python benchmarks/benchmark_table_io.py --rows 1000000 --operation validate --output validation.json
python benchmarks/benchmark_table_io.py --rows 100000 --runs 8 --operation merge --output merge.json
```

Run each table benchmark in a separate, otherwise idle process for both
revisions. Compare the result hashes as well as wall time, Python allocation
peak and process peak RSS. The reported wall time includes tracemalloc overhead;
RSS includes interpreter/import allocations. Do not compare a cold concurrent
test run with a warm isolated run as a speedup claim.

The Python 3.14 full-suite job shares its install with lint, typing, dependency
audit and distribution checks; the remaining supported Python versions have
fast lanes. A separate minimum-dependency job still runs the full suite. macOS
filesystem/CLI coverage runs for relevant paths, weekly and on manual dispatch.
Nightly real-tool coverage includes plain/gzip and single/paired private FASTQs,
default cleanup, merge, TMM and finalize. The wheel smoke test runs outside the
checkout with isolated Python, extracts bundled FASTA data, and runs yeast
selection on synthetic metadata. Network-free documentation tests execute the
tutorial's metadata edit and selection, private FASTQ metadata handoffs, and the
generated species-wise workspace guide; NCBI/taxonomy responses are fixtures.

TMM reference inputs and independent Decimal/edgeR oracles are described in
[`tests/reference/README.md`](tests/reference/README.md).

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

## Documentation and public Wiki

`.wiki/` is the canonical source for the public GitHub Wiki. Edit it in this
repository, not directly in the Wiki UI. `check_quality.py` includes
`check_docs.py`, which checks documented CLI commands, literal default-value
tables, internal links and headings without HTTP requests. Historical release
notes are exempt from current CLI/default checks. These checks do not execute
external tools or prove the biological validity of examples.

After the source changes pass tests and are committed and pushed to `master`,
update a separate Wiki checkout. Review any public-only edits before applying:

```bash
git clone https://github.com/kfuku52/amalgkit.wiki.git /tmp/amalgkit-wiki
# For an existing checkout, pull --ff-only first.
git -C /tmp/amalgkit-wiki log -1 --format=%H
python .github/scripts/sync_wiki.py --wiki-dir /tmp/amalgkit-wiki
```

The check exits nonzero if pages or publication provenance differ. Use the full
Wiki commit printed above as `REVIEWED_WIKI_COMMIT` after reviewing the source
and public differences:

```bash
python .github/scripts/sync_wiki.py --wiki-dir /tmp/amalgkit-wiki --apply --expect-head REVIEWED_WIKI_COMMIT
git -C /tmp/amalgkit-wiki diff
git -C /tmp/amalgkit-wiki add --all
git -C /tmp/amalgkit-wiki commit -m "Sync verified default-branch documentation"
git -C /tmp/amalgkit-wiki push origin master
python .github/scripts/sync_wiki.py --wiki-dir /tmp/amalgkit-wiki
```

The script refuses dirty checkouts, an unexpected origin/head, symbolic-link
targets, and public-only pages; it never deletes pages, commits, or pushes.
Review public-only pages manually rather than discarding them. Do not force
push if another editor updates the Wiki during publication. The generated
`_Footer.md` identifies the source version and commit on every published page.
Wiki publication uses the maintainer's existing Git permissions and is a
separate completion step; the main repository CI does not publish it or claim
that a successful source build means the public Wiki is synchronized.
