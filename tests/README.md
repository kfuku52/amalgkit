# Test suite

The default test layout follows production responsibilities. Large command
surfaces such as `getfastq`, `quant`, `select`, and utility helpers are split
across focused `test_<command>_<responsibility>.py` modules so failures can be
navigated without opening a monolithic test file.

Most tests are unmarked and belong to the fast unit-test lane. The registered
markers have narrower purposes:

- `integration`: a public workflow spanning multiple production components.
- `slow`: deliberate process waiting or real PDF rendering.
- `optional_dependency`: coverage requiring an optional project extra.
- `benchmark`: performance measurement kept out of the default correctness lane.

Environment setup, fast-lane and full delivery commands are maintained in
[CONTRIBUTING.md](../CONTRIBUTING.md). Run from the repository root with that
environment activated and its native-thread limits set.

## Choose checks for a change

Start with tests for the changed behavior, then include the file handoff or
consumer it affects. These are starting points, not substitutes for delivery
checks. Pass the listed paths/globs to `python -m pytest -q`; use `-k` only after
checking that it selects the intended tests.

| Changed surface | Focused checks and when to extend them |
| --- | --- |
| CLI arguments, help or examples | `tests/test_cli_parser.py tests/test_cli_help.py`; add `tests/test_documented_workflows.py` for example handoffs and run `python .github/scripts/check_docs.py` for syntax/default/link drift |
| A command's implementation | Its `tests/test_<command>*.py` files; include its workflow tests even though the fast lane excludes `integration` |
| Tables, identifiers, output validation or cleanup | `tests/test_output_contracts.py tests/test_workflow_contracts.py`; include affected producer/consumer tests such as `tests/test_quant_outputs.py tests/test_merge.py tests/test_util_metadata.py` |
| TMM, filtering or batch models | `tests/test_tmm_norm_factors.py tests/test_cstmm*.py`, relevant `tests/test_per_species*.py` / `tests/test_cross_species*.py` / `tests/test_batch_effect*.py`; preserve independent oracles in [reference/README.md](reference/README.md) and fixture provenance |
| External runner, download or resume | Relevant `tests/test_getfastq*.py`, `tests/test_quant*.py`, `tests/test_gsa*.py` and runner tests; real-tool coverage below is additional and needs installed binaries |
| Package metadata or bundled assets | `tests/test_dependency_floors.py tests/test_dataset.py`, distribution build and `twine check` from CONTRIBUTING.md; installed-wheel isolation is exercised by `.github/scripts/check_wheel.py` in CI |

A small, network-free workflow smoke test already exists:

```bash
python -m pytest -q tests/test_documented_workflows.py tests/test_workflow_contracts.py
```

It uses synthetic metadata/FASTQs, fixture network responses and pytest-managed
temporary directories. It exercises selection, private-input handoffs and
downstream output contracts; it does not validate external quantifiers or real
PDF rendering. Expect all selected tests to pass. Do not point examples at an
existing research workspace.

Do not append the fast-lane marker exclusion to a focused workflow or help
check: it would remove the very coverage being requested. Missing tools/extras
can yield skips even in a successful full suite. Use `-rs` to inspect their
reasons; install extras in a separate environment when needed. Nightly end-to-end
and benchmark commands are separate from this small fixture lane, and are not
prerequisites for ordinary documentation edits.

Warnings are errors unless `pyproject.toml` explicitly allows a known scientific
or fallback warning. Tests that introduce a new expected warning should assert
it locally or document a narrowly matched allow-list entry.

Integration tests use lightweight PDF placeholders when they are checking plot
orchestration and output naming. Dedicated `slow` tests retain real PDF
rendering coverage.

`test_quant_fragment_length.py` checks fragment sources, assumptions, validation,
provenance and safe reuse with mocked runners. The separate
`test_quant_fragment_length_integration.py` uses an installed kallisto (otherwise
skips) to compare wrapper output exactly with direct invocation and to check
mean/SD sensitivity on synthetic transcripts including short and shared sequences.
It tests numerical behavior, not a universally accurate prior for real libraries.

`test_documented_workflows.py` executes the Wiki's yeast metadata edit and
selection, private FASTQ metadata handoffs, the generated species-wise guide,
and the long-read CSTMM/filter/finalize chain. Network responses and taxonomy
lookups use small fixtures; Oarfish count/model tables are fixtures rather than
an invocation of the external quantifier. `test_doc_tools.py` exercises drift
detection and safe Wiki staging. Required documentation and scripts are included
in the source distribution so these tests also run from an unpacked sdist.

`test_real_tools_integration.py` runs actual fastp/SeqKit/fasterq-dump probes and
Oarfish on small local fixtures. It checks private source preservation and
resume for plain/gzip single/paired FASTQs, the post-filter mapping-rate
denominator, cleanup and Oarfish setting changes. These tests skip when tools
are absent locally; nightly CI requires the tools before running them.
