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

Run the fast lane while developing:

```bash
python -m pytest -q -n 2 -m "not integration and not slow and not optional_dependency"
```

Run every available test before merging:

```bash
python -m pytest -q -n 2 --cov=amalgkit --cov-branch --cov-fail-under=75
```

Warnings are errors unless `pyproject.toml` explicitly allows a known scientific
or fallback warning. Tests that introduce a new expected warning should assert
it locally or document a narrowly matched allow-list entry.

Integration tests use lightweight PDF placeholders when they are checking plot
orchestration and output naming. Dedicated `slow` tests retain real PDF
rendering coverage.

`test_documented_workflows.py` executes the Wiki's yeast metadata edit and
selection, private FASTQ metadata handoffs, the generated species-wise guide,
and the long-read CSTMM/filter/finalize chain. Network responses and taxonomy
lookups use small fixtures; Oarfish count/model tables are fixtures rather than
an invocation of the external quantifier. `test_doc_tools.py` exercises drift
detection and safe Wiki staging. Required documentation and scripts are included
in the source distribution so these tests also run from an unpacked sdist.
