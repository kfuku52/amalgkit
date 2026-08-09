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

Run the fast lane while developing:

```bash
python -m pytest -q -n 2 -m "not integration and not slow and not optional_dependency"
```

Run every available test before merging:

```bash
python -m pytest -q -n 2
```

Integration tests use lightweight PDF placeholders when they are checking plot
orchestration and output naming. Dedicated `slow` tests retain real PDF
rendering coverage.
