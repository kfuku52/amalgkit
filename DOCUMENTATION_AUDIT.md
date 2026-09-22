# Documentation audit — 2026-09-22

## Scope and baseline

Audited `master` at `bd7e041d078d863a313d3f070eff17b5c7a7ad81`
(version 0.16.92), initially clean and equal to fetched `origin/master`.
Followed `AGENTS.md`, CONTRIBUTING, architecture and the CLI-documentation and
push skills. The original documentation-only delivery incremented the patch to
0.16.93 as required for a push, without changing pipeline behavior or dependencies. Historical release notes
and previous changelog entries are preserved.

Prioritized installation, README quick starts, canonical `.wiki/` command guides,
yeast tutorial, generated workspace guide, parser/help, private-input handoff,
selection, quant/merge outputs, normalization and finalize interpretation.
There was no `docs/` directory. Read handlers and existing behavioral tests as
well as CLI definitions; the static checker alone was not treated as execution.

## Corrected documentation findings (A)

| Location and previous gap | Correction | Implementation and test evidence |
| --- | --- | --- |
| README and Installation: bare pip/mamba installation into the active environment | Isolated venv or named Conda environment; external executables still require separate installation | `pyproject.toml`; fresh editable installation and `pip check` below; Conda resolution was not executed |
| Installation: native default taxonomy location and environment precedence omitted | Document XDG native/legacy precedence and command download-directory cache | `download_utils.resolve_default_ncbi_taxonomy_data_dir`, `get_ncbi_taxonomy`, `resolve_download_dir`; taxonomy unit tests |
| Select: location of threshold/strategy settings unclear | Distinguish TSV parameters from CLI flags and environment; flag seed discrepancy separately as B | `cli_parser.py`, `SELECT_PARAMETER_DEFINITIONS`, `apply_select_config_parameters`; CLI/parser and documented workflow tests |
| Tutorial/getfastq: bare `getfastq_stats.tsv` path and undefined percentage denominator | Name `getfastq/<RUN>/getfastq_stats.tsv`, per-stage base denominator, missing values and seconds | `write_getfastq_stats`, `calculate_filtered_percent`; `test_merge.py::TestMergeFastpStatsIntoMetadata` |
| Quant: parent output directory and literal column names absent; Oarfish reuse restriction omitted | Name `quant/<RUN>/`, abundance columns/units, and redo requirement for changed or unknown Oarfish options | `quant_main`, `check_fragment_model_reuse`, `output_contracts.py`; `test_quant_workflow.py::test_oarfish_reuse_checks_resolved_settings_after_fastq_cleanup` |
| Merge: matrix axes and incomplete-run behavior omitted | Describe target/run axes, fractional counts, and missing-file omission; separately flag selection discrepancy | `collect_species_quant_outputs`, `write_species_merged_quant_tables`; merge tests and synthetic CLI run below |
| Finalize: species parent directory and expression/group-mean scale unspecified | Name `finalize/<Species>/`, selected expression scale, and link to existing linear tau definition | `_copy_species_tables`, `per_species_finalize_python.py`, `sample_group_mean`; documented long-read normalization workflow |

## Implementation findings (B), resolved in 0.16.94

The following reproductions describe the original audit, before the fixes.
Version 0.16.94 resolves all three: merge honors populated sampling flags;
explicit CLI seeds override rule-file seeds, with zero only as the final fallback;
and missing layout data skips the plot without a dtype assignment error.
`tests/test_documented_workflows.py` now checks saved seeds in regular and
species-wise selection and invokes the real merge CLI with stale unselected
quant output and sparse metadata. `tests/test_merge.py` covers legacy/blank,
case-insensitive and invalid flags. Original audit evidence is retained below.

### B1: Reselected-out runs can enter merge

The select guide promises downstream exclusion of `is_sampled != yes` rows.
`merge.collect_species_runs` instead filters only species and `exclusion=no`;
its local variable named `is_sampled` does not read that column. Existing tests
cover `exclusion=yes` rejection, but do not establish `is_sampled=no` rejection.

Reproduced using a temporary yeast workspace with metadata run `0001`,
`exclusion=no`, `is_sampled=no`, and valid synthetic quant output. Running
`amalgkit merge --out_dir ./demo` exited 0 and wrote the `0001` column into
`demo/merge/Saccharomyces_cerevisiae/Saccharomyces_cerevisiae_est_counts.tsv`.
Targets `0001` and `NA` remained distinct lexical IDs and counts were 2 and 3.
This can matter when quant results survive reselection. Intended selection
semantics were not redefined to legitimize this behavior; the guides now warn
about it. Fixing the producer/consumer contract requires a separate code change.

### B2: Rule-file random seed is silently shadowed

`SELECT_PARAMETER_DEFINITIONS` accepts `random_seed` as a nonnegative integer.
However, the parser supplies 0 and `apply_select_config_parameters` preserves
any non-None runtime value, including this implicit default.

Reproduced by adding an enabled `stage=parameter` row with
`parameter_name=random_seed`, `parameter_value=17` to yeast rules, then running
`amalgkit select --out_dir ./demo` on synthetic metadata. The selected metadata
reported `sampling_seed=0`. Repeating with `--random_seed 17` reported 17.
An explicit CLI value is a verified workaround. Whether rule-file seeds should
be supported or rejected is unresolved; changing seed semantics is outside this
documentation task. No test was added to enshrine the questionable behavior.

### B3: Missing library layout crashes merge plotting

An initial version of the B1 fixture omitted `lib_layout`, `total_bases` and
`spot_length`. Metadata loading/merge accepted it and staged abundance tables,
but `merge` exited 1 with `Invalid value ... for dtype 'float64'` while preparing
the library-layout plot. In `merge_plots.py`, `layout_df.loc[:, 'lib_layout']`
assigns strings into the all-missing numeric column (pandas 3.0.6).
`merge_main` validates only run, scientific name and exclusion as required
columns. The missing-layout branch appears intended to skip the plot, rather
than fail with a dtype error. No published merge directory was left behind.

Adding valid `lib_layout=single`, `total_bases=100000000`, `spot_length=50` to
the synthetic metadata allowed the B1 run to complete. This distinguishes
fixture completeness from proof that the missing-value path works. The contract
for absent plotting metadata needs a code-level resolution, not a new invented
mandatory-field requirement in documentation.

## Scientific questions (C) and coverage limits

No new scientific specification was established. Whether yeast genotype labels
are comparable across species, sampling thresholds suit real studies, or an
assumed fragment distribution is valid remains dataset-dependent. Existing
warnings and scientific defaults were retained. Fixture agreement verifies
file handoffs and stated transformations, not biological validity.

Not exhaustively audited: every function docstring, every BUSCO/MMseqs option,
provider fallback permutation, scheduler example, historical output migration,
or numerical derivation. No live NCBI/GSA research downloads, taxonomy download,
large reference builds, paid services or external-service writes were used for
validation. Remote HTTP links and Conda package availability were not checked;
local repository/Wiki targets and anchors were checked offline. Public Wiki
publication is separate from the source-repository push.

## Execution evidence

Fresh environment: `/tmp/amalgkit-doc-audit-20260922`, macOS arm64,
Python 3.14.7. Existing `.venv` was not modified. Installation used the documented
uv alternative with that interpreter, then `uv pip install --python
/tmp/amalgkit-doc-audit-20260922/bin/python -e '.[test,quality]' pip` (exit 0,
79 packages). This verifies local editable installation, not the literal remote
Git URL or Bioconda commands. `python -m pip check` exited 0.
All pytest/quality runs used `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
MKL_NUM_THREADS=1`.

- Literal README commands `amalgkit dataset --name init --out_dir ./work` and
  `amalgkit dataset --name yeast --out_dir ./demo` ran in a fresh temporary cwd,
  exited 0, and produced the expected template headers, two FASTAs and two BUSCO
  tables. No reference data were fetched.
- `/tmp/amalgkit-doc-reproduce.py` executed the synthetic select/merge cases
  above, checking saved metadata seeds, matrix columns, IDs/counts and output
  paths. Its first sparse-input run failed as B3; the complete fixture exited 0.
  The successful workspace was `amalgkit-doc-examples-tsrk_v0r` under the system
  temporary directory. These are reduced substitute inputs, not a live tutorial.
- `python -m pytest -q tests/test_cli_parser.py tests/test_cli_help.py
  tests/test_documented_workflows.py tests/test_doc_tools.py`: 115 passed,
  exit 0. The workflows extract documented commands, use small synthetic inputs,
  fixture network/taxonomy responses and lightweight plot placeholders. The
  private-input case checks metadata handoff, not a real external quantifier.
- `python .github/scripts/check_quality.py`: exit 0; Ruff lint/format and mypy
  passed; offline documentation check covered 40 files, 182 commands, 31 literal
  defaults and 136 links with zero errors. These are static checks.
- `python -m pytest -q -n 2 --cov=amalgkit --cov-branch --cov-fail-under=75 -rs`:
  exit 0, 2228 passed, 11 skipped, combined statement/branch coverage 82.82%.
  Skips: 3 require inmoose, 2 require kallisto, 4 require fastp/SeqKit/fasterq-dump,
  and 2 require Oarfish. Those real-tool/optional paths were not verified.
- `python -m build` and `python -m twine check dist/*`: exit 0; built
  `amalgkit-0.16.93.tar.gz` and `amalgkit-0.16.93-py3-none-any.whl`; all checked
  distributions passed. Existing ignored distributions were retained.
- `python -m amalgkit --version` and `python -m amalgkit help quant`: exit 0,
  version 0.16.93 and current quant help. No external tools or data retrieval.

Additional package checks (both exit 0): installed the built wheel into a second
fresh environment, then ran `python -I .github/scripts/check_wheel.py` using
that environment's interpreter and an absolute script path, with cwd `/tmp`.
It verified installed-package import, bundled FASTA extraction and synthetic
yeast selection outside the source checkout. Extracted the sdist in a temporary
directory and ran its `check_docs.py`: 40 files, zero errors.

The audit report is included in the source distribution so the command guides'
links remain valid in source-distribution documentation checks.

## Fix verification — 0.16.94

The follow-up fixes used the same isolated Python 3.14.7 environment and native
thread limits as the audit. Focused CLI/help/documentation/merge/selection checks
passed (198 tests). The complete delivery command
`python -m pytest -q -n 2 --cov=amalgkit --cov-branch --cov-fail-under=75 -rs`
passed with 2241 tests, 11 skips and 82.83% coverage. The skip reasons remain the
missing optional inmoose and external tools listed above. The real merge CLI
regression requires neither external tools nor remote data: it verifies the
published counts, lexical IDs, metadata and PDF output on synthetic inputs.
