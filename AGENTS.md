# Working in AMALGKIT

## Start here

Read [README.md](README.md) for the pipeline and [CONTRIBUTING.md](CONTRIBUTING.md)
for environment setup and delivery checks. Read [ARCHITECTURE.md](ARCHITECTURE.md)
before changing command boundaries or file contracts. Use the relevant local
`.wiki/amalgkit-<command>.md` for command behavior; `.wiki/` is the canonical
Wiki source. No remote Wiki lookup is needed to find these instructions.

The CLI enters through `amalgkit/cli_entry.py` / `amalgkit/__main__.py`, then
`main.py` and `cli_parser.py`; handlers dispatch to command modules. Shared
boundaries and their compatibility rules are mapped in ARCHITECTURE.md.

## Execute and verify

Run commands from the repository root with the environment in CONTRIBUTING.md
activated. `python -m amalgkit --version` and `python -m amalgkit help quant`
check the local CLI without fetching data.

- Use [tests/README.md](tests/README.md#choose-checks-for-a-change) to select
  component tests and consumer/workflow checks. The fast lane deliberately
  omits integration and real CLI subprocess tests; it is not full validation.
- `python .github/scripts/check_quality.py` runs repository lint, formatting
  and mypy on the configured boundaries, and offline documentation checks.
  Do not format or type-check the whole legacy codebase as an incidental edit.
- Run CONTRIBUTING.md's delivery checks before push. Report failed checks and
  skips, including missing external tools or optional extras; a skip is not
  evidence that that path works. CI adds platform/dependency lanes.
- For CLI/example changes, use
  [.agents/skills/verify-cli-docs/SKILL.md](.agents/skills/verify-cli-docs/SKILL.md).
  Push and performance work use the existing skills below.

## Preserve research and data contracts

Selection presets and CLI defaults are scientific choices, not environment
fixes. Do not adjust sample groups, thresholds, seeds, normalization or batch
models simply to make a run pass. Consult `.wiki/Metadata-and-normalization.md`
for original library sizes and kallisto/Oarfish length models, and
`.wiki/Batch-correction-models.md` for protected designs and failure policy.
Keep independent numerical references and their provenance; see
`tests/reference/README.md` before changing TMM fixtures or tolerances.

Preserve public command-module imports, lexical identifiers (including `0001`
and `NA`), output schemas and resume/provenance compatibility. ARCHITECTURE.md
defines the shared validators, rollback and private-FASTQ ownership rules.

Use pytest's temporary fixtures or a fresh temporary workspace for experiments.
Do not edit users' FASTQs, reference databases/indices, metadata, selection
rules, or analysis outputs to validate code. Bundled datasets and reference
fixtures are intentional source assets, not disposable outputs. Keep local
environments, caches, build products and benchmark output out of commits.

At completion review the diff and report changed behavior, executed checks,
static-only checks, unavailable coverage, and any remaining limitations.

<!-- BEGIN KF AGENT POLICY: source=https://github.com/kfuku52/kf-agent-policy; version=10; sha256=82e3c0eb467582a414d9a6b2feaaaf6f5c8ae330d30f2e3efbf8c303155d0e2e -->
# Common agent policy

Repository-specific instructions override these defaults.

- Follow the user's task scope within higher-priority instructions and execution
  permissions. Complete implementation through affected verification and a result
  report; a plan or investigation ends with its requested deliverable. Continue
  authorized work without repeated approval; identify actual blocking boundaries.
- Inspect the worktree and preserve unrelated changes. Refresh remote information
  when needed; do not merge, rebase, or switch branches merely to inspect it.
- Prefer the default branch when starting work without an established branch.
  Preserve an existing task branch; follow explicit user branch instructions.
  Never create or switch branches solely for a commit, push, release, or PR.
- Change or recommend branch protection only when explicitly asked. Honor explicit
  repository-specific direct-push exceptions; otherwise report a rejected push
  without bypassing protection or inventing a branch or PR.
- Unpublished implementation details may be redesigned; preserve existing public
  APIs, file formats, and saved-data compatibility unless a breaking change is
  authorized. Update affected producers, consumers, tests, examples, and docs.
- Fix verified root causes; do not hide failures with fallbacks or weaker checks.
  Document unavoidable workarounds and their removal conditions.
- Read relevant docs and run the repository's check entrypoint for the change and
  phase. Verify affected behavior; report checks run and omitted. Repeat or broaden
  successful checks only for new changes, failures, or unresolved concerns.
- For library metadata, require demonstrated incompatibility for exact pins or
  upper bounds; keep reproducibility locks separate.
- When editing READMEs, keep them concise with useful visuals inline; put extended
  guides in linked documentation.
- For GitHub push/release work, use `prepare-github-push` in `.agents/skills/`.
  Local-only commits need no version bump; GitHub pushes require one.
- For software performance work, use `benchmark-performance` in `.agents/skills/`.
  Performance claims require comparable measurements and equivalent output.
- For GitHub Actions edits, use `optimize-github-actions` in `.agents/skills/`.
  Preserve required coverage; never run untrusted PR code on self-hosted runners.
<!-- END KF AGENT POLICY -->
