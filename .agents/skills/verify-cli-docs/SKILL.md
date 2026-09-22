---
name: verify-cli-docs
description: Verify AMALGKIT CLI and documentation changes against the parser and small executable workflow examples. Use for argument, help, tutorial or generated-guide changes, not numerical validation or Wiki publication.
---

# Verify CLI and documentation together

Inputs: the proposed diff, affected command(s), and the user-visible behavior
that should remain valid. Work from the repository root in the activated
environment from [CONTRIBUTING.md](../../../CONTRIBUTING.md), including its
native-thread limits. This skill needs no external bioinformatics tools.

1. Inspect the affected parser/handler and local `.wiki/amalgkit-<command>.md`.
   For generated guides, inspect `amalgkit/dataset.py` as well. Use the local
   canonical pages, not a second copy from the public Wiki. Preserve defaults,
   scientific choices and compatibility unless their change is in scope.
2. Run `python .github/scripts/check_docs.py`. It checks CLI syntax, literal
   defaults and local links against current code without HTTP requests. Fix the
   actual mismatch; do not exempt current examples as historical release notes.
3. Run the existing behavioral checks:

   ```bash
   python -m pytest -q tests/test_cli_parser.py tests/test_cli_help.py tests/test_documented_workflows.py tests/test_doc_tools.py
   ```

   Do not add the fast-lane marker filter: it excludes the help subprocesses
   and executable documentation workflows. Tests use temporary workspaces and
   fixture network responses. Never run a tutorial against a user's data to
   validate its wording.
4. If the changed behavior is not exercised by those cases, select the affected
   command/consumer tests using [tests/README.md](../../../tests/README.md).
   Add a small behavioral regression only for uncovered behavior; a parser
   check alone cannot prove that a metadata path connects pipeline stages.
5. Run the quality/delivery checks required for the current phase in
   CONTRIBUTING.md. Do not duplicate their list here or treat the focused
   checks as proof of scientific validity or external-tool operation.

Output: a concise report of the affected command/example, commands executed,
exit status and test/skipped counts, plus untested behavior. Success requires
zero documentation errors and all selected tests passing. A nonzero exit or
an empty test selection is a failure, not successful validation. Report missing
environment dependencies explicitly; use the documented isolated setup rather
than weakening warnings, assertions or dependency constraints. Separate static
checks from executed workflows. This procedure does not publish the Wiki;
publication has its own instructions in CONTRIBUTING.md.
