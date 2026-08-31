# Security Policy

Security fixes are provided on the current default branch. Patch updates are
not published as GitHub Releases or Bioconda recipe updates; those channels
may lag behind the default branch. Older versions may not receive backports.

To update from the default branch:

```bash
pip install --upgrade git+https://github.com/kfuku52/amalgkit
amalgkit --version
```

Consult [CHANGELOG.md](CHANGELOG.md) and any published security advisory for
the affected and fixed versions. A latest-release badge alone does not establish
that an installation includes all default-branch fixes.

Please do not report suspected vulnerabilities in a public issue. Use the
repository's **Security** tab and **Report a vulnerability** to open a private
GitHub Security Advisory. Include affected commands, a minimal reproduction,
impact, and any suggested mitigation.

If private vulnerability reporting is unavailable, email `kfuku52@gmail.com`
with the subject `AMALGKIT security report`.

For ordinary correctness bugs without security impact, use the public issue
tracker.
