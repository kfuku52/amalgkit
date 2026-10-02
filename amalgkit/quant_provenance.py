"""Content identities for safe quantification reuse, independent of file paths."""

from __future__ import annotations

import hashlib
import os
from collections.abc import MutableMapping, Sequence
from typing import Any

from amalgkit.runtime_utils import validate_run_id

PROVENANCE_KEY = "amalgkit_quant_inputs"


def file_identity(path: str, cache: MutableMapping[Any, Any]) -> dict[str, Any]:
    """Hash each unchanged file once per invocation; detect concurrent writes."""

    def signature() -> tuple[int, ...]:
        stat = os.stat(path)
        return (stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns)

    before = signature()
    key = (os.path.realpath(path), before)
    if key not in cache:
        digest = hashlib.sha256()
        with open(path, "rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        if signature() != before:
            raise ValueError(f"Quant input changed while fingerprinting: {path}")
        cache[key] = {"size": before[2], "sha256": digest.hexdigest()}
    return dict(cache[key])


def build_quant_provenance(
    run_id: str, index: str, inputs: Sequence[str], cache: MutableMapping[Any, Any]
) -> dict[str, Any]:
    return {
        "schema_version": 1,
        "run": run_id,
        "index": file_identity(index, cache),
        "inputs": [{"name": os.path.basename(path), **file_identity(path, cache)} for path in inputs],
    }


def validate_quant_provenance(value: Any) -> str:
    if not isinstance(value, dict) or value.get("schema_version") != 1 or not isinstance(value.get("run"), str):
        return "Invalid quant input provenance."
    try:
        if validate_run_id(value["run"]) != value["run"]:
            return "Invalid quant run identity."
    except ValueError:
        return "Invalid quant run identity."
    inputs = value.get("inputs")
    if not isinstance(inputs, list) or len(inputs) not in {1, 2}:
        return "Invalid quant FASTQ provenance."
    names = []
    for entry in inputs:
        if not isinstance(entry, dict):
            return "Invalid quant FASTQ provenance."
        name = entry.get("name")
        if not isinstance(name, str) or os.path.basename(name) != name or "/" in name or "\\" in name:
            return "Invalid quant FASTQ name provenance."
        suffix = name.removeprefix(value["run"])
        if suffix == name or not suffix.startswith((".", "_1.", "_2.")) or not name.endswith((".fastq", ".fastq.gz")):
            return "Invalid quant FASTQ name provenance."
        names.append(name)
    if len(set(names)) != len(names):
        return "Duplicate quant FASTQ name provenance."
    for identity in [value.get("index"), *inputs]:
        if not isinstance(identity, dict):
            return "Invalid quant file identity."
        digest = identity.get("sha256")
        size = identity.get("size")
        if not isinstance(size, int) or isinstance(size, bool) or size < 0:
            return "Invalid quant file size provenance."
        if not isinstance(digest, str) or len(digest) != 64 or any(c not in "0123456789abcdef" for c in digest):
            return "Invalid quant content digest provenance."
    return ""
