"""Pure parsers and compatibility checks for fasterq-dump probes."""

from __future__ import annotations

import re


def parse_fasterq_dump_written_spots(stdout_txt: str, stderr_txt: str) -> int | None:
    combined = "\n".join([stdout_txt or "", stderr_txt or ""])
    matched = re.findall(
        r"^\s*spots\s+written\s*:\s*([0-9][0-9,]*)\s*$",
        combined,
        flags=re.IGNORECASE | re.MULTILINE,
    )
    if len(matched) == 0:
        return None
    return int(matched[-1].replace(",", ""))


def parse_fasterq_dump_written_reads(stdout_txt: str, stderr_txt: str) -> int | None:
    combined = "\n".join([stdout_txt or "", stderr_txt or ""])
    matched = re.findall(
        r"^\s*reads\s+written\s*:\s*([0-9][0-9,]*)\s*$",
        combined,
        flags=re.IGNORECASE | re.MULTILINE,
    )
    if len(matched) == 0:
        return None
    return int(matched[-1].replace(",", ""))


def detect_fasterq_spot_range_support(help_stdout_txt: str, help_stderr_txt: str = "") -> bool:
    """Detect sra-tools spot-range flags, tolerating empty legacy help."""
    combined = "\n".join([str(help_stdout_txt or ""), str(help_stderr_txt or "")])
    if combined.strip() == "":
        return True
    has_start = "-N|" in combined or "--minSpotId" in combined or "--min-spot-id" in combined.lower()
    has_end = "-X|" in combined or "--maxSpotId" in combined or "--max-spot-id" in combined.lower()
    return bool(has_start and has_end)


def parse_fasterq_dump_version(
    version_stdout_txt: str,
    version_stderr_txt: str = "",
) -> tuple[int, int, int] | None:
    combined = "\n".join([str(version_stdout_txt or ""), str(version_stderr_txt or "")])
    matched = re.search(r"\b([0-9]+)\.([0-9]+)(?:\.([0-9]+))?\b", combined)
    if matched is None:
        return None
    return int(matched.group(1)), int(matched.group(2)), int(matched.group(3) or 0)


def ensure_supported_fasterq_dump_version(
    version_stdout_txt: str,
    version_stderr_txt: str = "",
    executable_name: str = "fasterq-dump",
) -> tuple[int, int, int]:
    """Return the parsed version or raise when sra-tools is unsupported."""
    version_parts = parse_fasterq_dump_version(version_stdout_txt, version_stderr_txt)
    if version_parts is None:
        raise RuntimeError(
            f"Could not determine fasterq-dump version from {executable_name} output. "
            "sra-tools >= 3 is required for amalgkit getfastq."
        )
    if version_parts[0] < 3:
        version_txt = "{}.{}.{}".format(*version_parts)
        raise RuntimeError(
            f"Unsupported fasterq-dump version detected: {version_txt} from {executable_name}. "
            "sra-tools >= 3 is required for amalgkit getfastq."
        )
    return version_parts
