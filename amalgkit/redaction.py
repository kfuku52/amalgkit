"""Credential-safe text for human-facing diagnostics (never command execution)."""

from __future__ import annotations

import re
import urllib.parse
from typing import Any

_URL_TOKEN_PATTERN = re.compile(r"[A-Za-z][A-Za-z0-9+.-]*://[^\s\"'<>]+")


def _redact_url_token(url: str) -> str:
    scheme = url.split("://", 1)[0].casefold()
    try:
        parts = urllib.parse.urlsplit(url)
        host = parts.hostname
        port = parts.port
    except ValueError:
        return f"{scheme}://<redacted>"
    if host is None:
        return f"{scheme}://<redacted>"
    if ":" in host:
        host = f"[{host}]"
    netloc = host if port is None else f"{host}:{port}"
    return urllib.parse.urlunsplit((scheme, netloc, parts.path, "", ""))


def redact_url_for_logging(value: object) -> str:
    """Remove URL userinfo, query and fragment, including in exception chains.

    Malformed URL tokens fail closed. Non-URL text and the URL host/path remain
    available for diagnosis. Do not use this to transform machine-readable data.
    """
    return _URL_TOKEN_PATTERN.sub(lambda match: _redact_url_token(match.group(0)), str(value))


def redact_log_value(value: Any) -> Any:
    """Sanitize structured fields before JSON encoding, without mutating them."""
    if isinstance(value, str):
        return redact_url_for_logging(value)
    if isinstance(value, dict):
        return {redact_log_value(key): redact_log_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [redact_log_value(item) for item in value]
    return value
