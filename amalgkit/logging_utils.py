"""Structured logging configuration for CLI and workflow boundaries."""

from __future__ import annotations

import datetime
import json
import logging
import os
import sys
from typing import Any

from amalgkit.redaction import redact_log_value

LOGGER_NAME = "amalgkit"
DEFAULT_LOG_LEVEL = "INFO"


class JsonLogFormatter(logging.Formatter):
    """Render one stable JSON object per log record."""

    def format(self, record: logging.LogRecord) -> str:
        payload: dict[str, Any] = {
            "timestamp": datetime.datetime.fromtimestamp(
                record.created,
                tz=datetime.UTC,
            ).isoformat(),
            "level": record.levelname,
            "logger": record.name,
            "event": getattr(record, "event", record.getMessage()),
            "message": record.getMessage(),
        }
        for key in (
            "command",
            "run_id",
            "species",
            "task_count",
            "failure_count",
            "duration_seconds",
            "returncode",
        ):
            value = getattr(record, key, None)
            if value is not None:
                payload[key] = value
        if record.exc_info:
            payload["exception"] = self.formatException(record.exc_info)
        # Redact at the final boundary: exception causes and interpolated log
        # messages can contain raw commands even when the command field is safe.
        return json.dumps(redact_log_value(payload), ensure_ascii=False, sort_keys=True)


def get_logger(name: str | None = None) -> logging.Logger:
    """Return an AMALGKIT namespaced logger."""
    return logging.getLogger(LOGGER_NAME if name is None else f"{LOGGER_NAME}.{name}")


def configure_logging(level: str | None = None, log_file: str | None = None) -> logging.Logger:
    """Configure JSONL logging without changing third-party root loggers."""
    logger = get_logger()
    for handler in list(logger.handlers):
        logger.removeHandler(handler)
        handler.close()
    logger.propagate = False
    if level is None and not log_file:
        logger.setLevel(logging.CRITICAL + 1)
        logger.addHandler(logging.NullHandler())
        return logger
    resolved_level = getattr(logging, str(level or DEFAULT_LOG_LEVEL).upper(), logging.INFO)
    logger.setLevel(resolved_level)
    formatter = JsonLogFormatter()

    if level is not None:
        stream_handler = logging.StreamHandler(sys.stderr)
        stream_handler.setLevel(resolved_level)
        stream_handler.setFormatter(formatter)
        logger.addHandler(stream_handler)

    if log_file:
        resolved_path = os.path.realpath(os.path.expanduser(str(log_file)))
        parent = os.path.dirname(resolved_path)
        if parent:
            os.makedirs(parent, exist_ok=True)
        file_handler = logging.FileHandler(resolved_path, encoding="utf-8")
        file_handler.setLevel(resolved_level)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    return logger
