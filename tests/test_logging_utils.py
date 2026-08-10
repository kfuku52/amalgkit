import json
import subprocess
import sys
from pathlib import Path

from amalgkit.logging_utils import configure_logging, get_logger


REPO_ROOT = Path(__file__).resolve().parents[1]


def test_structured_file_logging_includes_workflow_context(tmp_path):
    log_path = tmp_path / "logs" / "amalgkit.jsonl"
    configure_logging(level="INFO", log_file=str(log_path))

    get_logger("test").info(
        "run_started",
        extra={"event": "run_started", "command": "quant", "run_id": "SRR001"},
    )

    payload = json.loads(log_path.read_text(encoding="utf-8").splitlines()[-1])
    assert payload["event"] == "run_started"
    assert payload["command"] == "quant"
    assert payload["run_id"] == "SRR001"
    assert payload["level"] == "INFO"


def test_cli_log_file_records_command_lifecycle_without_log_level(tmp_path):
    log_path = tmp_path / "cli.jsonl"

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "amalgkit",
            "--log_file",
            str(log_path),
            "dataset",
            "--list",
        ],
        cwd=REPO_ROOT,
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    assert result.returncode == 0, result.stderr
    events = [json.loads(line)["event"] for line in log_path.read_text(encoding="utf-8").splitlines()]
    assert events == ["command_start", "command_end"]
