"""Compute CLI — invoked by the Node orchestrator as:

    python -m compute.cli <job_name> <json_inputs>

Exactly one JSON document is written to stdout on success. All logging goes
to stderr. A non-zero exit code means the job failed; the orchestrator
surfaces the stderr tail to the client.
"""

from __future__ import annotations

import importlib
import json
import sys
import traceback
from typing import Any

# Explicit whitelist so callers can't trigger import of arbitrary modules.
JOBS: dict[str, str] = {
    "venn": "compute.jobs.venn",
}


def _die(msg: str, code: int = 2) -> None:
    print(msg, file=sys.stderr)
    sys.exit(code)


def main(argv: list[str] | None = None) -> None:
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 2:
        _die("usage: python -m compute.cli <job> <json_inputs>")

    job_name, raw_inputs = argv
    module_path = JOBS.get(job_name)
    if module_path is None:
        _die(f"unknown job: {job_name}")

    try:
        inputs: dict[str, Any] = json.loads(raw_inputs)
        if not isinstance(inputs, dict):
            raise ValueError("inputs must be a JSON object")
    except (ValueError, json.JSONDecodeError) as exc:
        _die(f"invalid JSON inputs: {exc}")

    try:
        mod = importlib.import_module(module_path)
        if not hasattr(mod, "run"):
            _die(f"job '{job_name}' has no `run(inputs)` function")
        result = mod.run(inputs)
    except Exception as exc:  # noqa: BLE001
        print(f"job '{job_name}' raised {type(exc).__name__}: {exc}", file=sys.stderr)
        traceback.print_exc(file=sys.stderr)
        sys.exit(1)

    # One JSON document, nothing else. No trailing newline issues — json.dump
    # writes raw; we append a single newline for readability.
    json.dump(result, sys.stdout, separators=(",", ":"), default=_default)
    sys.stdout.write("\n")


def _default(obj: Any) -> Any:
    # numpy/pandas scalars slip through json.dumps; funnel them through here.
    if hasattr(obj, "item"):
        return obj.item()
    raise TypeError(f"Object of type {type(obj).__name__} is not JSON serializable")


if __name__ == "__main__":
    main()
