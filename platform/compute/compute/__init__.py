"""Genomeblocks compute sidecar.

Each submodule under ``compute.jobs`` implements one disposable job. The
``compute.cli`` entrypoint dispatches a job by name and prints exactly one
JSON document to stdout — that is the contract Node relies on.
"""

__version__ = "0.1.0"
