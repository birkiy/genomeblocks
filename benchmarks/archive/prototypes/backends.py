"""Prototype: a backend registry for genomeblocks.

Each kind of work ("intervals", "bigwig", "motifs", ...) has named
implementations behind one small protocol. Callers ask for a kind; the
registry returns the user's choice, or lets an ``auto`` policy pick per call.

    from prototypes import backends as B
    B.register("intervals", "numpy", NumpyIntervals())
    with B.use(intervals="pyranges"):
        ...
    B.get("intervals", op="overlaps_any", n=len(a) + len(b))

Third-party packages could add backends through an entry point group
("genomeblocks.backends") without genomeblocks importing them eagerly.
"""
from __future__ import annotations

import contextlib
import contextvars
from typing import Callable, Dict, Protocol

_REGISTRY: Dict[str, Dict[str, object]] = {}
_CHOICE: contextvars.ContextVar = contextvars.ContextVar("gb_backends", default={})
_POLICY: Dict[str, Callable[..., str]] = {}


def register(kind: str, name: str, impl: object) -> None:
    _REGISTRY.setdefault(kind, {})[name] = impl


def set_policy(kind: str, fn: Callable[..., str]) -> None:
    """fn(op=..., n=...) -> backend name, used when the choice is 'auto'."""
    _POLICY[kind] = fn


def get(kind: str, **ctx):
    name = _CHOICE.get().get(kind, "auto")
    if name == "auto":
        name = _POLICY[kind](**ctx)
    return _REGISTRY[kind][name]


@contextlib.contextmanager
def use(**choices):
    token = _CHOICE.set({**_CHOICE.get(), **choices})
    try:
        yield
    finally:
        _CHOICE.reset(token)


class IntervalBackend(Protocol):
    name: str
    def overlaps_any(self, a, b): ...           # bool mask over a
    def merge(self, a): ...
