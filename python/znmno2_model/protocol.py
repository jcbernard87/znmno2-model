"""Cycling protocols.

A protocol string such as ``'cc I=100 Vmin=1.0; rest t=600; cc I=-100 Vmax=1.85; cv V=1.85 Imin=5'``
is parsed into a list of Step objects. Currents are specific currents in mA/g of active material,
positive for discharge. A GITT discharge is ``'cc I=50 t=720; rest t=1800'`` with ``cycles`` large
and ``end_on_cutoff`` true (the cutoff ends the protocol).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

from .params import Params

KINDS = ("cc", "cv", "rest")
_KEYS = {"cc": {"i", "t", "vmin", "vmax"}, "cv": {"v", "t", "imin"}, "rest": {"t"}}


@dataclass(frozen=True)
class Step:
    kind: str
    I: float = 0.0                  # cc: specific current [mA/g], positive = discharge
    V: float = 0.0                  # cv: held voltage [V]
    t: Optional[float] = None       # maximum duration [s]
    Vmin: float = 0.0
    Vmax: float = 0.0
    Imin: Optional[float] = None    # cv: end when |I| <= Imin [mA/g]


def parse(text: str, p: Params) -> list[Step]:
    text = (text or "").strip()
    if not text:
        raise ValueError("empty protocol")
    steps = []
    for n, part in enumerate(filter(None, (s.strip() for s in text.split(";"))), start=1):
        words = part.split()
        kind = words[0].lower()
        if kind not in KINDS:
            raise ValueError(f"step {n}: unknown step type {words[0]!r} (expected cc, cv or rest)")
        kv = {}
        for w in words[1:]:
            if "=" not in w:
                raise ValueError(f"step {n}: expected key=value, got {w!r}")
            k, v = w.split("=", 1)
            k = k.lower()
            if k not in _KEYS[kind]:
                raise ValueError(f"step {n}: {kind} does not take {k!r}")
            kv[k] = float(v.lower().replace("d", "e"))
        t = kv.get("t")
        if t is not None and t <= 0:
            raise ValueError(f"step {n}: t must be positive")
        if kind == "cc":
            if "i" not in kv:
                raise ValueError(f"step {n}: cc needs I=")
            steps.append(Step("cc", I=kv["i"], t=t, Vmin=kv.get("vmin", p.V_min), Vmax=kv.get("vmax", p.V_max)))
        elif kind == "cv":
            if "v" not in kv:
                raise ValueError(f"step {n}: cv needs V=")
            if t is None and "imin" not in kv:
                raise ValueError(f"step {n}: cv needs t= or Imin= to end")
            steps.append(Step("cv", V=kv["v"], t=t, Imin=kv.get("imin"), Vmin=p.V_min, Vmax=p.V_max))
        else:
            if t is None:
                raise ValueError(f"step {n}: rest needs t=")
            steps.append(Step("rest", t=t, Vmin=p.V_min, Vmax=p.V_max))
    return steps


def expand(p: Params) -> list[Step]:
    return parse(p.steps, p) * max(1, p.cycles)
