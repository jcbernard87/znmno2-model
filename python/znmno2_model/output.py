"""Configurable output (#15): the time-series columns to write, and profiles against position.

Profiles go to one long table: a header row, then one row per cell per profile time with t_h [h], x_um [um] and the
chosen variables of Model.profiles, each ES15.7 as in the time series. They are written on output rows, so they
never change the time steps: with every `profile_every`-th row, at the first row at or after each of
`profile_times`, and with the last row. A value that does not exist in a cell (the solids outside the cathode) is
NaN.
"""
from __future__ import annotations


def names(text: str) -> list:
    return [w for w in text.replace(",", " ").split() if w]


def times(text: str) -> list:
    return sorted(float(w) for w in names(text))


def choose(available, wanted: str, what: str) -> list:
    """The wanted names (case-insensitive) in the available order; all of them if `wanted` is empty."""
    if not wanted.strip():
        return list(available)
    low = {a.lower(): a for a in available}
    out = []
    for w in names(wanted):
        if w.lower() not in low:
            raise ValueError(f"unknown {what} '{w}'")
        out.append(low[w.lower()])
    return [a for a in available if a in out]


def write_profiles(path, var_names, profiles) -> None:
    """profiles: [(t [s], x [um], {variable: values})]."""
    with open(path, "w") as fh:
        fh.write(" ".join(f"{c:>15}" for c in ("t_h", "x_um", *var_names)) + "\n")
        for t, x, vals in profiles:
            for j in range(len(x)):
                fh.write(" ".join(f"{v:15.7E}" for v in (t / 3600.0, x[j], *(vals[n][j] for n in var_names))) + "\n")


class ProfileRecorder:
    """The profile schedule of a run: with every p.profile_every-th output row, at the first row at or after each of
    p.profile_times, and with the last row. `profile(state, free)` gives (x [um], {variable: values})."""

    def __init__(self, p, available, profile):
        self.on = bool(p.profile_file)
        self.every, self.profile = p.profile_every, profile
        self.times = times(p.profile_times) if self.on else []
        self.vars = choose(available, p.profile_vars, "profile variable") if self.on else []
        self.profiles, self.n_rows, self.last = [], 0, -1

    def row(self, t, state, free):
        if not self.on:
            return
        i = self.n_rows
        self.n_rows += 1
        due = self.every > 0 and i % self.every == 0
        while self.times and t >= self.times[0] * (1.0 - 1.0e-12):
            self.times.pop(0)
            due = True
        if due:
            self._take(t, state, free)
            self.last = i

    def finish(self, res, path, free_of):
        """The last row's profile, if not already taken; write the file."""
        if not self.on:
            return
        if self.last != self.n_rows - 1:
            self._take(res.end[0], res.final_state, free_of(res.final_state))
        res.profiles, res.profile_vars = self.profiles, tuple(self.vars)
        write_profiles(path, self.vars, self.profiles)

    def _take(self, t, state, free):
        x, vals = self.profile(state, free)
        self.profiles.append((t, x, {v: vals[v] for v in self.vars}))
