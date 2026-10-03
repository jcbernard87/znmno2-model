"""Reader for the shared input file (a Fortran namelist subset).

The Fortran and C++ ports will read the same file. Supported syntax:

    ! comment
    &group
      name = 1.0d-3, other = 22   ! comments after values
      mode = 'faithful'
    /

Values are integers, reals (e/E/d/D exponents), quoted strings, or .true./.false.
Names are case-insensitive. Arrays and repeat counts are not supported.
"""
from __future__ import annotations

import re
from dataclasses import fields

from .params import Params

GROUPS = {
    "run": ["mode", "data_dir"],             # read by the Fortran and C++ programs; Python runs mode = 'corrected'
    "cell": ["A_cell", "L_probe", "eps_probe", "tau_probe", "L_sep", "eps_sep", "tau_factor_sep",
             "bruggeman_sep", "L_cath", "eps_cath", "tau_factor_cath", "bruggeman_cath", "n_probe", "n_sep",
             "n_cath", "sigma"],
    "solids": ["vf_MnO2", "vf_ZMO", "vf_host", "vf_ZHS", "M_MnO2", "rho_MnO2", "rho_ZMO", "rho_host", "M_ZHS",
               "rho_ZHS", "r_MnO2", "r_ZMO", "r_host", "r_ZHS", "z_ZMO", "zmin", "zmax", "theta0", "mass_AM"],
    "electrolyte": ["c_ZnSO4", "c_MnSO4", "c_H2SO4", "D_Zn", "D_Mn", "D_SO4", "D_H"],
    "reactions": ["U1", "k1", "alpha1", "U2", "k2", "alpha2", "a_seed_R2", "k3", "alpha3", "V_at_zmin",
                  "V_at_zmax", "c_ref3", "k_an", "alpha_an", "logK_ZHS", "k_ZHS", "a_seed_ZHS"],
    "options": ["ph_mode", "pH_fixed", "species", "transport", "basis", "R1_on", "R2_on", "R3_on", "zhs",
                "zhs_nucleation", "ZnO_on", "ZnOH2_on", "logK_ZnO", "logK_ZnOH2", "k_ZnO", "k_ZnOH2", "M_ZnO",
                "rho_ZnO", "M_ZnOH2", "rho_ZnOH2", "r3_ocp", "r3_end_width", "U3_nernst", "logk_file", "equilibria_db",
                "D_OH",
                "D_HSO4", "D_complex"],
    "constants": ["R", "T", "F"],
    "protocol": ["steps", "cycles", "end_on_cutoff", "V_min", "V_max"],
    "numerics": ["dt", "newton_tol", "newton_max_iter"],
    "output": ["file", "write_interval"],
}

EXTRA = ("file", "mode", "data_dir")        # settings that are not model parameters

_TOKEN = re.compile(r"""\s*(?:'([^']*)'|"([^"]*)"|([^,\s/=]+))""")


def _strip_comment(line: str) -> str:
    out, quote = [], None
    for ch in line:
        if quote:
            if ch == quote:
                quote = None
        elif ch in "'\"":
            quote = ch
        elif ch == "!":
            break
        out.append(ch)
    return "".join(out)


def _value(text: str, quoted: bool):
    if quoted:
        return text
    low = text.lower()
    if low in (".true.", "t", ".t."):
        return True
    if low in (".false.", "f", ".f."):
        return False
    try:
        return int(text)
    except ValueError:
        return float(low.replace("d", "e"))


def parse(text: str) -> dict:
    """Parse namelist text into {group: {name: value}} with lower-case keys."""
    body = "\n".join(_strip_comment(l) for l in text.splitlines())
    out: dict = {}
    for m in re.finditer(r"&(\w+)((?:'[^']*'|\"[^\"]*\"|[^/'\"])*)/", body):   # a '/' inside quotes is text
        group, content = m.group(1).lower(), m.group(2)
        entries = {}
        for part in re.finditer(r"(\w+)\s*=\s*('[^']*'|\"[^\"]*\"|[^,\s/]+)", content):
            name, raw = part.group(1).lower(), part.group(2)
            quoted = raw[:1] in "'\""
            entries[name] = _value(raw[1:-1] if quoted else raw, quoted)
        out[group] = entries
    return out


def load(path) -> tuple[Params, dict]:
    """Read an input file; return (Params, extra) where extra holds non-model settings (output file)."""
    with open(path) as fh:
        data = parse(fh.read())
    known = {g: {n.lower(): n for n in names} for g, names in GROUPS.items()}
    values, extra = {}, {}
    for group, entries in data.items():
        if group not in known:
            raise ValueError(f"unknown namelist group &{group}")
        for key, val in entries.items():
            if key not in known[group]:
                raise ValueError(f"unknown name {key!r} in &{group}")
            name = known[group][key]
            (extra if name in EXTRA else values)[name] = val
    types = {f.name: f.type for f in fields(Params)}
    for name, val in list(values.items()):
        t = types[name]
        if t in ("int", int):
            values[name] = int(val)
        elif t in ("str", str):
            values[name] = str(val)
        elif t in ("bool", bool):
            values[name] = bool(val)
        else:
            values[name] = float(val)
    return Params(**values), extra
