"""gfortran list-directed output (`write(unit,*) ...`), reproduced byte for byte.

The original program writes Time_Voltage.txt with list-directed WRITE. gfortran formats each item
after a separator blank:
  REAL(8)       G25.17E3: F20.(17-k) plus 5 blanks when 0.1 <= |x| < 1e17 (k = decimal exponent),
                otherwise E25.17E3 (1P, three-digit exponent); NaN/Infinity right-justified in 25
  INTEGER(4)    I11
  CHARACTER     the text itself
"""
from __future__ import annotations

import math

_W, _D = 25, 17


def _real(x: float) -> str:
    if math.isnan(x):
        return "NaN".rjust(_W)
    if math.isinf(x):
        return ("Infinity" if x > 0 else "-Infinity").rjust(_W)
    if x == 0.0:
        body = ("-" if math.copysign(1.0, x) < 0 else "") + "0." + "0" * (_D - 1)
        return body.rjust(_W - 5) + " " * 5
    mant, exp = f"{abs(x):.{_D - 1}e}".split("e")      # 17 significant digits, correctly rounded
    e = int(exp)
    digits = mant.replace(".", "")
    sign = "-" if x < 0 else ""
    k = e + 1                                            # 10**(k-1) <= |x| < 10**k after rounding
    if 0 <= k <= _D:                                     # 0.1 <= |x| < 1e17: F editing
        if k == 0:
            body = "0." + digits
        else:
            body = digits[:k] + "." + digits[k:]
        return (sign + body).rjust(_W - 5) + " " * 5
    esign = "-" if e - 0 < 0 else "+"
    # 1P: one digit before the point; the exponent is that of the leading digit
    body = f"{digits[0]}.{digits[1:]}E{esign}{abs(e):03d}"
    return (sign + body).rjust(_W)


def _int(i: int) -> str:
    return f"{i:11d}"


def record(*items) -> str:
    """One list-directed record (without the newline)."""
    out = []
    for it in items:
        if isinstance(it, str):
            out.append(" " + it)
        elif isinstance(it, (bool,)):
            raise TypeError("logical output not supported")
        elif isinstance(it, int):
            out.append(" " + _int(it))
        else:
            out.append(" " + _real(float(it)))
    return "".join(out)
