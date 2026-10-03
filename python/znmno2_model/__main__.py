"""Command line: python -m znmno2_model [input.nml]"""
import sys
from pathlib import Path

from .namelist import load
from .simulate import run


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    path = Path(argv[0]) if argv else Path(__file__).resolve().parents[2] / "input" / "default.nml"
    p, extra = load(path)
    if extra.get("mode", "corrected") != "corrected":
        raise SystemExit("this program runs mode = 'corrected'; the faithful ports are znmno2_model.faithful "
                         "(Python) and the Fortran and C++ programs")
    out = extra.get("file", "znmno2_out.txt")
    r = run(p)
    r.write(out)
    print(f"exit '{r.exit_reason}' after {r.steps} steps, {r.rows[-1][3]:.1f} mAh/g; wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
