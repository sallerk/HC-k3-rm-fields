"""76.4: read [EK]'s auxiliary Igusa-Clebsch formulas (igusaD.txt of the arXiv source, D = 5, 8, 12, 13, 17) as
data, with the parser k3lib.parse_safe (no code from the source is run), and check them:
  - D = 5, 8: against the printed formulas ([EK] Cor. 15, Thm 17; k3lib.EK8);
  - all D: against [EK]'s example tables (read from hmspaper.tex), point by point in P(1,2,3,5).
Writes the parsed invariants and bad-locus factors to ../outputs/ek_ic_formulas.txt ('^' for powers, readable by
parse_safe), so later scripts do not read the source folder.
Usage: python -I -B c05_ek_igusa_parse.py SOURCE_FOLDER (the extracted arXiv source of [EK], arXiv:1209.3527v3)
"""
import importlib.util
import itertools
import re
import sys
from pathlib import Path

import sympy as sp

sys.stdout.reconfigure(encoding="utf-8")
HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("k3lib", HERE / "k3lib.py")
L = importlib.util.module_from_spec(spec)
spec.loader.exec_module(L)
if len(sys.argv) < 2:
    sys.exit("usage: c05_ek_igusa_parse.py SOURCE_FOLDER (the extracted arXiv source of [EK], arXiv:1209.3527v3)")
SRC = Path(sys.argv[1])
OUT = HERE.parent / "outputs" / "ek_ic_formulas.txt"
DS = (5, 8, 12, 13, 17)
NAMES = ["A1", "A", "B1", "B", "B2"]
FINAL = "-24*B1/A1,-12*A,96*(A/A1)*B1-36*B,-4*A1*B2"


def read_igusa(D):
    text = (SRC / str(D) / f"igusa{D}.txt").read_text(encoding="ascii")
    m = re.search(r"IG\((\w+),\s*(\w+),\s*A1,A,B1,B,B2\)\s*=(.*?)\[(.*?)\]", text, re.S)
    v1, v2, body, final = m.groups()
    assert re.sub(r"\s", "", final) == FINAL, final
    parts = [s for s in body.split(";") if s.strip()]
    assert [p.split("=", 1)[0].strip() for p in parts] == NAMES
    val = {p.split("=", 1)[0].strip(): L.parse_safe(p.split("=", 1)[1], {v1, v2}) for p in parts}
    A1, A, B1, B, B2 = (val[n] for n in NAMES)
    I = {2: -24 * B1 / A1, 4: -12 * A, 6: 96 * (A / A1) * B1 - 36 * B, 10: -4 * A1 * B2}
    I = {k: sp.factor(sp.cancel(v)) for k, v in I.items()}
    bad = re.search(r"badlocus\s*=\s*(.*?);", text, re.S)
    badf = [f for f, _ in sp.factor_list(L.parse_safe(bad.group(1), {v1, v2}))[1]] if bad else []
    return (sp.Symbol(v1), sp.Symbol(v2)), I, badf


def read_table(D, tex):
    sec = re.search(r"\\section\{Discriminant \$" + str(D) + r"\$\}(.*?)\\section\{", tex, re.S).group(1)
    tab = re.search(r"\\subsection\{Examples\}(.*?)\\end\{tabular\}", sec, re.S).group(1)
    rows = re.findall(r"\$\((-?\d+(?:/\d+)?),\s*(-?\d+(?:/\d+)?)\)\$\s*&\s*\$([^$]*)\$", tab)
    x = sp.Symbol('x')
    return [((sp.Rational(a), sp.Rational(b)), L.parse_safe(p, {"x"})) for a, b, p in rows], x


def wp_equal(I, J):
    """(I2:I4:I6:I10) = (J2:J4:J6:J10) in P(1,2,3,5)."""
    if any((I[k] == 0) != (J[k] == 0) for k in L.W):
        return False
    return all(sp.simplify(I[a] ** L.W[b] * J[b] ** L.W[a] - J[a] ** L.W[b] * I[b] ** L.W[a]) == 0
               for a, b in itertools.combinations(L.W, 2))


def main():
    L.ic_polys()
    tex = (SRC / "hmspaper.tex").read_text(encoding="utf-8", errors="strict")
    g, h = sp.symbols('g h')
    EK5 = {2: 6 * (4 * g + 1), 4: 9 * g ** 2, 6: 9 * (4 * h + 9 * g ** 3 + 2 * g ** 2), 10: 4 * h ** 2}
    lines = []
    for D in DS:
        (v1, v2), I, badf = read_igusa(D)
        print(f"== D = {D}, variables ({v1}, {v2})")
        for k in (2, 4, 6, 10):
            print(f"  I{k} = {I[k]}")
        print(f"  bad-locus factors: {badf}")
        if D == 5:
            print(f"  equal to [EK] Cor. 15: {all(sp.expand(I[k] - EK5[k]) == 0 for k in L.W)}")
        if D == 8:
            sub = {v1: L.r_, v2: L.s_}
            print(f"  equal to [EK] Thm 17 (k3lib.EK8): {all(sp.expand(I[k].subs(sub) - L.EK8[k]) == 0 for k in L.W)}")
        rows, x = read_table(D, tex)
        ok = []
        for (a, b), f in rows:
            J = L.ic(L.poly_coeffs(f, x))
            Ip = {k: sp.cancel(I[k].subs({v1: a, v2: b})) for k in L.W}
            ok.append(wp_equal(Ip, J))
        print(f"  table: {len(rows)} rows, {sum(ok)} match; failing rows: {[i for i, o in enumerate(ok) if not o]}")
        lines.append(f"{D}:vars:{v1},{v2}")
        lines += [f"{D}:I{k}:{sp.sstr(I[k]).replace('**', '^')}" for k in (2, 4, 6, 10)]
        lines += [f"{D}:bad:{sp.sstr(f).replace('**', '^')}" for f in badf]
    OUT.write_text("".join(s + "\n" for s in lines), encoding="utf-8")
    print(f"written {OUT.name}")


if __name__ == "__main__":
    main()
