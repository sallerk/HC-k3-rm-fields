"""76.4: the special curves left over by Theorems A13 and Prop. 6.6 (NOTE.md, Section 6) carry abelian surfaces of
Picard number >= 3. Fibre types (orders of c4, c6, Delta) of [EK]'s elliptic K3 families:
  - D = 13 (13.txt, three parameters f, g, h of weights 1, 2, 1, before f = 1 is set): the exceptional curve X_b of
    c06 is the locus f = 0 (coordinate s = g/h^2); checked here by the weighted limit of the invariants;
  - D = 12 ([EK]'s printed equation): the lines f = +-1, and the limit e = c f^3, f -> oo (the exceptional curve X3).
The D = 13 formulas are typed from 13.txt (data only). Reads ../outputs/ek_ic_formulas.txt (from c05).
Usage: python -I -B c07_ek13_fzero.py
"""
import importlib.util
import sys
from pathlib import Path

import sympy as sp

sys.stdout.reconfigure(encoding="utf-8")
HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("k3lib", HERE / "k3lib.py")
L = importlib.util.module_from_spec(spec)
spec.loader.exec_module(L)
CACHE = HERE.parent / "outputs" / "ek_ic_formulas.txt"
W = L.W
t, X, c = sp.symbols('t X c')


def ek13(fv, g, h):
    a = 4 * g + fv ** 2
    e = 2 * fv - h
    x0 = e * h
    b = -(3 * x0 ** 4 + 4 * h ** 2 * x0 ** 3 + (-6 * h ** 4 - 8 * a * h ** 2) * x0 ** 2 + (16 * a * h ** 4 - 12 * h ** 6) * x0
          - 5 * h ** 8 + 24 * a * h ** 6 - 16 * a ** 2 * h ** 4) / (64 * h ** 4)
    r = -(x0 ** 3 + 3 * h ** 2 * x0 ** 2 + (3 * h ** 4 - 4 * a * h ** 2) * x0 + h ** 6 - 4 * a * h ** 4) / (16 * h ** 3)
    c1 = (3 * r * x0 ** 2 + (-6 * h ** 2 * r - 4 * b * h) * x0 + (4 * a * h ** 2 - h ** 4) * r) / (4 * h)
    cc = r ** 2 + c1 * t
    p = X ** 3 + a * t ** 2 * X ** 2 + b * t ** 3 * (t - 1) * X + cc * t ** 4 * (t - 1) ** 2
    xs = t ** 2 * (t - 1) * (x0 + h ** 2 * t)
    ys = t ** 2 * (t - 1) * (h ** 3 * t ** 3 + (3 * h * x0 - h ** 3) / 2 * t ** 2
                             + (3 * x0 ** 2 - 6 * h ** 2 * x0 - h ** 4 + 4 * a * h ** 2) / (8 * h) * t + r)
    return sp.expand(p), sp.expand(xs), sp.expand(ys)


def c4c6D(p):
    """For y^2 = X^3 + a2 X^2 + a4 X + a6."""
    P = sp.Poly(p, X)
    a2, a4, a6 = P.coeff_monomial(X ** 2), P.coeff_monomial(X), P.coeff_monomial(1)
    b2, b4, b6, b8 = 4 * a2, 2 * a4, 4 * a6, 4 * a2 * a6 - a4 ** 2
    c4 = b2 ** 2 - 24 * b4
    c6 = -b2 ** 3 + 36 * b2 * b4 - 216 * b6
    D = -b2 ** 2 * b8 - 8 * b4 ** 3 - 27 * b6 ** 2 + 9 * b2 * b4 * b6
    return [sp.numer(sp.together(q)) for q in (c4, c6, D)]


def kodaira(o4, o6, oD):
    """Kodaira type from the orders of c4, c6, Delta (characteristic 0, minimal model)."""
    if oD == 0:
        return "I0"
    if o4 == 0:
        return f"I{oD}"
    if o4 >= 2 and o6 >= 3 and oD >= 6 and (o4 == 2 or o6 == 3):
        return f"I*{oD - 6}"
    return {2: "II", 3: "III", 4: "IV", 8: "IV* (E6)", 9: "III* (E7)", 10: "II* (E8)"}.get(oD, "?")


def fibres(p):
    qs = c4c6D(p)
    out = []
    for pt in (0, 1):
        o = [min(m[0] for m in sp.Poly(sp.expand(q.subs(t, t + pt)), t).monoms()) for q in qs]
        out.append(f"t = {pt}: {o} {kodaira(*o)}")
    degs = [sp.Poly(q, t).degree() for q in qs]
    o = [8 - degs[0], 12 - degs[1], 24 - degs[2]]
    out.append(f"t = oo: {o} {kodaira(*o)}")
    rest = [(sp.Poly(fa, t).degree(), m) for fa, m in sp.factor_list(qs[2], t)[1]
            if sp.Poly(fa, t).degree() > 0 and fa not in (t, t - 1)]
    return "; ".join(out) + f"; other factors of Delta (degree, multiplicity): {rest}"


def load(D, names):
    rows = [ln.split(":", 2) for ln in CACHE.read_text(encoding="utf-8").splitlines() if ln.startswith(f"{D}:")]
    return {int(r_[1][1:]): L.parse_safe(r_[2], names) for r_ in rows if r_[1].startswith("I")}


def weighted_limit(I, sub, var):
    """Weighted limit of (I2:I4:I6:I10) under the substitution, as var -> oo."""
    lead = {}
    for k in W:
        nu, de = sp.fraction(sp.cancel(I[k].subs(sub, simultaneous=True)))
        lead[k] = (sp.Poly(nu, var).degree() - sp.Poly(de, var).degree(),
                   sp.factor(sp.Poly(nu, var).LC() / sp.Poly(de, var).LC()))
    mx = max(sp.Rational(lead[k][0], W[k]) for k in W)
    return [lead[k][1] if sp.Rational(lead[k][0], W[k]) == mx else sp.Integer(0) for k in W]


def same_point(P, Q):
    w = list(W.values())
    return all(sp.simplify(P[i] ** w[j] * Q[j] ** w[i] - Q[i] ** w[j] * P[j] ** w[i]) == 0
               for i in range(4) for j in range(i + 1, 4))


def main():
    g, h, f, e, s = sp.symbols('g h f e s')
    # D = 13
    for fv, label in ((sp.Integer(1), "f = 1"), (sp.Integer(0), "f = 0")):
        p, xs, ys = ek13(fv, g, h)
        print(f"== D = 13 family, {label}: section P on the curve: {sp.expand(ys ** 2 - p.subs(X, xs)) == 0}; {fibres(p)}")
    I13 = load(13, {"g", "h"})
    # (1 : g : h) = (1 : s F^2 : F) ~ (1/F : s : 1) in P(1,2,1)_{f,g,h}: the point s = g/h^2 of the locus f = 0, as F -> oo
    lim = weighted_limit(I13, {g: s * f ** 2, h: f}, f)
    yv = sp.Symbol('y')
    xb = [8 * (4 * yv + 3), 16 * yv ** 2, 32 * yv ** 2 * (6 * yv + 13), 64 * yv ** 4]
    cands = sp.solve(sp.numer(sp.together(lim[1] * xb[0] ** 2 - xb[1] * lim[0] ** 2)), yv)
    ok = [y0 for y0 in cands if same_point(lim, [q.subs(yv, y0) for q in xb])]
    print(f"   weighted limit of IG13 at (g, h) = (s F^2, F), F -> oo: {lim}; equals the map on X_b at y = {ok}")
    print("   Shioda-Tate at f = 0: rho >= 2 + 8 + 7 + 1 + rank MW >= 19 (P has infinite order: a II* fibre forces"
          " trivial torsion). So the abelian surfaces on X_b have Picard number >= 3.")
    # D = 12
    for fv in (f, sp.Integer(-1), sp.Integer(1)):
        p = sp.expand(X ** 3 + ((1 - fv ** 2) * (1 - t) + t) * t * X ** 2 + 2 * e * t ** 3 * (t - 1) * X
                      + e ** 2 * (t - 1) ** 2 * t ** 5)
        print(f"== D = 12 family at f = {fv}: {fibres(p)}")
    I12 = load(12, {"e", "f"})
    lim = weighted_limit(I12, {e: c * f ** 3}, f)
    yv = sp.Symbol('y')
    x3 = [-8 * (3 * yv + 2), 4 * (15 * yv + 1), -4 * (102 * yv ** 2 + 67 * yv + 6), -4 * yv ** 3]
    cands = sp.solve(sp.numer(sp.together(lim[1] * x3[0] ** 2 - x3[1] * lim[0] ** 2)), yv)
    ok = [y0 for y0 in cands if same_point(lim, [q.subs(yv, y0) for q in x3])]
    print(f"   weighted limit of IG12 along e = c f^3, f -> oo: {lim}; equals the map on X3 at y = {ok}")
    p = sp.expand(X ** 3 - t * (1 - t) * X ** 2 + c ** 2 * (t - 1) ** 2 * t ** 5)
    print(f"   limit family y^2 = x^3 - t(1-t) x^2 + c^2 (t-1)^2 t^5: {fibres(p)}")


if __name__ == "__main__":
    main()
