"""76.4: two parity conditions on the moduli map of a pencil with RM over the base (NOTE.md, Section 6), applied to
the special curves left by Theorem A13 (X_b, D = 13) and by Proposition 6.6 (the line f = -1, D = 12).
1. The degree-15 invariant. I15^2 = c R30(I2, I4, I6, I10), with R30 of weight 15. R30 is found up to a constant as
   the weight-15 form vanishing on the sextics g(z^2) (nullspace over Q), and checked: nullity 1, and
   R30(I(f1)) R30(I(f2)) is a non-zero rational square for random sextics f1, f2.
   For a pencil, I_{2k}(t) = nu^k E_{2k}(phi(t)) gives nu^15 R30(E(phi(t))) = c I15(t)^2, so
   ord_t nu + ord_t R30(E(phi(t))) is even for every t.
2. The lift to Y_-(D): [EK]'s double cover z^2 = b, restricted to the curve. The b's are typed from [EK]'s theorems
   for D = 12 and 13 (data).
Usage: python -I -B c08_parity.py
"""
import importlib.util
import math
import random
import sys
from pathlib import Path

import sympy as sp
from flint import fmpz_mat

sys.stdout.reconfigure(encoding="utf-8")
HERE = Path(__file__).resolve().parent


def _load(name):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


L = _load("k3lib")
C6 = _load("c06_ek_boundary")
W = L.W
J = sp.symbols('J2 J4 J6 J10')
MONS = [(a, b, c, d) for d in range(4) for c in range(6) for b in range(8) for a in range(16)
        if a + 2 * b + 3 * c + 5 * d == 15]


def mon_values(I):
    return [I[2] ** a * I[4] ** b * I[6] ** c * I[10] ** d for a, b, c, d in MONS]


def r30():
    rnd = random.Random(76)
    rows = []
    while len(rows) < len(MONS) + 15:
        p, q, r = (rnd.randint(-9, 9) for _ in range(3))
        I = L.ic([r, 0, q, 0, p, 0, 1])
        if I[10] != 0:
            rows.append([int(v) for v in mon_values(I)])
    X, nullity = fmpz_mat(rows).nullspace()
    assert nullity == 1, nullity
    co = [int(X[i, 0]) for i in range(len(MONS))]
    g = math.gcd(*co)
    co = [v // g for v in co]
    if co[MONS.index((15, 0, 0, 0))] < 0:
        co = [-v for v in co]
    R = sum(cf * J[0] ** a * J[1] ** b * J[2] ** c * J[3] ** d for cf, (a, b, c, d) in zip(co, MONS))
    return R, nullity, rnd


def at(R, I):
    return R.subs(dict(zip(J, (I[2], I[4], I[6], I[10]))), simultaneous=True)


def show_factors(expr, gens, label):
    fl = sp.factor_list(sp.expand(expr), *gens)
    print(f"   {label}: constant {fl[0]}; factors (factor, multiplicity):")
    for fa, mu in fl[1]:
        print(f"      {fa}  ^{mu}   (degree {sp.Poly(fa, *gens).total_degree()})")
    return fl


def main():
    R, nullity, rnd = r30()
    print(f"== R30: weight-15 forms vanishing on the sextics g(z^2): nullity {nullity}; "
          f"{sum(1 for _ in sp.Poly(R, *J).terms())} terms; coefficient of J2^15: {sp.Poly(R, *J).coeff_monomial(J[0] ** 15)}")
    vals = []
    for _ in range(5):
        co = [rnd.randint(-5, 5) for _ in range(6)] + [rnd.randint(1, 5)]
        vals.append(int(at(R, L.ic(co))))
    sq = [v != 0 and vals[0] * v > 0 and math.isqrt(vals[0] * v) ** 2 == vals[0] * v for v in vals]
    print(f"   random sextics: R30 != 0 and R30(f1) R30(fi) a square: {sq}; sign of R30(f1): {1 if vals[0] > 0 else -1}")
    print(f"   R30 at (1:0:0:0): {R.subs({J[0]: 1, J[1]: 0, J[2]: 0, J[3]: 0})}; at the theta point (1, 1/16, 1/64, 0): "
          f"{R.subs({J[0]: 1, J[1]: sp.Rational(1, 16), J[2]: sp.Rational(1, 64), J[3]: 0})}")
    # D = 12: the line f = -u in (e : f : u)
    (e, f), N, m, E = C6.homog(12)
    u = C6.u
    El = {k: sp.expand(E[k].subs(f, -u)) for k in W}
    print(f"== D = 12, the line f = -u: E restricted = {[sp.factor(El[k]) for k in W]}")
    b12 = (f - 1) * (f + 1) * (f ** 6 - f ** 4 - 18 * e * f ** 2 + 27 * e ** 2 + 16 * e)
    print(f"   [EK] branch function on f = -1: {sp.expand(b12.subs(f, -1))} (the line is a branch component)")
    show_factors(at(R, El), (e, u), "R30(E) on the line, a form in (e, u)")
    # D = 13: the curve X_b, map R(y0 : y1) of weights (1, 2, 3, 5) with y = y1/y0 = s = g/h^2 (c07)
    y0, y1 = sp.symbols('y0 y1')
    Rb = {2: 8 * (4 * y1 + 3 * y0), 4: 16 * y1 ** 2, 6: 32 * y1 ** 2 * (6 * y1 + 13 * y0), 10: 64 * y1 ** 4 * y0}
    print("== D = 13, the curve X_b = {f = 0} of P(1,2,1)_{f,g,h}, coordinate y = g/h^2")
    show_factors(at(R, Rb), (y0, y1), "R30 on X_b, a form in (y0, y1)")
    g, h, ff = sp.symbols('g h f')
    b13 = 108 * g * h ** 3 - (27 * g ** 2 + 468 * g - 4) * h ** 2 + 8 * (82 * g ** 2 + 71 * g - 2) * h - 16 * (2 * g - 1) ** 3
    wt = {g: 2, h: 1}
    terms = sp.Poly(b13, g, h).terms()
    top = max(wt[g] * i + wt[h] * j for (i, j), _ in terms)
    b13w = sum(cf * g ** i * h ** j * ff ** (top - 2 * i - j) for (i, j), cf in terms)
    print(f"   [EK] branch function, weighted (f, g, h) of weights (1, 2, 1), weight {top}: on f = 0: "
          f"{sp.factor(b13w.subs(ff, 0))}; with g = y h^2: {sp.factor(b13w.subs(ff, 0).subs(g, sp.Symbol('y') * h ** 2) / h ** top)}")
    # D = 12: the curve X3 (excluded by the order argument; for completeness)
    R3 = {2: -8 * (3 * y1 + 2 * y0), 4: 4 * (15 * y1 + y0) * y0, 6: -4 * (102 * y1 ** 2 + 67 * y1 * y0 + 6 * y0 ** 2) * y0,
          10: -4 * y1 ** 3 * y0 ** 2}
    print("== D = 12, the curve X3, coordinate y = y1/y0 = c (e = c f^3, f -> oo)")
    show_factors(at(R, R3), (y0, y1), "R30 on X3")
    print(f"   [EK] branch function along e = c f^3: leading term {sp.factor(sp.LT(sp.Poly(sp.expand(b12.subs(e, sp.Symbol('c') * f ** 3)), f)))}")
    # which curves of the bielliptic locus: the D6 locus z^6 + a z^3 + 1 and the D4 locus z^5 + z^3 + a z
    a = sp.Symbol('a')
    line = {k: El[k].subs(u, 1) for k in W}
    xbs = {k: Rb[k].subs({y0: 1, y1: sp.Symbol('y')}) for k in W}
    fams = {"D6: z^6 + a z^3 + 1": [1, 0, 0, a, 0, 0, 1], "D4: z^5 + z^3 + a z": [0, a, 0, 1, 0, 1, 0]}
    for name, co in fams.items():
        P = L.ic(co)
        for cname, cur, var in (("D = 12 line f = -1", line, e), ("D = 13 X_b", xbs, sp.Symbol('y'))):
            # (P2:P4:P6:P10) = (C2:C4:C6:C10) in P(1,2,3,5): eliminate the curve parameter from two ratio equations
            q4 = sp.numer(sp.together(P[4] * cur[2] ** 2 - cur[4] * P[2] ** 2))
            q6 = sp.numer(sp.together(P[6] * cur[2] ** 3 - cur[6] * P[2] ** 3))
            q10 = sp.numer(sp.together(P[10] * cur[2] ** 5 - cur[10] * P[2] ** 5))
            res = sp.factor(sp.resultant(q4, q6, var))
            same = all(sp.resultant(q4, qq, var) == 0 for qq in (q6, q10))
            print(f"   {name} vs {cname}: Res(q4, q6) = {res if len(str(res)) < 200 else str(res)[:200] + '...'}; curves equal: {same}")


if __name__ == "__main__":
    main()
