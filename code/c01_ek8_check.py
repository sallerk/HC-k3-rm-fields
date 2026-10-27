"""76.2(a): normalisation of Elkies-Kumar Theorem 17 (discriminant 8) against their own example table, the
theta point (z^3 - 1)^2, and the TTV n = 4 curves y^2 = (x + 2)(x^4 - 4x^2 + 2 + sigma).
Usage: python -I -B c01_ek8_check.py
"""
import importlib.util
import sys
from pathlib import Path

import sympy as sp

sys.stdout.reconfigure(encoding="utf-8")
spec = importlib.util.spec_from_file_location("k3lib", Path(__file__).resolve().parent / "k3lib.py")
L = importlib.util.module_from_spec(spec)
spec.loader.exec_module(L)
X, r, s = L.X, L.r_, L.s_

TABLE = [  # (r, s), sextic: the sixteen rows of [EK] Table 2 (Section 7.3), in the published version
    ((sp.Rational(2, 25), sp.Rational(83, 50)), X**6 - X**5 + 3*X**4 + X**3 + X**2 + 2*X + 1),
    ((sp.Rational(-34, 9), sp.Rational(50, 9)), X**6 - 2*X**5 - 2*X**4 - 4*X**3 + 2*X**2 + 4*X - 3),
    ((sp.Integer(-22), sp.Rational(59, 2)), X**6 + 2*X**5 - 3*X**4 + 5*X**3 + X**2 - X - 1),
    ((sp.Rational(2, 49), sp.Rational(58, 49)), -X**6 - 4*X**5 - 6*X**4 + 2*X**2 + 2*X - 1),
    ((sp.Integer(-52), sp.Rational(83, 2)), -4*X**6 + 2*X**5 + 5*X**4 - 7*X**3 - X + 1),
    ((sp.Rational(1, 8), sp.Rational(59, 32)), -X**6 - X**5 - 7*X**2 + 5*X - 4),
    ((sp.Rational(-4, 9), sp.Rational(23, 8)), -2*X**6 + 6*X**5 - X**4 - 7*X**3 - X**2 - 3*X),
    ((sp.Integer(80), sp.Rational(83, 2)), X**6 + 5*X**5 + 8*X**4 + 5*X**3 + 5*X**2 - 8*X + 4),
    ((sp.Rational(19, 2), sp.Integer(22)), -X**6 - 4*X**5 - 8*X**4 - 8*X**3 - 8*X**2 + 4*X - 4),
    ((sp.Integer(94), sp.Integer(-54)), -X**6 + 6*X**4 - 8*X**3 + 6*X**2 + 6*X - 9),
    ((sp.Rational(13, 8), sp.Rational(-2, 9)), -X**6 + 6*X**5 - 9*X**4 + 3*X**2 - 6*X - 2),
    ((sp.Rational(86, 9), sp.Rational(-13, 2)), X**6 - 7*X**4 - 7*X**3 - X**2 + 9*X + 9),
    ((sp.Rational(-3, 8), sp.Rational(9, 4)), -3*X**6 - 9*X**5 - 2*X**4 + 10*X**3 + X**2 + 3*X),
    ((sp.Rational(1, 14), sp.Rational(10, 7)), -X**6 - 2*X**5 - 7*X**4 + 4*X**3 - 3*X**2 + 10*X - 5),
    ((sp.Rational(1, 18), sp.Rational(31, 18)), -3*X**6 - 10*X**5 - 7*X**4 - 5*X**2 + 2*X - 1),
    ((sp.Rational(1, 32), sp.Rational(19, 16)), 4*X**6 - 7*X**5 - 3*X**4 - 2*X**3 + 10*X**2 + 5*X + 1),
]


def weighted_match(I, J):
    """mu with I_k = mu^w J_k for all k (I2, J2 != 0), or None."""
    mu = sp.nsimplify(I[2] / J[2])
    ok = all(sp.simplify(I[k] - mu ** L.W[k] * J[k]) == 0 for k in (2, 4, 6, 10))
    return mu if ok else None


def main():
    L.ic_polys()
    print("== [EK] Section 7.3 table: (r, s) vs Igusa-Clebsch invariants of the listed sextic")
    for (rv, sv), f in TABLE:
        I = L.ic(L.poly_coeffs(f, X))
        J = {k: L.EK8[k].subs({r: rv, s: sv}) for k in L.EK8}
        mu = weighted_match(I, J)
        z2 = L.EK8_Y.subs({r: rv, s: sv})
        pts = L.ek8_points(I)
        print(f"  (r,s)=({rv},{sv}): match={mu is not None} mu={mu} | 2Psi={z2} square in Q: {sp.sqrt(z2).is_rational}"
              f" | all EK points of the sextic: {pts}")

    print("== theta point f = (z^3 - 1)^2 (two P^1 meeting in 3 points)")
    I = L.ic(L.poly_coeffs((X ** 3 - 1) ** 2, X))
    print("  invariants:", I)
    print("  ratios I4/I2^2 =", sp.Rational(I[4], I[2] ** 2), " I6/I2^3 =", sp.Rational(I[6], I[2] ** 3))
    for pt in [(0, 0), (0, sp.Rational(4, 3))]:
        J = {k: L.EK8[k].subs({r: pt[0], s: pt[1]}) for k in L.EK8}
        print(f"  EK8 at (r,s)={pt}: {J}; I4/I2^2 = {sp.Rational(J[4], J[2] ** 2)}, I6/I2^3 = {sp.Rational(J[6], J[2] ** 3)}")
    print("  points of H8 (r = 0 line) with these ratios: solve on r = 0:")
    eq4 = sp.expand(I[4] * L.EK8[2] ** 2 - I[2] ** 2 * L.EK8[4]).subs(r, 0)
    eq6 = sp.expand(I[6] * L.EK8[2] ** 3 - I[2] ** 3 * L.EK8[6]).subs(r, 0)
    print("   s-roots of the I4 relation:", sp.solve(eq4, s), " of the I6 relation:", sp.solve(eq6, s))

    print("== 2-nodal point f = (z^2 - 1)^2 (z^2 - c)")
    c = sp.Symbol('c')
    I = L.ic(L.poly_coeffs((X ** 2 - 1) ** 2 * (X ** 2 - c), X))
    print("  I10 =", sp.factor(I[10]), "; I4/I2^2 =", sp.factor(I[4] / I[2] ** 2))

    print("== TTV n = 4: y^2 = (x + 2)(x^4 - 4x^2 + 2 + sigma), numeric sigma")
    for sig in [1, 2, 3, 5, -1, sp.Rational(1, 2)]:
        f = (X + 2) * (X ** 4 - 4 * X ** 2 + 2 + sig)
        I = L.ic(L.poly_coeffs(f, X))
        pts = L.ek8_points(I)
        info = []
        for (rv, sv) in pts:
            z2 = sp.factor(L.EK8_Y.subs({r: rv, s: sv}))
            info.append(((rv, sv), z2, sp.sqrt(z2)))
        print(f"  sigma={sig}: {info}")


if __name__ == "__main__":
    main()
