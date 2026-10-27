"""78.7d: the curve over U = 0 in the model X of Y_-(8) (NOTE.md, Section 5.4).
1. Is the image of U = 0 in M2 the D4 locus? Test: the point (R : S : U) = (rho : 1 : 0) against the D4 family
   x^5 + x^3 + t x (and x(x^4 + a x^2 + 1)) in P(1, 2, 3, 5): eliminate t; an identically vanishing resultant in
   rho with a common root for every rho means the whole line lies on the D4 locus.
2. If so: the parameter map rho -> a (with a^2 or a rational function of rho), so that the special points of the
   line (D4 -> larger groups, E10 = 0) can be read off.
Usage: python -I -B c17_d8_u0.py
"""
import importlib.util
import itertools
import sys
from pathlib import Path

import sympy as sp

sys.stdout.reconfigure(encoding="utf-8")
HERE = Path(__file__).resolve().parent


def _load(name):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


L = _load("k3lib")
R, S, U, rho, a = sp.symbols('R S U rho a')
r_, s_ = L.r_, L.s_
KEYS = (2, 4, 6, 10)
W = L.W


def homog(expr, k):
    return sp.expand(sp.expand(expr.subs({r_: R / U, s_: S / U}) * U ** k))


EH = {k: homog(L.EK8[k], W[k]) for k in KEYS}


def wpairs(x, y):
    return [sp.expand(x[i] ** W[j] * y[j] ** W[i] - x[j] ** W[i] * y[i] ** W[j]) for i, j in itertools.combinations(KEYS, 2)]


def main():
    E0 = {k: sp.expand(EH[k].subs({R: rho, S: 1, U: 0})) for k in KEYS}
    print("E on U = 0 at (rho : 1 : 0):", {k: sp.factor(v) for k, v in E0.items()})
    D4 = L.ic([0, 1, 0, a, 0, 1, 0])          # x (x^4 + a x^2 + 1)
    eqs = [e for e in wpairs(E0, D4) if e != 0]
    print(f"   {len(eqs)} weighted-projective equations in (rho, a)")
    g = sp.Poly(eqs[0], a)
    for e in eqs[1:]:
        g = sp.gcd(g, sp.Poly(e, a))
    print(f"   gcd in Q(rho)[a] of all equations: {sp.factor(g.as_expr())}")
    if g.degree() > 0:
        print("   => every point of U = 0 lies on the D4 locus: the curve over U = 0 is the D4 family.")
        sol = sp.solve(g.as_expr(), a)
        print(f"   a as a function of rho: {sol}")
    else:
        print("   => the general point of U = 0 is not on the D4 locus.")


if __name__ == "__main__":
    main()
