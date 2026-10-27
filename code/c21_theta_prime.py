"""78.8 (the LaTeX note, Section 5): the second point theta' of Gamma with 16 I4 = I2^2.
1. On Gamma = {I10 = 0, I4^3 = (3 I6 - I2 I4)^2}, the condition 16 I4 = I2^2 forces I2 != 0 and
   I6 / I2^3 in {1/64, 5/192}: the points theta = (1 : 1/16 : 1/64 : 0) and theta' = (1 : 1/16 : 5/192 : 0).
2. The points of the note where 16 I4 = I2^2 holds and the point is not theta are theta':
   (z^2 - 1)^2 (z^2 - c) at c^2 + 4c + 1 = 0; for D = 8, E(0 : s : 1) at s = 4/3; for D = 12, E(0 : f : 1) at
   f^2 = 3 (the line e = 0), and the restriction of the map to the exceptional curve Z1 at y = 3/2.
Usage: python -I -B c21_theta_prime.py
"""
import importlib.util
import sys
from pathlib import Path

import sympy as sp

sys.stdout.reconfigure(encoding="utf-8")
HERE = Path(__file__).resolve().parent


def load(name):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


L = load("k3lib")
B = load("c06_ek_boundary")
THETA = (1, sp.Rational(1, 16), sp.Rational(1, 64), 0)
THETA1 = (1, sp.Rational(1, 16), sp.Rational(5, 192), 0)


def normal(pt):
    """(1 : I4/I2^2 : I6/I2^3 : I10/I2^5) for a point (I2 : I4 : I6 : I10) of P(1,2,3,5) with I2 != 0."""
    I2, I4, I6, I10 = (sp.radsimp(sp.expand(v)) for v in pt)
    assert I2 != 0
    return tuple(sp.radsimp(sp.expand(v)) for v in (I2 / I2, I4 / I2 ** 2, I6 / I2 ** 3, I10 / I2 ** 5))


def z(v):
    return sp.simplify(sp.expand(v))


def name(pt):
    n = normal(pt)
    if all(sp.simplify(a - b) == 0 for a, b in zip(n, THETA)):
        return "theta"
    if all(sp.simplify(a - b) == 0 for a, b in zip(n, THETA1)):
        return "theta'"
    return f"other {n}"


def part1():
    print("== 1. Gamma and 16 I4 = I2^2")
    i2, i4, i6 = sp.symbols('I2 I4 I6')
    # I2 = 0 would give I4 = 0 and then 3 I6 = 0 on Gamma: not a point of P(1,2,3,5) with I10 = 0
    print(f"   I2 = 0: I4 = 0 and Gamma gives (3 I6)^2 = 0, so I6 = 0: no point")
    sol = sp.solve([sp.Eq(i4, i2 ** 2 / 16), sp.Eq(i4 ** 3, (3 * i6 - i2 * i4) ** 2)], [i4, i6], dict=True)
    print(f"   solutions: {[{k: sp.factor(v) for k, v in s.items()} for s in sol]}")
    print(f"   I6 / I2^3 = {sorted(set(sp.cancel(s[i6] / i2 ** 3) for s in sol))}")


def coeffs(f):
    P = sp.Poly(sp.expand(f), L.X)
    return [P.coeff_monomial(L.X ** i) for i in range(7)]


def part2():
    print("== 2. the points of the note")
    c = sp.Symbol('c')
    I = L.ic(coeffs((L.X ** 2 - 1) ** 2 * (L.X ** 2 - c)))
    for c0 in [sp.Integer(0)] + sp.solve(c ** 2 + 4 * c + 1, c):
        pt = [I[k].subs(c, c0) for k in L.W]
        print(f"   (z^2 - 1)^2 (z^2 - c), c = {c0}: 16 I4 - I2^2 = {z(16 * pt[1] - pt[0] ** 2)}; point {name(pt)}")
    R, S, U = sp.symbols('R S U')
    E8 = {k: sp.expand(U ** L.W[k] * L.EK8[k].subs({L.r_: R / U, L.s_: S / U}, simultaneous=True)) for k in L.EK8}
    for s0 in (sp.Integer(0), sp.Rational(4, 3)):
        pt = [E8[k].subs({R: 0, S: s0, U: 1}) for k in L.W]
        print(f"   D = 8, E(0 : {s0} : 1) = {pt}: 16 E4 - E2^2 = {16 * pt[1] - pt[0] ** 2}; point {name(pt)}")
    (e, f), N, m, E = B.homog(12)
    u = B.u
    for f0 in [sp.Integer(-1), sp.sqrt(3), -sp.sqrt(3)]:
        pt = [E[k].subs({e: 0, f: f0, u: 1}) for k in L.W]
        print(f"   D = 12, E(0 : {f0} : 1): 16 E4 - E2^2 = {z(16 * pt[1] - pt[0] ** 2)}; point {name(pt)}")
    y = sp.Symbol('y')
    Z1 = [8 * (y - 3), 4 * y ** 2, 8 * y ** 2 * (y - 4), 0]
    print(f"   D = 12, on Z1 (the restriction of the note): 16 I4 - I2^2 = {sp.factor(16 * Z1[1] - Z1[0] ** 2)}")
    pt = [v.subs(y, sp.Rational(3, 2)) if v != 0 else 0 for v in Z1]
    print(f"   D = 12, Z1 at y = 3/2: {pt}; point {name(pt)}")


def main():
    part1()
    part2()


if __name__ == "__main__":
    main()
