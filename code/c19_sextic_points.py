"""78.8 (the LaTeX note, Section 5): the invariant points of the singular semistable sextics used in the note.
1. (z^3 - 1)^2 (three double roots) has the point theta = (1 : 1/16 : 1/64 : 0); x^3 y^3 has (1 : 0 : 0 : 0).
2. (z^2 - 1)^2 (z^2 - c), c not in {0, 1} (every sextic with exactly two double roots is equivalent to one of
   these): I2, I4, I6 have no common zero; the point lies on Gamma = {I10 = 0, I4^3 = (3 I6 - I2 I4)^2} and is
   neither theta nor (1 : 0 : 0 : 0).
3. Exactly one double root: z^2 (z - 1)(z - a)(z - b), with the sixth root at infinity: the Gamma polynomial
   I4^3 - (3 I6 - I2 I4)^2 vanishes only where two of the roots 0, 1, a, b, infinity meet.
Usage: python -I -B c19_sextic_points.py
"""
import importlib.util
import sys
from pathlib import Path

import sympy as sp

sys.stdout.reconfigure(encoding="utf-8")
spec = importlib.util.spec_from_file_location("k3lib", Path(__file__).resolve().parent / "k3lib.py")
L = importlib.util.module_from_spec(spec)
spec.loader.exec_module(L)
X = L.X
c, a, b = sp.symbols('c a b')


def coeffs(f):
    P = sp.Poly(sp.expand(f), X)
    return [P.coeff_monomial(X ** i) for i in range(7)]


def gamma_poly(I):
    return sp.factor(I[4] ** 3 - (3 * I[6] - I[2] * I[4]) ** 2)


def main():
    print("== 1. theta and the triple-root point")
    T = L.ic(coeffs((X ** 3 - 1) ** 2))
    print(f"   (z^3 - 1)^2: {T}; I4/I2^2 = {sp.nsimplify(T[4] / T[2] ** 2)}, I6/I2^3 = {sp.nsimplify(T[6] / T[2] ** 3)}")
    # x^3 y^3 as the polynomial X^3 (a3 = 1, the other coefficients 0)
    U = L.ic([0, 0, 0, 1, 0, 0, 0])
    print(f"   x^3 y^3: {U}")
    print("== 2. two double roots: (z^2 - 1)^2 (z^2 - c)")
    I = L.ic(coeffs((X ** 2 - 1) ** 2 * (X ** 2 - c)))
    for k in (2, 4, 6, 10):
        print(f"   I{k} = {sp.factor(I[k])}")
    g = sp.gcd(sp.gcd(sp.Poly(I[2], c), sp.Poly(I[4], c)), sp.Poly(I[6], c))
    print(f"   gcd(I2, I4, I6) = {g.as_expr()} (common zeros only at c = 0, 1 allowed)")
    print(f"   Gamma polynomial: {gamma_poly(I)}")
    th4 = sp.factor(sp.numer(sp.together(I[4] * 16 - I[2] ** 2)))
    th6 = sp.factor(sp.numer(sp.together(I[6] * 64 - I[2] ** 3)))
    common = sp.gcd(sp.Poly(th4, c), sp.Poly(th6, c))
    print(f"   theta conditions 16 I4 - I2^2 = {th4}, 64 I6 - I2^3 = {th6}; common factor {common.as_expr()}")
    z46 = sp.gcd(sp.Poly(I[4], c), sp.Poly(I[6], c))
    print(f"   (1:0:0:0) condition I4 = I6 = 0: common factor {z46.as_expr()}")
    print("== 3. one double root: z^2 (z - 1)(z - a)(z - b), sixth root at infinity")
    J = L.ic(coeffs(X ** 2 * (X - 1) * (X - a) * (X - b)))
    print(f"   I10 = {sp.factor(J[10])}")
    print(f"   Gamma polynomial: {gamma_poly(J)}")


if __name__ == "__main__":
    main()
