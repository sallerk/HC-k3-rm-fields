"""76.2(c)/76.5(b): exact checks behind the D = 8 obstruction for smooth branch sextics, a consistency check on a
D = 8 pencil with singular branch curve, and the D = 5 calibration on the [K3M] example.

A. The Elkies-Kumar map E: P^2_{R:S:U} -> P(1,2,3,5) (Theorem 17, homogenised) is a morphism; E10 = -8 S^2 R^3; the
   line R = 0 maps onto the curve Gamma = {I10 = 0, I4^3 = (3 I6 - I2 I4)^2} of 2-nodal sextics, with (0:0:1) -> theta
   point (z^3 - 1)^2 and (0:1:0) -> the triple-root point (closed orbit x^3 y^3); the line S = 0 maps onto Gamma too.
B. The pencils y^2 = (z + 2c)(z^4 - 4c^2 z^2 + 2c^4 + S), c linear and S quintic in t (Tautz-Top-Verberkmoes n = 4,
   RM by Z[sqrt 2] over the base): the zeros of R along the pencil are the t with S + 2c^4 = 0, where the line of the
   pencil passes through a node of the (reducible) branch curve, and there ord_t I10 = 3.
C. [EK] Corollary 15 (D = 5) checked on its own table; along the smooth [K3M] example pencil the D = 5 coordinates
   (g, h) are rational in t, I4(t) is a square, and I10 = const * h-part^2: no factor of odd order appears.
Usage: python -I -B c03_d8_obstruction.py
"""
import importlib.util
import random
import sys
from pathlib import Path

import sympy as sp

sys.stdout.reconfigure(encoding="utf-8")
spec = importlib.util.spec_from_file_location("k3lib", Path(__file__).resolve().parent / "k3lib.py")
L = importlib.util.module_from_spec(spec)
spec.loader.exec_module(L)
X, r, s = L.X, L.r_, L.s_
R, S, U, t, c, sig, w = sp.symbols('R S U t c sigma w')


def hom(e, k):
    """Homogenise an EK8 polynomial of weight k (degree k in r, s) to a form of degree k in R, S, U."""
    return sp.expand(U ** k * e.subs({r: R / U, s: S / U}, simultaneous=True))


def part_a():
    print("== A. Elkies-Kumar D = 8 map")
    E = {k: hom(L.EK8[k], L.W[k]) for k in L.EK8}
    print("  E10 =", sp.factor(E[10]))
    # morphism: no common zero of E2, E4, E6, E10 on P^2 (check on R = 0 and S = 0, where E10 vanishes)
    for name, sub in (("R = 0", {R: 0}), ("S = 0", {S: 0})):
        es = [sp.factor(E[k].subs(sub)) for k in (2, 4, 6)]
        print(f"  on {name}: E2, E4, E6 =", es)
    gam = lambda I: sp.expand(I[4] ** 3 - (3 * I[6] - I[2] * I[4]) ** 2)
    for name, sub in (("R = 0", {R: 0}), ("S = 0", {S: 0})):
        I = {k: E[k].subs(sub) for k in (2, 4, 6)}
        print(f"  Gamma equation I4^3 - (3 I6 - I2 I4)^2 restricted to {name}:", gam(I))
    I2n = L.ic(L.poly_coeffs((X ** 2 - 1) ** 2 * (X ** 2 - c), X))
    print("  2-nodal family (z^2-1)^2 (z^2-c): I10 =", sp.factor(I2n[10]), "; Gamma equation:", sp.factor(gam(I2n)))
    th = L.ic(L.poly_coeffs((X ** 3 - 1) ** 2, X))
    print("  theta (z^3-1)^2 on Gamma:", gam(th) == 0, "; E at (0:0:1):", [E[k].subs({R: 0, S: 0, U: 1}) for k in (2, 4, 6, 10)],
          "; theta IC:", [th[k] for k in (2, 4, 6, 10)])
    x3 = L.ic([0, 0, 0, 1, 0, 0, 0])
    print("  closed triple-root orbit x^3 y^3: IC =", [x3[k] for k in (2, 4, 6, 10)],
          "; E at (0:1:0):", [E[k].subs({R: 0, S: 1, U: 0}) for k in (2, 4, 6, 10)])
    # which points of P^2 map to (1:0:0:0): E4 = E6 = E10 = 0
    sols = []
    for sub in ({R: 0}, {S: 0}):
        e4, e6 = sp.expand(E[4].subs(sub)), sp.expand(E[6].subs(sub))
        sols.append(sp.factor(sp.gcd(sp.Poly(e4, R, S, U), sp.Poly(e6, R, S, U)).as_expr()))
    print("  common factor of E4, E6 on R = 0 and on S = 0 (preimages of (1:0:0:0)):", sols)
    # TTV n = 4 identity, symbolic in sigma: IC((x+2)(x^4-4x^2+2+sigma)) is proportional in P(1,2,3,5) to
    # EK8(r, s) with r = -(sigma+2)/16, s = (sigma-2)/16
    I = L.ic(L.poly_coeffs((X + 2) * (X ** 4 - 4 * X ** 2 + 2 + sig), X))
    J = {k: L.EK8[k].subs({r: -(sig + 2) / 16, s: (sig - 2) / 16}) for k in L.EK8}
    lam = sp.factor(I[2] / J[2])
    ok = all(sp.simplify(I[k] - lam ** L.W[k] * J[k]) == 0 for k in (2, 4, 6, 10))
    print("  TTV n = 4 symbolic: IC = lam^w EK8(-(sigma+2)/16, (sigma-2)/16) with lam =", lam, ":", ok)
    z2 = sp.factor(L.EK8_Y.subs({r: -(sig + 2) / 16, s: (sig - 2) / 16}))
    print("  Y_-(8) branch function along the TTV line:", z2)


def part_b():
    print("== B. TTV-type pencils with singular branch curve (line + quintic)")
    rng = random.Random(76)
    for trial in range(3):
        cc = rng.randint(1, 5) + rng.randint(1, 5) * t
        SS = sum(rng.randint(-6, 6) * t ** k for k in range(6))
        f = sp.expand((X + 2 * cc) * (X ** 4 - 4 * cc ** 2 * X ** 2 + 2 * cc ** 4 + SS))
        I = L.ic(L.poly_coeffs(f, X))
        I10 = sp.factor(I[10])
        fl = sp.factor_list(sp.Poly(I[10], t))
        mults = {}
        for fac, m in fl[1]:
            mults[sp.factor(fac.as_expr())] = m
        Rz = sp.Poly(sp.expand(SS + 2 * cc ** 4), t)
        Sz = sp.Poly(sp.expand(SS - 2 * cc ** 4), t)
        mR = [m for fac, m in fl[1] if sp.rem(Rz, fac) == 0 or sp.gcd(Rz, fac).degree() > 0]
        mS = [m for fac, m in fl[1] if sp.gcd(Sz, fac).degree() > 0]
        print(f"  trial {trial}: c = {cc}, S = {SS}")
        print(f"    deg I10 = {sp.degree(I[10], t)}; multiplicities of the factors of I10 sharing roots with S + 2c^4: {mR},"
              f" with S - 2c^4: {mS}; all factor multiplicities: {sorted(m for _, m in fl[1])}")
        # nodes of the branch curve: line z = -2c meets the quintic where S + 2c^4 = 0
        quint = sp.expand(X ** 4 - 4 * cc ** 2 * X ** 2 + 2 * cc ** 4 + SS)
        on_line = sp.expand(quint.subs(X, -2 * cc) - (SS + 2 * cc ** 4))
        print("    quintic restricted to the line z = -2c minus (S + 2c^4):", on_line)


def part_c():
    print("== C. D = 5 calibration")
    g, h = sp.symbols('g h')
    EK5 = {2: 6 * (4 * g + 1), 4: 9 * g ** 2, 6: 9 * (4 * h + 9 * g ** 3 + 2 * g ** 2), 10: 4 * h ** 2}
    Q_ = sp.Rational
    table5 = [  # the sixteen rows of [EK] Table 1 (Section 6.3), in the published version
        ((Q_(-8, 3), Q_(47, 2)), -X**5 + X**4 - X**3 - X**2 + 2*X - 1),
        ((Q_(-37, 6), 67), X**6 + 2*X**5 + X**4 - 2*X**3 + 2*X**2 - 4*X + 1),
        ((0, Q_(27, 50)), 3*X**5 + 5*X**3 + 1),
        ((Q_(-1, 24), Q_(1, 100)), 4*X**6 + 4*X**5 + 5*X**4 - 5*X**2 - 2*X - 2),
        ((Q_(-25, 54), Q_(59, 81)), -X**6 - 2*X**4 - 6*X**3 - 5*X**2 - 6*X - 1),
        ((Q_(-2, 3), Q_(-14, 25)), -7*X**6 - 7*X**5 - 5*X**4 + 5*X**2 - X - 1),
        ((Q_(47, 54), Q_(71, 81)), 3*X**6 - 6*X**5 + 7*X**4 - 2*X**3 - 2*X**2 - 1),
        ((Q_(-1, 6), Q_(1, 25)), X**6 - 4*X**5 + 10*X**4 - 10*X**3 + 5*X**2 + 2*X - 3),
        ((Q_(-4, 3), Q_(16, 25)), -X**6 - X**5 + 5*X**2 - 7*X - 12),
        ((Q_(-4, 3), Q_(49, 22)), 7*X**5 + 5*X**4 + 3*X**3 - 9*X**2 - 14*X - 7),
        ((0, Q_(-54, 11)), 2*X**6 - 6*X**5 - 6*X**4 + 3*X**3 - 18*X**2 - 6*X - 2),
        ((Q_(11, 6), Q_(53, 11)), 9*X**6 - 14*X**5 + 13*X**4 - 2*X**3 - 22*X**2 - 8*X - 7),
        ((Q_(-5, 24), Q_(1, 64)), 6*X**6 - 2*X**5 - 15*X**4 - 16*X**3 - 25*X**2 - 8*X - 4),
        ((Q_(11, 6), Q_(-89, 25)), -X**6 + 2*X**5 - 5*X**4 + 30*X**3 - 10*X**2 + 8*X - 1),
        ((Q_(-2, 3), Q_(68, 11)), -2*X**6 - 2*X**5 - 11*X**4 - 29*X**3 - 31*X**2 - 26*X - 6),
        ((Q_(-2, 3), Q_(26, 25)), -2*X**6 + 36*X**5 + 5*X**4 + 35*X**3 - 10*X**2 - 21*X - 17),
    ]
    def ek5_points(I):
        """(g, h) over Q with IC = EK5(g, h) in P(1,2,3,5): g = u/(1-4u), u = +-2 sqrt(I4)/I2, h from I6, test I10."""
        out = []
        q = sp.sqrt(sp.Rational(I[4]))
        if not q.is_rational:
            return out
        for u in {2 * q / I[2], -2 * q / I[2]}:
            if 1 - 4 * u == 0:
                continue
            gv = u / (1 - 4 * u)
            mu = sp.Rational(I[2]) / (6 * (4 * gv + 1))
            hv = (sp.Rational(I[6]) / mu ** 3 / 9 - 9 * gv ** 3 - 2 * gv ** 2) / 4
            if sp.Rational(I[10]) == mu ** 5 * 4 * hv ** 2:
                out.append((gv, hv))
        return out

    for (gv, hv), f in table5:
        I = L.ic(L.poly_coeffs(f, X))
        J = {k: EK5[k].subs({g: gv, h: hv}) for k in EK5}
        mu = sp.nsimplify(I[2] / J[2])
        ok = all(sp.simplify(I[k] - mu ** L.W[k] * J[k]) == 0 for k in (2, 4, 6, 10))
        print(f"  [EK] Cor. 15 at listed (g,h) = ({gv},{hv}): match = {ok}, mu = {mu}; rational (g,h) of the sextic:"
              f" {ek5_points(I)}")
    # [K3M] example: 4F = 4x0^6 + 40x0^3x1^2x2 - 20x0^3x2^3 + 4x0x1^5 - 20x0x1^4x2 + 4x1^6 - 20x1^3x2^3 + 4x1x2^5 + 5x2^6
    f = 4 * t ** 6 + 40 * t ** 3 * X - 20 * t ** 3 * X ** 3 + 4 * t - 20 * t * X + 4 - 20 * X ** 3 + 4 * X ** 5 + 5 * X ** 6
    I = L.ic(L.poly_coeffs(f, X))
    print("  [K3M] example pencil: deg_t I2, I4, I6, I10 =", [sp.degree(I[k], t) for k in (2, 4, 6, 10)])
    cf, fl4 = sp.factor_list(sp.Poly(I[4], t))
    assert all(m % 2 == 0 for _, m in fl4) and sp.sqrt(cf).is_rational
    sq4 = sp.sqrt(cf) * sp.prod([fa.as_expr() ** (m // 2) for fa, m in fl4])
    print("  I4(t) =", sp.factor(I[4]), "| = (", sp.factor(sq4), ")^2:", sp.expand(sq4 ** 2 - I[4]) == 0)
    # g from I4/I2^2 = g^2/(4(4g+1)^2): u = g/(4g+1) = +-2 sqrt(I4)/I2, g = u/(1-4u)
    for sgn in (1, -1):
        u = sgn * 2 * sq4 / I[2]
        gt = sp.cancel(u / (1 - 4 * u))
        mu = sp.cancel(I[2] / (6 * (4 * gt + 1)))
        ht = sp.cancel((I[6] / mu ** 3 / 9 - 9 * gt ** 3 - 2 * gt ** 2) / 4)
        ok10 = sp.cancel(I[10] - mu ** 5 * 4 * ht ** 2) == 0
        print(f"  sign {sgn}: I10 = mu^5 4 h^2: {ok10}")
        if ok10:
            print("     g(t) =", sp.factor(gt))
            print("     h(t) =", sp.factor(ht))
            print("     mu(t) =", sp.factor(mu))
            num_h, den_h = sp.fraction(sp.factor(ht))
            print("     numerator of h: degree", sp.degree(num_h, t), "factor multiplicities",
                  sorted(m for _, m in sp.factor_list(sp.Poly(num_h, t))[1]))
    fl = sp.factor_list(sp.Poly(I[10], t))
    print("  I10(t) factor multiplicities:", sorted(m for _, m in fl[1]), "; total degree", sp.degree(I[10], t))


if __name__ == "__main__":
    L.ic_polys()
    if sys.argv[1:] == ["c"]:
        part_c()
    else:
        part_a()
        part_b()
        part_c()
