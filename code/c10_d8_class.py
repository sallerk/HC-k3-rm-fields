"""76.6: the class of the moving D = 8 moduli curves on [EK]'s rational model P^2_{m,n} of Y_-(8) (NOTE.md, Section 5).
[EK] Section 8 (analysis): r = (m^2 - 1)/(32 n^2 - 16), s = s1 - r + 5/4 + 1/(32 r), s1 = m(16r - 1)/(32 r),
z = n(16 r - 1) with z^2 = 2(16rs^2 + 32r^2 s - 40rs - s + 16r^3 + 24r^2 + 12r + 2). Formulas typed from [EK] (data).
1. The homogeneous map Psi: (m : n : w) -> (R : S : U), its degree and base points; the multiplicity of the linear
   system at each base point (order of vanishing of a general member), and the image of the line m = -w.
2. Image degrees of test curves: general lines and conics; cubics with a node at a+ = (-1, 1/sqrt2) through
   a- = (-1, -1/sqrt2) (random members, parametrised by the lines through the node). Whether the image passes
   through (0 : 1 : 0), and its intersections with R = 0 and S = 0.
Usage: python -I -B c10_d8_class.py
"""
import random
import sys

import sympy as sp

sys.stdout.reconfigure(encoding="utf-8")
m, n, w, t, u, v = sp.symbols('m n w t u v')
SQ = sp.sqrt(2)


def psi():
    r = (m ** 2 - 1) / (32 * n ** 2 - 16)
    s1 = m * (16 * r - 1) / (32 * r)
    s = sp.cancel(s1 - r + sp.Rational(5, 4) + 1 / (32 * r))
    z = n * (16 * r - 1)
    b8 = 2 * (16 * r * s ** 2 + 32 * r ** 2 * s - 40 * r * s - s + 16 * r ** 3 + 24 * r ** 2 + 12 * r + 2)
    assert sp.simplify(sp.cancel(z ** 2 - b8)) == 0
    den = sp.lcm(sp.fraction(sp.cancel(r))[1], sp.fraction(sp.cancel(s))[1])
    P = [sp.expand(sp.cancel(r * den)), sp.expand(sp.cancel(s * den)), sp.expand(den)]
    g = sp.gcd(sp.gcd(P[0], P[1]), P[2])
    P = [sp.expand(sp.cancel(p / g)) for p in P]
    d = max(sp.Poly(p, m, n).total_degree() for p in P)
    return [sp.expand(sum(c * m ** i * n ** j * w ** (d - i - j) for (i, j), c in sp.Poly(p, m, n).terms())) for p in P], d


def mult_at(H, pt):
    """Order of vanishing at an affine point (w = 1) of a general member of the linear system spanned by H."""
    a, b = pt
    out = []
    for p in H:
        q = sp.expand(p.subs({w: 1, m: m + a, n: n + b}))
        out.append(min(i + j for (i, j), c in sp.Poly(q, m, n).terms() if sp.simplify(c) != 0) if q != 0 else 99)
    return min(out), out


def image_degree(H, mt, nt, wt=1):
    img = [sp.expand(p.subs({m: mt, n: nt, w: wt})) for p in H]
    g = sp.gcd(sp.gcd(img[0], img[1]), img[2])
    red = [sp.cancel(q / g) for q in img]
    return max(sp.degree(q, t) for q in red), red


def main():
    H, d = psi()
    print(f"== Psi: degree {d}")
    for name, p in zip("RSU", H):
        print(f"   {name} = {sp.factor(p)}")
    for lab, pt in (("a+", (-1, 1 / SQ)), ("a-", (-1, -1 / SQ)), ("b+", (1, 1 / SQ)), ("b-", (1, -1 / SQ))):
        print(f"   multiplicity at {lab} = {pt}: {mult_at(H, pt)}")
    print(f"   image of the line m = -w: {[sp.factor(p.subs(m, -w)) for p in H]}")
    print(f"   near (1:0:0) (m = 1): orders in (n, w): "
          f"{[min(i + j for (i, j), c in sp.Poly(sp.expand(p.subs(m, 1)), n, w).terms()) for p in H]}")
    rnd = random.Random(76)
    for k in (1, 2):
        mt = sum(rnd.randint(-5, 5) * t ** i for i in range(k + 1))
        nt = sum(rnd.randint(-5, 5) * t ** i for i in range(k + 1))
        print(f"== general curve of degree {k}: image degree {image_degree(H, mt, nt)[0]}")
    # rescale n' = sqrt2 n: the map has rational coefficients and a+-, b+- = (-+1, +-1)
    n2 = sp.Symbol('n2')
    Hq = [sp.expand(p.subs(n, n2 / SQ)) for p in H]
    assert all(c.is_rational for p in Hq for c in sp.Poly(p, m, n2, w).coeffs())
    print(f"   with n' = sqrt2 n: R, S, U = {[sp.factor(p) for p in Hq]}")
    # the line at infinity
    print(f"   on w = 0: {[sp.factor(p.subs(w, 0)) for p in Hq]}")
    for lab, pt in (("c+", 1), ("c-", -1)):
        loc = [sp.expand(p.subs({m: 1, n2: n2 + pt})) for p in Hq]
        print(f"   multiplicity at {lab} = (1 : {pt} : 0) (n' coordinates): "
              f"{[min(i + j for (i, j), c in sp.Poly(q, n2, w).terms()) for q in loc]}")
    loc = [sp.expand(p.subs(m, 1)) for p in Hq]
    print(f"   at (1:0:0): R, S, U mod (w, n'^2): {[sp.rem(sp.rem(q, w, w), n2 ** 2, n2) for q in loc]} "
          f"(base scheme (w, n'^2), length 2)")
    eps = {"a+": 2, "a-": 2, "b+": 1, "b-": 1, "oo": 1, "oo'": 1, "c+": 1, "c-": 1}
    print(f"   A = 4H - sum eps_i E_i with eps = {eps}: A^2 = {16 - sum(e * e for e in eps.values())} (degree of Psi = 2); "
          f"A.L = {4 - 2 - 2} for L = (m = -w); A.W = {4 - 1 - 1 - 1 - 1} for W = (w = 0)")
    print("   -K - A = -H + E_a+ + E_a-, so for C' = kH - sum mu_i E_i with A.C' = 6:")
    print("   -K.C' - 1 = 5 - (C'.L~) - ... ; C'.L~ = k - mu_a+ - mu_a-, C'.W~ = k - mu_oo - mu_oo' - mu_c+ - mu_c-.")
    print("   A curve avoiding L~ and W~ (both contracted to (0:1:0)) has -K.C' - 1 = k - mu_b+ - mu_b- - 1 = A.C' - 1 = 5.")
    # cubics with a node at a+ = (-1, 1) through a- = (-1, -1): C = C2(X, Y) + C3(X, Y), X = m + 1, Y = n' - 1;
    # lines through the node, (X, Y) = u (1, t), u = -C2(1, t)/C3(1, t):
    # (m : n' : w) = (-C2 - C3 : -t C2 + C3 : C3)
    X, Y = sp.symbols('X Y')
    for trial in range(3):
        c2 = [rnd.randint(-4, 4) for _ in range(3)]
        c3 = [rnd.randint(-4, 4) for _ in range(3)]
        C2 = c2[0] * X ** 2 + c2[1] * X * Y + c2[2] * Y ** 2
        C3 = c3[0] * X ** 3 + c3[1] * X ** 2 * Y + c3[2] * X * Y ** 2 + sp.Rational(c2[2], 2) * Y ** 3
        assert (C2 + C3).subs({X: 0, Y: -2}) == 0
        c2t, c3t = C2.subs({X: 1, Y: t}), C3.subs({X: 1, Y: t})
        mt, nt, wt = sp.expand(-c2t - c3t), sp.expand(-t * c2t + c3t), sp.expand(c3t)
        img = [sp.expand(p.subs({m: mt, n2: nt, w: wt})) for p in Hq]
        g = sp.gcd(sp.gcd(img[0], img[1]), img[2])
        Rt, St, Ut = [sp.expand(sp.cancel(q / g)) for q in img]
        dd = max(sp.degree(q, t) for q in (Rt, St, Ut))
        print(f"== nodal cubic {trial}: C2 = {C2}, C3 = {C3}; image degree {dd} (common factor of degree {sp.degree(g, t)})")
        print(f"   R(t): degree {sp.degree(Rt, t)}, factors (degree, multiplicity) {[(sp.degree(f_, t), e) for f_, e in sp.factor_list(Rt)[1]]}")
        print(f"   S(t): degree {sp.degree(St, t)}, factors (degree, multiplicity) {[(sp.degree(f_, t), e) for f_, e in sp.factor_list(St)[1]]}")
        print(f"   common zeros of R and U, i.e. passages through (0:1:0): degree {sp.degree(sp.gcd(Rt, Ut), t)}; "
              f"of R and S (theta corner (0:0:1)): degree {sp.degree(sp.gcd(Rt, St), t)}")


if __name__ == "__main__":
    main()
