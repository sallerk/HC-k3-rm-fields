"""78.7b (first step): the smooth special points on the lines L1: 16R = U and L2: 4R + 4S = 7U of P^2_{R:S:U}, and,
for comparison, on the TTV line 4R + 4S + U = 0 (NOTE.md, Section 5.4).
A smooth special point is a point of the line with E10 != 0 and R30(E) = 0, so the curve there is smooth with an
extra involution. For each one:
1. the order of R30(E) along the line (c13 factorisation, recomputed here);
2. whether the curve has a larger automorphism group, tested exactly: the point lies on the image of the D4 family
   x^5 + x^3 + t x, of the D6 family x^6 + x^3 + t, or equals one of x^5 - x (order 48), x^6 - 1 (order 24),
   x^5 - 1 (order 10), in P(1, 2, 3, 5). The test eliminates t by resultants from the weighted-projective equations
   E_i(s)^{w_j} I_j(t)^{w_i} = E_j(s)^{w_i} I_i(t)^{w_j}, keeps the common factors with the point's minimal
   polynomial, and confirms each hit with every pair (i, j): exactly (a gcd over Q) at a rational point, and
   numerically (40 digits) at other points; the negative answers come from the exact elimination.
Usage: python -I -B c15_d8_special.py
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
C8 = _load("c08_parity")
R, S, U, s, t = sp.symbols('R S U s t')
r_, s_ = L.r_, L.s_
KEYS = (2, 4, 6, 10)
W = L.W

# each line as (R, S, U) in terms of the affine parameter s (point at s = oo handled separately)
LINES = {
    "L1: 16R = U": {R: sp.Rational(1, 16), S: s, U: 1},
    "L2: 4R + 4S = 7U": {R: sp.Rational(7, 4) - s, S: s, U: 1},
    "TTV: 4R + 4S + U = 0": {R: -sp.Rational(1, 4) - s, S: s, U: 1},
}
# the point at infinity of each line, approached with S = 1, U = u -> 0
LINES_INF = {
    "L1: 16R = U": {R: s / 16, S: 1, U: s},
    "L2: 4R + 4S = 7U": {R: -1 + 7 * s / 4, S: 1, U: s},
    "TTV: 4R + 4S + U = 0": {R: -1 - s / 4, S: 1, U: s},
}


def homog(expr, k):
    return sp.expand(sp.expand(expr.subs({r_: R / U, s_: S / U}) * U ** k))


EH = {k: homog(L.EK8[k], W[k]) for k in KEYS}
R30, _, _ = C8.r30()
R30E = sp.expand(C8.at(R30, EH))

FAMILIES = {
    "D4 family x^5 + x^3 + t x": [0, t, 0, 1, 0, 1, 0],
    "D6 family x^6 + x^3 + t": [t, 0, 0, 1, 0, 0, 1],
}
ISOLATED = {
    "x^5 - x (|Aut| = 48)": [0, -1, 0, 0, 0, 1, 0],
    "x^6 - 1 (|Aut| = 24)": [-1, 0, 0, 0, 0, 0, 1],
    "x^5 - 1 (|Aut| = 10)": [-1, 0, 0, 0, 0, 1, 0],
}


def wpairs(a, b):
    """Weighted-projective equality conditions for P(1, 2, 3, 5)."""
    return [sp.expand(a[i] ** W[j] * b[j] ** W[i] - a[j] ** W[i] * b[i] ** W[j]) for i, j in itertools.combinations(KEYS, 2)]


def ic_w(coeffs):
    """Invariants of a sextic in the (1, 2, 3, 5) normalisation of [EK] (I_{2k} has weight k)."""
    return L.ic(coeffs)


def special_points(sub):
    """Irreducible factors p(s) of R30(E) on the line with E10 != 0 there, with their multiplicities."""
    r30 = sp.factor_list(sp.expand(R30E.subs(sub)), s)
    e10 = sp.expand(EH[10].subs(sub))
    out = []
    for fa, mu in r30[1]:
        if sp.Poly(fa, s).degree() == 0:
            continue
        if sp.rem(sp.Poly(e10, s), sp.Poly(fa, s)).is_zero:
            continue
        out.append((fa, mu))
    return out, r30


def test_family(p, Ep, coeffs):
    """Does the curve at the point p(s) = 0 (invariants Ep, polynomials in s) lie on the family's image?"""
    It = ic_w(coeffs)
    eqs = [e for e in wpairs(Ep, It) if e != 0]
    P = sp.Poly(p, s)
    hits = []
    # eliminate t between pairs of equations, keep common factors with p
    for e1, e2 in itertools.combinations(eqs, 2):
        if not (e1.has(t) and e2.has(t)):
            continue
        res = sp.resultant(e1, e2, t)
        if sp.expand(res) == 0:
            continue
        g = sp.gcd(sp.Poly(res, s), P)
        if g.degree() > 0:
            hits.append(g.as_expr())
            break
        return False, "resultant coprime to p"
    if not hits:
        return False, "no t-elimination possible"
    if P.degree() == 1:
        # rational point: exact confirmation, the gcd over Q of the equations in t
        s0 = sp.solve(p, s)[0]
        gt = None
        for e in eqs:
            if not e.has(t):
                if sp.simplify(e.subs(s, s0)) != 0:
                    return False, "an equation without t fails"
                continue
            q = sp.Poly(sp.expand(e.subs(s, s0)), t)
            if q.is_zero:
                continue
            gt = q if gt is None else sp.gcd(gt, q)
        if gt is None or gt.degree() == 0:
            return False, "no common root in t (exact)"
        return True, [f"exact: t a root of {sp.factor(gt.as_expr())}"]
    # confirm: over K = Q(s)/(p), the equations in t have a common root
    gens = [sp.Poly(e, t) for e in eqs if e.has(t)]
    # reduce coefficients modulo p and take the gcd in t over the number field via a generic numeric check
    roots = sp.Poly(p, s).nroots(n=40)
    ok = []
    for rt in roots:
        polys = [sp.Poly(sp.N(e.subs(s, rt), 40), t) for e in eqs if e.has(t)]
        cand = polys[0].nroots(n=40, maxsteps=200)
        found = None
        for c in cand:
            if all(abs(sp.N(q.as_expr().subs(t, c), 40)) < sp.Float(10) ** -20 * (1 + max(abs(cc) for cc in q.all_coeffs()))
                   for q in polys):
                found = c
                break
        ok.append(found)
    return all(v is not None for v in ok), ok


def test_isolated(p, Ep, coeffs):
    Ic = ic_w(coeffs)
    eqs = [e for e in wpairs(Ep, Ic)]
    P = sp.Poly(p, s)
    return all(sp.rem(sp.Poly(e, s), P).is_zero for e in eqs)


def main():
    for name, sub in LINES.items():
        print(f"== {name}")
        sp_pts, r30 = special_points(sub)
        print(f"   R30(E) on the line: constant {r30[0]}, factors {[(f, m) for f, m in r30[1]]}")
        Ep = {k: sp.expand(EH[k].subs(sub)) for k in KEYS}
        for p, mu in sp_pts:
            line = f"   point {p} = 0: ord R30 = {mu}"
            fam = []
            for fname, co in FAMILIES.items():
                hit, info = test_family(p, Ep, co)
                if hit:
                    if all(isinstance(v, str) for v in info):
                        fam.append(f"on the {fname} ({'; '.join(info)})")
                    else:
                        fam.append(f"on the {fname} (t = {[sp.N(v, 12) for v in info]}, numerical)")
            for iname, co in ISOLATED.items():
                if test_isolated(p, Ep, co):
                    fam.append(f"equal to {iname}")
            print(line + ("; " + "; ".join(fam) if fam else "; on neither the D4 nor the D6 family, none of the isolated curves: Aut = V4"))
        # the point at infinity
        sub_inf = LINES_INF[name]
        Ei = {k: sp.expand(EH[k].subs(sub_inf)) for k in KEYS}
        e10_0 = Ei[10].subs(s, 0)
        r30_inf = sp.expand(R30E.subs(sub_inf))
        ordinf = min(m[0] for m in sp.Poly(r30_inf, s).monoms()) if r30_inf != 0 else None
        print(f"   point at infinity (R:S:U) = ({sub_inf[R].subs(s, 0)} : 1 : 0): E10 = {e10_0}, ord R30 = {ordinf}")
        if e10_0 != 0:
            pt = {k: Ei[k].subs(s, 0) for k in KEYS}
            for iname, co in ISOLATED.items():
                Ic = ic_w(co)
                if all(sp.expand(e) == 0 for e in wpairs(pt, Ic)):
                    print(f"      the curve there is {iname}")
            for fname, co in FAMILIES.items():
                It = ic_w(co)
                eqs = [sp.expand(e) for e in wpairs(pt, It)]
                eqs = [e for e in eqs if e != 0]
                g = eqs[0]
                for e in eqs[1:]:
                    g = sp.gcd(g, e)
                print(f"      {fname}: common factor in t of the equations: {sp.factor(g)}")


if __name__ == "__main__":
    main()
