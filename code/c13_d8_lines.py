"""76.6: the curves of A-degree one on the model X of Y_-(8) (NOTE.md, Proposition 5.3(c)).
X -> P^2_{R:S:U} is the double cover branched along the quartic Q = U * B8, with B8 = [EK]'s branch cubic (k3lib.EK8_Y,
homogenised). A curve C0 with A.C0 = 1 maps birationally onto a line l with Q|_l a square (or l = {U = 0}).
1. The singular points of B8 and of Q, and their types (Hessian rank, cubic term).
2. The lines through each singular point of Q along which Q restricts to a square: the residual binary forms and
   their discriminants; the list of such lines.
3. The degree-15 invariant R30 (c08) on [EK]'s map E, homogenised: its factorisation, and its restriction to each
   line of 2 (identically zero means every curve over the line has an extra involution).
Usage: python -I -B c13_d8_lines.py
"""
import importlib.util
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
R, S, U, x, y, a = sp.symbols('R S U x y a')
r_, s_ = L.r_, L.s_


def homog(expr, k):
    """[EK]'s affine (r, s) expression of degree <= k as a form of degree k in (R, S, U)."""
    return sp.expand(sp.expand(expr.subs({r_: R / U, s_: S / U}) * U ** k))


def main():
    B8 = homog(L.EK8_Y, 3)
    Q = sp.expand(U * B8)
    print(f"B8 = {sp.factor(B8)}")
    # 1. singular points
    G = sp.groebner([B8, sp.diff(B8, R), sp.diff(B8, S), sp.diff(B8, U)], R, S, U, order='lex')
    print(f"1. singular locus of B8: Groebner basis {list(G.exprs)}")
    pc = {R: sp.Rational(1, 16), S: sp.Rational(27, 16), U: 1}
    print(f"   B8 and its gradient at (1/16 : 27/16 : 1): {[sp.expand(e.subs(pc)) for e in (B8, sp.diff(B8, R), sp.diff(B8, S), sp.diff(B8, U))]}")
    loc = sp.expand(B8.subs({U: 1, R: R + pc[R], S: S + pc[S]}))
    q2 = sum(c * R ** i * S ** j for (i, j), c in sp.Poly(loc, R, S).terms() if i + j == 2)
    q3 = sum(c * R ** i * S ** j for (i, j), c in sp.Poly(loc, R, S).terms() if i + j == 3)
    print(f"   local equation at the singular point: quadratic part {sp.factor(q2)}, cubic part {sp.factor(q3)}")
    lin = sp.factor_list(q2)[1][0][0]
    print(f"   quadratic part = const * ({lin})^2; cubic part on the line {lin} = 0 (non-zero means a cusp, A2): "
          f"{sp.simplify(q3.subs(S, sp.solve(lin, S)[0]))}")
    print(f"   Q at (0:1:0): B8 = {B8.subs({R: 0, S: 1, U: 0})}, U = 0; B8 on U = 0: {sp.factor(B8.subs(U, 0))}")
    # 2. lines through the singular points of Q along which Q is a square
    print("2. lines through the singular points of Q:")
    pencils = {
        "(0:1:0), R = a U": {R: a * y, S: x, U: y},
        "(1:-1:0), R + S = a U": {R: x, S: -x + a * y, U: y},
        "cusp (1/16 : 27/16 : 1), lines through it and (x : y) on the line at infinity": None,
    }
    for lab, sub in pencils.items():
        if sub is None:
            continue
        Ql = sp.factor(sp.expand(Q.subs(sub)))
        print(f"   {lab}: Q|_l = {Ql}")
    # residual forms: Q|_l = y^2 * (binary quadratic) for lines through (0:1:0) and (1:-1:0)
    for lab, sub in (("R = a U", {R: a * y, S: x, U: y}), ("R + S = a U", {R: x, S: -x + a * y, U: y})):
        res = sp.cancel(sp.expand(Q.subs(sub)) / y ** 2)
        disc = sp.factor(sp.discriminant(sp.expand(res.subs(y, 1)), x))
        print(f"   {lab}: Q|_l / y^2 = {sp.factor(res)}; discriminant in x: {disc}")
    # lines through the cusp: parametrise by the direction (p : q); Q|_l as a binary form in (lambda : mu)
    p, q, lam, mu = sp.symbols('p q lam mu')
    sub = {R: lam * pc[R] + mu * p, S: lam * pc[S] + mu * q, U: lam}
    Ql = sp.expand(Q.subs(sub))
    Ql2 = sp.cancel(Ql / mu ** 2)
    print(f"   lines through the cusp, direction (p : q : 0): Q|_l = mu^2 * ({sp.factor(Ql2)})")
    disc = sp.factor(sp.discriminant(sp.expand(Ql2.subs(mu, 1)), lam))
    print(f"   discriminant of the residual quadratic in (lam : mu): {disc}")
    lines = {"U = 0": U, "16R = U": 16 * R - U, "4R + 4S = 7U": 4 * R + 4 * S - 7 * U, "4R + 4S + U = 0 (TTV)": 4 * R + 4 * S + U}
    print("   squares check:")
    for lab, l in lines.items():
        if lab == "U = 0":
            print(f"     {lab}: a component of Q (branch line)")
            continue
        sol = sp.solve(l, R)[0]
        Ql = sp.factor(sp.expand(Q.subs(R, sol)))
        sq = all(e % 2 == 0 for _, e in sp.factor_list(sp.expand(Q.subs(R, sol)), S, U)[1])
        print(f"     {lab}: Q|_l = {Ql}; square up to a constant: {sq}")
    # 3. R30 on E
    R30, _, _ = C8.r30()
    E = {2: homog(L.EK8[2], 1), 4: homog(L.EK8[4], 2), 6: homog(L.EK8[6], 3), 10: homog(L.EK8[10], 5)}
    val = sp.expand(C8.at(R30, E))
    fl = sp.factor_list(val, R, S, U)
    print(f"3. R30(E(R, S, U)): constant {fl[0]}; factors: " + "; ".join(f"[{f}]^{e}" for f, e in fl[1]))
    for lab, l in lines.items():
        if lab == "U = 0":
            res = sp.expand(val.subs(U, 0))
        else:
            res = sp.expand(val.subs(R, sp.solve(l, R)[0]))
        print(f"   R30(E) on {lab}: {'identically zero' if res == 0 else 'not identically zero'}"
              + ("" if res == 0 else f"; factors {[(sp.factor(f), e) for f, e in sp.factor_list(res)[1]]}"))
    for lab, l in lines.items():
        if lab == "U = 0":
            res = sp.expand(E[10].subs(U, 0))
        else:
            res = sp.expand(E[10].subs(R, sp.solve(l, R)[0]))
        print(f"   E10 on {lab}: {sp.factor(res)}")


if __name__ == "__main__":
    main()
