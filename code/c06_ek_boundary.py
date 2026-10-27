"""76.4: the boundary of [EK]'s planes for D = 12, 13, 17, and the checks behind Theorems A13 and A17 (NOTE.md, Section 6).
Reads ../outputs/ek_ic_formulas.txt (written by c05).

1. Weighted homogenisation: N = the least polynomial with N^w I_{2k} polynomial; E_{2k} = N^w I_{2k} homogenised in
   (v1 : v2 : u) to degree m*w (w = 1, 2, 3, 5 for I2, I4, I6, I10). Base points; factors of E10.
2. Resolution of the base points by point blow-ups. Along every exceptional divisor X: the orders o_k of E_k,
   mu = min o_k / w, the excess o10 - 5 mu (= ord_t I10 at a transverse crossing of a general point, for a
   semistable fibre) and the restricted map; where the strict transforms of the special lines meet X.
3. The forbidden lines (g = 0 for D = 13, h = 0 for D = 17): at every point, ord_t I10 against Prop. 2.2.
Usage: python -I -B c06_ek_boundary.py
"""
import importlib.util
import sys
from fractions import Fraction
from pathlib import Path

import sympy as sp

sys.stdout.reconfigure(encoding="utf-8")
HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("k3lib", HERE / "k3lib.py")
L = importlib.util.module_from_spec(spec)
spec.loader.exec_module(L)
CACHE = HERE.parent / "outputs" / "ek_ic_formulas.txt"
W = L.W
u, x, y = sp.symbols('u x y')
THETA = (1, sp.Rational(1, 16), sp.Rational(1, 64), 0)


def load(D):
    rows = [ln.split(":", 2) for ln in CACHE.read_text(encoding="utf-8").splitlines() if ln.startswith(f"{D}:")]
    v = next(r[2] for r in rows if r[1] == "vars").split(",")
    I = {int(r[1][1:]): L.parse_safe(r[2], set(v)) for r in rows if r[1].startswith("I")}
    return sp.symbols(v), I


def homog(D):
    (a, b), I = load(D)
    N = sp.Integer(1)
    for k in W:
        for fac, e in sp.factor_list(sp.fraction(sp.cancel(I[k]))[1])[1]:
            need = -(-e // W[k])
            have = dict(sp.factor_list(N)[1]).get(fac, 0) if N != 1 else 0
            N *= fac ** max(0, need - have)
    P = {k: sp.expand(sp.cancel(N ** W[k] * I[k])) for k in W}
    m = max(-(-sp.Poly(P[k], a, b).total_degree() // W[k]) for k in W)
    E = {k: sp.expand(sum(c * a ** i * b ** j * u ** (m * W[k] - i - j) for (i, j), c in sp.Poly(P[k], a, b).terms()))
         for k in W}
    return (a, b), N, m, E


def wclass(pt):
    """Type of a point (I2:I4:I6:I10) of P(1,2,3,5) with I10 = 0."""
    I2, I4, I6, I10 = pt
    if I10 != 0:
        return "M2"
    if I4 == 0 and I6 == 0:
        return "(1:0:0:0)"
    if sp.simplify(I4 ** 3 - (3 * I6 - I2 * I4) ** 2) != 0:
        return "I10=0, off Gamma"
    if I2 != 0 and sp.simplify(I4 / I2 ** 2 - THETA[1]) == 0 and sp.simplify(I6 / I2 ** 3 - THETA[2]) == 0:
        return "theta"
    return "Gamma"


def ordv(p, v):
    P = sp.Poly(p, v)
    return min(mo[0] for mo in P.monoms())


def strip(p, divs):
    """p / (x^ox y^oy) for the coordinate divisors in divs, and the orders."""
    ox = ordv(p, x) if 'x' in divs else 0
    oy = ordv(p, y) if 'y' in divs else 0
    return sp.expand(sp.cancel(p / (x ** ox * y ** oy))), ox, oy


def is_base(Ek, divs):
    """True iff no normalised coordinate is non-zero at the origin (divs: {'x': mu, 'y': mu})."""
    for k in W:
        G, ox, oy = strip(Ek[k], divs)
        ax = ox - W[k] * divs['x'] if 'x' in divs else 0
        ay = oy - W[k] * divs['y'] if 'y' in divs else 0
        if ax == 0 and ay == 0 and G.subs({x: 0, y: 0}) != 0:
            return False
    return True


RESULT = []


def blowup(Ek, divs, lines, label, depth, maxdepth=6):
    """Blow up the origin. lines: name -> polynomial (pull-back of a special line, strict transform tracked)."""
    for chart in ("A", "B"):
        if chart == "A":
            sub, ex, par = {x: x, y: x * y}, x, y
        else:
            sub, ex, par = {x: x * y, y: y}, y, x
        P = {k: sp.expand(Ek[k].subs(sub, simultaneous=True)) for k in W}
        o = {k: ordv(P[k], ex) for k in W}
        mu = min(Fraction(o[k], W[k]) for k in W)
        name = label + chart
        # divisors through points of the new divisor: the new one, plus in chart A the old 'y', in chart B the old 'x'
        nd = {('x' if chart == "A" else 'y'): mu}
        old = 'y' if chart == "A" else 'x'
        if old in divs:
            nd[old] = divs[old]
        if chart == "A":
            Gr, ay = {}, {}
            for k in W:
                G, ox, oy = strip(P[k], nd)
                if ox - W[k] * mu == 0:
                    Gr[k] = sp.expand(G.subs(x, 0))
                    ay[k] = (oy - W[k] * nd['y']) if 'y' in nd else 0
            restr = {k: (Gr[k] * y ** ay[k] if k in Gr else sp.Integer(0)) for k in W}
            gen = [sp.factor(restr[k]) for k in W]
            yv = sp.Rational(7, 3)
            kind = wclass(tuple(restr[k].subs(y, yv) for k in W))
            ref = next(k for k in W if k in Gr)
            const = all(sp.simplify(sp.diff(restr[k] ** W[ref] / restr[ref] ** W[k], y)) == 0 for k in W)
            meets = {}
            for nm, Lp in lines.items():
                Lq = sp.expand(Lp.subs(sub, simultaneous=True))
                Ls = sp.expand(sp.cancel(Lq / ex ** ordv(Lq, ex))).subs(x, 0)
                meets[nm] = list(sp.Poly(Ls, y).ground_roots()) if sp.Poly(Ls, y).degree() > 0 else []
            RESULT.append((name, depth, [o[k] for k in W], mu, o[10] - 5 * mu, gen, kind, const))
            print(f"{'  ' * depth}[{name}] o = {[o[k] for k in W]}, mu = {mu}, excess = {o[10] - 5 * mu}, "
                  f"generic image: {kind}{', constant' if const else ', non-constant'}; restricted map (chart A) = {gen}")
            if any(meets.values()):
                print(f"{'  ' * depth}    strict transforms of the special lines meet it at {par} = {meets}")
            cands = None
            for k in Gr:
                cands = Gr[k] if cands is None else sp.gcd(cands, Gr[k])
            roots = set(sp.Poly(cands, y).ground_roots()) if sp.Poly(cands, y).degree() > 0 else set()
            if 'y' in nd:
                roots.add(sp.Integer(0))
            for c in sorted(roots, key=lambda r: (not r.is_rational, str(r))):
                Q = {k: sp.expand(P[k].subs(y, y + c)) for k in W}
                d2 = dict(nd) if c == 0 else {'x': mu}
                if not is_base(Q, d2):
                    continue
                if not c.is_rational:
                    print(f"{'  ' * depth}    base point at {par} = {c} (irrational, not followed)")
                    continue
                if depth + 1 >= maxdepth:
                    print(f"{'  ' * depth}    base point at {par} = {c} (depth limit)")
                    continue
                print(f"{'  ' * depth}    base point at {par} = {c}: blowing up")
                nl = {nm: sp.expand(sp.cancel(Lp.subs(sub, simultaneous=True) / ex ** ordv(sp.expand(Lp.subs(sub, simultaneous=True)), ex))).subs(y, y + c)
                      for nm, Lp in lines.items()}
                blowup(Q, d2, nl, name + f"({c})", depth + 1, maxdepth)
        else:
            Q = P
            if is_base(Q, nd):
                if depth + 1 >= maxdepth:
                    print(f"{'  ' * depth}    base point at the far end ({par} = 0) (depth limit)")
                    continue
                print(f"{'  ' * depth}    base point at the far end ({par} = 0 in chart B): blowing up")
                nl = {nm: sp.expand(sp.cancel(Lp.subs(sub, simultaneous=True) / ex ** ordv(sp.expand(Lp.subs(sub, simultaneous=True)), ex)))
                      for nm, Lp in lines.items()}
                blowup(Q, nd, nl, name + "(0)", depth + 1, maxdepth)


CHARTS = {  # D: list of (label, substitution putting the base point at the origin, special lines in (v1, v2, u))
    13: [("(1:0:0)", lambda a, b: {a: 1, b: x, u: y})],
    17: [("(-1/2,0)", lambda a, b: {a: (x - 1) / 2, b: y, u: 1}), ("(0:1:0)", lambda a, b: {a: x, b: 1, u: y})],
    12: [("(0,1)", lambda a, b: {a: x, b: 1 + y, u: 1}), ("(1:0:0)", lambda a, b: {a: 1, b: x, u: y})],
}


def special_lines(D, a, b):
    return {13: {"g=0": a, "h=0": b, "u=0": u}, 17: {"h=0": b, "u=0": u, "2g+u=0": 2 * a + u},
            12: {"e=0": a, "u=0": u, "f=u": b - u}}[D]


def main():
    for D in (13, 17, 12):
        (a, b), N, m, E = homog(D)
        print(f"== D = {D}: variables ({a}, {b}); N = {N}; m = {m}")
        print(f"   E10 = {sp.factor(E[10])}")
        G = sp.groebner([E[k].subs(u, 1) for k in W], a, b, order='lex')
        print(f"   affine base locus (lex Groebner basis): {G.exprs}")
        for pt in ({a: 1, b: y}, {a: y, b: 1}):
            print(f"   at infinity, ({pt[a]} : {pt[b]} : 0): E = {[sp.factor(E[k].subs(u, 0).subs(pt)) for k in W]}")
        for lab, sub in CHARTS[D]:
            s = sub(a, b)
            Ek = {k: sp.expand(E[k].subs(s, simultaneous=True)) for k in W}
            lines = {nm: sp.expand(Lp.subs(s, simultaneous=True)) for nm, Lp in special_lines(D, a, b).items()}
            lines = {nm: Lp for nm, Lp in lines.items() if Lp.subs({x: 0, y: 0}) == 0 and Lp != 0}
            print(f"   -- resolution of the base point {lab} (special lines through it: {list(lines)})")
            blowup(Ek, {}, lines, "", 2)
    forbidden_lines()


def forbidden_lines():
    print("== forbidden lines")
    # D = 13, g = 0
    (g, h), N, m, E = homog(13)
    print(f"  D = 13: E4 on g = 0: {sp.factor(E[4].subs(g, 0))} (non-zero for u != 0)")
    print(f"          E at (0:1:0): {[E[k].subs({g: 0, h: 1, u: 0}) for k in W]}; E10 = {sp.factor(E[10])}")
    # D = 17, h = 0
    (g, h), N, m, E = homog(17)
    print(f"  D = 17: E4 on h = 0: {sp.factor(E[4].subs(h, 0))}")
    print(f"          E at (g,h) = (-1, 0): {[E[k].subs({g: -1, h: 0, u: 1}) for k in W]}")
    print(f"          E at (1:0:0): {[E[k].subs({g: 1, h: 0, u: 0}) for k in W]} (theta point: {wclass(tuple(E[k].subs({g: 1, h: 0, u: 0}) for k in W))})")
    xs = {g: (x - 1) / 2, h: y, u: 1}
    Ex = {k: sp.expand(E[k].subs(xs, simultaneous=True)) for k in W}
    mults = {k: min(sum(mo) for mo in sp.Poly(Ex[k], x, y).monoms()) for k in W}
    print(f"          at p1 = (-1/2, 0), with x = 2g+1, y = h: multiplicities of E2, E4, E6, E10 = {[mults[k] for k in W]}")
    print(f"          E4 / x^4 at p1 = {sp.cancel(Ex[4] / x ** 4).subs({x: 0, y: 0})}; E10 = {sp.factor(Ex[10])}")
    print(f"          E6 / x^4: multiplicity {min(sum(mo) for mo in sp.Poly(sp.cancel(Ex[6] / x ** 4), x, y).monoms())}")
    # theta points and (1:0:0:0) points on h = 0 (affine)
    t = sp.Symbol('t')
    img = [sp.factor(E[k].subs({h: 0, u: 1}).subs(g, t)) for k in W]
    thet = sp.factor(sp.numer(sp.together(img[1] / img[0] ** 2 - THETA[1])))
    print(f"          image of h = 0: {img}; theta condition I4/I2^2 = 1/16 <=> {thet} = 0")
    print(f"  D = 17: elements of norm -3 in O_17: none, since (17/3) = {sp.jacobi_symbol(17, 3)} (3 inert)")
    # D = 12: the three lines e = 0, u = 0, f = u
    (e, f), N, m, E = homog(12)
    img = [sp.factor(E[k].subs({e: 0, u: 1})) for k in W]
    thet = sp.factor(sp.numer(sp.together(img[1] / img[0] ** 2 - THETA[1])))
    thet6 = sp.factor(sp.numer(sp.together(img[2] / img[0] ** 3 - THETA[2])))
    print(f"  D = 12: image of e = 0: {img}")
    print(f"          theta on e = 0 <=> {thet} = 0 and {thet6} = 0; (1:0:0:0) on e = 0 at f = 0")
    print(f"          image of f = u (affine e): {[sp.factor(E[k].subs({f: 1, u: 1})) for k in W]}")
    print(f"          image of u = 0: (e : f : 0) -> {[sp.factor(E[k].subs({u: 0})) for k in W]}")
    print(f"          units of O_12 = Z[sqrt3]: fundamental unit 2 + sqrt3 has norm {sp.expand((2 + sp.sqrt(3)) * (2 - sp.sqrt(3)))}")


if __name__ == "__main__":
    main()
