"""76.6: effective dimension of the [TTV] n = 4 double planes of NOTE.md Check 3.3,
    F = (x2 + 2c) (x1 (x2^4 - 4 c^2 x2^2 + 2 c^4) + S(x0, x1)),  c = c0 x1 + c1 x0,  S a quintic form,
i.e. f_t(z) = (z + 2c)(z^4 - 4c^2 z^2 + 2c^4 + S) with c linear and S quintic in t (P0 = (0:0:1) on the quintic).
Effective dimension = rank [dF/dparams | x_i dF/dx_j] - rank [x_i dF/dx_j] over Q at random rational parameters
(the second block spans the tangent space of the PGL3 x scalars orbit). For comparison the same number for all
double planes branched along a line plus a quintic (expected 14).
Also (NOTE.md, Theorem C, the real multiplication): symbolic checks of the genus-4 cover
Y: v^2 = x(x^8 + s x^4 + 1) -> T: eta^2 = (xi + 2)(D4(xi) + s), D4(xi) = xi^4 - 4 xi^2 + 2, xi = x + 1/x,
eta = v (1 + 1/x)/x^2, the automorphisms g = (i x, z8 v) and iota = (1/x, v/x^5) of Y, iota g iota = g^-1, the
action on the differentials x^j dx/v, and the substitution z = c xi turning the fibres of F into twists of T.
Usage: python -I -B c11_ttv_rank.py [seed]
"""
import random
import sys

import sympy as sp
from flint import fmpq_mat, fmpq

sys.stdout.reconfigure(encoding="utf-8")
x0, x1, x2 = sp.symbols('x0 x1 x2')
MONS = [x0 ** i * x1 ** j * x2 ** (6 - i - j) for i in range(7) for j in range(7 - i)]


def coeff_vec(G):
    P = sp.Poly(sp.expand(G), x0, x1, x2)
    d = dict(P.terms())
    return [d.get(sp.Poly(mo, x0, x1, x2).monoms()[0], 0) for mo in MONS]


def rank(cols):
    M = fmpq_mat(len(cols[0]), len(cols), [fmpq(int(sp.fraction(sp.Rational(c))[0]), int(sp.fraction(sp.Rational(c))[1]))
                                          for row in zip(*cols) for c in row])
    return M.rank()


def eff_dim(F, params, vals):
    Fv = F.subs(vals)
    Jp = [coeff_vec(sp.diff(F, p).subs(vals)) for p in params]
    Jg = [coeff_vec(xi * sp.diff(Fv, xj)) for xi in (x0, x1, x2) for xj in (x0, x1, x2)]
    rg = rank(Jg)
    return rank(Jp + Jg) - rg, rg


def cover_check():
    x, v, s, xi, c, St, z = sp.symbols('x v s xi c S z')
    z8 = sp.exp(sp.I * sp.pi / 4)
    Yeq = x * (x ** 8 + s * x ** 4 + 1)
    D4 = xi ** 4 - 4 * xi ** 2 + 2
    print("== the genus-4 cover (all checks should be True)")
    print("   D4(x + 1/x) = x^4 + x^-4:", sp.simplify(D4.subs(xi, x + 1 / x) - x ** 4 - x ** -4) == 0)
    g = (sp.I * x, z8 * v)
    io = (1 / x, v / x ** 5)
    print("   g preserves Y:", sp.simplify(g[1] ** 2 - Yeq.subs(x, g[0]) - z8 ** 2 * (v ** 2 - Yeq)) == 0)
    print("   iota preserves Y:", sp.simplify(io[1] ** 2 - Yeq.subs(x, io[0]) - (v ** 2 - Yeq) / x ** 10) == 0)

    def comp(p, q):  # p after q
        return tuple(sp.simplify(e.subs({x: q[0], v: q[1]}, simultaneous=True)) for e in p)
    ginv = (-sp.I * x, v / z8)
    print("   iota g iota = g^-1:", all(sp.simplify(a - b) == 0 for a, b in zip(comp(io, comp(g, io)), ginv)))
    print("   g^4 = (x, -v), so g has order 8:", comp(g, comp(g, comp(g, g))) == (x, -v))
    xi_, eta_ = x + 1 / x, v * (1 + 1 / x) / x ** 2
    print("   (xi, eta) is iota-invariant:", sp.simplify(xi_.subs(x, io[0]) - xi_) == 0
          and sp.simplify(eta_.subs({x: io[0], v: io[1]}, simultaneous=True) - eta_) == 0)
    print("   eta^2 = (xi + 2)(D4(xi) + s) on Y:",
          sp.simplify(Yeq * (1 + 1 / x) ** 2 / x ** 4 - (xi_ + 2) * (D4.subs(xi, xi_) + s)) == 0)
    print("   fixed points of iota: x = 1/x and v = v/x^5 at x = 1 (both points), x = -1 only if v = 0 (s = -2)")
    for j in range(4):
        lam = sp.I ** (j + 1) / z8
        ok = sp.simplify(sp.expand_complex((g[0] ** j * sp.diff(g[0], x) / g[1] - lam * x ** j / v) * v / x ** j)) == 0
        print(f"   g^* (x^{j} dx/v) = z8^{2 * j + 1} x^{j} dx/v: {ok and sp.simplify(sp.expand_complex(lam - z8 ** (2 * j + 1))) == 0}; "
              f"(g + g^-1)^* eigenvalue {sp.nsimplify(sp.simplify(sp.expand_complex(lam + 1 / lam)))}")
    print("   iota^* (x^j dx/v) = - x^(3-j) dx/v:",
          all(sp.simplify(io[0] ** j * sp.diff(io[0], x) / io[1] + x ** (3 - j) / v) == 0 for j in range(4)))
    print("   so the iota-invariant differentials are w0 - w3 (eigenvalue sqrt2) and w1 - w2 (eigenvalue -sqrt2)")
    f = (z + 2 * c) * (z ** 4 - 4 * c ** 2 * z ** 2 + 2 * c ** 4 + St)
    print("   f(c xi) = c^5 (xi + 2)(D4(xi) + S/c^4):",
          sp.simplify(f.subs(z, c * xi) - c ** 5 * (xi + 2) * (D4 + St / c ** 4)) == 0)


def main():
    cover_check()
    rnd = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 11)
    c0, c1 = sp.symbols('c0 c1')
    s = sp.symbols('s0:6')
    ch = c0 * x1 + c1 * x0
    Sh = sum(s[i] * x0 ** i * x1 ** (5 - i) for i in range(6))
    F = (x2 + 2 * ch) * (x1 * (x2 ** 4 - 4 * ch ** 2 * x2 ** 2 + 2 * ch ** 4) + Sh)
    params = [c0, c1, *s]
    for trial in range(3):
        vals = {p: sp.Rational(rnd.randint(-20, 20), rnd.randint(1, 7)) for p in params}
        d, rg = eff_dim(F, params, vals)
        print(f"TTV family, trial {trial}: effective dimension {d} (8 parameters; orbit tangent rank {rg})")
    # all double planes branched along a line plus a quintic
    lc = sp.symbols('l0:3')
    qc = sp.symbols('q0:21')
    Q5 = [x0 ** i * x1 ** j * x2 ** (5 - i - j) for i in range(6) for j in range(6 - i)]
    G = (lc[0] * x0 + lc[1] * x1 + lc[2] * x2) * sum(q * mo for q, mo in zip(qc, Q5))
    params = [*lc, *qc]
    vals = {p: sp.Rational(rnd.randint(-20, 20), rnd.randint(1, 7)) for p in params}
    d, rg = eff_dim(G, params, vals)
    print(f"line + quintic: effective dimension {d} (expected 14 = 20 - 6)")


if __name__ == "__main__":
    main()
