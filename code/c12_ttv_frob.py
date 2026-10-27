"""76.6: Frobenius on the TTV K3 surfaces of NOTE.md Check 3.3 and c11, to bound their Picard number.
S: w^2 = F(x0, x1, x2) = (x2 + 2c)(x1 (x2^4 - 4c^2 x2^2 + 2c^4) + S(x0, x1)), c = c0 x1 + c1 x0, S a quintic form;
B = line + quintic G with 5 nodes, P0 = (0:0:1) on G. H^2 of the resolution contains the 6 classes H, E_1..E_5;
on the 16-dimensional rest, tr(Frob_q) = sum over P in P^2(F_q) of chi(F(P)) = sum over (t, z) in F_q^2 of
chi(f_t(z)), f_t(z) = F(t, 1, z) (the line x1 = 0 contributes 0). With z = c x and sigma = S/c^4:
chi(f_t(z)) = chi(c) chi((x + 2)(D(x) + sigma)), D(x) = x^4 - 4x^2 + 2, so the sum is
sum_t chi(c(t)) a(sigma(t)) with a(sigma) = sum_x chi(x + 2) chi(D(x) + sigma), a correlation over (F_q, +)
computed by an FFT over (Z/p)^k. Brute force for small q as a check.
Output: the characteristic polynomial P(T) = det(1 - Frob T) on the 16-dimensional part, completed by the
functional equation, the absolute values of its roots / q, and its cyclotomic factors. rho(S over F_p-bar) <=
6 + (degree of the cyclotomic part); rho over Q-bar <= rho over F_p-bar. The sign eps of the functional equation
is kept only if it is consistent with all computed traces (for kmax = 8 this means c_8 = eps c_8). Also the
factorisation of the normalised polynomial over Q(sqrt2).
Usage: python -I -B c12_ttv_frob.py p kmax [c0 c1 s0 ... s5]
"""
import sys
from fractions import Fraction

import numpy as np
import sympy as sp

sys.stdout.reconfigure(encoding="utf-8")


def poly_mulmod(a, b, m, p):
    """a, b, m: coefficient lists (low to high) over F_p, m monic of degree k."""
    k = len(m) - 1
    res = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                res[i + j] = (res[i + j] + x * y) % p
    for d in range(len(res) - 1, k - 1, -1):
        c = res[d]
        if c:
            for j in range(k + 1):
                res[d - k + j] = (res[d - k + j] - c * m[j]) % p
    return (res + [0] * k)[:k]


def poly_powmod(a, e, m, p):
    r = [1] + [0] * (len(m) - 2)
    while e:
        if e & 1:
            r = poly_mulmod(r, a, m, p)
        a = poly_mulmod(a, a, m, p)
        e >>= 1
    return r


def primitive_poly(p, k):
    q = p ** k
    fac = sp.factorint(q - 1)
    one = [1] + [0] * (k - 1)
    X = ([0, 1] + [0] * k)[:k] if k > 1 else None
    import itertools
    for tail in itertools.product(range(p), repeat=k):
        m = list(tail) + [1]
        if m[0] == 0:
            continue
        if k == 1:
            g = (-m[0]) % p          # X = -m0 in F_p; need a generator of F_p^*
            if all(pow(g, (q - 1) // r, p) != 1 for r in fac):
                return m
            continue
        if poly_powmod(X, q - 1, m, p) != one:
            continue
        if all(poly_powmod(X, (q - 1) // r, m, p) != one for r in fac):
            return m
    raise ValueError("no primitive polynomial")


class Field:
    def __init__(self, p, k):
        self.p, self.k, self.q = p, k, p ** k
        q = self.q
        m = primitive_poly(p, k)
        self.pw = [p ** j for j in range(k)]
        exp = np.empty(q - 1, dtype=np.int64)
        if k == 1:
            g = (-m[0]) % p
            v = 1
            for i in range(q - 1):
                exp[i] = v
                v = v * g % p
        else:
            red = [[(-t * m[j]) % p for j in range(k)] for t in range(p)]
            cur = [1] + [0] * (k - 1)
            pw = self.pw
            for i in range(q - 1):
                exp[i] = sum(d * w for d, w in zip(cur, pw))
                top = cur[-1]
                cur = [0] + cur[:-1]
                if top:
                    r = red[top]
                    cur = [(a + b) % p for a, b in zip(cur, r)]
        log = np.full(q, -1, dtype=np.int64)
        log[exp] = np.arange(q - 1)
        assert (log[1:] >= 0).all()
        self.exp, self.log = exp, log

    def const(self, c):
        return np.int64(c % self.p)

    def add(self, a, b):
        a = np.asarray(a, dtype=np.int64)
        b = np.asarray(b, dtype=np.int64)
        out = np.zeros(np.broadcast(a, b).shape, dtype=np.int64)
        for w in self.pw:
            out += ((a // w + b // w) % self.p) * w
        return out

    def neg(self, a):
        a = np.asarray(a, dtype=np.int64)
        out = np.zeros(a.shape, dtype=np.int64)
        for w in self.pw:
            out += ((-(a // w)) % self.p) * w
        return out

    def mul(self, a, b):
        a, b = np.broadcast_arrays(np.asarray(a, dtype=np.int64), np.asarray(b, dtype=np.int64))
        out = np.zeros(a.shape, dtype=np.int64)
        nz = (a != 0) & (b != 0)
        out[nz] = self.exp[(self.log[a[nz]] + self.log[b[nz]]) % (self.q - 1)]
        return out

    def power(self, a, e):
        a = np.asarray(a, dtype=np.int64)
        out = np.zeros(a.shape, dtype=np.int64)
        nz = a != 0
        out[nz] = self.exp[(self.log[a[nz]] * e) % (self.q - 1)]
        if e == 0:
            out[:] = 1
        return out

    def chi(self, a):
        a = np.asarray(a, dtype=np.int64)
        out = np.zeros(a.shape, dtype=np.int64)
        nz = a != 0
        out[nz] = 1 - 2 * (self.log[a[nz]] % 2)
        return out


def poly_eval(Fq, coeffs, x):
    """Horner, coeffs low to high (integers)."""
    acc = np.full(np.shape(x), Fq.const(coeffs[-1]), dtype=np.int64)
    for c in reversed(coeffs[:-1]):
        acc = Fq.add(Fq.mul(acc, x), Fq.const(c))
    return acc


def trace_fft(Fq, c0, c1, S):
    q, p, k = Fq.q, Fq.p, Fq.k
    xs = np.arange(q, dtype=np.int64)
    x2 = Fq.mul(xs, xs)
    D = Fq.add(Fq.add(Fq.mul(x2, x2), Fq.mul(Fq.const(-4), x2)), Fq.const(2))
    h = np.bincount(D, weights=Fq.chi(Fq.add(xs, Fq.const(2))).astype(np.float64), minlength=q)
    ht = h[Fq.neg(xs)]
    g = Fq.chi(xs).astype(np.float64)
    shape = (p,) * k
    a = np.fft.ifftn(np.fft.fftn(ht.reshape(shape)) * np.fft.fftn(g.reshape(shape))).real.reshape(q)
    ai = np.rint(a).astype(np.int64)
    assert np.abs(a - ai).max() < 1e-6
    del a, h, ht, g, D, x2
    ct = Fq.add(Fq.const(c0), Fq.mul(Fq.const(c1), xs))
    St = poly_eval(Fq, S, xs)
    nz = ct != 0
    inv4 = Fq.power(ct[nz], -4 % (q - 1))
    sig = Fq.mul(St[nz], inv4)
    tot = int((Fq.chi(ct[nz]) * ai[sig]).sum())
    # the fibre with c(t) = 0: f = z (z^4 + S(t))
    for tz in xs[~nz]:
        Sz = St[tz]
        z = xs
        val = Fq.mul(z, Fq.add(Fq.power(z, 4), Sz))
        tot += int(Fq.chi(val).sum())
    return tot


def trace_brute(Fq, c0, c1, S):
    q = Fq.q
    xs = np.arange(q, dtype=np.int64)
    tot = 0
    for t in range(q):
        c = int(Fq.add(Fq.const(c0), Fq.mul(Fq.const(c1), np.int64(t))))
        St = int(poly_eval(Fq, S, np.int64(t)))
        z = xs
        lin = Fq.add(z, Fq.mul(Fq.const(2), np.int64(c)))
        c2 = Fq.mul(np.int64(c), np.int64(c))
        quart = Fq.add(Fq.add(Fq.power(z, 4), Fq.mul(Fq.mul(Fq.const(-4), c2), Fq.mul(z, z))),
                       Fq.add(Fq.mul(Fq.const(2), Fq.mul(c2, c2)), np.int64(St)))
        tot += int(Fq.chi(Fq.mul(lin, quart)).sum())
    return tot


def good_reduction(p, c0, c1, S):
    x0, x1, x2 = sp.symbols('x0 x1 x2')
    ch = c0 * x1 + c1 * x0
    Sh = sum(S[i] * x0 ** i * x1 ** (5 - i) for i in range(6))
    G = sp.expand(x1 * (x2 ** 4 - 4 * ch ** 2 * x2 ** 2 + 2 * ch ** 4) + Sh)
    eqs = [G] + [sp.diff(G, v) for v in (x0, x1, x2)]
    for chart in (x0, x1, x2):
        E = [sp.expand(e.subs(chart, 1)) for e in eqs]
        vs = [v for v in (x0, x1, x2) if v != chart]
        Gb = sp.groebner(E, *vs, modulus=p, order='grevlex')
        if list(Gb.exprs) != [1]:
            return False, f"quintic singular mod {p} (chart {chart} = 1)"
    g = sp.Poly(sp.expand(G.subs(x2, -2 * ch)), x0, x1)
    gt = sp.Poly(sp.expand(g.as_expr().subs(x1, 1)), x0)
    if gt.degree() < 5 and sp.Poly(sp.expand(g.as_expr().subs(x0, 1)), x1).degree() < 5:
        pass
    d = sp.discriminant(sp.Poly(g.as_expr().subs(x1, 1), x0)) if gt.degree() == 5 else None
    if gt.degree() != 5 or d % p == 0 or sp.Poly(gt, modulus=p).degree() != 5:
        return False, f"line meets the quintic badly mod {p} (degree {gt.degree()}, disc {d})"
    return True, "ok"


def charpoly(traces, q, deg=16):
    """traces[k-1] = tr Frob^k, k = 1..K. Power series of exp(-sum T_k x^k / k), then the functional equation
    c_{deg-i} = eps q^{deg-2i} c_i."""
    K = len(traces)
    c = [Fraction(1)] + [Fraction(0)] * K
    # Newton: k c_k = -sum_{i=1..k} T_i c_{k-i}
    for kk in range(1, K + 1):
        c[kk] = -sum(Fraction(traces[i - 1]) * c[kk - i] for i in range(1, kk + 1)) / kk
    out = {}
    for eps in (1, -1):
        full = [None] * (deg + 1)
        ok = True
        for i in range(deg + 1):
            j = deg - i
            if i <= K:
                full[i] = c[i]
            if j <= K:
                val = eps * Fraction(q) ** (deg - 2 * j) * c[j]
                if full[i] is not None and full[i] != val:
                    ok = False
                full[i] = val
        if ok and all(v is not None for v in full):
            out[eps] = full
    return out


def analyse(full, q):
    T = sp.Symbol('T')
    # R(T) = sum_j c_j q^{-j} T^{deg-j}: roots alpha_i / q
    deg = len(full) - 1
    R = sum(sp.Rational(full[j].numerator, full[j].denominator) / sp.Integer(q) ** j * T ** (deg - j)
            for j in range(deg + 1))
    R = sp.Poly(R, T)
    roots = np.roots([float(cf) for cf in R.all_coeffs()])
    absr = sorted(abs(r) for r in roots)
    fl = sp.factor_list(R.as_expr())
    cyc = 0
    facs = []
    for f, e in fl[1]:
        fp = sp.Poly(f, T)
        is_cyc = False
        mono = [sp.Rational(cf) for cf in fp.monic().all_coeffs()]
        for n in range(1, 200):
            if sp.totient(n) == fp.degree() and \
                    [sp.Rational(cf) for cf in sp.Poly(sp.cyclotomic_poly(n, T), T).all_coeffs()] == mono:
                is_cyc = n
                break
        if is_cyc:
            cyc += fp.degree() * e
        facs.append((fp.degree(), e, is_cyc))
    # over Q(sqrt2): with RM by Q(sqrt2) defined over F_p, R splits into two conjugate octics when p = +-1 mod 8
    facs2 = [(sp.degree(f, T), e) for f, e in sp.factor_list(R.as_expr(), T, extension=sp.sqrt(2))[1]]
    return absr, facs, cyc, R, facs2


def main():
    p = int(sys.argv[1])
    kmax = int(sys.argv[2])
    if len(sys.argv) > 3:
        c0, c1, *S = [int(v) for v in sys.argv[3:11]]
    else:
        c0, c1, S = 1, 2, [3, 1, -2, 1, 5, -1]
    print(f"TTV member c = {c0} + {c1} t, S = {S}, p = {p}")
    ok, msg = good_reduction(p, c0, c1, S)
    print(f"good reduction at {p}: {ok} ({msg})")
    if not ok:
        return
    traces = []
    for k in range(1, kmax + 1):
        Fq = Field(p, k)
        tr = trace_fft(Fq, c0, c1, S)
        if Fq.q <= 125:
            tb = trace_brute(Fq, c0, c1, S)
            assert tb == tr, (k, tb, tr)
            print(f"   q = {Fq.q}: brute force agrees ({tb})")
        traces.append(tr)
        print(f"   k = {k}, q = {Fq.q}: trace on the 16-dimensional part = {tr}", flush=True)
    cps = charpoly(traces, p)
    for eps, full in cps.items():
        absr, facs, cyc, R, facs2 = analyse(full, p)
        print(f"== eps = {eps}: |roots/q| in [{absr[0]:.6f}, {absr[-1]:.6f}]")
        print(f"   factors of R (degree, multiplicity, cyclotomic index or False): {facs}")
        print(f"   factors of R over Q(sqrt2) (degree, multiplicity): {facs2}")
        print(f"   cyclotomic degree {cyc}: rho(S over F_{p}-bar) <= 6 + {cyc} = {6 + cyc}")
        print(f"   R(T) = {R.as_expr()}")


if __name__ == "__main__":
    main()
