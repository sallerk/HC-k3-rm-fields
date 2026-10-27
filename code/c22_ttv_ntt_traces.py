"""86 (after referee report R4, item 11): the traces of Step 7 of the note's Theorem 4.1 and of its
Remark 4.2 by an exact number-theoretic transform, in place of the floating-point FFT of c12.
For q = p^k the sum over (t, z) in F_q^2 of chi(f_t(z)) is the sum over the t with c(t) != 0 of
chi(c(t)) A(sigma(t)), plus the fibres with c(t) = 0, where A(s) = sum_x chi(x + 2) chi(D(x) + s),
D(x) = x^4 - 4x^2 + 2, is a correlation on the additive group (Z/p)^k of F_q (as in c12, trace_fft). Here A is
computed by the k-dimensional discrete Fourier transform of size p in each axis over the prime field Z/P,
P = 2^31 - 1, which contains the p-th roots of unity for p = 3, 7 (p | P - 1). This is exact: A(s) is an integer
with |A(s)| <= q < P/2, so it is the centred residue modulo P.
Members as in c20: member 1 (c = 1 + 2t, P = 3 + t - 2t^2 + t^3 + 5t^4 - t^5) at p = 3, k <= 12, and member 2
(c = 5 + 3t, P = 5 - t + 6t^2 + 5t^3 + 5t^4 + 4t^5) at p = 7, k <= 8. Checks: c12's trace_brute (exact counting)
for q <= 3^8 and q <= 7^4, the whole function A against a direct sum for q <= 729, and c12's FFT values. Then the
characteristic polynomial from the exact traces alone (c12.charpoly, c12.analyse).
Usage: python -I -B c22_ttv_ntt_traces.py
"""
import importlib.util
import sys
import time
from pathlib import Path

import numpy as np
import sympy as sp

sys.stdout.reconfigure(encoding="utf-8")
spec = importlib.util.spec_from_file_location("c12", Path(__file__).resolve().parent / "c12_ttv_frob.py")
C12 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(C12)

P = 2 ** 31 - 1

# name, p, kmax, largest q for exact counting, c0, c1, coefficients of P(t)
MEMBERS = [
    ("member 1", 3, 12, 3 ** 8, 1, 2, [3, 1, -2, 1, 5, -1]),
    ("member 2", 7, 8, 7 ** 4, 5, 3, [5, -1, 6, 5, 5, 4]),
]


def root_of_unity(p):
    """A primitive p-th root of unity modulo P (p prime, p | P - 1)."""
    assert (P - 1) % p == 0
    for h in range(2, P):
        w = pow(h, (P - 1) // p, P)
        if w != 1:
            return w
    raise ValueError("no root of unity")


def dft(arr, w, p, k):
    """The DFT of size p in each of the k axes over Z/P; arr has shape (p,)*k and entries in [0, P).
    Each product is below 2^62 and is reduced modulo P before it is added, so int64 does not overflow."""
    W = [[pow(w, i * j, P) for j in range(p)] for i in range(p)]
    for axis in range(k):
        a = np.moveaxis(arr, axis, 0)
        out = np.empty_like(a)
        for i in range(p):
            acc = np.zeros(a.shape[1:], dtype=np.int64)
            for j in range(p):
                acc = (acc + W[i][j] * a[j] % P) % P
            out[i] = acc
        arr = np.moveaxis(out, 0, axis)
    return arr


def correlation_ntt(Fq, w, winv):
    """A(s) = sum_x chi(x + 2) chi(D(x) + s) for all s in F_q, exactly."""
    q, p, k = Fq.q, Fq.p, Fq.k
    xs = np.arange(q, dtype=np.int64)
    x2 = Fq.mul(xs, xs)
    D = Fq.add(Fq.add(Fq.mul(x2, x2), Fq.mul(Fq.const(-4), x2)), Fq.const(2))
    chi2 = Fq.chi(Fq.add(xs, Fq.const(2)))
    h = (np.bincount(D[chi2 == 1], minlength=q) - np.bincount(D[chi2 == -1], minlength=q)).astype(np.int64)
    ht = h[Fq.neg(xs)]
    g = Fq.chi(xs)
    del x2, D, chi2, h
    shape = (p,) * k
    Ht = dft((ht % P).reshape(shape), w, p, k)
    G = dft((g % P).reshape(shape), w, p, k)
    prod = Ht * G % P
    del Ht, G, ht, g
    a = dft(prod, winv, p, k).reshape(q) * pow(q, -1, P) % P
    a = np.where(a > P // 2, a - P, a)
    assert np.abs(a).max() <= q
    return a


def check_correlation(Fq, a):
    """The whole function A against the direct sums (small q only)."""
    xs = np.arange(Fq.q, dtype=np.int64)
    x2 = Fq.mul(xs, xs)
    D = Fq.add(Fq.add(Fq.mul(x2, x2), Fq.mul(Fq.const(-4), x2)), Fq.const(2))
    chi2 = Fq.chi(Fq.add(xs, Fq.const(2)))
    for s in range(Fq.q):
        if int((chi2 * Fq.chi(Fq.add(D, np.int64(s)))).sum()) != a[s]:
            return False
    return True


def trace_ntt(Fq, c0, c1, S, a):
    """As c12.trace_fft, with the exact correlation a."""
    q = Fq.q
    xs = np.arange(q, dtype=np.int64)
    ct = Fq.add(Fq.const(c0), Fq.mul(Fq.const(c1), xs))
    St = C12.poly_eval(Fq, S, xs)
    nz = ct != 0
    inv4 = Fq.power(ct[nz], -4 % (q - 1))
    sig = Fq.mul(St[nz], inv4)
    tot = int((Fq.chi(ct[nz]) * a[sig]).sum())
    # the fibre with c(t) = 0: f = z (z^4 + P(t))
    for tz in xs[~nz]:
        val = Fq.mul(xs, Fq.add(Fq.power(xs, 4), St[tz]))
        tot += int(Fq.chi(val).sum())
    return tot


def main():
    print(f"P = 2^31 - 1 = {P}, prime: {sp.isprime(P)}; P mod 3 = {P % 3}, P mod 7 = {P % 7}")
    assert sp.isprime(P)
    for name, p, kmax, qbrute, c0, c1, S in MEMBERS:
        w = root_of_unity(p)
        winv = pow(w, -1, P)
        print(f"== {name}: c = {c0} + {c1} t, P coefficients {S}, p = {p}; "
              f"primitive {p}-th root of unity modulo P: {w}")
        ok, msg = C12.good_reduction(p, c0, c1, S)
        print(f"   good reduction at {p}: {ok} ({msg})")
        traces = []
        for k in range(1, kmax + 1):
            Fq = C12.Field(p, k)
            t0 = time.time()
            a = correlation_ntt(Fq, w, winv)
            tn = trace_ntt(Fq, c0, c1, S, a)
            t1 = time.time()
            line = f"   k = {k}, q = {Fq.q}: NTT sum = {tn} ({t1 - t0:.1f} s)"
            if Fq.q <= 729:
                okA = check_correlation(Fq, a)
                assert okA
                line += f"; A equals the direct sums: {okA}"
            if Fq.q <= qbrute:
                tb = C12.trace_brute(Fq, c0, c1, S)
                assert tb == tn, (k, tb, tn)
                line += f"; exact count = {tb}, equal: {tb == tn}"
            tf = C12.trace_fft(Fq, c0, c1, S)
            line += f"; FFT value = {tf}, equal: {tf == tn}"
            print(line, flush=True)
            del a
            traces.append(tn)
        print(f"== {name}: the characteristic polynomial from the NTT traces k = 1..{kmax} only")
        cps = C12.charpoly(traces, p)
        print(f"   signs eps consistent with all {kmax} traces: {sorted(cps)}")
        for eps, full in cps.items():
            absr, facs, cyc, R, facs2 = C12.analyse(full, p)
            print(f"   eps = {eps}: c_1..c_8 = {[str(v) for v in full[1:9]]}")
            print(f"   factors of R (degree, multiplicity, cyclotomic index or False): {facs}")
            print(f"   factors of R over Q(sqrt2) (degree, multiplicity): {facs2}")
            print(f"   cyclotomic degree {cyc}: rho(S over F_{p}-bar) <= 6 + {cyc} = {6 + cyc}")
            print(f"   R(T) = {R.as_expr()}")


if __name__ == "__main__":
    main()
