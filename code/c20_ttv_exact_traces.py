"""78.8 (after the blind referee R1, point M6): the traces of Step 7 of the note's Theorem 4.1 by exact counting.
For the first member (c = 1 + 2t, P = 3 + t - 2t^2 + t^3 + 5t^4 - t^5, p = 3), the sums
    sum over (t, z) in F_q^2 of chi(f_t(z)),  q = 3^k, k = 1..8,
are computed with c12's trace_brute (exact integer arithmetic in F_q, no floating point) and compared with c12's
FFT values; the characteristic polynomial and the bound rho <= 6 are then derived from the exact traces alone.
As a further check of the FFT, the second member (c = 5 + 3t, P = 5 - t + 6t^2 + 5t^3 + 5t^4 + 4t^5, p = 7) is
compared for k = 1..4 (q <= 2401).
Usage: python -I -B c20_ttv_exact_traces.py
"""
import importlib.util
import sys
import time
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")
spec = importlib.util.spec_from_file_location("c12", Path(__file__).resolve().parent / "c12_ttv_frob.py")
C12 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(C12)

MEMBERS = [
    ("member 1", 3, 8, 1, 2, [3, 1, -2, 1, 5, -1]),
    ("member 2", 7, 4, 5, 3, [5, -1, 6, 5, 5, 4]),
]


def main():
    exact_member1 = None
    for name, p, kmax, c0, c1, S in MEMBERS:
        print(f"== {name}: c = {c0} + {c1} t, P coefficients {S}, p = {p}")
        ok, msg = C12.good_reduction(p, c0, c1, S)
        print(f"   good reduction at {p}: {ok} ({msg})")
        exact = []
        for k in range(1, kmax + 1):
            Fq = C12.Field(p, k)
            t0 = time.time()
            tb = C12.trace_brute(Fq, c0, c1, S)
            t1 = time.time()
            tf = C12.trace_fft(Fq, c0, c1, S)
            print(f"   k = {k}, q = {Fq.q}: exact sum = {tb}, FFT value = {tf}, equal: {tb == tf}"
                  f" (exact count {t1 - t0:.1f} s)", flush=True)
            exact.append(tb)
        if name == "member 1":
            exact_member1 = exact
    print("== member 1: the characteristic polynomial from the exact traces k = 1..8 only")
    cps = C12.charpoly(exact_member1, 3)
    for eps, full in cps.items():
        absr, facs, cyc, R, facs2 = C12.analyse(full, 3)
        print(f"   eps = {eps}: |roots/q| in [{absr[0]:.6f}, {absr[-1]:.6f}]")
        print(f"   factors of R (degree, multiplicity, cyclotomic index or False): {facs}")
        print(f"   factors of R over Q(sqrt2) (degree, multiplicity): {facs2}")
        print(f"   cyclotomic degree {cyc}: rho(S over F_3-bar) <= 6 + {cyc} = {6 + cyc}")
        print(f"   R(T) = {R.as_expr()}")


if __name__ == "__main__":
    main()
