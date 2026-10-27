"""78.7c-d: checks for the cover counts of NOTE.md Section 5.4 (draft DRAFT_78_7.md).
1. The D4 family C_a: y^2 = x(x^4 + a x^2 + 1) and its quotient by tau: (x, y) -> (1/x, y/x^3): with
   u = x + 1/x and v = y (x^3 + 1)/x^3, the point (u, v/(u - 1)) lies on E1: Y^2 = (u + 2)(u^2 + a - 2) (checked
   modulo the equation of C_a); tau fixes u and v; j(E1) as a rational function of a, and its degree.
2. Enumeration of ramification profiles, for each case, of the bound N <= 2e - 2 - sum over constrained points of
   sum (k_i - 1), with at most 6 - e unstable fibres:
   - a point with local monodromy of order m on C0/iota: the parts k with m not dividing k are unstable fibres;
   - U = 0: over rho = -1 the odd parts are unstable; over rho = 0 the parts are free; the number of odd parts over
     {rho = -1, rho = 0} is 2b with b in {1, 2} and 12 (3 - b) >= 3e (p_g of the elliptic surface over B' is 1,
     chi = 3 - b, deg j = 3e).
   Calibration: the TTV line (m = 4 at sigma = oo) must give N <= 5 as in Proposition 5.4.
Usage: python -I -B c18_d8_counts.py
"""
import sys

import sympy as sp

sys.stdout.reconfigure(encoding="utf-8")


def partitions(n, maxpart=None):
    if maxpart is None:
        maxpart = n
    if n == 0:
        yield ()
        return
    for k in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - k, k):
            yield (k,) + rest


def cost(p):
    return sum(k - 1 for k in p)


def check_quotient():
    x, y, a, u = sp.symbols('x y a u')
    f = x * (x ** 4 + a * x ** 2 + 1)
    uu = x + 1 / x
    v = y * (x ** 3 + 1) / x ** 3
    Y = v / (uu - 1)
    lhs = sp.together(Y ** 2 - (uu + 2) * (uu ** 2 + a - 2))
    num = sp.numer(lhs)
    rem = sp.rem(sp.Poly(sp.expand(num), y), sp.Poly(y ** 2 - f, y))
    ok_eq = sp.simplify(rem.as_expr()) == 0
    tau = {x: 1 / x, y: y / x ** 3}
    ok_u = sp.simplify(uu.subs(tau, simultaneous=True) - uu) == 0
    ok_v = sp.simplify(v.subs(tau, simultaneous=True) - v) == 0
    # j-invariant of Y^2 = (u + 2)(u^2 + a - 2)
    cub = sp.Poly(sp.expand((u + 2) * (u ** 2 + a - 2)), u)
    c3, c2, c1, c0 = cub.all_coeffs()
    # depressed: u = w - c2/3
    p_ = sp.simplify(c1 - c2 ** 2 / 3)
    q_ = sp.simplify(c0 - c1 * c2 / 3 + 2 * c2 ** 3 / 27)
    j = sp.factor(sp.simplify(1728 * 4 * p_ ** 3 / (4 * p_ ** 3 + 27 * q_ ** 2)))
    num, den = sp.fraction(j)
    deg = max(sp.Poly(num, a).degree(), sp.Poly(den, a).degree())
    print(f"1. E1 = C_a/tau: (u, v/(u-1)) on Y^2 = (u+2)(u^2+a-2) modulo C_a: {ok_eq}; u, v tau-invariant: {ok_u}, {ok_v}")
    print(f"   j(E1) = {j}; degree in a: {deg}")
    return deg


def bound_simple(e, points):
    """points: list of m (monodromy orders); returns the best N and the profiles attaining it."""
    best, arg = None, None

    def rec(i, chosen, nbad, c):
        nonlocal best, arg
        if nbad > 6 - e:
            return
        if i == len(points):
            n = 2 * e - 2 - c
            if best is None or n > best:
                best, arg = n, tuple(chosen)
            return
        for p in partitions(e):
            bad = sum(1 for k in p if k % points[i] != 0)
            rec(i + 1, chosen + [p], nbad + bad, c + cost(p))

    rec(0, [], 0, 0)
    return best, arg


def bound_u0(e, degj):
    best, arg = None, None
    for p1 in partitions(e):            # over rho = -1
        odd1 = sum(1 for k in p1 if k % 2)
        if odd1 > 6 - e:
            continue
        for p0 in partitions(e):        # over rho = 0
            odd0 = sum(1 for k in p0 if k % 2)
            tot = odd1 + odd0
            if tot % 2 or tot == 0:
                continue
            b = tot // 2
            chi = 3 - b
            if chi < 1 or 12 * chi < degj * e:
                continue
            n = 2 * e - 2 - cost(p1) - cost(p0)
            if best is None or n > best:
                best, arg = n, (p1, p0, b)
    return best, arg


def main():
    degj = check_quotient()
    cases = {
        "TTV (calibration; sigma = oo, m = 4)": [4],
        "L1: 16R = U (cusp, m = 3)": [3],
        "L2: 4R + 4S = 7U (cusp m = 3, (1:-1:0) m = 4)": [3, 4],
    }
    print("2. N <= 2e - 2 - sum(k - 1), at most 6 - e unstable fibres:")
    for name, pts in cases.items():
        res = {e: bound_simple(e, pts) for e in range(1, 7)}
        top = max(v[0] for v in res.values() if v[0] is not None)
        print(f"   {name}: " + "; ".join(f"e={e}: N<={v[0]} {v[1]}" for e, v in res.items() if v[0] is not None)
              + f"  => N <= {top}")
    res = {e: bound_u0(e, degj) for e in range(1, 7)}
    top = max(v[0] for v in res.values() if v[0] is not None)
    print("   U = 0 (D4 locus): " + "; ".join(f"e={e}: N<={v[0]} {v[1]}" if v[0] is not None else f"e={e}: none"
                                          for e, v in res.items()) + f"  => N <= {top}")


if __name__ == "__main__":
    main()
