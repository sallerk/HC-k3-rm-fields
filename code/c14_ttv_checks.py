"""NOTE.md, Theorem C, proof steps 5 and 7 (after referee R2, points M3 and M9).
Step 5 (upper bound for the moduli): g = (x0, x1, x2) -> (a x0 + b x1, d x1, lam x2) gives F o g = d lam^5 F'' exactly,
where F'' has the shape of F with c'' = c o g / lam and P'' = P o g / (d lam^4); scalars act trivially, and the
orbits of this group in the 8-dimensional (c, P)-space have dimension 3 (rank of the orbit map at the identity, exact,
at random rational points; by lower semicontinuity rank 3 at one point gives orbit dimension 3 on a dense open set).
Step 7 (the sign eps): from the Frobenius traces for k <= 8 in the c12 outputs, Newton's identities give c_1, ..., c_8
of det(T - Frob | V); the functional equation at i = 8 reads c_8 = eps c_8, so c_8 != 0 forces eps = +1. The traces
for k >= 9 in the same output are then compared with the resulting polynomial (a consistency check), and the
normalised R(T) is tested for evenness and factored over Q and Q(sqrt 2).
Usage: python -I -B c14_ttv_checks.py [seed]
"""
import random
import re
import sys
from fractions import Fraction
from pathlib import Path

import sympy as sp

sys.stdout.reconfigure(encoding="utf-8")
OUT = Path(__file__).resolve().parent.parent / "outputs"
x0, x1, x2, a, b, d, lam, T = sp.symbols('x0 x1 x2 a b d lam T')
c0, c1 = sp.symbols('c0 c1')
s = sp.symbols('s0:6')


def F_of(cl, P):
    return sp.expand((x2 + 2 * cl) * (x1 * (x2 ** 4 - 4 * cl ** 2 * x2 ** 2 + 2 * cl ** 4) + P))


def step5(seed):
    cl = c0 * x1 + c1 * x0
    P = sum(s[i] * x0 ** i * x1 ** (5 - i) for i in range(6))
    g = {x0: a * x0 + b * x1, x1: d * x1, x2: lam * x2}
    cpp = sp.expand(cl.subs(g, simultaneous=True) / lam)
    Ppp = sp.expand(P.subs(g, simultaneous=True) / (d * lam ** 4))
    Fg = sp.expand(F_of(cl, P).subs(g, simultaneous=True))
    print("F o g - d lam^5 F'' == 0:", sp.expand(Fg - d * lam ** 5 * F_of(cpp, Ppp)) == 0)
    new = [sp.Poly(cpp, x0, x1).coeff_monomial(x1), sp.Poly(cpp, x0, x1).coeff_monomial(x0)]
    new += [sp.Poly(Ppp, x0, x1).coeff_monomial(x0 ** i * x1 ** (5 - i)) for i in range(6)]
    mu = sp.symbols('mu')
    scal = {a: mu, b: 0, d: mu, lam: mu}
    print("scalars act trivially:", all(sp.simplify(e.subs(scal) - p) == 0 for e, p in zip(new, [c0, c1, *s])))
    J = sp.Matrix([[sp.diff(e, v) for v in (a, b, d, lam)] for e in new]).subs({a: 1, b: 0, d: 1, lam: 1})
    rnd = random.Random(seed)
    for trial in range(3):
        vals = {p: sp.Rational(rnd.randint(-20, 20), rnd.randint(1, 7)) for p in (c0, c1, *s)}
        print(f"   trial {trial}: rank of the orbit map at the identity = {J.subs(vals).rank()}"
              " (4 group parameters, scalars trivial)")


def traces_from(fname):
    tr = {}
    for line in (OUT / fname).read_text(encoding="utf-8").splitlines():
        m = re.search(r"k = (\d+), q = (\d+): trace on the 16-dimensional part = (-?\d+)", line)
        if m:
            tr[int(m.group(1))] = (int(m.group(2)), int(m.group(3)))
    return tr


def step7(fname, p):
    tr = traces_from(fname)
    q = p
    c = [Fraction(1)]                       # det(T - Frob | V) = sum c_i T^(16 - i)
    for k in range(1, 9):
        c.append(-sum(Fraction(tr[i][1]) * c[k - i] for i in range(1, k + 1)) / k)
    print(f"{fname}: c_8 = {c[8]}  (c_8 = eps c_8, so eps = +1: {c[8] != 0})")
    full = c + [Fraction(0)] * 8
    for i in range(9, 17):                  # c_i = eps q^(2i - 16) c_(16 - i), eps = +1
        full[i] = Fraction(q) ** (2 * i - 16) * c[16 - i]
    # power sums of the roots from the full polynomial (Newton), compared with the remaining traces
    ps = {}
    for k in range(1, 13):
        ps[k] = -k * (full[k] if k <= 16 else 0) - sum(full[i] * ps[k - i] for i in range(1, k))
    extra = [k for k in sorted(tr) if k > 8]
    print("   traces for k > 8 agree with the polynomial:", all(ps[k] == tr[k][1] for k in extra), "k =", extra)
    R = sp.expand(sum(sp.Rational(full[i].numerator, full[i].denominator) * (q * T) ** (16 - i)
                      for i in range(17)) / q ** 16)
    print("   R(T) =", R)
    print("   R even:", sp.expand(R - R.subs(T, -T)) == 0)
    print("   factors over Q:", [sp.degree(f, T) for f, _ in sp.factor_list(R)[1]],
          "; over Q(sqrt 2):", [sp.degree(f, T) for f, _ in sp.factor_list(R, extension=sp.sqrt(2))[1]])


if __name__ == "__main__":
    step5(int(sys.argv[1]) if len(sys.argv) > 1 else 5)
    step7("c12_ttv_frob_p3.txt", 3)
    step7("c12_ttv_frob_m2_p7.txt", 7)
