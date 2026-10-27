"""76.2 sanity check (the local orders are standard; this only guards against slips): for a smooth sextic B and
P0 = (0:0:1) not on B, with the line {x0 = 0} of the pencil of lines through P0 bitangent / tritangent / flex-tangent
to B, the order at t = 0 of I10 (= disc) of f_t(z) = F(t, 1, z) is 2 / 3 / 2. B smooth is certified by the
vanishing of the degree-13 part of the Jacobian ring mod two primes (as in the [K3M] repository, b14_example.py).
Usage: python -I -B c04_local_orders.py
"""
import importlib.util
import itertools
import random
import sys
from pathlib import Path

import flint
import sympy as sp

sys.stdout.reconfigure(encoding="utf-8")
spec = importlib.util.spec_from_file_location("k3lib", Path(__file__).resolve().parent / "k3lib.py")
L = importlib.util.module_from_spec(spec)
spec.loader.exec_module(L)
x0, x1, x2 = xs = sp.symbols('x0 x1 x2')
t = sp.Symbol('t')


def mons(n):
    return [e for e in itertools.product(range(n + 1), repeat=3) if sum(e) == n]


def smooth_mod(F, p):
    """True iff the Jacobian ideal of F contains all forms of degree 13 mod p (=> F smooth mod p, hence over Q)."""
    Fi = sp.Poly(sp.expand(F), *xs)
    ctx = flint.nmod_mpoly_ctx.get(('x0', 'x1', 'x2'), modulus=p)
    g = ctx.from_dict({tuple(m): int(c) % p for m, c in Fi.as_dict().items()})
    m13, m8 = mons(13), mons(8)
    idx = {m: i for i, m in enumerate(m13)}
    X = ctx.gens()
    rows = []
    for dv in (g.derivative(0), g.derivative(1), g.derivative(2)):
        for m in m8:
            pr = dv * X[0] ** m[0] * X[1] ** m[1] * X[2] ** m[2]
            row = [0] * len(m13)
            for ex, cf in pr.to_dict().items():
                row[idx[tuple(int(u) for u in ex)]] = int(cf)
            rows.append(row)
    return flint.nmod_mat(len(rows), len(m13), [v for r_ in rows for v in r_], p).rank() == len(m13)


def rand_form(rng, d):
    return sum(rng.randint(-3, 3) * x0 ** a * x1 ** b * x2 ** c for a, b, c in mons(d))


def main():
    rng = random.Random(2024)
    L.ic_polys()
    cases = {
        "bitangent (2,2,1,1)": lambda: x0 * rand_form(rng, 5) + rand_form_line(rng, 2) ** 2 * rand_form_line(rng, 2),
        "tritangent (2,2,2)": lambda: x0 * rand_form(rng, 5) + rand_form_line(rng, 3) ** 2,
        "flex (3,1,1,1)": lambda: x0 * rand_form(rng, 5) + (x2 - x1) ** 3 * rand_form_line(rng, 3),
    }
    for name, make in cases.items():
        for _ in range(20):
            F = sp.expand(make())
            if sp.Poly(F, *xs).coeff_monomial(x2 ** 6) == 0 or not (smooth_mod(F, 1000003) and smooth_mod(F, 32003)):
                continue
            f = sp.expand(F.subs({x0: t, x1: 1}).subs(x2, L.X))
            I10 = sp.Poly(L.ic(L.poly_coeffs(f, L.X))[10], t)
            ordt = min(m[0] for m in I10.monoms()) if not I10.is_zero else None
            f0 = sp.factor(f.subs(t, 0))
            print(f"  {name}: B smooth, P0 not on B; f_0 = {f0}; ord_t=0 I10 = {ordt}")
            break


def rand_form_line(rng, d):
    """binary form in x1, x2 of degree d (its restriction to x0 = 0 is itself)."""
    return sum(rng.randint(-3, 3) * x1 ** (d - k) * x2 ** k for k in range(d + 1)) + x2 ** d


if __name__ == "__main__":
    main()
