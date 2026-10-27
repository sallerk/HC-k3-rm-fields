"""78.7b: checks of the local models of M2 at curves with extra automorphisms, used for the local monodromy of the
moduli curves over L1, L2 and the TTV line (NOTE.md, Section 5.4).
Model: near [C*], M2 = Def(C*)/Gbar with Gbar = Aut(C*)/<iota>, and I15 vanishes to order one on each fixed
hyperplane of an involution of Gbar, so R30 = c I15^2 is, up to a unit, the square of the product of the linear forms
of these hyperplanes:
- V4 (Gbar = C2): one hyperplane, R30 ~ l^2: order 2 along a generic smooth family through C*;
- D6 (Gbar = S3 acting on Def = 1 + std): three planes, R30 ~ (l1 l2 l3)^2: order 6;
- |Aut| = 48 (Gbar = S4 acting on Def = the reflection representation of A3): six planes, R30 ~ Vandermonde^2:
  order 12.
A smooth family f* + u g (g a random integral sextic) has a holomorphic lift to Def, so ord_u R30 is the order of
the model function on a generic line through 0. The script computes ord_{u=0} R30(I(f* + u g)) exactly for several
random g and each C*.
Usage: python -I -B c16_local_models.py
"""
import importlib.util
import random
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
u = sp.Symbol('u')
R30, _, _ = C8.r30()

CURVES = {
    "V4: x^6 + 2 x^4 + 3 x^2 + 1 (x -> -x)": [1, 0, 3, 0, 2, 0, 1],
    "D6: x^6 + x^3 + 1/2 (the cusp of B8)": [sp.Rational(1, 2), 0, 0, 1, 0, 0, 1],
    "|Aut| = 48: x^5 - x": [0, -1, 0, 0, 0, 1, 0],
}


def order_r30(base, g):
    co = [sp.sympify(b) + u * c for b, c in zip(base, g)]
    I = L.ic(co)
    val = sp.Poly(sp.expand(C8.at(R30, I)), u)
    if val.is_zero:
        return None
    return min(m[0] for m in val.monoms())


def main():
    rnd = random.Random(787)
    print("ord_{u=0} R30(I(f* + u g)) for random integral sextics g (model prediction: V4 2, D6 6, order-48 12):")
    for name, base in CURVES.items():
        I0 = L.ic(base)
        assert I0[10] != 0
        ords = []
        for _ in range(4):
            g = [rnd.randint(-5, 5) for _ in range(7)]
            ords.append(order_r30(base, g))
        print(f"   {name}: {ords}")


if __name__ == "__main__":
    main()
