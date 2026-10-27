"""76.5: the degree-15 invariant on [EK]'s Humbert planes. For D = 5, 8, 12, 13, 17: factor R30(E(v1, v2)) (R30 from
c08, E = [EK]'s Igusa-Clebsch map times N^w, so that it is polynomial) and compare its odd-multiplicity part with
[EK]'s branch function b of Y_-(D) -> H_D (typed from [EK]'s theorems; data). Prints the factorisation, b, and
whether (odd part of R30(E)) / b is a constant times boundary factors.
Usage: python -I -B c09_r30_branch.py
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


C8 = _load("c08_parity")
C6 = C8.C6
W = C8.W


def branch(D, a, b):
    """[EK]'s branch functions (theorems for D = 5, 8, 12, 13, 17), in [EK]'s plane coordinates (a, b)."""
    if D == 5:
        g, h = a, b
        return 2 * (6250 * h ** 2 - 4500 * g ** 2 * h - 1350 * g * h - 108 * h - 972 * g ** 5 - 324 * g ** 4 - 27 * g ** 3)
    if D == 8:
        r, s = a, b
        return 2 * (16 * r * s ** 2 + 32 * r ** 2 * s - 40 * r * s - s + 16 * r ** 3 + 24 * r ** 2 + 12 * r + 2)
    if D == 12:
        e, f = a, b
        return (f - 1) * (f + 1) * (f ** 6 - f ** 4 - 18 * e * f ** 2 + 27 * e ** 2 + 16 * e)
    if D == 13:
        g, h = a, b
        return 108 * g * h ** 3 - (27 * g ** 2 + 468 * g - 4) * h ** 2 + 8 * (82 * g ** 2 + 71 * g - 2) * h - 16 * (2 * g - 1) ** 3
    if D == 17:
        g, h = a, b
        return (-256 * h ** 3 + (192 * g ** 2 + 464 * g + 185) * h ** 2 - 2 * (2 * g + 1) * (12 * g ** 3 - 65 * g ** 2 - 54 * g - 9) * h
                + (g + 1) ** 4 * (2 * g + 1) ** 2)


def main():
    R, _, _ = C8.r30()
    for D in (5, 8, 12, 13, 17):
        (a, b), I = C6.load(D)
        N = sp.Integer(1)
        for k in W:
            for fac, e in sp.factor_list(sp.fraction(sp.cancel(I[k]))[1])[1]:
                need = -(-e // W[k])
                have = dict(sp.factor_list(N)[1]).get(fac, 0) if N != 1 else 0
                N *= fac ** max(0, need - have)
        P = {k: sp.expand(sp.cancel(N ** W[k] * I[k])) for k in W}
        val = sp.expand(C8.at(R, P))
        fl = sp.factor_list(val, a, b)
        odd = [fa for fa, mu in fl[1] if mu % 2]
        bb = sp.factor_list(sp.expand(branch(D, a, b)), a, b)
        bfac = [fa for fa, mu in bb[1] if mu % 2]
        print(f"== D = {D} ({a}, {b}); N = {N}")
        print(f"   R30(E): constant {fl[0]}; factors: " + "; ".join(f"[{fa}]^{mu}" for fa, mu in fl[1]))
        print(f"   odd part of R30(E): {odd}")
        print(f"   odd part of [EK]'s branch function: {bfac}")
        ex = [fa for fa in odd if not any(sp.simplify(fa / g_) .is_number for g_ in bfac)]
        miss = [g_ for g_ in bfac if not any(sp.simplify(fa / g_).is_number for fa in odd)]
        print(f"   odd factors of R30(E) not in the branch function: {ex}; branch factors not odd in R30(E): {miss}")


if __name__ == "__main__":
    main()
