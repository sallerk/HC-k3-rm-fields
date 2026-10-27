"""Shared helpers for the k3-rm-fields computations (loaded by the scripts with importlib from this folder).

Igusa-Clebsch invariants I2, I4, I6, I10 of a binary sextic sum a_i z^i, as isobaric polynomials in a0..a6,
obtained by exact interpolation over Q from the root formulas (Mestre's sums over 15, 10, 60 terms, times
lc^{2,4,6,10}), as in the code of the [K3M] repository (b16_ic_interp_q.py); cached in ../outputs/ic_polys.txt.
Normalisation: I2 = 6 a3^2 - 16 a2 a4 + 40 a1 a5 - 240 a0 a6, I10 = disc (checked when the cache is built).
"""
import itertools
import random
import re
from pathlib import Path

import sympy as sp
from flint import fmpq_mat

HERE = Path(__file__).resolve().parent
CACHE = HERE.parent / "outputs" / "ic_polys.txt"
R6 = range(6)
X = sp.Symbol('X')
A = sp.symbols('a0:7')


def _matchings(s):
    if not s:
        yield []
        return
    for j in range(1, len(s)):
        for m in _matchings(s[1:j] + s[j + 1:]):
            yield [(s[0], s[j])] + m


M15 = list(_matchings(list(R6)))
SPLITS = [(T, tuple(i for i in R6 if i not in T)) for T in itertools.combinations(R6, 3) if 0 in T]
I6T = [(T1, T2, P) for T1, T2 in SPLITS for P in itertools.permutations(T2)]
assert (len(M15), len(SPLITS), len(I6T)) == (15, 10, 60)


def _d2(r, i, j):
    return (r[i] - r[j]) ** 2


def _tri(r, T):
    return _d2(r, T[0], T[1]) * _d2(r, T[1], T[2]) * _d2(r, T[2], T[0])


def ic_roots(lc, r):
    """Invariants from the leading coefficient and the six roots (any ring)."""
    i2 = sum(_d2(r, *m[0]) * _d2(r, *m[1]) * _d2(r, *m[2]) for m in M15)
    i4 = sum(_tri(r, T1) * _tri(r, T2) for T1, T2 in SPLITS)
    i6 = sum(_tri(r, T1) * _tri(r, T2) * _d2(r, T1[0], P[0]) * _d2(r, T1[1], P[1]) * _d2(r, T1[2], P[2])
             for T1, T2, P in I6T)
    i10 = 1
    for i, j in itertools.combinations(R6, 2):
        i10 *= _d2(r, i, j)
    return {2: lc ** 2 * i2, 4: lc ** 4 * i4, 6: lc ** 6 * i6, 10: lc ** 10 * i10}


def _coeffs(lc, r):
    return [int(c) for c in sp.Poly(lc * sp.prod([X - ri for ri in r]), X).all_coeffs()[::-1]]


def _interpolate(k, rng):
    ms = [c for c in itertools.combinations_with_replacement(range(7), k) if sum(c) == 3 * k]
    n = len(ms) + 25
    samples = []
    for _ in range(n + 5):
        lc, r = rng.choice([1, 2, 3, -1, -2, 5]), [rng.randint(-9, 9) for _ in R6]
        samples.append((_coeffs(lc, r), ic_roots(lc, r)[k]))
    M = fmpq_mat(n, len(ms) + 1)
    for i, (cs, val) in enumerate(samples[:n]):
        for j, m in enumerate(ms):
            M[i, j] = int(sp.prod([cs[q] for q in m]))
        M[i, len(ms)] = int(val)
    Rr, rk = M.rref()
    assert rk == len(ms)
    sol = [None] * len(ms)
    for i in range(rk):
        piv = next(j for j in range(len(ms) + 1) if Rr[i, j] != 0)
        assert piv < len(ms)
        sol[piv] = sp.Rational(int(Rr[i, len(ms)].p), int(Rr[i, len(ms)].q))
    poly = sp.expand(sum(sol[j] * sp.prod([A[q] for q in ms[j]]) for j in range(len(ms))))
    held = all(poly.subs(dict(zip(A, cs))) == val for cs, val in samples[n:])
    assert held
    return poly


def ic_polys():
    """{2: I2, 4: I4, 6: I6, 10: I10} as sympy polynomials in a0..a6."""
    if CACHE.exists():
        loc = {str(a): a for a in A}
        out = {}
        for line in CACHE.read_text(encoding="utf-8").splitlines():
            k, expr = line.split(":", 1)
            out[int(k)] = sp.sympify(expr, locals=loc)
        return out
    rng = random.Random(20261011)
    out = {k: _interpolate(k, rng) for k in (2, 4, 6, 10)}
    gen = sum(A[i] * X ** i for i in range(7))
    assert sp.expand(sp.discriminant(gen, X) - out[10]) == 0
    assert sp.expand(out[2] - (6 * A[3] ** 2 - 16 * A[2] * A[4] + 40 * A[1] * A[5] - 240 * A[0] * A[6])) == 0
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    CACHE.write_text("".join(f"{k}:{sp.sstr(out[k])}\n" for k in (2, 4, 6, 10)), encoding="utf-8")
    return out


_POLYS = None


def ic(coeffs):
    """Invariants of sum coeffs[i] z^i (coeffs of length <= 7; entries numbers or sympy expressions)."""
    global _POLYS
    if _POLYS is None:
        _POLYS = {k: sp.Poly(v, *A) for k, v in ic_polys().items()}
    c = list(coeffs) + [0] * (7 - len(coeffs))
    sub = dict(zip(A, c))
    return {k: sp.expand(P.as_expr().subs(sub)) for k, P in _POLYS.items()}


_TOK = re.compile(r"\s*(?:(\d+)|([A-Za-z][A-Za-z0-9]*)|([-+*/^()]))")


def parse_safe(text, names):
    """Rational expression in the given symbol names from untrusted text, by a recursive-descent parser over
    integers, the names, + - * / ^ and parentheses (implicit product as in '2x'); anything else is rejected.
    Nothing is evaluated as code."""
    toks, pos = [], 0
    text = text.strip()
    while pos < len(text):
        m = _TOK.match(text, pos)
        if not m or m.end() == pos:
            raise ValueError(f"bad character at {pos}: {text[pos:pos + 10]!r}")
        if m.group(2) is not None and m.group(2) not in names:
            raise ValueError(f"unknown name {m.group(2)!r}")
        toks.append(m.group(1) or m.group(2) or m.group(3))
        pos = m.end()
    syms = {n: sp.Symbol(n) for n in names}
    i = 0

    def peek():
        return toks[i] if i < len(toks) else None

    def take(t=None):
        nonlocal i
        tok = peek()
        if t is not None and tok != t:
            raise ValueError(f"expected {t!r}, got {tok!r}")
        i += 1
        return tok

    def expr():
        v = term()
        while peek() in ('+', '-'):
            v = v + term() if take() == '+' else v - term()
        return v

    def term():
        v = factor()
        while peek() is not None and (peek() in ('*', '/') or peek() == '(' or peek().isalnum()):
            op = take() if peek() in ('*', '/') else '*'
            v = v * factor() if op == '*' else v / factor()
        return v

    def factor():
        if peek() in ('+', '-'):
            return factor() if take() == '+' else -factor()
        b = atom()
        if peek() == '^':
            take()
            return b ** factor()
        return b

    def atom():
        tok = take()
        if tok == '(':
            v = expr()
            take(')')
            return v
        if tok is not None and tok.isdigit():
            return sp.Integer(int(tok))
        if tok in syms:
            return syms[tok]
        raise ValueError(f"unexpected token {tok!r}")

    v = expr()
    if i != len(toks):
        raise ValueError(f"trailing tokens {toks[i:]}")
    return v


def poly_coeffs(expr, var):
    """Coefficients a0..a6 of a polynomial expression in var of degree <= 6."""
    p = sp.Poly(sp.expand(expr), var)
    c = p.all_coeffs()[::-1]
    assert len(c) <= 7
    return c + [0] * (7 - len(c))


# Elkies-Kumar, discriminant 8 (Theorem 17): H8 = P^2_{r,s}, weights (1, 2, 3, 5) for (I2, I4, I6, I10)
r_, s_ = sp.symbols('r s')
EK8 = {2: -4 * (3 * s_ + 8 * r_ - 2),
       4: 4 * (9 * r_ * s_ + 4 * r_ ** 2 + 4 * r_ + 1),
       6: -4 * (36 * r_ * s_ ** 2 + 94 * r_ ** 2 * s_ - 35 * r_ * s_ + 4 * s_ + 48 * r_ ** 3 + 40 * r_ ** 2 + 4 * r_ - 2),
       10: -8 * s_ ** 2 * r_ ** 3}
EK8_Y = 2 * (16 * r_ * s_ ** 2 + 32 * r_ ** 2 * s_ - 40 * r_ * s_ - s_ + 16 * r_ ** 3 + 24 * r_ ** 2 + 12 * r_ + 2)
W = {2: 1, 4: 2, 6: 3, 10: 5}


def ek8_points(I):
    """All RATIONAL (r, s) with (I2:I4:I6:I10) = EK8(r, s) in P(1,2,3,5), for rational invariants with I2 != 0
    (and EK8[2](r, s) != 0). Eliminates s from I4 E2^2 = I2^2 E4 and I6 E2^3 = I2^3 E6 by a resultant, takes the
    rational roots, then tests I10 E2^5 = I2^5 E10."""
    e4 = sp.expand(I[4] * EK8[2] ** 2 - I[2] ** 2 * EK8[4])
    e6 = sp.expand(I[6] * EK8[2] ** 3 - I[2] ** 3 * EK8[6])
    e10 = sp.expand(I[10] * EK8[2] ** 5 - I[2] ** 5 * EK8[10])
    res = sp.Poly(sp.resultant(e4, e6, s_), r_)
    out = []
    rroots = set(sp.Poly(res, r_).ground_roots().keys()) if not res.is_zero else set()
    for rr in rroots:
        if not rr.is_rational:
            continue
        for ss in sp.Poly(e4.subs(r_, rr), s_).ground_roots().keys():
            if not ss.is_rational:
                continue
            if EK8[2].subs({r_: rr, s_: ss}) == 0:
                continue
            if e6.subs({r_: rr, s_: ss}) == 0 and e10.subs({r_: rr, s_: ss}) == 0:
                out.append((rr, ss))
    return sorted(set(out))
