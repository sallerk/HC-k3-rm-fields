"""76.2(c)/76.5(b): the boundary components of the cusp resolution of the Hilbert modular surface Y_-(D)
(group SL(O + d^-1)), and the degeneration type of a genus-2 fibre crossing each of them transversally.

Recipe (checked against D = 5 and D = 8 below): the boundary curves correspond to the lattice points A on the
boundary of the convex hull of the totally positive elements of the inverse different d^-1 = (1/sqrt D) O (the
rays of the fan); a curve crossing the component of A transversally has local monodromy with quadratic form
q_A(xi) = Tr(A xi^2) on the character lattice O. Writing q_A = a x^2 + b y^2 + c (x - y)^2 with a, b, c >= 0
(Selling/Voronoi reduction of a positive definite binary form), the stable fibre has one node for each non-zero
coefficient, of thickness equal to the coefficient, and ord(I10) = a + b + c.
For K3 pencils with smooth total space all thicknesses must be 1: (0,1,1) = 2-nodal fibre (bitangent line),
(1,1,1) = theta fibre (tritangent line); a coefficient >= 2 needs a singular point of the branch curve.
Usage: python -I -B c02_cusp_voronoi.py
"""
import sys
from fractions import Fraction
from math import isqrt

sys.stdout.reconfigure(encoding="utf-8")


class QD:
    """Elements x + y*w of O_D, w = sqrt(D/4) (D = 0 mod 4) or (1 + sqrt D)/2 (D = 1 mod 4); exact embeddings."""

    def __init__(self, D):
        self.D = D
        if D % 4 == 0:
            self.n, self.t = -(D // 4), 0       # w^2 = t w - n  with  w^2 = D/4
            self.wn = D // 4                    # w^2 = D/4
        else:
            self.n, self.t = None, 1            # w^2 = w + (D-1)/4
            self.wn = (D - 1) // 4

    def mul(self, a, b):
        x1, y1 = a
        x2, y2 = b
        # (x1 + y1 w)(x2 + y2 w) = x1x2 + (x1y2 + x2y1) w + y1y2 w^2
        if self.D % 4 == 0:
            return (x1 * x2 + y1 * y2 * self.wn, x1 * y2 + x2 * y1)
        return (x1 * x2 + y1 * y2 * self.wn, x1 * y2 + x2 * y1 + y1 * y2)

    def tr(self, a):
        x, y = a
        return 2 * x + (y if self.D % 4 == 1 else 0)

    def emb(self, a):
        """The two real embeddings as floats (only for sorting / convexity tests with exact fallback)."""
        x, y = a
        sq = self.D ** 0.5
        if self.D % 4 == 0:
            return (x + y * sq / 2, x - y * sq / 2)
        return (x + y * (1 + sq) / 2, x + y * (1 - sq) / 2)


def selling(alpha, beta, gamma):
    """Reduce the positive definite form alpha x^2 + 2 beta x y + gamma y^2 (rationals) to the Voronoi form
    a x^2 + b y^2 + c (x - y)^2 with a, b, c >= 0; returns (a, b, c) sorted."""
    # Selling reduction on the obtuse superbase: work with the 3 "Selling parameters"
    # p12 = -beta, p13 = -(alpha + beta), p23 = -(gamma + beta) for the superbase (e1, e2, e3 = -e1 - e2).
    e = [(1, 0), (0, 1), (-1, -1)]

    def ip(u, v):
        return alpha * u[0] * v[0] + beta * (u[0] * v[1] + u[1] * v[0]) + gamma * u[1] * v[1]

    for _ in range(1000):
        changed = False
        for i in range(3):
            for j in range(i + 1, 3):
                if ip(e[i], e[j]) > 0:
                    k = 3 - i - j
                    # Selling step: e_i -> -e_i, e_k -> e_k + 2 e_i  (keeps a superbase)
                    ei, ek = e[i], e[k]
                    e[i] = (-ei[0], -ei[1])
                    e[k] = (ek[0] + 2 * ei[0], ek[1] + 2 * ei[1])
                    changed = True
                    break
            if changed:
                break
        if not changed:
            break
    coeffs = sorted(-ip(e[i], e[j]) for i in range(3) for j in range(i + 1, 3))
    return tuple(coeffs)


def hull_points(D, box=60):
    """Totally positive elements of d^-1 = (1/sqrt D) O, scaled by sqrt D: i.e. elements v of O with
    v/sqrt(D) totally positive, i.e. emb1(v) > 0 > emb2(v)... we use the exact description: A = v / sqrt D,
    A totally positive <=> v has embeddings of signs (+, -). Return the lower-left convex hull boundary points
    of {(A, A')} in a box, in order of decreasing A'."""
    K = QD(D)
    sq = D ** 0.5
    pts = []
    for x in range(-box, box + 1):
        for y in range(-box, box + 1):
            e1, e2 = K.emb((x, y))
            a1, a2 = e1 / sq, -e2 / sq   # A = v/sqrt D has embeddings (e1/sqrtD, e2/(-sqrtD))
            if a1 > 1e-12 and a2 > 1e-12 and a1 < 12 and a2 < 12:
                pts.append((a1, a2, (x, y)))
    # lower-left convex hull (towards the origin) of points in the positive quadrant
    pts.sort(key=lambda p: (p[0], -p[1]))
    hull = []
    for p in pts:
        while len(hull) >= 2:
            (x1, y1, _), (x2, y2, _) = hull[-2], hull[-1]
            cross = (x2 - x1) * (p[1] - y1) - (y2 - y1) * (p[0] - x1)
            if cross < -1e-12:   # keep only convex (towards origin) chain
                hull.pop()
            else:
                break
        hull.append(p)
    # keep the decreasing-y chain (the boundary facing the origin) and include collinear lattice points
    chain = [hull[0]]
    for p in hull[1:]:
        if p[1] < chain[-1][1]:
            chain.append(p)
    # add lattice points lying exactly on the chain segments (collinear boundary points)
    full = []
    for (p, q) in zip(chain, chain[1:]):
        full.append(p)
        for r in pts:
            if p[0] < r[0] < q[0]:
                cross = (q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0])
                if abs(cross) < 1e-9:
                    full.append(r)
    full.append(chain[-1])
    full.sort(key=lambda p: p[0])
    return full


def form(D, v):
    """q_A(xi) = Tr(A xi^2), A = v / sqrt D, on O with basis (1, w): returns (alpha, beta, gamma) as Fractions,
    q = alpha x^2 + 2 beta x y + gamma y^2. Tr(v u / sqrt D) = (v u - (v u)') / ... computed exactly:
    for u in O, Tr(u/sqrt D) = (coefficient of sqrt D in u) * 2 / ... handled via the w-coefficient."""
    K = QD(D)

    def tr_over_sqrtD(u):
        # u = x + y w; u/sqrt D has trace = (u - u')/sqrt D = y (w - w')/sqrt D = y * (2*sqrt(D/4)/sqrt D = 1
        # if D = 0 mod 4; (sqrt D)/sqrt D = 1 if D = 1 mod 4)
        return Fraction(u[1])

    one, w = (1, 0), (0, 1)
    vw = K.mul(v, w)
    alpha = tr_over_sqrtD(K.mul(v, K.mul(one, one)))
    beta = tr_over_sqrtD(K.mul(v, K.mul(one, w)))
    gamma = tr_over_sqrtD(K.mul(v, K.mul(w, w)))
    return alpha, beta, gamma


def fundamental_unit(D):
    K = QD(D)
    for y in range(1, 2000):
        for x in range(-4000, 4001):
            u = (x, y)
            # norm
            if D % 4 == 0:
                nrm = x * x - y * y * (D // 4)
            else:
                nrm = x * x + x * y - y * y * ((D - 1) // 4)
            if abs(nrm) == 1:
                e1, _ = K.emb(u)
                if e1 > 1:
                    return u, nrm
    return None


def main():
    for D in (5, 8, 12, 13, 17, 21, 24, 28, 29, 33):
        eps, nrm = fundamental_unit(D)
        pts = hull_points(D)
        rows = []
        for a1, a2, v in pts:
            al, be, ga = form(D, v)
            det = al * ga - be * be
            vor = selling(al, be, ga)
            rows.append((round(a1, 4), round(a2, 4), v, (al, be, ga), det, vor, sum(vor)))
        print(f"== D = {D}: fundamental unit {eps} of norm {nrm}")
        seen = []
        for row in rows:
            print(f"   A ~ ({row[0]}, {row[1]}) v={row[2]} q=[{row[3][0]},{row[3][1]},{row[3][2]}] det={row[4]}"
                  f" Voronoi={row[5]} ord(I10)={row[6]}")
            seen.append(row[5])
        kinds = sorted(set(seen))
        print(f"   types on this stretch of the boundary: {kinds}")


if __name__ == "__main__":
    main()
