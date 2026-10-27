# HC-k3-rm-fields

Code and data supporting the note *An explicit five-dimensional family of K3 surfaces with real multiplication by
Q(√2), and the limits of the genus-two pencil method* (K. Saller, 2026), DOI:
[10.5281/zenodo.23710067](https://doi.org/10.5281/zenodo.23710067).

**Subject.** On a K3 double plane w² = F(x0, x1, x2), the lines through a point P0 cut out a pencil of curves of
genus two. If their Jacobians have real multiplication by a real quadratic order over the base of the pencil, the
real multiplication is induced on the K3 surface by an algebraic cycle. The note studies this mechanism for real
quadratic fields other than Q(√5). For Q(√2) it gives a five-dimensional family of double planes, branched along a
line and a quintic, whose pencils consist of quadratic twists of curves of Tautz, Top and Verberkmoes; √2 acts
through an algebraic cycle on every member, and the very general member has Picard number six, so that the Hodge
conjecture holds for its powers S^n, n ≤ 7, and its Hilbert schemes S^[n], n ≤ 35. The note then bounds the
mechanism: for a smooth branch curve, a centre off the curve and the maximal order, it gives real multiplication
only by Q(√5) among the discriminants 5, 8, 12, 13 and 17; for the maximal order Z[√2], branch curves with
simple singularities and a centre that is not a singular point of the curve, it gives at most six moduli, while
the period domain has dimension eight.

**Role of the computations.** Section 8 of the note lists sixteen computations, (C1)–(C16); several steps of the
proofs rest on them. Among them are the Frobenius traces that bound the Picard number in Theorem 4.1 ((C9) and
(C11)), the Igusa–Clebsch maps of Elkies and Kumar for the discriminants 12, 13 and 17, read from the auxiliary
files of the arXiv source of their paper ((C4) and (C5)), and the special points of Lemma 6.7 ((C12)). The
arithmetic of these items is exact.
[`README.txt`](README.txt) gives the requirements, the commands, the script behind each item, how the labels in
the scripts translate into the numbering of the note, and the result of the rerun made when this repository was
assembled.

**Status:** not peer reviewed.

**Disclosure:** the arguments, the computations and the code were produced with the assistance of an AI system
(Claude, by Anthropic).

## Licences

- **Code** (`*.py`, `*.sh`): MIT License, see [`LICENSE`](LICENSE).
- **Everything else** (README files and computed outputs): Creative Commons Attribution 4.0 International
  (CC BY 4.0), see [`LICENSE-CC-BY-4.0`](LICENSE-CC-BY-4.0).
- No third-party files are included. Some scripts contain formulas and table entries typed from the paper of
  Elkies and Kumar, and `outputs/ek_ic_formulas.txt` holds their Igusa–Clebsch maps as read by `c05` from the
  arXiv source of that paper (arXiv:1209.3527v3), which is not included. These are cited data; the licences
  cover the author's own work.
