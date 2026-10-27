Supporting code and data for the note
"An explicit five-dimensional family of K3 surfaces with real multiplication by Q(sqrt2), and the limits of the
genus-two pencil method"
================================================================================================================

Status: accompanies the note (K. Saller, 2026), DOI 10.5281/zenodo.23710067; not peer reviewed.

Setting. S is the minimal resolution of a double plane w^2 = F(x0, x1, x2) branched along a sextic B with at most
simple singularities, P0 = (0:0:1), and the pencil is that of the lines through P0, with fibres
C_t: w^2 = f_t(z) = F(t, 1, z). I2, I4, I6, I10 are the Igusa-Clebsch invariants of binary sextics in the
normalisation of [EK] (N. D. Elkies and A. Kumar, K3 surfaces and equations for Hilbert modular surfaces, Algebra
Number Theory 8 (2014); arXiv:1209.3527v3), with I10 the discriminant, and R30 = R^2 is the square of the
invariant R of degree 15, a form of weight 15 in (I2, I4, I6, I10). For D = 5, 8, 12, 13, 17, [EK] give the
Hilbert modular surface Y_-(D) as a double cover of a plane and the Igusa-Clebsch invariants as functions on that
plane; E = (E2 : E4 : E6 : E10) is this map, homogenised. Gamma is the curve of the invariant points of the sextics
with two double roots, theta its point (z^3 - 1)^2. TTV are Tautz, Top and Verberkmoes; the curves
T_sigma: eta^2 = (xi + 2)(xi^4 - 4 xi^2 + 2 + sigma) are their curves for n = 4.

Role of the computations. Section 8 of the note lists the items (C1)-(C16); Section 3 below says which script
makes each of them. The note marks the statements that use them as computed and cites the item. The arithmetic is
exact, over Z, Q, number fields and finite fields, with three exceptions that the note does not use: c12 computes
the Frobenius traces with a floating-point FFT (the note's traces come from the exact transform of c22, which also
compares them with c12's values) and the absolute values of the roots of the characteristic polynomial in floating
point; c15 would confirm a hit at an irrational point numerically (40 digits), but in its output every hit is at a
rational point and is confirmed exactly; and c02, which the note does not use, prints floating-point
approximations. Where the note says "random", the choices come from seeded generators (fixed seeds, Section 2).

Not checked by any script: the proofs (the scripts work at the points, lines, curves and members named in the
note); the cited results, among them those of [EK] beyond the formulas and tables that c01, c03 and c05 check, of
[TTV], [CQ02] (the automorphism groups of curves of genus two and the normal forms of Theorems 8-10), [K3M] and
[Pow].

Labels. The headers of the scripts name statements of the working notes NOTE.md (not part of this repository),
whose numbering changed while the note was written, and steps of the working log (for example "76.4" or "78.7b");
Section 4 translates them into the numbering of the note. The output files name no statements.

Licence: the code (*.py, *.sh) is under the MIT licence; the outputs and this text are under CC BY 4.0 (see the
LICENSE files of the repository). No third-party files are included; see README.md for the formulas of [EK] that
the scripts and outputs/ek_ic_formulas.txt contain.


1. Requirements
---------------
Python 3 with SymPy (every script except c02), python-flint (k3lib.py and c11; k3lib.py is loaded by c01, c03,
c04, c05, c06, c07, c08, c09, c13, c15, c16, c17, c19 and c21, directly or through c06 or c08) and NumPy (c12, and
c20 and c22, which load c12). compare_outputs.py uses the standard library only. Tested versions: Python 3.12.9,
SymPy 1.13.1, python-flint 0.9.0 and NumPy 2.1.3, on Windows 10; these are the versions that Section 8 of the note
names.

c05 needs the arXiv source of [EK], version 3: the file https://arxiv.org/e-print/1209.3527v3, a gzipped tar
archive, extracted into a folder. c05 reads hmspaper.tex and the files D/igusaD.txt for D = 5, 8, 12, 13, 17 of
that folder as data, with a parser that accepts only arithmetic expressions (no code from the source is run).


2. Layout and how to run
------------------------
code/      c01_ek8_check.py ... c22_ttv_ntt_traces.py, and k3lib.py (helpers: the Igusa-Clebsch invariants, the
           formulas of [EK] for D = 8, the safe parser)
outputs/   the reference outputs: NAME.txt, the screen output of each run (five runs of c12, see below); and two
           files that the scripts write themselves: ic_polys.txt (I2, I4, I6, I10 as polynomials in the
           coefficients of a sextic, written by k3lib.py on first use, which checks their normalisation) and
           ek_ic_formulas.txt (the Igusa-Clebsch maps of [EK] for D = 5, 8, 12, 13, 17, written by c05 and read
           by c06 and c07, and through c06 by c08, c09 and c21)
run_all.sh, compare_outputs.py

The scripts load k3lib.py, c06, c08 or c12 from their own folder, and read and write the files of ../outputs
relative to that folder. run_all.sh therefore copies code/ to OUT/code and runs the scripts there, so that the
reference outputs in outputs/ are not touched.

From the top folder of this repository:
    sh run_all.sh rerun EK_SOURCE_FOLDER
    python3 -I compare_outputs.py outputs rerun/outputs
run_all.sh runs every script once, in the order c01, ..., c22, each as python3 -I -B SCRIPT [ARGUMENTS] with its
screen output in OUT/outputs/NAME.txt, and prints a line "NAME exit N" for each run. The arguments are: the
source folder for c05; for c12 the prime, the largest k and, for the second member of Remark 4.2, its
coefficients (c0 c1 s0 ... s5):
    c12_ttv_frob_p3      3 12          (first member, p = 3: Step 7 of Theorem 4.1)
    c12_ttv_frob_p5      5 9           (first member, p = 5; not used by the note)
    c12_ttv_frob_p7      7 8           (first member, p = 7: bad reduction; not used by the note)
    c12_ttv_frob_m2_p3   3 12 5 3 5 -1 6 5 5 4   (second member, p = 3; not used by the note)
    c12_ttv_frob_m2_p7   7 8 5 3 5 -1 6 5 5 4    (second member, p = 7: Remark 4.2)
and none for the other scripts; c11 and c14 then use their default seeds 11 and 5. Without EK_SOURCE_FOLDER, c05
is skipped and the reference ek_ic_formulas.txt is copied to OUT/outputs instead. The shell script calls python3;
set the environment variable PYTHON to use another interpreter name (for example PYTHON=python).
compare_outputs.py compares each .txt file with the reference, line by line; it ignores the line ends (the
reference outputs have CRLF line ends, from Windows) and the running times printed by c20 and c22. Its expected
last line is "ALL EQUAL". Running times: about half an hour for the whole run on the machine of Section 5, most
of it in c22, two of the runs of c12 and c09 (Section 5).


3. The scripts and the items of Section 8 of the note
-----------------------------------------------------
Statement numbers are those of the note.

(C1)  c01: D = 8, [EK] Theorem 17 against the sixteen rows of [EK] Table 2 (each listed point (r, s) is the only
      rational point of the plane with the invariants of the listed sextic, and the double cover defining
      Y_-(8) has a rational point over it); the invariants of (z^3 - 1)^2 and E on the line r = 0
      (16 E4 = E2^2 only for s = 0, 4/3, and 64 E6 = E2^3 only at s = 0 among these); the points of the plane
      over T_sigma for six values of sigma.
(C2)  c03: part A, E on the lines R = 0 and S = 0 and the equation of Gamma on R = 0; part B, the curves T_sigma
      lie over (r, s) = (-(sigma + 2)/16, (sigma - 2)/16), as an identity in sigma; part C, D = 5: [EK]
      Corollary 15 against the sixteen rows of [EK] Table 1, and the pencil of Remark 5.7 (the degrees of the
      I_2k(t), the lift (g(t), h(t)), the factor 2000(t^3 + 1)^2, h_1 of degree 15 with simple roots).
(C3)  c04: the order of I10 at a line through P0 of contact type (2,2,1,1), (2,2,2) and (3,1,1,1) with a random
      smooth sextic B (P0 not on B): 2, 3 and 2; smoothness certified modulo two primes (degree-13 part of the
      Jacobian ring).
(C4)  c05: the Igusa-Clebsch maps of [EK] for D = 5, 8, 12, 13, 17 from the arXiv source; against [EK]
      Corollary 15 and Theorem 17 for D = 5, 8, and against the five tables of examples. c09: the factors of odd
      multiplicity of R30(E) against the branch functions of Y_-(D) in [EK], including the extra factor f - 1
      for D = 12.
(C5)  c06: for D = 12, 13, 17 (Lemma 5.9): N, m and the forms E_2k; the base points; the blow-ups that resolve
      E, with the orders o_2k along the exceptional curves and the restricted maps; the images of the special
      lines; the values of epsilon.
(C6)  c08: R30 as the unique form of weight 15, up to a scalar, vanishing on the sextics g(z^2); its restrictions
      to the line f = -u (D = 12) and to R_b(X_b) (D = 13), and the identification of their images
      (Lemma 5.13).
(C7)  c10: Lemma 6.1, the map Psi, its base points with multiplicities, its degree, the class A with A^2 = 2, and
      the curves contracted by Psi.
(C8)  c11: Theorem 4.1, the identification of C_t with the twist of T_sigma(t) by c(t); the cover Y_sigma of
      T_sigma of genus four, its differentials and the action of gamma^* + (gamma^-1)^* by +-sqrt2; the rank
      computation of Step 5 (effective dimension 5) at random rational parameters.
(C9)  c22: the sums Tr(Frob_q | V) of Step 7 of Theorem 4.1 (first member, q = 3^k, k <= 12) and of Remark 4.2
      (second member, q = 7^k, k <= 8), by a number-theoretic transform modulo 2^31 - 1, compared with direct
      counting for q <= 3^8 and q <= 7^4 (c12's exact counting routine) and with c12's floating-point FFT values.
      c20: the same exact counts, for the first member for k <= 8 and the second member for k <= 4, compared with
      the FFT values. c12: the original computation with the FFT, five runs (Section 2).
(C10) c13: Lemma 6.1 and Proposition 6.4(d): the singular points of Q = U B8 and their types, the binary forms of
      the proof, the lines along which Q restricts to a square, and R30(E) and E10 on these lines.
(C11) c14: Steps 5 and 7 of Theorem 4.1: the identity F o g = d lambda^5 F'' and the rank of the orbit map at
      three random rational points; from the traces in c12_ttv_frob_p3.txt and c12_ttv_frob_m2_p7.txt (equal to
      the exact traces of c22) the coefficients of the characteristic polynomial by Newton's identities, the sign
      epsilon, the comparison with the traces for k >= 9, and the factorisations over Q and Q(sqrt2). c22 derives
      the same polynomials from its exact traces alone.
(C12) c15: Lemma 6.7 and the paragraph after it: on the lines 16R = U, 4R + 4S = 7U and 4R + 4S + U = 0, the
      points with E10 != 0 and R30(E) = 0 and the order of R30(E) there; for each point, whether the curve lies on
      the image of the families y^2 = x^5 + x^3 + t x, y^2 = x^6 + x^3 + t, or is one of x^5 - x, x^6 - 1,
      x^5 - 1 (eliminating t by resultants, and reducing modulo the minimal polynomial of the point).
(C13) c16: Step 3 of the proof of Lemma 6.6: ord_u R30(f* + u g) for random integral sextics g, for
      f* = x^6 + x^3 + 1/2 and f* = x^5 - x, and for a curve with automorphism group V4 as a control.
(C14) c17: Lemma 6.9(a), the curve over the point (kappa : 1 : 0) of U = 0 (the script calls the parameter rho).
(C15) c18: Lemma 6.9(c) and Propositions 6.8 and 6.10: the curve E_tau (the script's E1) and j(E_tau), and the
      enumeration of the ramification profiles and of the resulting bounds; as a calibration, the same
      enumeration for the line of the curves T_sigma (Proposition 6.5).
(C16) c19: Lemma 5.1, the invariants of x^3 y^3, (z^3 - 1)^2, (z^2 - 1)^2 (z^2 - c) and
      z^2 (z - 1)(z - a)(z - b). c21: the two points of Gamma with 16 I4 = I2^2, and that the points of Lemmas
      5.1(c), 5.3 and 5.9(iii) with 16 I4 = I2^2 other than theta are theta'.

Not used by the note: c02 (the cusp resolutions of Y_-(D) and the degeneration types of genus-two fibres crossing
them, D <= 33), c07 (the fibre types of the elliptic K3 families of [EK] over the special curves of D = 12 and 13),
and the c12 runs p5, p7 and m2_p3.


4. Labels in the scripts
------------------------
The headers name statements of NOTE.md, in the numbering of the day each script was written. In the numbering of
the note:
  c06  "Theorems A13 and A17"                          Theorems 5.11 (D = 13) and 5.10 (D = 17)
       "Prop. 2.2"                                     Proposition 3.2 and Table 1
  c07, c08  "Theorem A13", "Prop./Proposition 6.6"     Theorem 5.11, Proposition 5.12 (the line f = -1, D = 12)
       "NOTE.md, Section 6"                            Section 5 (discriminants 12, 13 and 17)
  c10  "NOTE.md, Section 5"                            Section 6 (Lemma 6.1)
  c11, c12  "Check 3.3", "Theorem C"                   the family of Theorem 4.1
  c13  "Proposition 5.3(c)"                            Proposition 6.4(d)
  c14  "Theorem C, proof steps 5 and 7"                Theorem 4.1, Steps 5 and 7
  c15, c16, c17  "NOTE.md, Section 5.4"                Lemma 6.7, Lemma 6.6 (Step 3), Lemma 6.9(a)
  c18  "Section 5.4", "Proposition 5.4"                Lemma 6.9(c), Propositions 6.8 and 6.10; Proposition 6.5
  c19, c21  "the LaTeX note, Section 5"                Section 5 (Lemma 5.1 and (C16))
  c20, c22  "Theorem 4.1", "Remark 4.2"                the same numbers
"the [K3M] repository" (c04, k3lib.py) is the repository of the earlier note on Q(sqrt5), cited in the note as
[K3M]; "referee R1, R2, R4" are reviews of drafts of the note; DRAFT_78_7.md (c18) is the working draft of what
became NOTE.md, Section 5.4, also not part of this repository.


5. Reference outputs and the rerun of 27 October 2026
-----------------------------------------------------
Reference outputs. The files in outputs/ were made in the working tree between 24 and 26 October 2026, each by the
last version of its script; the platform and the versions of those runs are not recorded in them.

Changes against the working tree. c05_ek_igusa_parse.py: the source folder is now a required argument; in the
working tree it defaulted to a folder on the author's machine (lines 7 and 22-24). c22_ttv_ntt_traces.py: one
phrase of the docstring (line 1). run_all.sh, compare_outputs.py, the README files and the licences are new. All
other files are unchanged.

The rerun. When this repository was assembled, on 27 October 2026, run_all.sh was run from its files on the
author's machine (Windows 10, the sh of Git for Windows, and the versions listed in Section 1), with
PYTHON=python, the extracted arXiv source of [EK], version 3, as EK_SOURCE_FOLDER, and a folder outside the
repository as OUT. It took 32 minutes: c22 about 11, the run c12_ttv_frob_m2_p7 about 8.5, c09 about 6.5, the run
c12_ttv_frob_p5 about 2, the runs c12_ttv_frob_p3 and c12_ttv_frob_m2_p3 about 1 each, and every other run less
than half a minute. All 26 runs ended with exit status 0, and compare_outputs.py found all 28 files equal to the
reference outputs ("ALL EQUAL"). 26 of them, among them ek_ic_formulas.txt (written by c05 from the arXiv source)
and ic_polys.txt, are equal byte for byte; c20_ttv_exact_traces.txt and c22_ttv_ntt_traces.txt differ only in the
running times they print. The reference outputs were made on the same machine, so the comparison shows only that
two runs there agree.
