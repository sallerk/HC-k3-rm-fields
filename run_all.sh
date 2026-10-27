#!/bin/sh
# run_all.sh -- rerun the computations of this repository.
# Usage:   sh run_all.sh [OUT [EK]]
#   OUT  the folder for the rerun (default: rerun). The scripts are copied to OUT/code and write their results
#        into OUT/outputs, so the reference outputs in outputs/ are not touched.
#   EK   the folder with the extracted arXiv source of Elkies-Kumar, arXiv:1209.3527v3 (needed by c05 only).
#        Without it, c05 is skipped and the reference file outputs/ek_ic_formulas.txt, which c06 and c07 read
#        (and through c06 also c08, c09 and c21), is copied to OUT/outputs instead.
# Then:    python3 -I compare_outputs.py outputs OUT/outputs
# The interpreter is python3, or the value of the environment variable PYTHON (for example PYTHON=python).
# The scripts run one after the other (see README.txt, Sections 2 and 5, for the running times). A script that
# fails does not stop the others; each line "NAME exit N" reports the exit status of one run.
OUT=${1:-rerun}
EK=${2:-}
PY=${PYTHON:-python3}
T=$(dirname "$0")
C="$OUT/code"
O="$OUT/outputs"
mkdir -p "$C" "$O" || exit 1
cp "$T"/code/*.py "$C"/ || exit 1
run() {
  out=$1; script=$2; shift 2
  "$PY" -I -B "$C/$script.py" "$@" > "$O/$out.txt"
  echo "$out exit $?"
}
run c01_ek8_check c01_ek8_check
run c02_cusp_voronoi c02_cusp_voronoi
run c03_d8_obstruction c03_d8_obstruction
run c04_local_orders c04_local_orders
if [ -n "$EK" ]; then
  run c05_ek_igusa_parse c05_ek_igusa_parse "$EK"
else
  cp "$T/outputs/ek_ic_formulas.txt" "$O/" || exit 1
  echo "c05_ek_igusa_parse skipped (no source folder given); the reference ek_ic_formulas.txt is used"
fi
run c06_ek_boundary c06_ek_boundary
run c07_ek13_fzero c07_ek13_fzero
run c08_parity c08_parity
run c09_r30_branch c09_r30_branch
run c10_d8_class c10_d8_class
run c11_ttv_rank c11_ttv_rank
run c12_ttv_frob_p3 c12_ttv_frob 3 12
run c12_ttv_frob_p5 c12_ttv_frob 5 9
run c12_ttv_frob_p7 c12_ttv_frob 7 8
run c12_ttv_frob_m2_p3 c12_ttv_frob 3 12 5 3 5 -1 6 5 5 4
run c12_ttv_frob_m2_p7 c12_ttv_frob 7 8 5 3 5 -1 6 5 5 4
run c13_d8_lines c13_d8_lines
run c14_ttv_checks c14_ttv_checks
run c15_d8_special c15_d8_special
run c16_local_models c16_local_models
run c17_d8_u0 c17_d8_u0
run c18_d8_counts c18_d8_counts
run c19_sextic_points c19_sextic_points
run c20_ttv_exact_traces c20_ttv_exact_traces
run c21_theta_prime c21_theta_prime
run c22_ttv_ntt_traces c22_ttv_ntt_traces
