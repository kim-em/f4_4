#!/usr/bin/env bash
# The invariant forms G of C(F4,4) and the checks of the paper's Section "Unitarity".
# usage: run_forms.sh [WORK] [JOBS]      (WORK default ../../work/unitary; table: $F4_ARCHIVE)
# Residues of G at 15 primes (the 2nd to 16th of lin.find_primes(52, 40), as in logs/unitary/f4_forms_runs.log) and 8 roots, about 1 minute
# with 2 jobs; exact G over k from the first 13 primes (2 held out); then positivity, square classes
# and the isometry check of the whole table at two (prime, root) pairs.
set -eo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
WORK=$(mkdir -p "${1:-$HERE/../../work/unitary}" && cd "${1:-$HERE/../../work/unitary}" && pwd)
JOBS=${2:-2}
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
P="1073740201 1073739577 1073739473 1073739421 1073737549 1073737393 1073737081 1073735573 1073735417 1073733233 1073733077 1073732609 1073732557 1073731829 1073731777"
mkdir -p "$WORK/G"
for p in $P; do echo $p; done | xargs -P "$JOBS" -I{} python3 "$HERE/f4_forms.py" "$WORK/G" {} 1 3 5 7 9 11 25 23
python3 "$HERE/f4_forms_recon.py" "$WORK/G" 13 "$WORK/G_exact.pkl" $P
python3 "$HERE/f4_forms_analyse.py" "$WORK/G_exact.pkl" 40
python3 "$HERE/f4_classes.py" "$WORK/G_exact.pkl" 300
python3 "$HERE/f4_isometry_check.py" "$WORK/G_exact.pkl" 1073741189 1 "$JOBS"
python3 "$HERE/f4_isometry_check.py" "$WORK/G_exact.pkl" 1073726317 7 "$JOBS"
