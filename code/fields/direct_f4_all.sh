#!/usr/bin/env bash
# The whole direct pipeline for C(F4,4) in the normal-form gauge (paper, Sections 4 and 6).
# usage: direct_f4_all.sh WORK [JOBS]        (NORECON=1: residues only)
# Memory: one direct_f4.py run peaks at about 2-3 GB; f4_full_recon with 4 processes at a few GB.
#   WORK/P.txt: primes used for the reconstruction (6 roots each), WORK/Q.txt: held-out primes (12 real
#   roots); copied from data/P.txt and data/Q.txt (the primes of the archived tables) if absent.
# Steps: direct_f4.py at every (p, w) -> f4_full_recon.py (F over k, 6 roots per prime)
#        -> direct_R_recon.py (R, Rinv over Q(zeta_13), 12 roots per prime).
# Outputs: WORK/bin/F4_4_F_part{1,2}.bin and WORK/bin/F4_4_R_direct.txt.  f4_full_recon.py reads the
# block layout from $F4_ARCHIVE/F4_4_F_part{1,2}.bin (index only).  Disk: about 0.45 GB per (p, w).
set -eo pipefail
WORK=$(mkdir -p "$1" && cd "$1" && pwd); JOBS=${2:-2}
HERE=$(cd "$(dirname "$0")" && pwd)
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
[ -f "$WORK/P.txt" ] || cp "$HERE/data/P.txt" "$WORK/P.txt"
[ -f "$WORK/Q.txt" ] || cp "$HERE/data/Q.txt" "$WORK/Q.txt"
( for p in $(cat "$WORK/P.txt"); do for w in 1 3 5 7 9 11; do echo $p $w; done; done
  for p in $(cat "$WORK/Q.txt"); do for w in 1 3 5 7 9 11 15 17 19 21 23 25; do echo $p $w; done; done ) |
while read p w; do [ -f "$WORK/F_${p}_$w.npy" ] || echo "$p $w"; done |
  xargs -P "$JOBS" -L 1 sh -c "python3 '$HERE/direct_f4.py' \$0 \$1 '$WORK' --nproc 1 >> '$WORK/direct.log' 2>&1"
[ -n "$NORECON" ] && exit 0
python3 "$HERE/f4_full_recon.py" "$WORK" "$WORK/bin" --nproc 4 | tee "$WORK/recon.log"
python3 "$HERE/direct_R_recon.py" "$WORK" "$WORK/bin/F4_4_R_direct.txt" | tee "$WORK/recon_R.log"
