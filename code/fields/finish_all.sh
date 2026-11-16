#!/usr/bin/env bash
# The long verification runs on the archived normal-form table, and the comparison of the tables
# reconstructed by the direct pipeline (direct_f4_all.sh WORK, run before) with the archived ones.
# Writes LOGDIR/finish_all.log (default ../../logs/fields) and one log per run.  Total about 20 hours
# with 2 workers per verifier; memory at most about 16 GB.
# usage: finish_all.sh WORK [LOGDIR]       archive: $F4_ARCHIVE (default ../../archive)
# (The run recorded in deposit/logs/finish_all.log was started after verify_hex_bin.py and
# direct_f4_all.sh with NORECON=1, and waited for them; here the steps run one after the other.)
set -eo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
W=$(cd "$1" && pwd)
LOGS=$(mkdir -p "${2:-$HERE/../../logs/fields}" && cd "${2:-$HERE/../../logs/fields}" && pwd)
A=${F4_ARCHIVE:-$HERE/../../archive}
B="$A/F4_4_F_part1.bin $A/F4_4_F_part2.bin"
R=$A/F4_4_R_normalform.txt
F4=$HERE/../f4
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
LOG=$LOGS/finish_all.log
echo "[$(date +%T)] started" >> $LOG
python3 $F4/verify_hex_bin.py $R $B --jobs 2 > $LOGS/verify_hex_descended.log 2>&1
echo "[$(date +%T)] hexagons: $(grep CONCLUSION $LOGS/verify_hex_descended.log | tr '\n' ' ')" >> $LOG
python3 $F4/pent_exact_subset.py 0000,0001,1000,0004,2000 $B --jobs 2 > $LOGS/pent_exact_subset5_descended.log 2>&1
python3 $F4/pent_exact_subset.py 0000,0001,0002,0010,1000,0004,2000 $B --jobs 2 > $LOGS/pent_exact_subset7_descended.log 2>&1
python3 $F4/pent_exact_subset.py maxdim=50 $B --jobs 2 > $LOGS/pent_exact_maxdim50_descended.log 2>&1
echo "[$(date +%T)] subsets: $(grep -h CONCLUSION $LOGS/pent_exact_subset5_descended.log $LOGS/pent_exact_subset7_descended.log $LOGS/pent_exact_maxdim50_descended.log | tr '\n' ' ')" >> $LOG
python3 $F4/pent_reduced.py $B --jobs 2 > $LOGS/pent_reduced_descended.log 2>&1
echo "[$(date +%T)] pentagons with fourth object V: $(grep CONCLUSION $LOGS/pent_reduced_descended.log)" >> $LOG
python3 $HERE/f4_full_recon.py $W $W/bin --nproc 4 > $LOGS/direct_f4_recon.log 2>&1
for i in 1 2; do
  if cmp -s $W/bin/F4_4_F_part$i.bin $A/F4_4_F_part$i.bin; then r=identical; else r=DIFFERENT; fi
  echo "[$(date +%T)] direct F table part $i vs archive: $r" >> $LOG
done
python3 $HERE/direct_R_recon.py $W $W/bin/F4_4_R_direct.txt > $LOGS/direct_R_recon.log 2>&1
if cmp -s $W/bin/F4_4_R_direct.txt $R; then r=identical; else r=DIFFERENT; fi
echo "[$(date +%T)] direct R, Rinv vs archive F4_4_R_normalform.txt: $r" >> $LOG
for pw in "1073741189 1" "1073740201 3" "1073726161 25"; do
  python3 $HERE/f4_nf_checks.py braid $W $pw >> $LOGS/f4_nf_braid_checks.log 2>&1
done
echo "[$(date +%T)] twists/S: $(grep -c True $LOGS/f4_nf_braid_checks.log) True, $(grep -c False $LOGS/f4_nf_braid_checks.log) False" >> $LOG
echo "[$(date +%T)] done" >> $LOG
