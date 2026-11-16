#!/usr/bin/env bash
# Exact gauges for the two sigma_5-semilinear lifts (pi and pi^3) on the normal-form table
# ($F4_ARCHIVE/F4_4_F_part{1,2}.bin), written to OUTDIR/h_pi1_a5.pkl and OUTDIR/h_pi3_a5.pkl.
# usage: run_gauges.sh [OUTDIR [LOGDIR]]      defaults: data, ../../logs/f4-forms
# Two processes, about 8.5 minutes and 0.65 GB each.
set -e
cd "$(dirname "$0")"
OUT=${1:-data}; LOGS=${2:-../../logs/f4-forms}
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
PLAN=../fields/data/F4_4_canon.plan.pkl
PI1="0001:1000,1000:2000,2000:0004,0004:0001,0010:0100,0100:0020,0020:0003,0003:0010,1001:1002,1002:1010,1010:0012,0012:1001"
PI3="0001:0004,0004:2000,2000:1000,1000:0001,0010:0003,0003:0020,0020:0100,0100:0010,1001:0012,0012:1010,1010:1002,1002:1001"
mkdir -p "$OUT" "$LOGS"
python3 relabel_exact_nf.py $PLAN 5 $PI1 14 11 "$OUT/h_pi1_a5.pkl" > "$LOGS/relabel_exact_pi1_a5.log" 2>&1 &
python3 relabel_exact_nf.py $PLAN 5 $PI3 14 11 "$OUT/h_pi3_a5.pkl" > "$LOGS/relabel_exact_pi3_a5.log" 2>&1 &
wait
