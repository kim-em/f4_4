#!/usr/bin/env bash
# Runs the documented commands from a fresh clone; one log per step in logs/clean-clone/.
set -u
cd "$1"
L=logs/clean-clone; mkdir -p $L work
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1
A=$F4_ARCHIVE; R=$A/F4_4_R_normalform.txt; B="$A/F4_4_F_part1.bin $A/F4_4_F_part2.bin"
step() { local name=$1; shift; echo "## $*" > $L/$name.log; ( command time -f "%e s, %M kB (peak RSS)" bash -c "$*" ) >> $L/$name.log 2>&1; echo "exit status $?" >> $L/$name.log; echo "[$(date +%T)] $name: $(tail -2 $L/$name.log | head -1)"; }
case "${2:-all}" in
a|all)
step a1_s_nondegenerate python3 code/checks/s_nondegenerate.py deposit/F4_4_invariants.txt
step a2_r_invariants python3 code/checks/r_invariants.py deposit/F4_4_R_normalform.txt.xz
step a3_autoeq_checks python3 code/checks/autoeq_checks.py deposit/F4_4_R_normalform.txt.xz
step a4_nf_galois python3 code/checks/nf_galois.py deposit/F4_4_H_normalform.txt.xz
step a5_exact_checks "python3 code/checks/exact_checks.py deposit/F4_4_R_normalform.txt.xz $B work/inv.txt && cmp work/inv.txt deposit/F4_4_invariants.txt && echo 'work/inv.txt identical to deposit/F4_4_invariants.txt'"
step a6_nf_witness python3 code/checks/nf_witness.py $B
step a7_pseudo_unitary python3 code/checks/pseudo_unitary.py $R $B
step a8_pent_exact_subset5 python3 code/f4/pent_exact_subset.py 0000,0001,1000,0004,2000 $B --jobs 2
step a9_fsym_k_VVVV python3 code/fields/fsym_k.py $B 0001,0001,0001,0001
step a10_nf_unit_check python3 code/fields/f4_nf_checks.py unit $B
;;&
b|all)
step b1_direct_f4 python3 code/fields/direct_f4.py 1073741189 1 work/direct --nproc 4
step b2_compare_direct python3 code/fields/compare_direct.py work/direct 1073741189 1
step b3_nf_braid_check python3 code/fields/f4_nf_checks.py braid work/direct 1073741189 1
step b4_basic python3 code/f4/basic.py
step b5_cells_check "cd code/f4 && python3 cells.py F4 4 1073740201 3 ../../work/c.pkl 0,0,0,1 && python3 check_cells.py ../../work/c.pkl"
rm -rf work/direct work/c.pkl
;;&
r|all)
step r1_verify_hex_first_prime python3 code/f4/verify_hex_bin.py $R $B --jobs 4 --limit-primes 1
step r2_pent_reduced_maxdim100_first_prime python3 code/f4/pent_reduced.py $B --jobs 4 --maxdim 100 --limit-primes 1
;;&
o|all)
step o1_fusion_auts python3 code/autoeq/fusion_auts.py F4 4
step o2_stab_check python3 code/autoeq/stab_check.py
step o3_braid_transport python3 code/autoeq/braid_transport_f4.py 1073741189
step o4_unitary_run_forms code/unitary/run_forms.sh work/unitary 2
I2="0001:2000,2000:0001,1000:0004,0004:1000,0010:0020,0020:0010,0003:0100,0100:0003,1001:1010,1010:1001,0012:1002,1002:0012"
step o5_f4_relabel_exact_pi2 "cd code/fields && python3 f4_relabel_exact.py data/F4_4_canon.plan.pkl 1 $I2 17 14"
step o6_run_gauges "code/f4-forms/run_gauges.sh $PWD/work/f4forms $PWD/$L/run_gauges && cmp work/f4forms/h_pi1_a5.pkl code/f4-forms/data/h_pi1_a5.pkl && cmp work/f4forms/h_pi3_a5.pkl code/f4-forms/data/h_pi3_a5.pkl && echo 'gauges identical to code/f4-forms/data'"
step o7_descent_u "cd code/f4-forms && python3 descent_u.py data/h_pi1_a5.pkl data/h_pi3_a5.pkl ../../deposit/F4_4_R_normalform.txt.xz ../../work/f4forms && cmp ../../work/f4forms/u_values.txt data/u_values.txt && echo 'u_values.txt identical to data/u_values.txt'"
step o8_brauer "cd code/f4-forms && python3 brauer.py data/u_values.pkl ../../work/f4forms/brauer.gp && cmp ../../work/f4forms/brauer.gp brauer.gp && echo 'brauer.gp identical' && gp -q brauer.gp"
step o9_endo "cd code/f4-forms && python3 uvals_gp.py data/u_values.pkl ../../work/f4forms/uvals.gp && cmp ../../work/f4forms/uvals.gp uvals.gp && echo 'uvals.gp identical' && gp -q endo.gp"
step o10_kforms "cd code/f4-forms && gp -q kforms.gp"
step o11_fusion_facts "cd code/f4-forms && python3 fusion_facts.py ../../deposit/F4_4_R_normalform.txt.xz"
;;&
m|all)
step m1_compare_so13 "cd code/metaplectic && python3 compare_so13.py ../../work/gauge.txt && cmp ../../work/gauge.txt data/gauge_so13_ours_to_AFT.txt && echo 'gauge identical to data/gauge_so13_ours_to_AFT.txt'"
step m2_compare_so13_check "cd code/metaplectic && python3 compare_so13.py --check data/gauge_so13_ours_to_AFT.txt"
step m3_check_formulas_p123 "cd code/metaplectic && python3 check_formulas.py 1 2 3 --jobs 4"
step m4_gauge_covariance "cd code/metaplectic && python3 gauge_covariance.py"
step m5_verify_fr_so13 python3 code/f4/verify_fr.py code/metaplectic/data/B6_2_FR_exact.txt --jobs 4
;;
esac
