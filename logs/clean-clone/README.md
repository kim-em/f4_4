# Commands rerun from a fresh clone

The repository was cloned into a new empty directory and the commands below were run there by
`run.sh` (one log per command, named as in the table), in

    nix-shell --pure -p "python314.withPackages (ps: with ps; [numpy scipy sympy mpmath])" pari time

(Python 3.14.8, numpy 2.5.3, scipy 1.18.1, sympy 1.14.0, mpmath 1.4.1, PARI/GP 2.17.3, no other
Python packages, `PYTHONPATH` unset; `environment.log`), with the archive files given only through
`F4_ARCHIVE` and at most 4 processes. Each log starts with the command (archive paths written as
`$F4_ARCHIVE`) and ends with the wall time and peak memory measured by GNU time and the exit status.
All exit statuses are 0. The paper was built with `make -C paper` (`paper_build.log`: 33 pages, no
errors, undefined references or overfull boxes).

| log | command | result | time, memory |
|---|---|---|---|
| `a1_s_nondegenerate` | `code/checks/s_nondegenerate.py` | `S~^2 = D^2 Id`: CONCLUSION as in the paper | 1.5 s, 31 MB |
| `a2_r_invariants` | `code/checks/r_invariants.py` | as `deposit/logs/reproduction_normal_form.log` | 12 s, 60 MB |
| `a3_autoeq_checks` | `code/checks/autoeq_checks.py` | as above | 0.8 s, 61 MB |
| `a4_nf_galois` | `code/checks/nf_galois.py` | stabiliser {1, 25, 27, 51}; identical to `logs/checks/nf_galois.log` | 14 s, 0.30 GB |
| `a5_exact_checks` | `code/checks/exact_checks.py` | all checks True; invariants file identical to the deposited one | 50 s, 0.53 GB |
| `a6_nf_witness` | `code/checks/nf_witness.py` | identical to `logs/checks/nf_witness.log` | 3 s, 0.44 GB |
| `a7_pseudo_unitary` | `code/checks/pseudo_unitary.py` | identical to `logs/checks/pseudo_unitary.log` | 23 s, 0.68 GB |
| `a8_pent_exact_subset5` | `code/f4/pent_exact_subset.py 0000,0001,1000,0004,2000 ... --jobs 2` | CONCLUSION: pentagon for all 1585 five-tuples | 95 s, 0.86 GB |
| `a9_fsym_k_VVVV` | `code/fields/fsym_k.py ... 0001,0001,0001,0001` | the block of the paper's worked example | 1 s, 0.44 GB |
| `a10_nf_unit_check` | `code/fields/f4_nf_checks.py unit` | 6859 unit-leg blocks, 0 entries differ from the identity | 2 s, 0.44 GB |
| `b1_direct_f4` | `code/fields/direct_f4.py 1073741189 1 work/direct --nproc 4` | plan rule reproduced; 58740 blocks, 54881076 entries | 297 s, 1.10 GB |
| `b2_compare_direct` | `code/fields/compare_direct.py work/direct 1073741189 1` | 54881076 entries, 0 differ from the archived table | 24 s, 0.92 GB |
| `b3_nf_braid_check` | `code/fields/f4_nf_checks.py braid work/direct 1073741189 1` | dimensions, twists, monodromy, S: all True | 2 s, 0.50 GB |
| `b4_basic` | `code/f4/basic.py` | basic data of C(F4,4) | 1 s, 33 MB |
| `b5_cells_check` | `cells.py F4 4 1073740201 3 ...; check_cells.py` | failures: 0 | 6 s, 0.22 GB |
| `r1_verify_hex_first_prime` | `code/f4/verify_hex_bin.py R B --jobs 4 --limit-primes 1` | all 58740 four-tuples at the first prime, 24 roots: 0 failures (partial run) | 237 s, 4.1 GB |
| `r2_pent_reduced_maxdim100_first_prime` | `code/f4/pent_reduced.py B --jobs 4 --maxdim 100 --limit-primes 1` | invertibility of all 58740 F-matrices; 40632 five-tuples at the first prime: 0 failures (partial run) | 204 s, 2.0 GB |
| `o1_fusion_auts` | `code/autoeq/fusion_auts.py F4 4` | identical to `logs/autoeq/fusion_auts_F4_4.log` | 1 s, 31 MB |
| `o2_stab_check` | `code/autoeq/stab_check.py` | identical to `logs/autoeq/stab_check.log` | 0.02 s, 15 MB |
| `o3_braid_transport` | `code/autoeq/braid_transport_f4.py 1073741189` | same statements as `logs/autoeq/braid_transport_f4_1073741189.log` (automorphisms in another order) | 16 s, 0.82 GB |
| `o4_unitary_run_forms` | `code/unitary/run_forms.sh work/unitary 2` | identical to `logs/unitary/` (forms, reconstruction, positivity, classes, isometry check) | 244 s, 0.49 GB |
| `o5_f4_relabel_exact_pi2` | `code/fields/f4_relabel_exact.py ... 1 PI2 17 14` | CONCLUSION: h . F^sigma_1 = pi.F, as in `logs/checks/nf_relabel_exact.log` | 540 s, 0.63 GB |
| `o6_run_gauges` | `code/f4-forms/run_gauges.sh` (to `work/`) | both gauges identical to `code/f4-forms/data/h_pi{1,3}_a5.pkl`; logs in `run_gauges/` | 494 s, 0.62 GB per process |
| `o7_descent_u` | `descent_u.py` | `u_values.txt` identical; log identical to `logs/f4-forms/descent_u.log` | 17 s, 0.15 GB |
| `o8_brauer` | `brauer.py`, `gp -q brauer.gp` | `brauer.gp` identical; output identical to `logs/f4-forms/brauer.log` | 20 s, 68 MB |
| `o9_endo` | `uvals_gp.py`, `gp -q endo.gp` | `uvals.gp` identical; output identical to `logs/f4-forms/endo.log` | 0.1 s |
| `o10_kforms` | `gp -q kforms.gp` | identical to `logs/f4-forms/kforms.log` | 21 s, 1.1 GB |
| `o11_fusion_facts` | `fusion_facts.py` | identical to `logs/f4-forms/fusion_facts.log` | 0.1 s |
| `m1_compare_so13` | `code/metaplectic/compare_so13.py` | gauge file identical to `code/metaplectic/data/gauge_so13_ours_to_AFT.txt` | 17 s |
| `m2_compare_so13_check` | `compare_so13.py --check` | CONCLUSION: agreement | 4 s |
| `m3_check_formulas_p123` | `check_formulas.py 1 2 3 --jobs 4` | ALL OK; identical to `logs/metaplectic/check_formulas_p123.log` | 24 s |
| `m4_gauge_covariance` | `gauge_covariance.py` | identical to `logs/metaplectic/gauge_covariance.log` | 6 s |
| `m5_verify_fr_so13` | `code/f4/verify_fr.py code/metaplectic/data/B6_2_FR_exact.txt --jobs 4` | pentagon, hexagons, inverse: hold exactly; identical to `logs/f4/verify_fr_B6_2.log` | 114 s |

"Identical" means identical apart from run times, worker counts and paths.

Not rerun (their logs are in `../fields/`, `../f4-forms/`, `../metaplectic/` and `../../deposit/logs/`):
the complete runs of `verify_hex_bin.py` (8593 s with 2 workers) and `pent_reduced.py` (39,892 s),
`pent_exact_subset.py` on the 7-object subset (1192 s) and with `maxdim=50` (21,721 s), the pipeline at
all 78 (prime, root) pairs with the reconstructions (`direct_f4_all.sh`, about 9 hours) and
`finish_all.sh`, the second exact relabelling (`f4_relabel_exact.py` for `pi^3`, about 9 minutes; the
same computation as `pi^3` in `run_gauges.sh`, which was rerun), `check_formulas.py` for
`p = 4, 5, 6`, and the SageMath runs of `code/metaplectic` (SageMath is not a requirement).
