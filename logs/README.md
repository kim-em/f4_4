# Logs

Logs of the runs reported in the paper, by the directory of `../code/` whose programs produced them.
The deposit (`../deposit/logs/`) contains a selection of the same logs under the names cited in the
paper. Paths of the machine on which the runs were made are written as `ARCHIVE/` (the archive
directory), `WORK/` (a scratch directory) and `OUT/` (an output directory).

| directory | content |
|---|---|
| `fields/` | `finish_all.log` (the long verifications and the byte-for-byte comparison of the direct pipeline with the archive); `verify_hex_descended.log`, `pent_reduced_descended.log`, `pent_exact_*_descended.log` (the verifiers on the archived table); `f4_full_recon.log` (the reconstruction of the archived table), `direct_f4_recon.log`, `direct_R_recon.log` (the reconstruction from the direct pipeline); `f4_nf_braid_checks.log`, `f4_nf_unit_check.log` |
| `checks/` | `nf_exact_checks.log`, `nf_galois.log`, `nf_witness.log`, `nf_relabel_exact.log` (`../code/fields/f4_relabel_exact.py`), `pseudo_unitary.log`; `reproduction_normal_form.log` (one pipeline run and the short checks, from a clean copy); `worked_example_descended.log`, `minpolys.log` (short computations for the paper's examples) |
| `autoeq/` | `fusion_auts_F4_4.log`, `stab_check.log`, `braid_transport_f4_<p>.log` |
| `unitary/` | `f4_forms_runs.log`, `f4_forms_recon.log`, `f4_forms_analyse.log`, `f4_classes.log`, `f4_isometry_check.log` |
| `f4-forms/` | `relabel_exact_pi{1,3}_a5.log` (`run_gauges.sh`), `descent_u.log`, `brauer.log`, `endo.log`, `kforms.log`, `fusion_facts.log` |
| `f4/` | `verify_fr_B6_2.log` (the exact verification of the C(so(13),2) data of `../code/metaplectic/data/`) |
| `metaplectic/` | the C(so(13),2) comparison (see `../code/metaplectic/README.md`) |
| `clean-clone/` | the commands of `../README.md` rerun from a fresh clone of this repository (see its `README.md`) |
