# Programs

Run every command from the top directory of the repository unless stated otherwise. The archive
files (`F4_4_F_part1.bin`, `F4_4_F_part2.bin`, `F4_4_R_normalform.txt`) are read from
`$F4_ARCHIVE` (default `archive/`); scratch files go to `work/` (or `$F4_WORK`), which git ignores.
Section numbers refer to the paper (`paper/f4-fsymbols.tex`); conventions are those of its
Section 3. Logs of all runs are in `../logs/<directory>/`.

## `f4/`: the construction modulo p, and the verifiers

| file | content |
|---|---|
| `rootsys.py` | root data, Weyl and quantum dimensions, Freudenthal, Kac-Walton fusion |
| `basic.py` | the basic facts of C(F4,4) (`python3 code/f4/basic.py`; writes `work/basic_F4_4.pkl`) |
| `lin.py` | linear algebra modulo `p` |
| `qg.py` | `U_q(g)` modules over `F_p`, highest weight vectors in `M (x) V`, `U^-` maps |
| `cells.py`, `check_cells.py` | the cells (braiding, cup, cap, split, merge on two-step paths); their relations at every base |
| `pathmodel.py`, `objects.py` | path spaces and generator matrices; spectral paths, idempotents `p_S`, auxiliary bases |
| `hfamily.py`, `run_h.py` | F-symbols with third leg `V` |
| `recursion.py`, `run_split.py` | the pentagon recursion, split by first leg |
| `rmatrix.py`, `run_r.py` | R-symbols |
| `fsym_bin.py` | the binary format of the F-table: writer and reader |
| `verify_hex_bin.py` | exact proof of both hexagons and `R Rinv = 1` (paper, Section 7.2) |
| `pent_reduced.py` | exact proof of the pentagons with fourth object `V` and of invertibility (Section 7.4) |
| `pent_exact_subset.py` | exact proof of the pentagon on subsets of objects, or `maxdim=M` |
| `verify_fr.py` | standalone exact verifier for text data with F, R, Rinv (used for C(so(13),2), `metaplectic/`); arithmetic helpers of the verifiers |
| `pentagon-reduction.md` | the reduction lemma (pentagons with fourth object `V` suffice) and the matrix layout of the pentagon used by all checkers |

The construction scripts are driven by `fields/direct_f4.py`; by hand (`p`, root `z^w`):

    cd code/f4
    python3 cells.py F4 4 1073741189 1 ../../work/c.pkl 0,0,0,1
    python3 check_cells.py ../../work/c.pkl
    python3 run_h.py ../../work/c.pkl ../../work/H.pkl 4

The verifiers take the binary files as arguments, `--jobs J` workers, and (`verify_hex_bin.py`,
`pent_reduced.py`) a scratch directory `--work DIR` (default `work/verify_hex`, `work/pent_reduced`;
one file of 2.6 GB per prime, deleted after use). For partial runs: `pent_exact_subset.py` with a
small subset; `pent_reduced.py --maxdim M --limit-primes K`; `verify_hex_bin.py --limit-primes K`.
A partial run ends with `CONCLUSION: no failures, but not all primes were run`.

## `fields/`: normal form, pipeline, reconstruction

| file | content |
|---|---|
| `f4_canon.py`, `torus.py` | the normal form (Section 4.7) modulo `p`; `f4_canon.PARTS` is the archived table |
| `data/F4_4_canon.plan.pkl` | the stored index choices of the normal form (asserted by `direct_f4.py` at every prime) |
| `direct_f4.py` | the pipeline at one (prime, root): cells, H family, normal form, recursion, export, R |
| `compare_direct.py` | compares one run of `direct_f4.py` with the archived table reduced at the same root |
| `direct_f4_all.sh` | `direct_f4.py` at all 78 (prime, root) pairs, then the two reconstructions |
| `data/P.txt`, `data/Q.txt` | the primes used for the reconstruction and the held-out primes |
| `f4_full_recon.py` | exact reconstruction of F over `k = Q(zeta_13)^+`, written in the binary format |
| `direct_R_recon.py` | exact reconstruction of R, Rinv over `Q(zeta_13)` (header from `data/F4_4_R_header.txt`) |
| `recon.py`, `cyc.py` | CRT, rational reconstruction; exact arithmetic in `Q(zeta_52)` |
| `finish_all.sh` | the long verifications on the archive and the byte-for-byte comparison of the reconstruction with it |
| `f4_nf_checks.py` | unit legs exactly; dimensions, twists, monodromy and S modulo `p` from the direct residues |
| `fsym_k.py` | reader of the table in the basis `eta_j` of `k`; text output of blocks |
| `f4_relabel.py`, `f4_relabel_exact.py` | relabellings by fusion ring automorphisms modulo `p`; exact gauges (Section 8) |
| `fusion_aut.py` | fusion ring automorphisms (library) |

    python3 code/fields/direct_f4.py 1073741189 1 work/direct --nproc 4      # 305 s, 1.1 GB
    python3 code/fields/compare_direct.py work/direct 1073741189 1           # 25 s
    code/fields/direct_f4_all.sh work/all 2                                  # about 9 h; 35 GB of residues
    code/fields/finish_all.sh work/all                                       # about 20 h
    python3 code/fields/f4_nf_checks.py unit archive/F4_4_F_part1.bin archive/F4_4_F_part2.bin

`f4_relabel_exact.py` (the run of `deposit/logs/nf_relabel_exact.log`, about 9 minutes each):

    cd code/fields
    I2="0001:2000,2000:0001,1000:0004,0004:1000,0010:0020,0020:0010,0003:0100,0100:0003,1001:1010,1010:1001,0012:1002,1002:0012"
    P3="0001:0004,0004:2000,2000:1000,1000:0001,0010:0003,0003:0020,0020:0100,0100:0010,1001:0012,0012:1010,1010:1002,1002:1001"
    python3 f4_relabel_exact.py data/F4_4_canon.plan.pkl 1 $I2 17 14
    python3 f4_relabel_exact.py data/F4_4_canon.plan.pkl 5 $P3 17 14

`f4_full_recon.py` and `f4_nf_checks.py braid` read the block layout (keys, row and column labels,
not the entries) from the archived table, so `direct_f4_all.sh` and `finish_all.sh` need the archive.

## `checks/`: the exact checks of the paper

| file | content |
|---|---|
| `exact_checks.py` | objects, fusion graph, twists, S (Kac-Peterson), balancing, unit legs, `F[X,X,X,X]_11 d_X = 1`; writes the invariants file |
| `s_nondegenerate.py` | `S~^2 = D^2 Id` exactly |
| `r_invariants.py`, `nf_galois.py` | invariants of R; Galois action on the F-symbols with third leg `V` |
| `pseudo_unitary.py` | `d_X` is a character; `(F^-1)_11 = F_11` modulo two primes |
| `autoeq_checks.py`, `nf_witness.py` | automorphisms, Mueger centre, transport of `R[X,X,Z]`; monoidal obstructions |

Commands: `../README.md`, "Verifying the data".

## `autoeq/`: autoequivalences and braidings (Section 8)

    python3 code/autoeq/fusion_auts.py F4 4                         # fusion ring automorphisms; 1 s
    python3 code/autoeq/stab_check.py                              # stabiliser of the normal form; 1 s
    python3 code/autoeq/braid_transport_f4.py 1073741189            # transported R modulo p; about 1 min per prime

`braid_transport_f4.py` reads the archived F-table and the R-file (`deposit/F4_4_R_normalform.txt.xz`,
or `$F4_RFILE`). The logs in `../logs/autoeq/` were produced from the tables in the auxiliary bases of
the construction; all quantities printed are gauge invariant.

## `unitary/`: the invariant forms (Section 7.7)

    code/unitary/run_forms.sh [WORK] [JOBS]

computes the forms `G` modulo 15 primes (`f4_forms.py`), reconstructs them over `k`
(`f4_forms_recon.py`), and runs `f4_forms_analyse.py` (positivity at 256 bits), `f4_classes.py`
(square classes at 300 primes, `kprimes.py`) and `f4_isometry_check.py` (the whole table at two
(prime, root) pairs).

## `f4-forms/`: forms over the cubic field (Section 8.5)

Run in `code/f4-forms` (PARI/GP needed):

    ./run_gauges.sh                         # relabel_exact_nf.py for pi and pi^3: data/h_pi{1,3}_a5.pkl; 8.5 min, 2 processes
    python3 descent_u.py data/h_pi1_a5.pkl data/h_pi3_a5.pkl ../../deposit/F4_4_R_normalform.txt.xz data   # 17 s
    python3 brauer.py data/u_values.pkl brauer.gp && gp -q brauer.gp                                         # 21 s
    python3 uvals_gp.py data/u_values.pkl uvals.gp && gp -q endo.gp                                          # 1 s
    gp -q kforms.gp                                                                                          # 23 s
    python3 fusion_facts.py ../../deposit/F4_4_R_normalform.txt.xz

`data/` holds the outputs of the first two commands (the gauges and the scalars `u_X`), so the
later commands run without the archive; `brauer.gp` and `uvals.gp` are outputs of `brauer.py` and
`uvals_gp.py`.

## `metaplectic/`: C(so(13),2) and the closed formulas of Ardonne, Finch and Titsworth

See `metaplectic/README.md`.
