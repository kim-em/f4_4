# An explicit presentation of the modular tensor category C(F4,4)

This repository contains the paper "An explicit presentation of the modular tensor category
C(F_4,4)" (`paper/`), the data deposited with it (`deposit/`), the programs that construct and
verify the data (`code/`) and the logs of the runs reported in the paper (`logs/`). The complete
F-table (13 GB) is in a separate data archive, **[DATA ARCHIVE DOI: to be supplied]**
(`archive/README.md`). Public URL: **[REPOSITORY URL: to be supplied]**.

The data are the F-symbols of C(F4,4) over `k = Q(zeta_13)^+` and its R-symbols over `Q(zeta_13)`,
in a normal-form gauge. The paper proves from the data files alone that they satisfy the pentagon
and hexagon equations, and that the category they define is a modular tensor category with the
dimensions, twists and S-matrix of C(F4,4).

Code, data and paper source are released under the Apache License 2.0 (`LICENSE`).

## Layout

    paper/          the paper (f4-fsymbols.tex, f4-fsymbols.bib; `make -C paper` builds paper/build/f4-fsymbols.pdf)
    deposit/        the files deposited with the paper: F with third leg V, R, invariants, logs, README with checksums
    archive/        README.md: the archive (two binary files of the F-table, the R-file), checksums; put the files here
    code/f4/        construction modulo p from the quantum group; the exact verifiers
    code/fields/    normal form, the pipeline at one (prime, root) and at all of them, exact reconstruction
    code/checks/    the exact checks of the paper
    code/autoeq/    autoequivalences and braidings
    code/unitary/   invariant forms, positivity, square classes
    code/f4-forms/  forms over the cubic field (Python and PARI/GP)
    code/metaplectic/  C(so(13),2) against the closed formulas of Ardonne, Finch and Titsworth
    logs/           logs of the runs, by directory; logs/clean-clone/ the commands below rerun from a fresh clone

`code/README.md` lists every program.

## Requirements

Python 3 with numpy, scipy, sympy and mpmath (`requirements.txt`; the versions used were Python
3.14.8, numpy 2.5.3, scipy 1.18.1, sympy 1.14.0, mpmath 1.4.1), and PARI/GP (2.17.3) for the `.gp`
scripts of `code/f4-forms`. With nix:

    nix-shell -p "python314.withPackages (ps: with ps; [numpy scipy sympy mpmath])" pari

Set `OPENBLAS_NUM_THREADS=1` when running several processes. All commands below are run from the
top directory of the repository. The archive files are read from `$F4_ARCHIVE` (default
`archive/`); scratch files go to `work/` (or `$F4_WORK`).

## Verifying the data in minutes

Without the archive (about 30 s in total):

    python3 code/checks/s_nondegenerate.py deposit/F4_4_invariants.txt        # S~^2 = D^2 Id
    python3 code/checks/r_invariants.py deposit/F4_4_R_normalform.txt.xz     # invariants of R
    python3 code/checks/autoeq_checks.py deposit/F4_4_R_normalform.txt.xz    # braidings, Mueger centre
    python3 code/checks/nf_galois.py deposit/F4_4_H_normalform.txt.xz        # minimality of k

With the archive (a few minutes, except the last two):

    A=${F4_ARCHIVE:-archive}; R=$A/F4_4_R_normalform.txt
    B="$A/F4_4_F_part1.bin $A/F4_4_F_part2.bin"
    python3 code/checks/exact_checks.py deposit/F4_4_R_normalform.txt.xz $B work/inv.txt   # 50 s, 0.5 GB
    cmp work/inv.txt deposit/F4_4_invariants.txt
    python3 code/checks/nf_witness.py $B                                    # 3 s
    python3 code/checks/pseudo_unitary.py $R $B                             # 30 s, 0.7 GB
    python3 code/f4/pent_exact_subset.py 0000,0001,1000,0004,2000 $B --jobs 2   # 85 s
    python3 code/f4/verify_hex_bin.py $R $B --jobs 2                        # 2.4 h: both hexagons
    python3 code/f4/pent_reduced.py $B --jobs 2                             # 11 h: the pentagon

The verifiers end with a line starting `CONCLUSION:`. Partial runs of the last two:
`verify_hex_bin.py ... --limit-primes 1` (all four-tuples, first prime only) and
`pent_reduced.py ... --maxdim 100 --limit-primes 1`.

## Rerunning one prime

    python3 code/fields/direct_f4.py 1073741189 1 work/direct --nproc 4    # about 5 min, 1.1 GB
    python3 code/fields/compare_direct.py work/direct 1073741189 1          # 25 s

The first command constructs the whole F-table and R modulo `p = 1073741189` at `q = z` from the
quantum group and brings it into the normal form; the second reduces the archived exact table at the
same root and compares: all 54,881,076 entries agree.

## Regenerating everything

    code/fields/direct_f4_all.sh work/all 2     # 78 (prime, root) pairs, about 9 h; then F and R over the fields
    code/fields/finish_all.sh work/all          # the long verifications; comparison with the archive (about 20 h)

The reconstructed files are identical byte for byte to the archived ones
(`deposit/logs/finish_all.log`). The other computations of the paper (autoequivalences, unitarity,
forms over the cubic field, the so(13)_2 comparison) are listed in `code/README.md`.

## Building the paper

    make -C paper          # latexmk; needs a TeX distribution with amsart, booktabs, enumitem, hyperref
