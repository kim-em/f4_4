#!/usr/bin/env python3
"""Check  F Gram_R F^T = Gram_L  for every block of the normal-form C(F4,4) table, modulo p at a
root, with the exact forms G (reconstructed from the H family only).  Since the F-symbols outside
the H family were not used to compute G, this tests the whole table against the formula.

usage: f4_isometry_check.py G_exact.pkl p w [JOBS]      table: $F4_ARCHIVE/F4_4_F_part{1,2}.bin
"""
import os, sys, time
sys.dont_write_bytecode = True
import numpy as np
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(1, os.path.join(HERE, "..", "f4"))
from f4_forms_analyse import load
from lin import Fp, primitive_root_of_unity
import fsym_bin

ARCHIVE = os.environ.get("F4_ARCHIVE", os.path.join(HERE, "..", "..", "archive"))
PARTS = [os.path.join(ARCHIVE, f"F4_4_F_part{i}.bin") for i in (1, 2)]


def gram(labels, legs, Gp, p):
    """block-diagonal Gram matrix on labels (Z, i, j) with spaces legs(Z) = (s1, s2)."""
    n = len(labels)
    M = np.zeros((n, n), dtype=np.int64)
    idx = {l: a for a, l in enumerate(labels)}
    for a, (Z, i, j) in enumerate(labels):
        s1, s2 = legs(Z)
        g1, g2 = Gp[s1], Gp[s2]
        for ii in range(g1.shape[0]):
            for jj in range(g2.shape[0]):
                M[a, idx[(Z, ii, jj)]] = g1[i, ii] * g2[j, jj] % p
    return M


def work(args):
    part, keys, p, w, Gp = args
    r = fsym_bin.Reader(part)
    F = Fp(p)
    z = primitive_root_of_unity(p, 52, w); zi = pow(z, -1, p)
    cv = [1] + [(pow(z, k, p) + pow(zi, k, p)) % p for k in range(1, 12)]
    bad = 0
    for key in keys:
        X, Y, W, T = key
        rows, cols, _, _, _ = r.block(key)
        M = r.residues(key, p, cv)
        GL = gram(rows, lambda Z: ((X, Y, Z), (Z, W, T)), Gp, p)
        GR = gram(cols, lambda Z: ((Y, W, Z), (X, Z, T)), Gp, p)
        if not (F.mm(F.mm(M, GR), M.T.copy()) % p == GL % p).all():
            bad += 1
    return len(keys), bad


def main(argv):
    Gx = load(argv[1])
    p, w = int(argv[2]), int(argv[3])
    jobs = int(argv[4]) if len(argv) > 4 else 2
    z = primitive_root_of_unity(p, 52, w); zi = pow(z, -1, p)
    eta = [(pow(z, 4 * j, p) + pow(zi, 4 * j, p)) % p for j in range(1, 7)]
    Gp = {t: np.array([[sum(c.numerator * e * pow(c.denominator, -1, p) for c, e in zip(x, eta)) % p for x in row] for row in M], dtype=np.int64)
          for t, M in Gx.items()}
    tasks = []
    for part in PARTS:
        keys = [tuple(k) for k in fsym_bin.Reader(part).keys()]
        for i in range(0, len(keys), 500):
            tasks.append((part, keys[i:i + 500], p, w, Gp))
    t0 = time.time()
    tot = bad = 0
    with Pool(jobs) as pool:
        for n, b in pool.imap_unordered(work, tasks):
            tot += n; bad += b
    print(f"p={p} w={w}: {tot} blocks checked, F Gram_R F^T = Gram_L fails in {bad}  ({time.time() - t0:.0f} s, {jobs} workers)")


if __name__ == "__main__":
    main(sys.argv)
