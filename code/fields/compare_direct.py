#!/usr/bin/env python3
"""Compare the residues of one run of direct_f4.py with the archived normal-form F-table reduced at the
same (prime, root).

Reads WORK/F_<p>_<w>.npy and WORK/F_<p>_<w>.npy.keys.pkl (direct_f4.py, blocks in sorted key order,
each block row-major) and the binary table F1.bin F2.bin (default $F4_ARCHIVE/F4_4_F_part{1,2}.bin),
reduces every block of the table modulo p at zeta -> z^w (fsym_bin.Reader.residues) and counts the
entries that differ.

usage: compare_direct.py WORK p w [F1.bin F2.bin]
"""
import os, pickle, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(1, os.path.join(HERE, "..", "f4"))
import fsym_bin
import f4_canon as fc


def main(argv):
    work, p, w = argv[1], int(argv[2]), int(argv[3])
    parts = argv[4:6] if len(argv) > 5 else fc.PARTS
    A = np.load(os.path.join(work, f"F_{p}_{w}.npy"), mmap_mode="r")
    keys = pickle.load(open(os.path.join(work, f"F_{p}_{w}.npy.keys.pkl"), "rb"))
    rs = [fsym_bin.Reader(x) for x in parts]
    where = {k: r for r in rs for k in r.keys()}
    cv = fc.cvals(p, w)
    off, diff, tot = 0, 0, 0
    for ent in keys:
        key = tuple(ent)
        M = where[key].residues(key, p, cv) % p
        n = M.size
        B = np.asarray(A[off:off + n]).reshape(M.shape) % p
        diff += int((B != M).sum()); tot += n; off += n
    assert off == A.size and len(keys) == len(where)
    print(f"direct residues at (p={p}, w={w}) vs the archived normal-form table reduced there: {tot} entries, {diff} differ")
    return 1 if diff else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
