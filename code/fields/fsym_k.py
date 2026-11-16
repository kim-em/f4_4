#!/usr/bin/env python3
"""Reader and per-block text extractor for the descended C(F4,4) F-table over k = Q(zeta_13)^+.

The files F4_4_F_part{1,2}.bin of the archive are in the binary format of ../f4/fsym_bin.py
(documented there: JSON header with per-block keys, row and column labels and offsets; per block a
common denominator and the numerators of the 12 coefficients in the real basis
1, c_1, ..., c_11 of Q(zeta_52)^+, c_k = zeta^k + zeta^-k).  Every entry lies in k; this module
converts the real-basis coefficients exactly to the integral basis eta_j = zeta13^j + zeta13^-j
(j = 1..6, zeta13 = zeta^4) of k, and asserts that the entry lies in k.
Part 1 holds the blocks with first leg 0011, 0012, 1010, part 2 the others.

As a library:   from fsym_k import KReader;  r = KReader([part1, part2]);  r.block_eta(key)
                -> rows, cols, matrix of 6-tuples of Fractions (eta coefficients)
                r.reader(key) is the underlying fsym_bin.Reader (for residues modulo p etc.).
As a script:    fsym_k.py F1.bin [F2.bin] X,Y,W,T [...]        text lines
                    F X Y W T | Z i j | Z' k l | e1 .. e6     (zero entries omitted)
                fsym_k.py F1.bin [F2.bin] --all > table.txt    (the whole table, large)
                fsym_k.py F1.bin [F2.bin] --real X,Y,W,T       the raw real-basis lines a0 .. a11
"""
import os, sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(1, os.path.join(HERE, "..", "f4"))
import fsym_bin
import cyc

# real-basis coordinates of eta_1..eta_6, and a left inverse on the 6-dim subspace they span
_E = [[Fraction(x) for x in cyc.to_real(e)] for e in cyc.ETA_BASIS]


def real_to_eta(v):
    """v: 12 real-basis coefficients of an element of k -> its 6 eta coefficients (exact)."""
    s = cyc._solve([tuple(row) + (Fraction(0),) * 12 for row in _E], tuple(Fraction(x) for x in v) + (Fraction(0),) * 12)
    assert s is not None, "entry not in Q(zeta_13)^+"
    return s


class KReader:
    def __init__(self, paths):
        self.readers = [fsym_bin.Reader(p) for p in paths]
        self.where = {k: r for r in self.readers for k in r.keys()}

    def keys(self):
        return sorted(self.where)

    def reader(self, key):
        return self.where[tuple(key)]

    def block_eta(self, key):
        rows, cols, M = self.where[tuple(key)].fractions(tuple(key))
        return rows, cols, [[real_to_eta(v) if any(v) else [Fraction(0)] * 6 for v in row] for row in M]


def main(argv):
    bins = [a for a in argv[1:] if a.endswith(".bin")]
    r = KReader(bins)
    raw = "--real" in argv
    keys = r.keys() if "--all" in argv else [tuple(a.split(",")) for a in argv[1:] if not a.endswith(".bin") and not a.startswith("--")]
    fmt = lambda x: str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"
    for key in keys:
        if raw:
            rows, cols, M = r.where[key].fractions(key)
        else:
            rows, cols, M = r.block_eta(key)
        for a, row in zip(rows, M):
            for b, v in zip(cols, row):
                if any(v):
                    print(f"F {' '.join(key)} | {a[0]} {a[1]} {a[2]} | {b[0]} {b[1]} {b[2]} | " + " ".join(fmt(x) for x in v))


if __name__ == "__main__":
    main(sys.argv)
