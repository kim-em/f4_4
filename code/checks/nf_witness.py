#!/usr/bin/env python3
"""Exact monoidal obstruction for relabellings, on an F-table (binary fsym_bin format).

For the blocks B1 = F[X,X,X,X] and B2 = F[X,X,X,0002], X in the orbit {0001, 1000, 2000, 0004},
all multiplicity spaces are one dimensional, so every gauge acts by F -> D1 F D2 with D1, D2
diagonal and I_ab = F_ab (F^-1)_ba is an invariant of the gauge class.  For each fusion ring
automorphism pi (as in autoeq_checks.py) and each a in (Z/52)^*, test exactly in Q(zeta_52)
whether I(F[pi V, ...])[pi a, pi b] = sigma_a(I(F[V, ...])[a, b]) for both blocks; this is
necessary for pi.F to be gauge equivalent to sigma_a.F.

usage: nf_witness.py F1.bin F2.bin
"""
import sys, os, math
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from exact_checks import fsym_bin, ZP, NB, COND, mul, add, scal, real_to_zeta, is_zero

ONE = [Fr(1)] + [Fr(0)] * (NB - 1)


def inv(x):
    """inverse in Q(zeta) by solving the 24 x 24 multiplication system."""
    cols = [mul(x, ZP[j]) for j in range(NB)]           # column j = x * zeta^j
    A = [[Fr(cols[j][i]) for j in range(NB)] + [ONE[i]] for i in range(NB)]
    n = NB
    for c in range(n):
        p = next(r for r in range(c, n) if A[r][c] != 0)
        A[c], A[p] = A[p], A[c]
        pv = A[c][c]
        A[c] = [v / pv for v in A[c]]
        for r in range(n):
            if r != c and A[r][c] != 0:
                f = A[r][c]
                A[r] = [u - f * v for u, v in zip(A[r], A[c])]
    return [A[i][n] for i in range(n)]


def matinv(M):
    n = len(M)
    A = [list(row) + [ONE if i == j else [Fr(0)] * NB for j in range(n)] for i, row in enumerate(M)]
    for c in range(n):
        p = next(r for r in range(c, n) if not is_zero(A[r][c]))
        A[c], A[p] = A[p], A[c]
        iv = inv(A[c][c])
        A[c] = [mul(v, iv) for v in A[c]]
        for r in range(n):
            if r != c and not is_zero(A[r][c]):
                f = A[r][c]
                A[r] = [add(u, scal(-1, mul(f, v))) for u, v in zip(A[r], A[c])]
    return [row[n:] for row in A]


def sigma(a, v):
    out = [Fr(0)] * NB
    for j, c in enumerate(v):
        if c:
            out = add(out, scal(c, ZP[(a * j) % COND]))
    return out


readers = [fsym_bin.Reader(b) for b in sys.argv[1:3]]
where = {k: r for r in readers for k in r.keys()}


def invariants(key):
    rows, cols, M = where[key].fractions(key)
    Z = [[real_to_zeta(list(v)) for v in row] for row in M]
    Zi = matinv(Z)
    I = {}
    for a, ra in enumerate(rows):
        for b, cb in enumerate(cols):
            I[(tuple(ra), tuple(cb))] = mul(Z[a][b], Zi[b][a])
    return I


orbit = ["0001", "1000", "2000", "0004"]
pis = {"id": dict(zip(orbit, orbit)),
       "(0001 1000 2000 0004)": dict(zip(orbit, ["1000", "2000", "0004", "0001"])),
       "(0001 0004 2000 1000)": dict(zip(orbit, ["0004", "0001", "1000", "2000"])),
       "(0001 2000)(1000 0004)": dict(zip(orbit, ["2000", "0004", "0001", "1000"]))}
# full permutations (labels of rows/columns), from autoeq_checks.py
full = {
    "id": {},
    "(0001 1000 2000 0004)": {"0001": "1000", "1000": "2000", "2000": "0004", "0004": "0001", "0010": "0100", "0100": "0020",
                              "0020": "0003", "0003": "0010", "1001": "1002", "1002": "1010", "1010": "0012", "0012": "1001"},
    "(0001 0004 2000 1000)": {"0001": "0004", "0004": "2000", "2000": "1000", "1000": "0001", "0010": "0003", "0003": "0020",
                              "0020": "0100", "0100": "0010", "1001": "0012", "0012": "1010", "1010": "1002", "1002": "1001"},
    "(0001 2000)(1000 0004)": {"0001": "2000", "2000": "0001", "1000": "0004", "0004": "1000", "0010": "0020", "0020": "0010",
                               "0003": "0100", "0100": "0003", "1001": "1010", "1010": "1001", "0012": "1002", "1002": "0012"}}
UNITS = [a for a in range(1, COND) if math.gcd(a, COND) == 1]
I = {}
for X in orbit:
    I[(X, X, X, X)] = invariants((X, X, X, X))
    I[(X, X, X, "0002")] = invariants((X, X, X, "0002"))
V = "0001"
for name, p in full.items():
    pm = lambda x: p.get(x, x)
    ok = []
    for a in UNITS:
        good = True
        for key in [(V, V, V, V), (V, V, V, "0002")]:
            pk = tuple(pm(x) for x in key)
            for (ra, cb), val in I[key].items():
                pra = (pm(ra[0]),) + ra[1:]
                pcb = (pm(cb[0]),) + cb[1:]
                if I[pk][(pra, pcb)] != sigma(a, val):
                    good = False
                    break
            if not good:
                break
        if good:
            ok.append(a)
    print(f"pi = {name}: I(F[pi V,...]) = sigma_a(I(F[V,...])) on both blocks exactly for a in {ok}")
