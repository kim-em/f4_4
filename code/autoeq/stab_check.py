#!/usr/bin/env python3
"""Inputs of the proof that the stabiliser of the F4 associator in the unit-normalised gauge group
is exactly the torus T = {g^{XY}_Z = u_X u_Y / u_Z} (paper, Section "Autoequivalences and braidings").

Reads the normal-form plan of ../fields/f4_canon.py and checks:
 (a) step 1: the D-rows of the normalised entries have rank r and r unit pivots, so the subsystem
     map Z^{edge vars} -> Z^P is onto and its kernel torus K is connected of dimension #vars - r;
 (b) the image of T in the edge variables has the same dimension (rank of the matrix of the
     characters u_A u_V / u_C over Q), so T|edges = K;
 (c) steps 2 and 3 cover every multiplicity space not covered by step 1 (each gets a gauge
     determined by already-fixed data), so the normalising gauge is unique up to T.
usage: stab_check.py
"""
import os, sys, pickle
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
plan = pickle.load(open(os.path.join(HERE, "..", "fields", "data", "F4_4_canon.plan.pkl"), "rb"))
N, obj = plan["N"], plan["obj"]
UNIT, V = "0000", "0001"
ev = plan["edge_vars"]
print(f"(a) step 1: {len(ev)} edge variables, rank {plan['rank']}, unit pivots {len(plan['torus'])}, residual {plan['residual']}")
nonunit = [x for x in obj if x != UNIT]
idx = {x: i for i, x in enumerate(nonunit)}
rows = []
for (A, B, C) in ev:
    r = [Fraction(0)] * len(nonunit)
    for x, s in ((A, 1), (B, 1), (C, -1)):
        if x != UNIT:
            r[idx[x]] += s
    rows.append(r)
def rank(M):
    M = [list(r) for r in M]; rk = 0; ncol = len(M[0])
    for c in range(ncol):
        piv = next((i for i in range(rk, len(M)) if M[i][c] != 0), None)
        if piv is None: continue
        M[rk], M[piv] = M[piv], M[rk]
        for i in range(len(M)):
            if i != rk and M[i][c] != 0:
                f = M[i][c] / M[rk][c]; M[i] = [a - f * b for a, b in zip(M[i], M[rk])]
        rk += 1
    return rk
print(f"(b) dim of T restricted to the edge variables: {rank(rows)} (T has dim {len(nonunit)})")
loops = [t for t in N if t[1] == V and t[0] != UNIT and N[t] > 1]
piv = plan["piv"]
covered = set(ev) | set(loops) | set(piv)
missing = [t for t in N if UNIT not in (t[0], t[1]) and t not in covered]
print(f"(c) multiplicity spaces without unit leg: {sum(1 for t in N if UNIT not in (t[0], t[1]))}; "
      f"step 1: {len(ev)}, step 2 (loops): {loops}, step 3: {len(piv)}; not covered: {missing}")
