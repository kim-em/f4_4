#!/usr/bin/env python3
"""Exact nondegeneracy of the S-matrix computed from the exported R-symbols.

Reads deposit/F4_4_invariants.txt (written by exact_checks.py from the R-file and the F-table): the elements
D*S~_{XY} (S~_{XY} = sum_Z d_Z Tr(R[(Y,X,Z)] R[(X,Y,Z)])), D*d_X and D, all in Q(zeta_52).
Checks exactly that S~ S~ = (sum_X d_X^2) * Id (all objects are self-dual, so the charge
conjugation matrix is the identity) and that sum_X d_X^2 != 0; hence S~ is invertible.

usage: s_nondegenerate.py deposit/F4_4_invariants.txt
"""
import sys, os
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from exact_checks import mul, add, scal, is_zero, NB

St, Nx, DEN, obj = {}, {}, None, []
for line in open(sys.argv[1]):
    if line.startswith("#"):
        continue
    head, _, coeffs = line.partition("|")
    t = head.split()
    v = [Fr(x) for x in coeffs.split()] if coeffs else None
    if t[0] == "OBJ":
        obj.append(t[1])
    elif t[0] == "DEN":
        DEN = v
    elif t[0] == "DIM":
        Nx[t[1]] = v
    elif t[0] == "S":
        St[(t[1], t[2])] = v
# (D S~)(D S~) = D^2 S~^2 ; want S~^2 = (sum d^2) Id  <=>  St St = (sum Nx^2) Id   (both sides times D^2)
tot = [Fr(0)] * NB
for X in obj:
    tot = add(tot, mul(Nx[X], Nx[X]))
ok = True
for X in obj:
    for Y in obj:
        acc = [Fr(0)] * NB
        for Z in obj:
            acc = add(acc, mul(St[(X, Z)], St[(Z, Y)]))
        want = tot if X == Y else [Fr(0)] * NB
        if not is_zero(add(acc, scal(-1, want))):
            ok = False
print(f"{len(obj)} objects; S~ S~ = (sum_X d_X^2) Id exactly: {ok}; sum_X d_X^2 nonzero: {not is_zero(tot)}")
print("CONCLUSION: the S-matrix computed from the R-symbols is invertible" if ok and not is_zero(tot) else "FAILED")
