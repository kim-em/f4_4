#!/usr/bin/env python3
"""Galois action on the normal-form F-symbols with third leg V (eta basis of k = Q(zeta_13)^+).

sigma_a (a in (Z/52)^*) acts on k through zeta_13 -> zeta_13^{a mod 13}: eta_j -> eta_{a j mod 13}
(eta_{13-j} = eta_j).  Reports, for each a, whether some entry of the file is moved.

usage: nf_galois.py F4_4_H_normalform.txt.xz
"""
import sys, lzma, math
from fractions import Fraction as Fr
vals = []
for line in lzma.open(sys.argv[1], "rt"):
    if line.startswith("#") or not line.strip():
        continue
    vals.append([Fr(x) for x in line.split("|")[-1].split()])
def sig(a, v):
    out = [Fr(0)] * 6
    for j, c in enumerate(v, start=1):
        m = (a * j) % 13
        m = min(m, 13 - m)
        out[m - 1] += c
    return out
print(f"{len(vals)} nonzero entries")
for a in [a for a in range(1, 52) if math.gcd(a, 52) == 1]:
    moved = any(sig(a, v) != v for v in vals)
    print(f"sigma_{a}: moves some entry: {moved}")
