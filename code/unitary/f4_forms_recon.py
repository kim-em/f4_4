#!/usr/bin/env python3
"""Exact reconstruction of the forms G (residues from f4_forms.py) over k = Q(zeta_13)^+.

usage: f4_forms_recon.py DIR NUSE OUT.pkl p1 p2 ...   (first NUSE primes used, the rest held out)
Residues at the 6 roots w = 1,3,5,7,9,11 per prime (one per embedding of k); the partner roots 25w
are compared where present (tau invariance).
"""
import os, pickle, sys, time
sys.dont_write_bytecode = True
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "fields"))
sys.path.insert(1, os.path.join(HERE, "..", "f4"))
from recon import coeffs, ratrec, partner, KREP

d, nuse, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
primes = [int(x) for x in sys.argv[4:]]
use, held = primes[:nuse], primes[nuse:]
ks = [4 * j for j in range(1, 7)]                     # eta_j = c_{4j}
t0 = time.time()
data = {p: {w: np.load(os.path.join(d, f"G_{p}_{w}.npz"))["G"].astype(object) for w in KREP} for p in primes}
tau = []
for p in primes:
    for w in KREP:
        f = os.path.join(d, f"G_{p}_{partner(w)}.npz")
        if os.path.exists(f):
            tau.append(bool((np.load(f)["G"] == data[p][w]).all()))
print(f"tau invariance (G at w equals G at 25w): {sum(tau)}/{len(tau)} checks pass")
cs = {p: coeffs(data[p], p, KREP, ks) for p in primes}
n = cs[primes[0]].shape[0]
M = 1
acc = np.zeros(cs[primes[0]].shape, dtype=object)
for p in use:
    inv = pow(M % p, -1, p)
    acc = acc + M * (((cs[p] - acc) % p) * inv % p)
    M *= p
acc = acc % M
rec, fails, maxbits = [], 0, 0
for i in range(n):
    row = []
    for j in range(6):
        a = int(acc[i, j])
        r = (0, 1) if a == 0 else ratrec(a, M)
        if r is None:
            fails += 1
        else:
            maxbits = max(maxbits, abs(r[0]).bit_length(), r[1].bit_length())
        row.append(r)
    rec.append(row)
bad = 0
if not fails:
    for p in held:
        for i in range(n):
            for j in range(6):
                a, b = rec[i][j]
                if (a - b * int(cs[p][i, j])) % p:
                    bad += 1
print(f"G: {n} entries x 6 eta coefficients; modulus {M.bit_length()} bits ({len(use)} primes); "
      f"failures {fails}; largest num/den {maxbits} bits; mismatches at {len(held)} held-out primes: {bad} ({time.time() - t0:.0f} s)")
if not fails:
    lay = pickle.load(open(os.path.join(d, "layout.pkl"), "rb"))
    pickle.dump(dict(layout=lay, G_eta=rec, used=use, held=held), open(out, "wb"))
    print("wrote", out)
