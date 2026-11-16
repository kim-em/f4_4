#!/usr/bin/env python3
"""Pseudo-unitarity of the category defined by the normal-form data (paper, Section "Modularity").

  1. d_X (quantum dimensions, exact in Q(zeta_52) as Nx/D, exact_checks.py) form a character of the
     fusion ring: sum_Z N_{XY}^Z d_Z = d_X d_Y exactly, and d_X > 0 at zeta = exp(2 pi i/52).
  2. For each X the entries F11 = F[X,X,X,X][(1,0,0),(1,0,0)] and G11 = (F^{-1})[(1,0,0),(1,0,0)]:
     F11 d_X = 1 exactly (exact_checks.py); G11 = F11 modulo two primes at all 12 real roots (G11 by
     solving x F = e_(1,0,0) modulo p).  The paper proves G11 = F11 in general (the zigzag argument of
     Proposition 3.3); this is a check.  Hence |X|^2 = 1/(F11 G11) = d_X^2 and the global dimension is
     sum d_X^2 = FPdim.

usage: pseudo_unitary.py R.txt F1.bin F2.bin
"""
import sys, os, math, cmath
from fractions import Fraction as Fr
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from exact_checks import read_R, RootSystem, ZP, NB, COND, mul, add, scal, is_zero, zeta_minus_inv, ONE, fsym_bin

obj, unit, N, R, Ri = read_R(sys.argv[1])
n = lambda a, b, c: N.get((a, b, c), 0)
Rs = RootSystem("F4")
pos = Rs.positive_roots
pair = lambda lam: [int(Rs.ip(tuple(l + 1 for l in lam), a)) for a in pos]
wt = {o: tuple(int(c) for c in o) for o in obj}
Nx = {}
for o in obj:
    v = ONE
    for a in pair(wt[o]):
        v = mul(v, zeta_minus_inv(a))
    Nx[o] = v
D = ONE
for b in pair((0, 0, 0, 0)):
    D = mul(D, zeta_minus_inv(b))
# 1. character: sum_Z N d_Z = d_X d_Y  <=>  D * sum_Z N Nx[Z] = Nx[X] Nx[Y]
ok = all(is_zero(add(mul(D, [sum(n(X, Y, Z) * Nx[Z][j] for Z in obj) for j in range(NB)]), scal(-1, mul(Nx[X], Nx[Y]))))
         for X in obj for Y in obj)
emb = lambda x: sum(complex(c) * cmath.exp(2j * math.pi * j / COND) for j, c in enumerate(x))
dnum = {o: emb(Nx[o]) / emb(D) for o in obj}
print(f"1. sum_Z N_XY^Z d_Z = d_X d_Y exactly for all 256 pairs: {ok}; d_X > 0 at zeta = e^(2 pi i/52): "
      f"{all(abs(v.imag) < 1e-9 and v.real > 0 for v in dnum.values())}; sum d_X^2 = {sum(v.real**2 for v in dnum.values()):.10f}")

# 2. G11 = F11 modulo primes
def inv_row_modp(M, p):
    """solution x of x M = e_0 modulo p (numpy, int64; p < 2^31)."""
    nn = M.shape[0]
    A = np.concatenate([(M.T % p).astype(np.int64), np.eye(nn, 1, dtype=np.int64)], axis=1)
    for c in range(nn):
        r = c + int(np.nonzero(A[c:, c] % p)[0][0])
        if r != c:
            A[[c, r]] = A[[r, c]]
        A[c] = A[c] * pow(int(A[c, c]), -1, p) % p
        f = A[:, c].copy(); f[c] = 0
        nz = np.nonzero(f)[0]
        if len(nz):
            A[nz] = (A[nz] - (f[nz, None] * A[c][None, :]) % p) % p
    return [int(v) for v in A[:, nn]]

readers = [fsym_bin.Reader(b) for b in sys.argv[2:4]]
where = {k: r for r in readers for k in r.keys()}
def cvals(p, w):
    from math import gcd
    g = 2
    while True:
        if all(pow(g, (p - 1) // q, p) != 1 for q in {2, 13} | {f for f in _factors(p - 1)}):
            break
        g += 1
    z = pow(g, (p - 1) // 52, p)
    r = pow(z, w, p); ri = pow(r, -1, p)
    return [1] + [(pow(r, k, p) + pow(ri, k, p)) % p for k in range(1, 12)]
def _factors(m):
    out, d = set(), 2
    while d * d <= m:
        while m % d == 0:
            out.add(d); m //= d
        d += 1
    if m > 1:
        out.add(m)
    return out
primes = [1073741189, 1073740201]
okG = True
for p in primes:
    for w in [1, 3, 5, 7, 9, 11, 15, 17, 19, 21, 23, 25]:
        cv = cvals(p, w)
        for X in obj:
            key = (X, X, X, X)
            rows, cols, _, _, _ = where[key].block(key)
            M = where[key].residues(key, p, cv)
            i0 = rows.index((unit, 0, 0)); j0 = cols.index((unit, 0, 0))
            # (F^{-1})[(1,0,0) as a column label of F, (1,0,0) as a row label of F]: solve x F = e_{j0}? use permutation
            # G = F^{-1}; G[j0, i0] = entry of row j0 of G at column i0; row j0 of G solves x F = e_{j0}
            perm = [j0] + [j for j in range(M.shape[1]) if j != j0]
            Mp = M[:, perm]                     # now column j0 is first
            x = inv_row_modp(Mp.T.T if False else Mp, p)  # first row of Mp^{-1}: x Mp = e_0 -> row j0 of M^{-1} (rows of Mp^{-1} permuted)
            G11 = x[i0] % p
            F11 = int(M[i0, j0]) % p
            if G11 != F11:
                okG = False
print(f"2. (F^-1)[(1,0,0),(1,0,0)] = F[(1,0,0),(1,0,0)] for F = F[X,X,X,X], all 16 X, modulo {primes} at all 12 real roots: {okG}")
