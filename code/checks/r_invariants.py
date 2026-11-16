#!/usr/bin/env python3
"""Gauge-invariant information about the field of the R-symbols of C(F4,4).

For X = Y the matrix R[(X,X,Z)] is an endomorphism of Hom(X (x) X, Z), so its characteristic
polynomial is gauge invariant.  We compute it exactly (Faddeev-LeVerrier) for all 100 + ... blocks
and report whether its coefficients lie in Q(zeta^2) = Q(zeta_13) (only even powers of zeta in the
power basis, since Phi_52(x) = Phi_13(-x^2)), and in which cases the one-dimensional invariants
R[(X,X,Z)] are q^e with e odd.  Also reports whether the R-symbols in the exported gauge involve
odd powers of zeta.

usage: r_invariants.py R.txt
"""
import sys, os
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from exact_checks import read_R, mat, matmul, trace, add, scal, mul, ZERO, ONE, NB, ZP, is_zero


def charpoly(A):
    """Faddeev-LeVerrier: coefficients c_0..c_n of det(x - A) = sum c_k x^k (c_n = 1)."""
    n = len(A)
    I = [[ONE if i == j else ZERO for j in range(n)] for i in range(n)]
    M = [[ZERO] * n for _ in range(n)]
    c = [None] * (n + 1)
    c[n] = ONE
    for k in range(1, n + 1):
        # M_k = A M_{k-1} + c_{n-k+1} I
        AM = matmul(A, M) if k > 1 else [[ZERO] * n for _ in range(n)]
        M = [[add(AM[i][j], c[n - k + 1] if i == j else ZERO) for j in range(n)] for i in range(n)]
        t = trace(matmul(A, M))
        c[n - k] = scal(Fr(-1, k), t)
    return c


def even_only(x):
    return all(x[j] == 0 for j in range(1, NB, 2))


obj, unit, N, R, Ri = read_R(sys.argv[1])
n = lambda a, b, c: N.get((a, b, c), 0)
odd_present = sum(1 for key in R for v in R[key].values() if not even_only(v))
print(f"R entries (nonzero) with an odd power of zeta in the exported gauge: {odd_present} of "
      f"{sum(len(v) for v in R.values())}")
one_dim, one_odd, multi, multi_even = 0, [], 0, 0
for X in obj:
    for Z in obj:
        m = n(X, X, Z)
        if not m:
            continue
        A = mat(R[(X, X, Z)], m)
        if m == 1:
            one_dim += 1
            if not even_only(A[0][0]):
                one_odd.append((X, Z))
        else:
            multi += 1
            cp = charpoly(A)
            if all(even_only(cc) for cc in cp):
                multi_even += 1
            else:
                print("  characteristic polynomial with odd powers:", X, Z, m)
print(f"one-dimensional R[(X,X,Z)]: {one_dim}; with odd powers of zeta: {len(one_odd)} {one_odd}")
print(f"R[(X,X,Z)] of size >= 2: {multi}; characteristic polynomial in Q(zeta^2) = Q(zeta_13): {multi_even} of {multi}")
