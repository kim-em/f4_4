#!/usr/bin/env python3
"""Exact reconstruction of the F4 canonical H family (over k = Q(zeta_13)^+) and of the gauge
g (over Q(zeta_52)^+) from the residues written by f4_canon.py run.

usage: recon.py DIR NUSE OUT.pkl [--text F_TEXT]
  uses the first NUSE primes of DIR (sorted as in primes.txt order given by the files present),
  checks every reconstructed value against all remaining primes at all 12 real roots, and checks
  tau-invariance of F_can (residue at w equals residue at 25w) at every prime.
"""
import glob, os, pickle, sys, time
import numpy as np

sys.path.insert(1, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "f4"))
from lin import primitive_root_of_unity

REAL = [1, 3, 5, 7, 9, 11, 15, 17, 19, 21, 23, 25]
KREP = [1, 3, 5, 7, 9, 11]


def partner(w):
    r = 25 * w % 52
    return r if r < 26 else 52 - r


def basis_matrix(p, ws, ks):
    z = primitive_root_of_unity(p, 52, 1)
    zi = pow(z, -1, p)
    return [[(pow(z, k * w, p) + pow(zi, k * w, p)) % p if k else 1 for k in ks] for w in ws]


def inv_mat(A, p):
    n = len(A)
    M = [list(r) + [int(i == j) for j in range(n)] for i, r in enumerate(A)]
    for c in range(n):
        piv = next(r for r in range(c, n) if M[r][c] % p)
        M[c], M[piv] = M[piv], M[c]
        iv = pow(M[c][c], -1, p)
        M[c] = [x * iv % p for x in M[c]]
        for r in range(n):
            if r != c and M[r][c]:
                f = M[r][c]
                M[r] = [(x - f * y) % p for x, y in zip(M[r], M[c])]
    return [r[n:] for r in M]


def coeffs(res, p, ws, ks):
    """res: dict w -> int64 array; returns array (len, len(ks)) of coefficients mod p."""
    Ai = np.array(inv_mat(basis_matrix(p, ws, ks), p), dtype=object)
    X = np.array([res[w] for w in ws], dtype=object)       # (nw, n)
    return (Ai.dot(X) % p).T                                 # (n, nk)


def ratrec(a, m):
    """rational reconstruction of a mod m with |num|, den <= sqrt(m/2); None if none."""
    import math
    bound = math.isqrt(m // 2)
    r0, r1, s0, s1 = m, a % m, 0, 1
    while r1 > bound:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        s0, s1 = s1, s0 - q * s1
    if s1 == 0 or abs(s1) > bound:
        return None
    if s1 < 0:
        r1, s1 = -r1, -s1
    from math import gcd
    if gcd(r1, s1) != 1:
        return None
    return (r1, s1)


def crt_rec(cs_by_p, use):
    """Garner CRT of the residue arrays cs_by_p[p] (p in use) and rational reconstruction of
    every entry: (rows of (num, den), failures, largest bit length, modulus)."""
    M = 1
    acc = np.zeros(cs_by_p[use[0]].shape, dtype=object)
    for p in use:
        c = cs_by_p[p]
        inv = pow(M % p, -1, p)
        acc = acc + M * (((c - acc) % p) * inv % p)
        M *= p
    acc = acc % M
    out = []
    fails = 0
    mb = 0
    for i in range(acc.shape[0]):
        row = []
        for j in range(acc.shape[1]):
            a = int(acc[i, j])
            if a == 0:
                row.append((0, 1)); continue
            r = ratrec(a, M)
            if r is None:
                fails += 1; row.append((0, 1)); continue
            mb = max(mb, abs(r[0]).bit_length(), r[1].bit_length())
            row.append(r)
        out.append(row)
    return out, fails, mb, M


def main(argv):
    d, nuse, out = argv[1], int(argv[2]), argv[3]
    primes = [int(x) for x in open(os.path.join(d, "..", "primes.txt")).read().split()]
    primes = [p for p in primes if os.path.exists(os.path.join(d, f"canon_{p}_1.npz"))]
    use, held = primes[:nuse], primes[nuse:]
    t0 = time.time()
    data = {p: {w: np.load(os.path.join(d, f"canon_{p}_{w}.npz")) for w in REAL} for p in primes}
    # tau invariance at every prime
    tau_ok = all((data[p][w]["F"] == data[p][partner(w)]["F"]).all() for p in primes for w in REAL)
    print(f"tau-invariance of F_can (residue at w = residue at 25w) at all {len(primes)} primes and 12 real roots: {tau_ok}")
    rawdiff = sum(int((data[p][1]["H"] != data[p][25]["H"]).sum()) for p in primes[:1])
    print(f"(for comparison, the input H family differs between w=1 and w=25 in {rawdiff} entries at the first prime)")
    result = {}
    for name, ws, ks in (("F", KREP, [4 * j for j in range(1, 7)]), ("g", REAL, list(range(12)))):
        cs = {p: coeffs({w: data[p][w][name].astype(object) for w in ws}, p, ws, ks) for p in primes}
        n = cs[primes[0]].shape[0]
        M = 1
        acc = np.zeros(cs[primes[0]].shape, dtype=object)
        for p in use:                                         # CRT
            c = cs[p]
            inv = pow(M % p, -1, p)
            acc = acc + M * (((c - acc) % p) * inv % p)
            M *= p
        acc = acc % M
        rec = []
        fails = 0
        maxbits = 0
        for i in range(n):
            row = []
            for j in range(len(ks)):
                a = int(acc[i, j])
                if a == 0:
                    row.append((0, 1)); continue
                r = ratrec(a, M)
                if r is None:
                    fails += 1; row.append(None); continue
                row.append(r)
                maxbits = max(maxbits, abs(r[0]).bit_length(), r[1].bit_length())
            rec.append(row)
        bad = 0
        if not fails:
            for p in held:
                c = cs[p]
                for i in range(n):
                    for j in range(len(ks)):
                        a, b = rec[i][j]
                        if (a - b * int(c[i, j])) % p:
                            bad += 1
        print(f"[{name}] {n} entries x {len(ks)} coefficients; modulus {M.bit_length()} bits from {len(use)} primes; "
              f"reconstruction failures {fails}; largest num/den {maxbits} bits; "
              f"mismatches at {len(held)} held-out primes: {bad}  ({time.time() - t0:.0f} s)")
        result[name] = rec
    layout = pickle.load(open(os.path.join(d, "layout.pkl"), "rb"))
    with open(out, "wb") as f:
        pickle.dump(dict(layout=layout, F_eta=result["F"], g_real=result["g"], primes_used=use, held=held), f)
    print("wrote", out)


if __name__ == "__main__":
    main(sys.argv)
