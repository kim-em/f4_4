#!/usr/bin/env python3
"""Exact verification of both hexagon equations and of R Rinv = 1 for C(F4,4), from the exported
data: F in the binary format of fsym_bin.py (export convention, real basis) and R, Rinv in the
text format of the R-file F4_4_R_normalform.txt (zeta basis; format in its header).

The identities are those of verify_fr.py:
  (H1)  P1 F[(Y,X,W,T)] P2 = F[(X,Y,W,T)] P3 F[(Y,W,X,T)]   (P1, P2, P3 built from Rinv)
  (H2)  Q5 F[(X,W,Y,T)] Q4 = F[(W,X,Y,T)] Q6 F[(X,Y,W,T)]   (Q4, Q5, Q6 built from R)
For each four-tuple the identity, multiplied by the denominators of the blocks involved
(the common denominator D of each F-block, and for each R- or Rinv-factor the lcm of the
(per block: lcm of the denominators of the R- or Rinv-blocks assembled into one factor)
becomes an identity in Z[zeta] with coefficients bounded by
   B = n^2 c^2 (H_A1 H_A2 H_A3 s_B1 s_B2 s_B3 + H_B1 H_B2 H_B3 s_A1 s_A2 s_A3),
n the block size, c the growth constant of verify_fr.py, H the coefficient bounds of the scaled
factors in the zeta basis and s their scale factors.  The four-tuple is checked at all 24 roots
of the first k primes p = 1 mod 52 (below 2^31, not dividing any denominator), k least with
product > 2B, so each identity is proved exactly.

usage: verify_hex_bin.py R.txt F1.bin [F2.bin ...] [--work DIR] [--jobs J] [--max-primes K] [--limit-primes K]
  --work DIR: scratch directory for one memory-mapped coefficient file per prime (about 2.6 GB,
  deleted after each prime); default work/verify_hex in the repository, or $F4_WORK/verify_hex.
  --limit-primes K: run only the first K primes (a partial check; the CONCLUSION line then says so).
"""
import sys, os, math, time, itertools, pickle
import multiprocessing as mp
import numpy as np
from fractions import Fraction

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
import fsym_bin
from verify_fr import (Data as _TextData, ZPOW, CONDUCTOR, NB, growth_constant, primes_1_mod_52,
                       primitive_roots_52, assemble, zeros)


def opt(argv, name, default):
    return argv[argv.index(name) + 1] if name in argv else default


argv = sys.argv
rfile = argv[1]
bins = [a for a in argv[2:] if a.endswith(".bin")]
work = opt(argv, "--work", os.path.join(os.environ.get("F4_WORK", os.path.join(here, "..", "..", "work")), "verify_hex"))
jobs = int(opt(argv, "--jobs", 12))
maxp = int(opt(argv, "--max-primes", 60))
os.makedirs(work, exist_ok=True)

t0 = time.time()
RD = _TextData(rfile)                       # R, Rinv (no F lines)
obj, N = RD.obj, RD.N
n = lambda X, Y, Z: N.get((X, Y, Z), 0)
readers = [fsym_bin.Reader(b) for b in bins]
where = {}
for r in readers:
    for k in r.keys():
        where[k] = r
blocks = RD.blocks                          # label orders (A rows, B columns) from the fusion rules
assert set(where) == set(blocks), "binary files do not contain exactly the F-blocks"

# conversion real basis -> zeta basis: growth factor of coefficient sizes
T = [[1 if j == 0 else 0 for j in range(NB)]] + [[ZPOW[k][j] + ZPOW[CONDUCTOR - k][j] for j in range(NB)] for k in range(1, 12)]
CONV = max(sum(abs(T[k][j]) for k in range(12)) for j in range(NB))
c = growth_constant()

# block metadata: D and a bound for the scaled coefficients in the zeta basis
meta = {}
for key, r in where.items():
    rows, cols, D, signs, limbs = r.block(key)
    assert [tuple(x) for x in rows] == [(a, i, j) for a, i, j in blocks[key][0]] and \
           [tuple(x) for x in cols] == [(a, i, j) for a, i, j in blocks[key][1]], ("label order differs", key)
    nl = limbs.shape[2]
    top = int(limbs[:, :, nl - 1].max()) if limbs.size else 0
    Hb = (top + 1) << (32 * (nl - 1))
    meta[key] = (D, CONV * Hb)
# per-block denominators and scaled heights of R and Rinv (zeta basis)
def rblocks(src):
    out = {}
    for (key, ij), cc in src.items():
        out.setdefault(key, []).append(cc)
    res = {}
    for key, lst in out.items():
        D = 1
        for cc in lst:
            for x in cc:
                D = math.lcm(D, x.denominator)
        res[key] = (D, max(abs(x.numerator) * (D // x.denominator) for cc in lst for x in cc))
    return res
RB, RIB = rblocks(RD.R), rblocks(RD.Ri)


def rmat(kind, keys):
    """scale and height of a block-diagonal matrix assembled from R (or Rinv) blocks."""
    tab = RB if kind == "R" else RIB
    keys = [k for k in keys if k in tab]
    D = 1
    for k in keys:
        D = math.lcm(D, tab[k][0])
    return D, max((tab[k][1] * (D // tab[k][0]) for k in keys), default=0)


print(f"read R data and {len(where)} F-blocks from {len(bins)} binary files ({time.time()-t0:.0f}s); "
      f"largest R/Rinv block denominator {max(v[0].bit_length() for v in list(RB.values()) + list(RIB.values()))} bits, "
      f"largest F block denominator {max(v[0].bit_length() for v in meta.values())} bits, c = {c}, conversion factor {CONV}",
      flush=True)


def bound(A, B, m):
    """A, B: lists of (scale, height) of the three factors on each side."""
    return m * m * c * c * (math.prod(h for _, h in A) * math.prod(s for s, _ in B) +
                            math.prod(h for _, h in B) * math.prod(s for s, _ in A))


inst = []
for X, Y, W, Tt in itertools.product(obj, repeat=4):
    if (X, Y, W, Tt) not in blocks:
        continue
    m = len(blocks[(X, Y, W, Tt)][0])
    O = obj
    P1 = rmat("Rinv", [(X, Y, Z) for Z in O]); P2 = rmat("Rinv", [(X, W, Z) for Z in O])
    P3 = rmat("Rinv", [(X, Z, Tt) for Z in O if n(Y, W, Z)])
    Q5 = rmat("R", [(X, W, Z) for Z in O]); Q4 = rmat("R", [(Y, W, Z) for Z in O])
    Q6 = rmat("R", [(Z, W, Tt) for Z in O if n(X, Y, Z)])
    Fb = lambda k: meta[k]
    b1 = bound([P1, Fb((Y, X, W, Tt)), P2], [Fb((X, Y, W, Tt)), P3, Fb((Y, W, X, Tt))], m)
    b2 = bound([Q5, Fb((X, W, Y, Tt)), Q4], [Fb((W, X, Y, Tt)), Q6, Fb((X, Y, W, Tt))], m)
    inst.append(((X, Y, W, Tt), 2 * max(b1, b2)))
dens = [meta[k][0] for k in meta] + [v[0] for v in RB.values()] + [v[0] for v in RIB.values()]
primes = []
for p in primes_1_mod_52(4 * maxp):
    if all(d % p for d in dens):
        primes.append(p)
    if len(primes) == maxp:
        break
prods = [1]
for p in primes:
    prods.append(prods[-1] * p)
need = {}
for key, b2x in inst:
    k = next((k for k in range(1, len(prods)) if prods[k] > b2x), None)
    assert k is not None, ("not enough primes", key, b2x.bit_length())
    need[key] = k
kmax = max(need.values())
print(f"{len(inst)} four-tuples; 2 x bound up to {max(b for _, b in inst).bit_length()} bits; "
      f"primes needed per four-tuple: max {kmax}, histogram "
      f"{ {k: list(need.values()).count(k) for k in sorted(set(need.values()))} }", flush=True)
full = int(opt(argv, "--limit-primes", kmax)) >= kmax
kmax = min(kmax, int(opt(argv, "--limit-primes", kmax)))

order = sorted(where)
offs = {}
o = 0
for k in order:
    offs[k] = o
    o += len(blocks[k][0]) * len(blocks[k][1])
TOTAL = o


_LOCAL = {}


def local_reader(key):
    """Readers must not be shared between forked processes (shared file offsets)."""
    pid = os.getpid()
    if _LOCAL.get("pid") != pid:
        _LOCAL.clear()
        _LOCAL["pid"] = pid
        _LOCAL["r"] = {id(r): fsym_bin.Reader(r.f.name) for r in readers}
    return _LOCAL["r"][id(where[key])]


def decode_job(args):
    p, keys, path = args
    mm = np.memmap(path, dtype=np.int32, mode="r+", shape=(TOTAL, 12))
    base = (1 << 32) % p
    for key in keys:
        rows, cols, D, signs, limbs = local_reader(key).block(key)
        acc = np.zeros(limbs.shape[:2], dtype=np.int64)
        for i in range(limbs.shape[2] - 1, -1, -1):
            acc = (acc * base + (limbs[:, :, i].astype(np.int64) % p)) % p
        acc = acc * (signs.astype(np.int64) % p) % p * pow(D, -1, p) % p
        mm[offs[key]:offs[key] + acc.shape[0]] = acc.astype(np.int32)
    mm.flush()
    return len(keys)


def root_job(args):
    p, w, path, keys = args
    mm = np.memmap(path, dtype=np.int32, mode="r", shape=(TOTAL, 12))
    z = primitive_roots_52(p)[0]
    # primitive_roots_52 lists z^u for units u in increasing order; take the w-th unit power of z
    g = pow(z, 1, p)
    r = pow(g, w, p)
    ri = pow(r, -1, p)
    cv = np.array([1] + [(pow(r, k, p) + pow(ri, k, p)) % p for k in range(1, 12)], dtype=np.int64)
    pw = [pow(r, k, p) for k in range(NB)]
    cache = {}

    def Fm(key):
        if key not in cache:
            a = mm[offs[key]:offs[key] + len(blocks[key][0]) * len(blocks[key][1])].astype(np.int64)
            v = np.zeros(a.shape[0], dtype=np.int64)
            for k in range(12):
                v = (v + a[:, k] * cv[k]) % p
            cache[key] = v.reshape(len(blocks[key][0]), len(blocks[key][1]))
        return cache[key]

    def ev(cc):
        s = 0
        for k, x in enumerate(cc):
            if x:
                s += x.numerator * pow(x.denominator, -1, p) % p * pw[k]
        return s % p
    R, Ri = {}, {}
    for tgt, src in ((R, RD.R), (Ri, RD.Ri)):
        for key in N:
            tgt[key] = np.zeros((N[key], N[key]), dtype=np.int64)
        for (key, (i, j)), cc in src.items():
            tgt[key][i, j] = ev(cc)
    from verify_fr import mm as mmod, P_CURRENT
    P_CURRENT[0] = p
    bad = []
    inv_ok = all(np.array_equal(mmod(R[k], Ri[k], p), np.eye(N[k], dtype=np.int64)) for k in N)
    for (X, Y, W, Tt) in keys:
        A_xyw, B_xyw = blocks[(X, Y, W, Tt)]
        A_yxw, B_yxw = blocks[(Y, X, W, Tt)]
        A_ywx, B_ywx = blocks[(Y, W, X, Tt)]
        A_xwy, B_xwy = blocks[(X, W, Y, Tt)]
        A_wxy, B_wxy = blocks[(W, X, Y, Tt)]
        P1 = assemble(A_xyw, A_yxw, [(((Z, i, j), (Z, i2, j)), Ri[(X, Y, Z)][i][i2]) for (Z, i, j) in A_xyw for i2 in range(n(Y, X, Z))])
        P2 = assemble(B_yxw, B_ywx, [(((Z, k, l), (Z, k2, l)), Ri[(X, W, Z)][k][k2]) for (Z, k, l) in B_yxw for k2 in range(n(W, X, Z))])
        P3 = assemble(B_xyw, A_ywx, [(((Z, k, l), (Z, k, l2)), Ri[(X, Z, Tt)][l][l2]) for (Z, k, l) in B_xyw for l2 in range(n(Z, X, Tt))])
        ok1 = np.array_equal(mmod(mmod(P1, Fm((Y, X, W, Tt)), p), P2, p), mmod(mmod(Fm((X, Y, W, Tt)), P3, p), Fm((Y, W, X, Tt)), p))
        Q5 = assemble(A_wxy, A_xwy, [(((Z, i2, j), (Z, i, j)), R[(X, W, Z)][i2][i]) for (Z, i2, j) in A_wxy for i in range(n(X, W, Z))])
        Q4 = assemble(B_xwy, B_xyw, [(((Z, k2, l), (Z, k, l)), R[(Y, W, Z)][k2][k]) for (Z, k2, l) in B_xwy for k in range(n(Y, W, Z))])
        Q6 = assemble(B_wxy, A_xyw, [(((Z, i, j2), (Z, i, j)), R[(Z, W, Tt)][j2][j]) for (Z, i, j2) in B_wxy for j in range(n(Z, W, Tt))])
        ok2 = np.array_equal(mmod(mmod(Q5, Fm((X, W, Y, Tt)), p), Q4, p), mmod(mmod(Fm((W, X, Y, Tt)), Q6, p), Fm((X, Y, W, Tt)), p))
        if not (ok1 and ok2):
            bad.append((X, Y, W, Tt, ok1, ok2))
    return p, w, inv_ok, bad, len(keys)


if __name__ == "__main__":
    units = [u for u in range(1, CONDUCTOR) if math.gcd(u, CONDUCTOR) == 1]
    failures, invfail = [], 0
    passed_full = {}
    for idx, p in enumerate(primes[:kmax]):
        t1 = time.time()
        keys_needed = [key for key, _ in inst if need[key] > idx]
        path = f"{work}/coef_{p}.int32"
        np.memmap(path, dtype=np.int32, mode="w+", shape=(TOTAL, 12)).flush()
        chunks = [order[i::jobs * 4] for i in range(jobs * 4)]
        with mp.get_context("fork").Pool(jobs) as pool:
            pool.map(decode_job, [(p, ch, path) for ch in chunks])
            res = pool.map(root_job, [(p, u, path, keys_needed) for u in units], chunksize=1)
        os.remove(path)
        nb = sum(len(b) for _, _, _, b, _ in res)
        ninv = sum(not ok for _, _, ok, _, _ in res)
        failures += [b for _, _, _, bb, _ in res for b in bb]
        invfail += ninv
        print(f"prime {idx + 1}/{kmax} p={p}: {len(keys_needed)} four-tuples at 24 roots, failures {nb}, "
              f"R Rinv = 1 failures {ninv} ({time.time()-t1:.0f}s)", flush=True)
    print(f"total time {time.time()-t0:.0f}s")
    if failures or invfail:
        print("CONCLUSION: FAILURES", failures[:5])
        sys.exit(1)
    print(f"CONCLUSION: both hexagon equations hold exactly in K for all {len(inst)} four-tuples, and "
          f"R Rinv = 1 holds modulo {kmax} primes at all roots" if full else
          f"CONCLUSION: no failures, but not all primes were run (first {kmax} only)")
