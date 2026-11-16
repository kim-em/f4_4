#!/usr/bin/env python3
"""Exact proof of the pentagon equation for all five-tuples (A, B, C, V, E), V = 0001, for the
exported C(F4,4) F-symbols (binary format of fsym_bin.py, export convention, real basis).
With the reduction lemma of pentagon-reduction.md this proves the pentagon for all five-tuples.

Reads only the exported files; shares no code with the producers (it uses fsym_bin.Reader and
the arithmetic helpers of verify_fr.py).

Equation, layout of pentagon-reduction.md: for simple A, B, C, D, E
   M1 M2 = N1 N2 N3,  M1 = (+)_{P,i1} F[(P,C,D,E)],  M2 = (+)_{R,k} F[(A,B,R,E)],
   N1 = (+)_{Q,i3} F[(A,B,C,Q)],  N2 = (+)_{T,a} F[(A,T,D,E)],  N3 = (+)_{S,n} F[(B,C,D,S)].
Height bound per five-tuple.  The scalar equation in row r and column c involves one block of
each of M1, M2, N1, N3 (determined by r or c) and the blocks of N2.  Multiply it by
D_M1 D_M2 D_N1 L_N2 D_N3 (block denominators; L_N2 the lcm over the N2-blocks): it becomes an
identity in Z[zeta] with coefficients at most
   B = m c H_M1 H_M2 D_N1 L_N2 D_N3 + m^2 c^2 H_N1 H_N2 H_N3 D_M1 D_M2,
with m the dimension, H the coefficient bounds of the scaled blocks in the zeta basis (maxima
over the blocks of each factor), c the growth constant modulo Phi_52.  The five-tuple is
checked (full matrices) at the first k primes p = 1 mod 52 below 2^31 not dividing any
denominator, k least with product > 2B.  Since F is real, checking at the 12 roots
zeta -> z^w, w = 1, 3, ..., 25 (w != 13) covers all 24 (the residue at z^-w of a real element
equals that at z^w).  Also checks that every F-matrix is invertible (modulo one prime).

usage: pent_reduced.py F1.bin [F2.bin ...] [--fourth 0001] [--jobs J] [--maxdim M] [--limit-primes K] [--work DIR]
  --work DIR: scratch directory for one memory-mapped coefficient file per prime (about 2.6 GB,
  deleted after each prime); default work/pent_reduced in the repository, or $F4_WORK/pent_reduced.
"""
import sys, os, math, time, itertools
import multiprocessing as mp
import numpy as np

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
import fsym_bin
from verify_fr import ZPOW, CONDUCTOR, NB, growth_constant, primes_1_mod_52, primitive_roots_52, mm

argv = sys.argv
bins = [a for a in argv[1:] if a.endswith(".bin")]
opt = lambda name, d: argv[argv.index(name) + 1] if name in argv else d
fourth = opt("--fourth", "0001")
jobs = int(opt("--jobs", 24))
maxdim = int(opt("--maxdim", 10 ** 9))
work = opt("--work", os.path.join(os.environ.get("F4_WORK", os.path.join(here, "..", "..", "work")), "pent_reduced"))
os.makedirs(work, exist_ok=True)
t0 = time.time()
readers = [fsym_bin.Reader(b) for b in bins]
where = {k: r for r in readers for k in r.keys()}
objs = sorted({k[0] for k in where})
unit = "0000"
N = {}
for key in where:
    X, Y, W, T = key
    if W == unit:
        N[(X, Y, T)] = len(where[key].index[key]["rows"])
n = lambda a, b, c: N.get((a, b, c), 0)
T = [[1 if j == 0 else 0 for j in range(NB)]] + [[ZPOW[k][j] + ZPOW[CONDUCTOR - k][j] for j in range(NB)] for k in range(1, 12)]
CONV = max(sum(abs(T[k][j]) for k in range(12)) for j in range(NB))
cgr = growth_constant()

# block metadata (denominator, coefficient bound of D*F in the zeta basis), read once
meta = {}
for r in readers:
    for b in r.meta["blocks"]:
        e = len(b["rows"]) * len(b["cols"])
        r.f.seek(r.base + b["offset"])
        D = int.from_bytes(r.f.read(4 * b["dlimbs"]), "little")
        r.f.seek(r.base + b["offset"] + 4 * b["dlimbs"] + 12 * e)
        nl = b["nlimbs"]
        limbs = np.frombuffer(r.f.read(4 * e * 12 * nl), dtype=np.uint32).reshape(e, 12, nl)
        top = int(limbs[:, :, nl - 1].max()) if e else 0
        meta[tuple(b["key"])] = (D, CONV * ((top + 1) << (32 * (nl - 1))))
print(f"read {len(meta)} F-blocks from {len(bins)} files ({time.time()-t0:.0f}s); c = {cgr}, conversion factor {CONV}", flush=True)


def blk_keys(A, B, C, D, E):
    M1 = [(P, C, D, E) for P in objs if n(A, B, P) and (P, C, D, E) in where]
    M2 = [(A, B, R, E) for R in objs if n(C, D, R) and (A, B, R, E) in where]
    N1 = [(A, B, C, Q) for Q in objs if n(Q, D, E) and (A, B, C, Q) in where]
    N2 = [(A, Tt, D, E) for Tt in objs if n(B, C, Tt) and (A, Tt, D, E) in where]
    N3 = [(B, C, D, S) for S in objs if n(A, S, E) and (B, C, D, S) in where]
    return M1, M2, N1, N2, N3


inst = []
for A, B, C, E in itertools.product(objs, repeat=4):
    D = fourth
    m = sum(n(A, B, P) * n(P, C, Q) * n(Q, D, E) for P in objs for Q in objs)
    if not m or m > maxdim:
        continue
    M1, M2, N1, N2, N3 = blk_keys(A, B, C, D, E)
    mx = lambda ks: (max(meta[k][0] for k in ks), max(meta[k][1] for k in ks))
    dM1, hM1 = mx(M1); dM2, hM2 = mx(M2); dN1, hN1 = mx(N1); dN3, hN3 = mx(N3)
    L2 = 1
    for k in N2:
        L2 = math.lcm(L2, meta[k][0])
    hN2 = max(meta[k][1] * (L2 // meta[k][0]) for k in N2)
    b = m * cgr * hM1 * hM2 * dN1 * L2 * dN3 + m * m * cgr * cgr * hN1 * hN2 * hN3 * dM1 * dM2
    inst.append(((A, B, C, D, E), m, 2 * b))
dens = [v[0] for v in meta.values()]
primes = []
for p in primes_1_mod_52(1000):
    if all(d % p for d in dens):
        primes.append(p)
    if len(primes) >= 200:
        break
prods = [1]
for p in primes:
    prods.append(prods[-1] * p)
need = {t: next(k for k in range(1, len(prods)) if prods[k] > b) for t, m, b in inst}
kmax = max(need.values())
kmax = min(kmax, int(opt("--limit-primes", kmax)))
order = sorted(where)
offs, o = {}, 0
for k in order:
    offs[k] = o
    o += len(where[k].index[k]["rows"]) * len(where[k].index[k]["cols"])
TOTAL = o
print(f"five-tuples with fourth object {fourth}: {len(inst)}, largest dimension {max(m for _, m, _ in inst)}, "
      f"sum of dim^3 {sum(m ** 3 for _, m, _ in inst):.3e}; 2 x bound up to {max(b for _, _, b in inst).bit_length()} bits; "
      f"primes needed: max {max(need.values())}, mean {sum(need.values()) / len(need):.2f}, histogram "
      f"{ {k: list(need.values()).count(k) for k in sorted(set(need.values()))} }", flush=True)

_LOC = {}


def rdr(key):
    if _LOC.get("pid") != os.getpid():
        _LOC.clear(); _LOC["pid"] = os.getpid()
        _LOC["r"] = {id(r): fsym_bin.Reader(r.f.name) for r in readers}
    return _LOC["r"][id(where[key])]


def decode_job(args):
    p, keys, path = args
    mmap = np.memmap(path, dtype=np.int32, mode="r+", shape=(TOTAL, 12))
    base = (1 << 32) % p
    for key in keys:
        rows, cols, D, signs, limbs = rdr(key).block(key)
        acc = np.zeros(limbs.shape[:2], dtype=np.int64)
        for i in range(limbs.shape[2] - 1, -1, -1):
            acc = (acc * base + (limbs[:, :, i].astype(np.int64) % p)) % p
        acc = acc * (signs.astype(np.int64) % p) % p * pow(D, -1, p) % p
        mmap[offs[key]:offs[key] + acc.shape[0]] = acc.astype(np.int32)
    mmap.flush()
    return len(keys)


def pent_job(args):
    p, w, path, tuples = args
    mmap = np.memmap(path, dtype=np.int32, mode="r", shape=(TOTAL, 12))
    z = primitive_roots_52(p)[0]
    r = pow(z, w, p); ri = pow(r, -1, p)
    cv = np.array([1] + [(pow(r, k, p) + pow(ri, k, p)) % p for k in range(1, 12)], dtype=np.int64)
    cache = {}

    def F(key):
        if key not in cache:
            b = where[key].index[key]
            nr, nc = len(b["rows"]), len(b["cols"])
            a = mmap[offs[key]:offs[key] + nr * nc].astype(np.int64)
            v = np.zeros(nr * nc, dtype=np.int64)
            for k in range(12):
                v = (v + a[:, k] * cv[k]) % p
            cache[key] = ([tuple(x) for x in b["rows"]], [tuple(x) for x in b["cols"]], v.reshape(nr, nc))
        return cache[key]
    bad, done = [], 0
    for (A, B, C, D, E) in tuples:
        O = objs
        I0 = [(P, i1, Q, i2, i3) for P in O for i1 in range(n(A, B, P)) for Q in O for i2 in range(n(P, C, Q)) for i3 in range(n(Q, D, E))]
        I1 = [(P, i1, R, k, l) for P in O for i1 in range(n(A, B, P)) for R in O for k in range(n(C, D, R)) for l in range(n(P, R, E))]
        I2 = [(R, k, S, m_, nn) for R in O for k in range(n(C, D, R)) for S in O for m_ in range(n(B, R, S)) for nn in range(n(A, S, E))]
        J1 = [(Tt, a, b, Q, i3) for Tt in O for a in range(n(B, C, Tt)) for Q in O for b in range(n(A, Tt, Q)) for i3 in range(n(Q, D, E))]
        J2 = [(Tt, a, S, cc, nn) for Tt in O for a in range(n(B, C, Tt)) for S in O for cc in range(n(Tt, D, S)) for nn in range(n(A, S, E))]
        ix = [{x: i for i, x in enumerate(L)} for L in (I0, I1, I2, J1, J2)]
        m = len(I0)

        def pieces(which):
            out = []
            if which == "M1":
                for P in O:
                    for i in range(n(A, B, P)):
                        key = (P, C, D, E)
                        if key in where:
                            rl, cl, M = F(key)
                            out.append(([ix[0][(P, i) + a] for a in rl], [ix[1][(P, i) + b] for b in cl], M))
            elif which == "M2":
                for R in O:
                    for k in range(n(C, D, R)):
                        key = (A, B, R, E)
                        if key in where:
                            rl, cl, M = F(key)
                            out.append(([ix[1][(a[0], a[1], R, k, a[2])] for a in rl], [ix[2][(R, k) + b] for b in cl], M))
            elif which == "N1":
                for Q in O:
                    for i in range(n(Q, D, E)):
                        key = (A, B, C, Q)
                        if key in where:
                            rl, cl, M = F(key)
                            out.append(([ix[0][(a[0], a[1], Q, a[2], i)] for a in rl], [ix[3][b + (Q, i)] for b in cl], M))
            elif which == "N2":
                for Tt in O:
                    for x in range(n(B, C, Tt)):
                        key = (A, Tt, D, E)
                        if key in where:
                            rl, cl, M = F(key)
                            out.append(([ix[3][(Tt, x, a[1], a[0], a[2])] for a in rl], [ix[4][(Tt, x) + b] for b in cl], M))
            else:
                for S in O:
                    for nn in range(n(A, S, E)):
                        key = (B, C, D, S)
                        if key in where:
                            rl, cl, M = F(key)
                            out.append(([ix[4][(a[0], a[1], S, a[2], nn)] for a in rl], [ix[2][(b[0], b[1], S, b[2], nn)] for b in cl], M))
            return out

        def place(pcs):
            X = np.zeros((m, m), dtype=np.int64)
            for ri, ci, M in pcs:
                X[np.ix_(ri, ci)] = M
            return X

        def times(X, pcs):
            out = np.zeros((m, m), dtype=np.int64)
            for ri, ci, M in pcs:
                sub = X[:, ri]
                out[:, ci] = (out[:, ci] + mm(sub, M, p)) % p
            return out
        lhs = times(place(pieces("M1")), pieces("M2"))
        rhs = times(times(place(pieces("N1")), pieces("N2")), pieces("N3"))
        if not np.array_equal(lhs, rhs):
            bad.append((A, B, C, D, E))
        done += 1
    return p, w, bad, done


if __name__ == "__main__":
    W12 = [w for w in range(1, 26) if w % 2 and w != 13]
    # invertibility of every F-matrix, modulo one prime < 2^30 (lin.Fp), implies invertibility over K
    from lin import Fp
    p0 = next(p for p in primes_1_mod_52(400, below=2 ** 30) if all(d % p for d in dens))
    F0 = Fp(p0)
    z0 = primitive_roots_52(p0)[0]
    cv0 = [1] + [(pow(z0, k, p0) + pow(z0, -k, p0)) % p0 for k in range(1, 12)]
    t1 = time.time()
    inv_ok = True
    for key in order:
        Mx = where[key].residues(key, p0, cv0)
        if F0.rank(Mx) < Mx.shape[0]:
            inv_ok = False
            print("singular modulo p0:", key)
    print(f"all {len(order)} F-matrices invertible modulo {p0} (hence over K): {inv_ok} ({time.time()-t1:.0f}s)", flush=True)
    failures = []
    for idx, p in enumerate(primes[:kmax]):
        t1 = time.time()
        tl = [t for t, m, b in sorted(inst, key=lambda x: -x[1]) if need[t] > idx]   # largest first
        path = os.path.join(work, f"pentred_coef_{p}.int32")
        np.memmap(path, dtype=np.int32, mode="w+", shape=(TOTAL, 12)).flush()
        nsplit = max(1, (2 * jobs) // len(W12))
        with mp.get_context("fork").Pool(jobs) as pool:
            pool.map(decode_job, [(p, order[i::jobs * 4], path) for i in range(jobs * 4)])
            tasks = [(p, w, path, tl[j::nsplit]) for w in W12 for j in range(nsplit)]
            res = pool.map(pent_job, tasks, chunksize=1)
        os.remove(path)
        fb = [b for _, _, bb, _ in res for b in bb]
        failures += fb
        print(f"prime {idx + 1}/{kmax} p={p}: {len(tl)} five-tuples at 12 real roots, failures {len(fb)} "
              f"({time.time()-t1:.0f}s)", flush=True)
    print(f"total time {time.time()-t0:.0f}s")
    full = kmax == max(need.values())
    if failures:
        print("CONCLUSION: FAILURES", failures[:5]); sys.exit(1)
    print(f"CONCLUSION: the pentagon equation holds exactly in K for all {len(inst)} five-tuples with fourth object {fourth}"
          if full else "CONCLUSION: no failures, but not all primes were run")
