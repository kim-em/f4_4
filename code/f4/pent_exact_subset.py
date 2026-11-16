"""Exact proof of the pentagon equation for the five-tuples with entries in a subset of objects,
for the exported F4 data (binary format, export convention, real basis).

For each five-tuple, M1 M2 = N1 N2 N3 (layout of pentagon-reduction.md) is multiplied by the
scales s of the five factors (s = lcm of the common denominators of the F-blocks assembled into
the factor) and becomes an identity in Z[zeta] with coefficients bounded by
   B = m c H_M1 H_M2 s_N1 s_N2 s_N3 + m^2 c^2 H_N1 H_N2 H_N3 s_M1 s_M2
(H: coefficient bounds of the scaled factors in the zeta basis, m the dimension, c as in
verify_fr.py).  The five-tuple is checked as full matrices at all 24 roots of the first k primes
p = 1 mod 52 below 2^31 not dividing any denominator, with product > 2B.

usage: pent_exact_subset.py SUBSET F1.bin [F2.bin ...] [--jobs J]
       SUBSET = comma-separated objects, or maxdim=M for all five-tuples with dim Hom(E, ABCD) <= M
"""
import sys, os, math, time, itertools
import multiprocessing as mp
import numpy as np

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
import fsym_bin
from verify_fr import ZPOW, CONDUCTOR, NB, growth_constant, primes_1_mod_52, primitive_roots_52, mm, assemble

subset = sys.argv[1].split(",") if not sys.argv[1].startswith("maxdim=") else None
maxdim = int(sys.argv[1].split("=")[1]) if subset is None else None
bins = [a for a in sys.argv[2:] if a.endswith(".bin")]
jobs = int(sys.argv[sys.argv.index("--jobs") + 1]) if "--jobs" in sys.argv else 16
t0 = time.time()
readers = [fsym_bin.Reader(b) for b in bins]
where = {k: r for r in readers for k in r.keys()}
obj = readers[0].meta.get("objects") or sorted({k[0] for k in where})
# fusion rules from the block labels: N(X,Y,Z) = number of (Z, i, .) row labels of F[(X,Y,1,Z)]
unit = "0000"
N = {}
for key in where:
    X, Y, W, T = key
    if W == unit:
        rows = where[key].index[key]["rows"]
        N[(X, Y, T)] = len(rows)
n = lambda a, b, cc: N.get((a, b, cc), 0)
objs = sorted({k[0] for k in where})
if subset is None:
    subset = objs
T = [[1 if j == 0 else 0 for j in range(NB)]] + [[ZPOW[k][j] + ZPOW[CONDUCTOR - k][j] for j in range(NB)] for k in range(1, 12)]
CONV = max(sum(abs(T[k][j]) for k in range(12)) for j in range(NB))
c = growth_constant()
meta = {}


def blockmeta(key):
    if key not in meta:
        rows, cols, D, signs, limbs = where[key].block(key)
        nl = limbs.shape[2]
        top = int(limbs[:, :, nl - 1].max()) if limbs.size else 0
        meta[key] = (D, CONV * ((top + 1) << (32 * (nl - 1))))
    return meta[key]


def factor(keys):
    keys = [k for k in keys if k in where]
    D = 1
    for k in keys:
        D = math.lcm(D, blockmeta(k)[0])
    return D, max((blockmeta(k)[1] * (D // blockmeta(k)[0]) for k in keys), default=0)


inst = []
for A, B, C, D_, E in itertools.product(subset, repeat=5):
    m = sum(n(A, B, P) * n(P, C, Q) * n(Q, D_, E) for P in objs for Q in objs)
    if not m or (maxdim is not None and m > maxdim):
        continue
    M1 = factor([(P, C, D_, E) for P in objs if n(A, B, P)])
    M2 = factor([(A, B, R, E) for R in objs if n(C, D_, R)])
    N1 = factor([(A, B, C, Q) for Q in objs if n(Q, D_, E)])
    N2 = factor([(A, Tt, D_, E) for Tt in objs if n(B, C, Tt)])
    N3 = factor([(B, C, D_, S) for S in objs if n(A, S, E)])
    b = m * c * M1[1] * M2[1] * N1[0] * N2[0] * N3[0] + m * m * c * c * N1[1] * N2[1] * N3[1] * M1[0] * M2[0]
    inst.append(((A, B, C, D_, E), m, 2 * b))
dens = [v[0] for v in meta.values()]
primes = [p for p in primes_1_mod_52(400) if all(d % p for d in dens)]
prods = [1]
for p in primes:
    prods.append(prods[-1] * p)
need = {t: next(k for k in range(1, len(prods)) if prods[k] > b) for t, m, b in inst}
kmax = max(need.values())
print(f"subset {subset if maxdim is None else 'all objects, dimension <= ' + str(maxdim)}: {len(inst)} five-tuples, largest dimension {max(m for _, m, _ in inst)}, "
      f"2 x bound up to {max(b for _, _, b in inst).bit_length()} bits, primes needed: max {kmax}, "
      f"mean {sum(need.values()) / len(need):.1f} ({time.time()-t0:.0f}s)", flush=True)
keys_used = sorted({k for k in meta})


def job(args):
    p, w, tuples = args
    rd = {id(r): fsym_bin.Reader(r.f.name) for r in readers}
    r = primitive_roots_52(p)[0]
    rr = pow(r, w, p); ri = pow(rr, -1, p)
    cv = [1] + [(pow(rr, k, p) + pow(ri, k, p)) % p for k in range(1, 12)]
    Fm = {}

    def F(key):
        if key not in Fm:
            Fm[key] = rd[id(where[key])].residues(key, p, cv)
        return Fm[key]

    def pieces(key, rmap, cmap):
        if key not in where:
            return []
        b = where[key].index[key]
        rows = [tuple(x) for x in b["rows"]]; cols = [tuple(x) for x in b["cols"]]
        M = F(key)
        return [((rmap(a), cmap(bb)), M[i, j]) for i, a in enumerate(rows) for j, bb in enumerate(cols)]
    from verify_fr import P_CURRENT
    P_CURRENT[0] = p
    bad = 0
    O = objs
    for (A, B, C, D, E) in tuples:
        cv5 = lambda L: L
        I0 = [(P, i1, Q, i2, i3) for P in O for i1 in range(n(A, B, P)) for Q in O for i2 in range(n(P, C, Q)) for i3 in range(n(Q, D, E))]
        I1 = [(P, i1, R, k, l) for P in O for i1 in range(n(A, B, P)) for R in O for k in range(n(C, D, R)) for l in range(n(P, R, E))]
        I2 = [(R, k, S, m_, nn) for R in O for k in range(n(C, D, R)) for S in O for m_ in range(n(B, R, S)) for nn in range(n(A, S, E))]
        J1 = [(Tt, a, b, Q, i3) for Tt in O for a in range(n(B, C, Tt)) for Q in O for b in range(n(A, Tt, Q)) for i3 in range(n(Q, D, E))]
        J2 = [(Tt, a, S, cc, nn) for Tt in O for a in range(n(B, C, Tt)) for S in O for cc in range(n(Tt, D, S)) for nn in range(n(A, S, E))]
        M1 = assemble(I0, I1, [x for P in O for i in range(n(A, B, P)) for x in pieces((P, C, D, E), lambda a, P=P, i=i: (P, i) + a, lambda b, P=P, i=i: (P, i) + b)])
        M2 = assemble(I1, I2, [x for R in O for k in range(n(C, D, R)) for x in pieces((A, B, R, E), lambda a, R=R, k=k: (a[0], a[1], R, k, a[2]), lambda b, R=R, k=k: (R, k) + b)])
        N1 = assemble(I0, J1, [x for Q in O for i in range(n(Q, D, E)) for x in pieces((A, B, C, Q), lambda a, Q=Q, i=i: (a[0], a[1], Q, a[2], i), lambda b, Q=Q, i=i: b + (Q, i))])
        N2 = assemble(J1, J2, [x for Tt in O for t in range(n(B, C, Tt)) for x in pieces((A, Tt, D, E), lambda a, Tt=Tt, t=t: (Tt, t, a[1], a[0], a[2]), lambda b, Tt=Tt, t=t: (Tt, t) + b)])
        N3 = assemble(J2, I2, [x for S in O for nn in range(n(A, S, E)) for x in pieces((B, C, D, S), lambda a, S=S, nn=nn: (a[0], a[1], S, a[2], nn), lambda b, S=S, nn=nn: (b[0], b[1], S, b[2], nn))])
        if not np.array_equal(mm(M1, M2, p), mm(mm(N1, N2, p), N3, p)):
            bad += 1
    return p, w, bad, len(tuples)


if __name__ == "__main__":
    units = [u for u in range(1, CONDUCTOR) if math.gcd(u, CONDUCTOR) == 1]
    tasks = []
    for idx, p in enumerate(primes[:kmax]):
        tl = [t for t, m, b in inst if need[t] > idx]
        tasks += [(p, w, tl) for w in units]
    with mp.get_context("fork").Pool(jobs) as pool:
        res = pool.map(job, tasks, chunksize=1)
    bad = sum(r[2] for r in res)
    checks = sum(r[3] for r in res)
    print(f"{len(tasks)} (prime, root) runs, {checks} five-tuple checks, failures {bad}; {time.time()-t0:.0f}s")
    what = f"from {subset}" if maxdim is None else f"with dim Hom(E, ABCD) <= {maxdim}"
    print("CONCLUSION:", f"the pentagon equation holds exactly in K for all {len(inst)} five-tuples {what}"
          if bad == 0 else "FAILURES")
    sys.exit(1 if bad else 0)
