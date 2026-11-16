#!/usr/bin/env python3
"""Exact proof of F^{sigma_a} ~ pi.F for C(F4,4) (pi a fusion ring automorphism).

With g the normal-form gauge of F and g_pi that of pi.F (f4_canon.py, same plan), and
NF(pi.F) = NF(F)^{sigma_a} (f4_relabel.py), the gauge h = g_pi^-1 sigma_a(g) satisfies
h . F^{sigma_a} = pi.F.  Here h is computed modulo p (h(w) = g_pi(w)^-1 g(a w), since
sigma_a(g) at the root z^w is g at z^{aw}), reconstructed over Q(zeta_52)^+, and the identity
      L_h F^{sigma_a}[key] = (pi.F)[key] R_h        for every H-block key = (X, Y, V, T)
is proved exactly by a height bound (Lemma "matrix bound" of the paper), F^{sigma_a}[key] at z^w being
F[key] at z^{aw} and (pi.F)[key] being F[pi^-1 key] relabelled.  By Lemma 4.1 of the paper (both
sides are pentagon solutions with strict unit and V-generated), equality on the H family gives
h . F^{sigma_a} = pi.F everywhere.

usage: f4_relabel_exact.py PLAN a PI NPR NUSE      PI = comma-separated X:Y pairs (the nontrivial part)
The table is f4_canon.PARTS ($F4_ARCHIVE/F4_4_F_part{1,2}.bin); PLAN is data/F4_4_canon.plan.pkl.
"""
import math, os, pickle, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(1, os.path.join(HERE, "..", "f4"))
import f4_canon as fc
from lin import Fp, find_primes
from recon import coeffs, ratrec, REAL
from f4_relabel import cls

C_GROWTH = 46


def relabelled_H(where, p, w, pi, H0):
    inv = {y: x for x, y in pi.items()}
    Vp = inv[fc.V]
    cv = fc.cvals(p, w)
    out = {}
    for key, r in where.items():
        if key[2] != Vp:
            continue
        rows, cols, _, _, _ = r.block(key)
        M = r.residues(key, p, cv)
        k2 = tuple(pi[x] for x in key)
        rows2 = [(pi[a[0]],) + tuple(a[1:]) for a in rows]
        cols2 = [(pi[b[0]],) + tuple(b[1:]) for b in cols]
        rows0, cols0, _ = H0[k2]
        ia = [rows2.index(a) for a in rows0]; ib = [cols2.index(b) for b in cols0]
        out[k2] = (rows0, cols0, M[np.ix_(ia, ib)], key, ia, ib)
    return out


def main(argv):
    plan = pickle.load(open(argv[1], "rb"))
    a = int(argv[2])
    pi = {o: o for o in plan["obj"]}
    for pr in argv[3].split(","):
        x, y = pr.split(":"); pi[x] = y
    npr, nuse = int(argv[4]), int(argv[5])
    where = fc.readers()
    primes = find_primes(52, npr + 40)
    t0 = time.time()
    gkeys = None
    hres = {}
    for p in primes[:npr]:
        F = Fp(p)
        gs, gps = {}, {}
        for w in REAL:
            H0 = fc.load_H(p, w, where)
            gs[w], _ = fc.canon(H0, plan, F)
            Hp = relabelled_H(where, p, w, pi, H0)
            gps[w], out = fc.canon({k: v[:3] for k, v in Hp.items()}, plan, F)
        if gkeys is None:
            gkeys = sorted(t for t in gs[1] if t[0] != "inv")
        hres[p] = {}
        for w in REAL:
            g2 = gs[cls(a * w)]
            hres[p][w] = np.concatenate([F.mm(gps[w][("inv",) + t], g2[t]).ravel() for t in gkeys])
        print(f"p={p}: h residues at 12 roots ({time.time() - t0:.0f} s)", flush=True)
    cs = {p: coeffs({w: hres[p][w].astype(object) for w in REAL}, p, REAL, list(range(12))) for p in hres}
    use, held = list(hres)[:nuse], list(hres)[nuse:]
    Mm = 1
    acc = np.zeros(cs[use[0]].shape, dtype=object)
    for p in use:
        acc = acc + Mm * (((cs[p] - acc) % p) * pow(Mm % p, -1, p) % p); Mm *= p
    acc %= Mm
    hr = [[(0, 1) if int(acc[i, j]) == 0 else ratrec(int(acc[i, j]), Mm) for j in range(12)] for i in range(acc.shape[0])]
    fails = sum(1 for r in hr for x in r if x is None)
    bad = sum(1 for p in held for i, r in enumerate(hr) for j, x in enumerate(r) if x and (x[0] - x[1] * int(cs[p][i, j])) % p)
    mb = max((max(abs(x[0]).bit_length(), x[1].bit_length()) for r in hr for x in r if x), default=0)
    print(f"h: {len(hr)} entries x 12; failures {fails}; mismatches at {len(held)} held-out primes {bad}; largest {mb} bits")
    if fails or bad:
        return
    # height bounds
    shapes = {}
    off = 0
    gpos = {}
    p0 = primes[0]
    H0 = fc.load_H(p0, 1, where)
    g0, _ = fc.canon(H0, plan, Fp(p0))
    for t in gkeys:
        sh = g0[t].shape; gpos[t] = (off, sh); off += sh[0] * sh[1]
    def hd(rows):
        D = 1
        for row in rows:
            for x, y in row:
                D = math.lcm(D, y)
        return D, max((2 * sum(abs(x) * (D // y) for x, y in row) for row in rows), default=0)
    hDH = {t: hd(hr[o:o + sh[0] * sh[1]]) for t, (o, sh) in gpos.items()}
    def blockDH(key):
        _, _, DF, signs, limbs = where[key].block(key)
        nl = limbs.shape[2]
        mx = 0
        for i in range(nl):
            pass
        mags = (limbs.astype(object) * np.array([1 << (32 * i) for i in range(nl)], dtype=object)).sum(axis=2)
        return DF, 2 * 12 * int(mags.max()) if mags.size else 0
    inv = {y: x for x, y in pi.items()}
    Hkeys = sorted(H0)
    bounds = {}
    for key in Hkeys:
        X, Y, W, T = key
        DS, HS = blockDH(key)
        DT, HT = blockDH(tuple(inv[x] for x in key))
        rows, cols, _ = H0[key]
        ts = {(X, Y, Z) for (Z, i, j) in rows} | {(Z, W, T) for (Z, i, j) in rows} | {(Y, W, Z) for (Z, k, l) in cols} | {(X, Z, T) for (Z, k, l) in cols}
        Dg = math.lcm(*[hDH[t][0] for t in ts]); Hg = max(hDH[t][1] * (Dg // hDH[t][0]) for t in ts)
        bounds[key] = len(rows) * C_GROWTH ** 2 * Hg ** 2 * (HS * DT + HT * DS)
    print(f"height bound: largest 2B {(2 * max(bounds.values())).bit_length()} bits ({time.time() - t0:.0f} s)", flush=True)
    proved = {k: False for k in Hkeys}; prod = {k: 1 for k in Hkeys}
    for q in primes:
        if all(proved.values()):
            break
        F = Fp(q)
        try:
            Ch = np.zeros((len(hr), 12), dtype=np.int64)
            for i, r in enumerate(hr):
                for j, (x, y) in enumerate(r):
                    if x:
                        assert y % q
                        Ch[i, j] = x % q * pow(y, -1, q) % q
            for w in REAL:
                cv = fc.cvals(q, w)
                hv = np.zeros(len(hr), dtype=np.int64)
                for j in range(12):
                    hv = (hv + Ch[:, j] * cv[j]) % q
                h = {t: hv[o:o + sh[0] * sh[1]].reshape(sh) for t, (o, sh) in gpos.items()}
                Hs = fc.load_H(q, cls(a * w), where)          # F^{sigma_a} at z^w
                Ht = relabelled_H(where, q, w, pi, fc.load_H(q, w, where))
                for key in Hkeys:
                    if proved[key]:
                        continue
                    rows, cols, Ms = Hs[key]
                    Mt = Ht[key][2]
                    X, Y, W, T = key
                    ri = {r_: i for i, r_ in enumerate(rows)}; ci = {c_: i for i, c_ in enumerate(cols)}
                    L = np.zeros((len(rows), len(rows)), dtype=np.int64)
                    for i_, (Z, i, j) in enumerate(rows):
                        g1, g2 = h[(X, Y, Z)], h[(Z, W, T)]
                        for aa in range(g1.shape[1]):
                            for bb in range(g2.shape[1]):
                                L[i_, ri[(Z, aa, bb)]] = g1[i, aa] * g2[j, bb] % q
                    R = np.zeros((len(cols), len(cols)), dtype=np.int64)
                    for b_, (Z, k, l) in enumerate(cols):
                        h1, h2 = h[(Y, W, Z)], h[(X, Z, T)]
                        for kk in range(h1.shape[1]):
                            for ll in range(h2.shape[1]):
                                R[b_, ci[(Z, kk, ll)]] = h1[k, kk] * h2[l, ll] % q
                    if not (F.mm(L, Ms) == F.mm(Mt, R)).all():
                        print(f"FAILURE q={q} w={w} {key}"); return
        except AssertionError:
            continue
        for k in Hkeys:
            if not proved[k]:
                prod[k] *= q; proved[k] = prod[k] > 2 * bounds[k]
        print(f"q={q}: identity at 12 real roots; proved {sum(proved.values())}/{len(Hkeys)} ({time.time() - t0:.0f} s)", flush=True)
    detok = all(Fp(q).rank(h[t]) == h[t].shape[0] for t in h)
    print(f"h invertible mod {q}: {detok}")
    print("CONCLUSION:", f"h . F^sigma_{a} = pi.F exactly on the H family" if all(proved.values()) and detok else "not proved")


if __name__ == "__main__":
    main(sys.argv)
