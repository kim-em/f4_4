#!/usr/bin/env python3
"""Field of moduli of C(F4,4) allowing relabellings: for every fusion ring automorphism pi and every
Galois class a of k = Q(zeta_13)^+ (a in {1,3,5,7,9,11} modulo +-1, +-25), test whether
F^{sigma_a} is gauge equivalent to pi.F, by comparing canonical forms modulo p:
    NF(pi.F) at the root w   versus   NF(F) at the root a w.
(pi.F)[(X,Y,W,T)] = F[(pi^-1 X, ...)] with rows and columns relabelled; its H family is the family
of F with third leg pi^-1(V), read from the full binary table.  The canonical form uses the same
plan (index choices) as f4_canon.py; the nonvanishing conditions are asserted.

usage: f4_relabel.py PLAN p
"""
import os, pickle, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import f4_canon as fc
from lin import Fp
from fusion_aut import fusion_automorphisms

REAL = [1, 3, 5, 7, 9, 11, 15, 17, 19, 21, 23, 25]


def cls(w):
    w %= 52
    return w if w < 26 else 52 - w


def main(argv):
    plan = pickle.load(open(argv[1], "rb"))
    p = int(argv[2])
    F = Fp(p)
    where = fc.readers()
    N = plan["N"]
    auts = fusion_automorphisms(plan["obj"], fc.UNIT, N)
    print(f"{len(auts)} fusion automorphisms")
    nf = {}
    for w in REAL:
        H = fc.load_H(p, w, where)
        _, out = fc.canon(H, plan, F)
        nf[w] = np.concatenate([out[k][2].ravel() for k in sorted(out)])
    for pi in auts:
        inv = {y: x for x, y in pi.items()}
        Vp = inv[fc.V]
        cv = fc.cvals(p, 1)
        H = {}
        for key, r in where.items():
            if key[2] != Vp:
                continue
            rows, cols, _, _, _ = r.block(key)
            M = r.residues(key, p, cv)
            key2 = tuple(pi[x] for x in key)
            rows2 = [(pi[a[0]],) + tuple(a[1:]) for a in rows]
            cols2 = [(pi[b[0]],) + tuple(b[1:]) for b in cols]
            H[key2] = (rows2, cols2, M)
        # reorder to the row/column order of the original H family
        H0 = fc.load_H(p, 1, where)
        Hs = {}
        for k2, (rows2, cols2, M) in H.items():
            rows0, cols0, _ = H0[k2]
            ia = [rows2.index(a) for a in rows0]; ib = [cols2.index(b) for b in cols0]
            Hs[k2] = (rows0, cols0, M[np.ix_(ia, ib)])
        try:
            _, out = fc.canon(Hs, plan, F)
        except AssertionError as e:
            print(f"pi {pi}: normal form conditions fail ({e})"); continue
        v = np.concatenate([out[k][2].ravel() for k in sorted(out)])
        match = sorted({cls(a) for a in REAL if (v == nf[cls(a)]).all()})
        print(f"pi {{{', '.join(f'{x}->{y}' for x, y in pi.items() if x != y)}}}: NF(pi.F) at w=1 equals NF(F) at the roots w in {match}", flush=True)


if __name__ == "__main__":
    main(sys.argv)
