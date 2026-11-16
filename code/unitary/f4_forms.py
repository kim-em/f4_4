#!/usr/bin/env python3
"""C(F4,4): the invariant forms G on all 2,302 multiplicity spaces, modulo p, from the H family.

Statement used (paper, Section "Unitarity"): a family of nondegenerate forms G = (G^{XY}_Z) makes every F-matrix
an isometry,  F Gram_R F^T = Gram_L  (Gram_L = sum_Z G^{XY}_Z (x) G^{ZW}_T on the rows (Z,i,j),
Gram_R = sum_Z' G^{YW}_Z' (x) G^{XZ'}_T on the columns (Z',k,l); F real),  iff
F^-T = (G^-1) . F  for the gauge action of ../fields/f4_canon.py.  So G is a gauge from F to
F^-T.  With NF the normal form of ../fields (a function of the gauge class), NF = g1.F = g2.F^-T
gives G = g1^-1 g2 (up to the stabiliser torus, which acts by scalars).  F^-T satisfies the same
pentagon, so by Lemma 4.1 of the paper ("H determines F") it suffices to use the H family.

Input: the normal-form table $F4_ARCHIVE/F4_4_F_part{1,2}.bin (f4_canon.PARTS; read only),
for which g1 = 1 (checked).  Output: residues of G at (p, w).

usage: f4_forms.py OUTDIR p w [w ...]
"""
import os, pickle, sys, time
sys.dont_write_bytecode = True
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
FIELDS = os.path.join(HERE, "..", "fields")
sys.path.insert(0, FIELDS)
sys.path.insert(1, os.path.join(HERE, "..", "f4"))
import f4_canon
from lin import Fp

PLAN = os.path.join(FIELDS, "data", "F4_4_canon.plan.pkl")


def forms_modp(H, plan, F):
    g1, out1 = f4_canon.canon(H, plan, F)
    f4_canon.check_canon(out1, plan, F)
    p = F.p
    ident = all((g1[t] % p == np.eye(g1[t].shape[0], dtype=np.int64)).all() for t in g1 if t[0] != "inv")
    Hi = {k: (r, c, F.inverse(M).T.copy()) for k, (r, c, M) in H.items()}
    g2, out2 = f4_canon.canon(Hi, plan, F)
    f4_canon.check_canon(out2, plan, F)
    same = all((out1[k][2] % p == out2[k][2] % p).all() for k in out1)
    G = {}
    for t in g1:
        if t[0] == "inv":
            continue
        G[t] = F.mm(g1[("inv",) + t], g2[t]) % p
    sym = all((G[t] == G[t].T).all() for t in G)
    return G, ident, same, sym


def main(argv):
    outdir, p = argv[1], int(argv[2])
    os.makedirs(outdir, exist_ok=True)
    plan = pickle.load(open(PLAN, "rb"))
    F = Fp(p)
    where = f4_canon.readers()
    for w in map(int, argv[3:]):
        t0 = time.time()
        H = f4_canon.load_H(p, w, where)
        t1 = time.time()
        G, ident, same, sym = forms_modp(H, plan, F)
        keys = sorted(G)
        np.savez(os.path.join(outdir, f"G_{p}_{w}.npz"), G=np.concatenate([G[t].ravel() for t in keys]))
        lay = os.path.join(outdir, "layout.pkl")
        if not os.path.exists(lay):
            pickle.dump(dict(keys=keys, shapes=[G[t].shape for t in keys]), open(lay, "wb"))
        print(f"[forms] p={p} w={w}: g1 = 1: {ident}; NF(F^-T) = NF(F): {same}; G symmetric: {sym}; "
              f"{len(keys)} spaces; load {t1 - t0:.1f} s, forms {time.time() - t1:.1f} s", flush=True)


if __name__ == "__main__":
    main(sys.argv)
