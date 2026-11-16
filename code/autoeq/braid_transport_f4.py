#!/usr/bin/env python3
"""Transport of the braiding of C(F4,4) along the fusion ring automorphisms, modulo p.

For a pentagon solution F with braiding R, the canonical gauge g = g(F) of ../fields/f4_canon.py
(unique up to the torus T = {u_X u_Y / u_Z}, which fixes every R) gives the gauge invariant
    R_can(F, R)[(X,Y,Z)] = g^{YX}_Z R[(X,Y,Z)] (g^{XY}_Z)^-1 .
Two braided structures (F, R), (F', R') with NF(F) = NF(F') are braided equivalent by an
equivalence that is the identity on objects iff R_can(F,R) = R_can(F',R').
For every fusion ring automorphism pi we compute, at the root z^1 of p,
    A_pi = R_can(pi.F, pi.R)        (pi.F)[key] = F[pi^-1 key],  (pi.R)[(X,Y,Z)] = R[pi^-1 (X,Y,Z)]
and compare it with
    B_w   = R_can(F, R)     at z^w   and   Brev_w = R_can(F, R^rev) at z^w,
    R^rev[(X,Y,Z)] = R[(Y,X,Z)]^-1   (the reverse braiding c^rev_{X,Y} = c_{Y,X}^-1),
for all 24 roots w.  NF(pi.F)(z^1) = NF(F)(z^w) is tested at the same time.  Reads only the
exact data (the binary F-table f4_canon.PARTS, i.e. $F4_ARCHIVE/F4_4_F_part{1,2}.bin, and the R-file,
default ../../deposit/F4_4_R_normalform.txt.xz, or $F4_RFILE) and the plan of ../fields.
All quantities compared are gauge invariant, so any gauge of (F, R) gives the same output.

usage: braid_transport_f4.py p [p ...]
"""
import os, sys, pickle, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "fields")); sys.path.insert(1, os.path.join(HERE, "..", "f4"))
import f4_canon as fc
from f4_relabel_exact import relabelled_H
from lin import Fp, primitive_root_of_unity
from verify_fr import Data
from fusion_aut import fusion_automorphisms

PLAN = os.path.join(HERE, "..", "fields", "data", "F4_4_canon.plan.pkl")
UN = [k for k in range(1, 52) if k % 2 and k % 13]
cls = lambda w: (w % 52) if (w % 52) < 26 else 52 - (w % 52)


def cyc(pi, obj):
    seen, out = set(), []
    for x in obj:
        if x in seen or pi[x] == x:
            continue
        c = [x]; seen.add(x); y = pi[x]
        while y != x:
            c.append(y); seen.add(y); y = pi[y]
        out.append("(" + " ".join(c) + ")")
    return " ".join(out) or "id"


def main(argv):
    plan = pickle.load(open(PLAN, "rb"))
    d = Data(os.environ.get("F4_RFILE", os.path.join(HERE, "..", "..", "deposit", "F4_4_R_normalform.txt.xz")))
    where = fc.readers()
    auts = fusion_automorphisms(plan["obj"], fc.UNIT, plan["N"])
    trip = sorted(d.N)
    t0 = time.time()
    for p in map(int, argv[1:]):
        F = Fp(p)
        def Rmod(w):
            z = primitive_root_of_unity(p, 52, w)
            zp = [pow(z, k, p) for k in range(24)]
            out = {}
            for t in trip:
                m = d.N[t]
                out[t] = np.array([[sum(int(c.numerator) % p * pow(c.denominator, -1, p) * zp[k]
                                        for k, c in enumerate(d.R.get((t, (j, i)), [])) if c) % p
                                    for i in range(m)] for j in range(m)], dtype=np.int64)
            return out
        def rcan(g, R):
            return np.concatenate([F.mm(F.mm(g[(Y, X, Z)], R[(X, Y, Z)]), g[("inv", X, Y, Z)]).ravel() for (X, Y, Z) in trip])
        gw, nf, B, Brev = {}, {}, {}, {}
        Rw = {}
        for w in sorted({cls(u) for u in UN}):
            H = fc.load_H(p, w, where)
            gw[w], out = fc.canon(H, plan, F)
            fc.check_canon(out, plan, F)
            nf[w] = np.concatenate([out[k][2].ravel() for k in sorted(out)])
        H1 = fc.load_H(p, 1, where)
        for w in UN:
            R = Rmod(w); Rw[w] = R
            Rr = {(X, Y, Z): F.inverse(R[(Y, X, Z)]) for (X, Y, Z) in trip}
            B[w] = rcan(gw[cls(w)], R); Brev[w] = rcan(gw[cls(w)], Rr)
        # sanity: R^rev is not R (else the comparison is empty) and R_can is fixed by sigma_27
        print(f"p={p}: R_can(F,R) = R_can(F,R^rev) at z^1: {(B[1] == Brev[1]).all()}; "
              f"roots w with R_can(F,R)(z^w) = R_can(F,R)(z^1): {[w for w in UN if (B[w] == B[1]).all()]}; "
              f"roots w with R_can(F,R^rev)(z^w) = R_can(F,R)(z^1): {[w for w in UN if (Brev[w] == B[1]).all()]}  ({time.time()-t0:.0f} s)", flush=True)
        for pi in auts:
            inv = {y: x for x, y in pi.items()}
            Hp = relabelled_H(where, p, 1, pi, H1)
            gp, out = fc.canon({k: v[:3] for k, v in Hp.items()}, plan, F)
            fc.check_canon(out, plan, F)
            v = np.concatenate([out[k][2].ravel() for k in sorted(out)])
            nfm = sorted(w for w in UN if (v == nf[cls(w)]).all())
            Rpi = {(X, Y, Z): Rw[1][(inv[X], inv[Y], inv[Z])] for (X, Y, Z) in trip}
            A = rcan(gp, Rpi)
            mR = [w for w in UN if (A == B[w]).all()]
            mRr = [w for w in UN if (A == Brev[w]).all()]
            print(f"p={p} pi = {cyc(pi, plan['obj'])}:\n   NF(pi.F)(z^1) = NF(F)(z^w) for w in {nfm}\n"
                  f"   R_can(pi.F, pi.R)(z^1) = R_can(F, R)(z^w) for w in {mR}\n"
                  f"   R_can(pi.F, pi.R)(z^1) = R_can(F, R^rev)(z^w) for w in {mRr}   ({time.time()-t0:.0f} s)", flush=True)


if __name__ == "__main__":
    main(sys.argv)
