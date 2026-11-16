#!/usr/bin/env python3
"""Reconstruct R and Rinv of C(F4,4) in the normal-form gauge over Q(zeta_13) from the direct
residues WORK/R_<p>_<w>.npz (direct_f4.py), never touching another gauge.

Q(zeta_13) has 12 embeddings; the roots w = 1,3,5,7,9,11 and 52 - w = 51,...,41 are distinct
modulo 13, so they are one root per orbit of sigma_27 (the other 12 roots w + 26 give the same
values, which is checked at the held-out primes wherever both are present).  Basis zeta13^j = zeta^{4j},
j = 0..11.  Output in the text format of the R-file F4_4_R_normalform.txt (zeta_52 basis); the header lines
(conventions and fusion rules N) are taken from data/F4_4_R_header.txt.  The primes are read from
WORK/P.txt and WORK/Q.txt.

usage: direct_R_recon.py WORK OUT.txt
"""
import os, pickle, sys
from fractions import Fraction
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(1, os.path.join(HERE, "..", "f4"))
import cyc
from lin import primitive_root_of_unity
from recon import inv_mat
from recon import crt_rec

W12 = [1, 3, 5, 7, 9, 11, 41, 43, 45, 47, 49, 51]


def main(argv):
    work, out = argv[1], argv[2]
    P = [int(x) for x in open(os.path.join(work, "P.txt")).read().split()]
    Q = [int(x) for x in open(os.path.join(work, "Q.txt")).read().split()]
    tk = None
    cs = {}
    allres = {}
    for p in P + Q:
        res = {}
        for w in range(1, 52, 2):
            f = os.path.join(work, f"R_{p}_{w}.npz")
            if os.path.exists(f):
                z = np.load(f)
                res[w] = np.concatenate([z["R"], z["Rinv"]]).astype(object)
                if tk is None:
                    tk = pickle.load(open(os.path.join(work, f"R_{p}_{w}.keys.pkl"), "rb"))
        allres[p] = res
        ws = [w for w in W12 if w in res]
        if len(ws) < 12:
            ws = []
            for c in range(1, 13):                    # any root in each class mod 13
                cand = [w for w in res if w % 13 == c]
                if cand:
                    ws.append(cand[0])
        if len(ws) < 12:
            continue
        z = primitive_root_of_unity(p, 52, 1)
        Vi = np.array(inv_mat([[pow(z, 4 * j * w, p) for j in range(12)] for w in ws], p), dtype=object)
        cs[p] = (Vi.dot(np.array([res[w] for w in ws])) % p).T
    use = [p for p in P if p in cs]
    rec, fails, mb, M = crt_rec(cs, use)
    # check every available root of every held-out prime (and sigma_27 pairs w, w+26 at all primes)
    bad = 0
    for p in Q:
        z = primitive_root_of_unity(p, 52, 1)
        for w, v in allres[p].items():
            b = [pow(z, 4 * j * w, p) for j in range(12)]
            for i, row in enumerate(rec):
                val = sum(x % p * pow(y, -1, p) * b[j] for j, (x, y) in enumerate(row) if x) % p
                if val != int(v[i]) % p:
                    bad += 1
    print(f"[R direct] {len(rec)} entries (R and Rinv) x 12 from {len(use)} primes: failures {fails}, largest {mb} bits; "
          f"mismatches at the roots of the held-out primes {Q}: {bad} "
          f"(roots checked: {sorted({w for p in Q for w in allres[p]})})")
    if fails or bad:
        return
    src = os.path.join(HERE, "data", "F4_4_R_header.txt")
    head = [l for l in open(src) if l.startswith("#")]
    N = {}
    for l in head:
        t = l[1:].split()
        if t[:1] == ["N"] and len(t) == 5:
            N[(t[1], t[2], t[3])] = int(t[4])
    def zeta(row):
        v = cyc.ZERO
        for j, (x, y) in enumerate(row):
            if x:
                v = cyc.add(v, cyc.scal(Fraction(x, y), cyc.zpow(4 * j)))
        return v
    pos, off = {}, 0
    for t in tk:
        pos[t] = off; off += N[t] ** 2
    half = off
    # The header text below is part of the archived R-file (kept verbatim so that the output is
    # identical byte for byte).
    with open(out, "w") as f:
        f.write("# C(F4,4): R-symbols and inverse R-symbols in the normal-form gauge of ../fields (the gauge of\n")
        f.write("# F4_4_F_descended_part{1,2}.bin).  All entries lie in Q(zeta_13) (zeta_13 = zeta^4) but are written, like\n")
        f.write("# ../f4/data/F4_4_R_exact.txt, in the power basis of zeta = exp(2 pi i/52).  Format and conventions:\n")
        for l in head[1:]:
            f.write(l)
        for t in sorted(N):
            m = N[t]
            for j in range(m):
                for i in range(m):
                    v = zeta(rec[pos[t] + j * m + i])
                    if not cyc.is_zero(v):
                        f.write("R " + " ".join(t) + f" | {j} {i} | " + " ".join(str(x) for x in v) + "\n")
        for t in sorted(N):
            m = N[t]
            for i in range(m):
                for j in range(m):
                    v = zeta(rec[half + pos[t] + i * m + j])
                    if not cyc.is_zero(v):
                        f.write("Rinv " + " ".join(t) + f" | {i} {j} | " + " ".join(str(x) for x in v) + "\n")
    print("wrote", out)


if __name__ == "__main__":
    main(sys.argv)
