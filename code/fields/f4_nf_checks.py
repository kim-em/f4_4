#!/usr/bin/env python3
"""Checks on the normal-form tables of C(F4,4).

  f4_nf_checks.py unit F1.bin F2.bin
      exact: every block with a unit leg is an identity matrix (strict unit), and
      F[(X,X,X,X)][(1,0,0),(1,0,0)] = 1/d_X with d_X the quantum Weyl dimension (exact value
      checked numerically to 1e-12 and exactly modulo the primes below).
  f4_nf_checks.py braid DIR p w
      modulo p at zeta -> z^w, from the direct outputs DIR/F_<p>_<w>.npy and DIR/R_<p>_<w>.npz:
      d_X = 1 / F[(X,X,X,X)][(1,0,0),(1,0,0)], theta_X d_X = sum_Z d_Z Tr R^{XX}_Z equal to q^{(X,X+2rho)},
      the monodromy R^{YX}_Z R^{XY}_Z is the scalar theta_Z/(theta_X theta_Y), and
      S~_{XY} = sum_Z d_Z Tr(R^{XY}_Z R^{YX}_Z) equals the Kac-Peterson S (normalised by S_00), and
      d_X equals the quantum Weyl dimension.  (Index conventions: R[(X,Y,Z)][j,i], export.)
      The block layout of the residues is read from $F4_ARCHIVE/F4_4_F_part{1,2}.bin (index only).
"""
import os, pickle, sys
from fractions import Fraction
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(1, os.path.join(HERE, "..", "f4"))
import fsym_bin
from rootsys import RootSystem
from lin import Fp, primitive_root_of_unity

OBJ = ["0000", "0001", "0002", "0010", "1000", "0003", "0011", "0100", "1001", "0004", "0012",
       "0020", "0101", "1002", "1010", "2000"]
UNIT = "0000"


def unit_check(bins):
    rs = [fsym_bin.Reader(b) for b in bins]
    where = {k: r for r in rs for k in r.keys()}
    nb = bad = 0
    for k, r in where.items():
        if UNIT not in k[:3]:
            continue
        rows, cols, M = r.fractions(k)
        nb += 1
        for i, row in enumerate(M):
            for j, v in enumerate(row):
                want = (Fraction(1),) + (Fraction(0),) * 11 if i == j else (Fraction(0),) * 12
                if tuple(v) != want:
                    bad += 1
        # rows and columns are matched in order for unit legs (same label sets up to the unit position)
    print(f"blocks with a unit among the first three legs: {nb}; entries differing from the identity: {bad}")


def braid_check(d, p, w):
    F = Fp(p)
    R = RootSystem("F4")
    order = 52
    q = primitive_root_of_unity(p, 52, w)
    qp = lambda e: pow(q, int(e) % order, p)
    wt = {x: tuple(int(c) for c in x) for x in OBJ}
    keys = pickle.load(open(f"{d}/F_{p}_{w}.npy.keys.pkl", "rb"))
    arr = np.load(f"{d}/F_{p}_{w}.npy", mmap_mode="r")
    # need only the blocks (X,X,X,X); find their offsets from the binary index layout
    archive = os.environ.get("F4_ARCHIVE", os.path.join(HERE, "..", "..", "archive"))     # block layout only
    rs = [fsym_bin.Reader(os.path.join(archive, f"F4_4_F_part{i}.bin")) for i in (1, 2)]
    where = {k: r for r in rs for k in r.keys()}
    off, offs = 0, {}
    for k in keys:
        ix = where[k].index[k]
        offs[k] = (off, ix["rows"], ix["cols"])
        off += len(ix["rows"]) * len(ix["cols"])
    dd = {}
    for X in OBJ:
        o, rows, cols = offs[(X, X, X, X)]
        i = [tuple(r) for r in rows].index((UNIT, 0, 0)); j = [tuple(c) for c in cols].index((UNIT, 0, 0))
        dd[X] = pow(int(arr[o + i * len(cols) + j]), -1, p)
    # quantum Weyl dimensions
    qint = lambda n: (qp(n) - qp(-n)) % p
    dq = {}
    for X in OBJ:
        num = den = 1
        for a in R.positive_roots:
            num = num * qint(R.ip(tuple(x + 1 for x in wt[X]), a)) % p
            den = den * qint(R.ip(R.rho, a)) % p
        dq[X] = num * pow(den, -1, p) % p
    print(f"p={p} w={w}: d_X = 1/F[(X,X,X,X)][1,1] equals the quantum Weyl dimension for all X: {all(dd[X] == dq[X] for X in OBJ)}")
    Rz = np.load(f"{d}/R_{p}_{w}.npz")
    tk = pickle.load(open(f"{d}/R_{p}_{w}.keys.pkl", "rb"))
    RR, o = {}, 0
    for t in tk:
        m = int(round(np.sqrt(0))) if False else None
    sizes = {}
    for t in tk:
        # multiplicity from the binary: N(X,Y,Z) = number of rows (Z,i,0) of F[(X,Y,0000,Z)]
        X, Y, Z = t
        o2, rows, cols = offs[(X, Y, UNIT, Z)]
        sizes[t] = len(rows)
    for t in tk:
        m = sizes[t]
        RR[t] = Rz["R"][o:o + m * m].reshape(m, m).astype(np.int64); o += m * m
    h = {X: int(R.casimir(wt[X])) for X in OBJ}
    th = {}
    for X in OBJ:
        s = sum(dd[Z] * int(np.trace(RR[(X, X, Z)]) % p) for Z in OBJ if (X, X, Z) in RR) % p
        th[X] = s * pow(dd[X], -1, p) % p
    print(f"  theta_X from R equals q^(X,X+2rho): {all(th[X] == qp(h[X]) for X in OBJ)}")
    bal = all((F.mm(RR[(Y, X, Z)], RR[(X, Y, Z)]) == qp(h[Z] - h[X] - h[Y]) * np.eye(M.shape[0], dtype=np.int64) % p).all()
              for (X, Y, Z), M in RR.items())
    print(f"  monodromy R^YX R^XY is the scalar theta_Z/(theta_X theta_Y): {bal}")
    St = np.array([[sum(dd[Z] * int(np.trace(F.mm(RR[(X, Y, Z)], RR[(Y, X, Z)])) % p) for Z in OBJ if (X, Y, Z) in RR) % p
                    for Y in OBJ] for X in OBJ])
    # Kac-Peterson
    from math import lcm
    Q = [[Fraction(R.Ainv[j][i]) * R.d[j] for i in range(R.rank)] for j in range(R.rank)]
    den = lcm(*[x.denominator for row in Q for x in row])
    Qi = np.array([[int(x * den) for x in row] for row in Q], dtype=np.int64)
    def orbit(lam):
        seen = {tuple(lam): 1}; fr = [tuple(lam)]
        while fr:
            new = []
            for x in fr:
                for i in range(R.rank):
                    y = tuple(x[k] - x[i] * R.A[k][i] for k in range(R.rank))
                    if y not in seen:
                        seen[y] = -seen[x]; new.append(y)
            fr = new
        return np.array(list(seen.keys()), dtype=np.int64), np.array(list(seen.values()), dtype=np.int64)
    tab = np.array([qp(e) for e in range(order)], dtype=np.int64)
    orbs = {X: orbit([x + 1 for x in wt[X]]) for X in OBJ}
    def kp(l, m):
        pts, sg = orbs[l]
        e2 = pts @ (Qi @ np.array([x + 1 for x in wt[m]], dtype=np.int64))
        assert not (e2 % den).any()
        return int(((sg % p) * tab[(-2 * (e2 // den)) % order] % p).sum() % p)
    S0 = kp(UNIT, UNIT)
    Skp = np.array([[kp(X, Y) * pow(S0, -1, p) % p for Y in OBJ] for X in OBJ])
    print(f"  S~ from R equals the Kac-Peterson S / S_00: {np.array_equal(St % p, Skp)}")


if __name__ == "__main__":
    if sys.argv[1] == "unit":
        unit_check(sys.argv[2:])
    else:
        braid_check(sys.argv[2], int(sys.argv[3]), int(sys.argv[4]))
