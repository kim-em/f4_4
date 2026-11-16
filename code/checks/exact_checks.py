#!/usr/bin/env python3
"""Exact checks on the exported C(F4,4) data used in the paper (read-only on the data).

  1. Basic data of the 16 objects: Weyl dimension, quantum dimension d_X, conformal weight,
     depth in the fusion graph of V = 0001, valuation of d_X at the primes above 13 (parity
     count, confirmed by the 13-adic valuation of the norm of d_X).
  2. Fusion graph of V (from the N lines of the R file).
  3. From the exact R-symbols (export convention), exactly in K = Q(zeta), zeta = exp(2 pi i/52):
       theta_X d_X = sum_Z d_Z Tr R[(X,X,Z)]          versus q^{(X, X + 2 rho)}, q = zeta;
       S~_{XY} = sum_Z d_Z Tr(R[(Y,X,Z)] R[(X,Y,Z)])  versus the Kac-Peterson ratio
                 sum_w eps(w) q^{-2(w(X+rho), Y+rho)} / sum_w eps(w) q^{-2(w rho, rho)};
     and the monodromy R[(Y,X,Z)] R[(X,Y,Z)] = theta_Z/(theta_X theta_Y) (scalar matrix).
  4. From the exact F-table (binary): F[(X,1,X,1)] = 1 and F[(X,X,X,X)][(1,0,0),(1,0,0)] d_X = 1
     exactly, for all 16 X; unit-leg F-matrices are identities.
  5. Writes the gauge-invariant summary file (dimensions, twists, S) for the deposit.

usage: exact_checks.py R.txt F1.bin F2.bin OUT_SUMMARY.txt
"""
import sys, os, itertools, math, time
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
F4DIR = os.environ.get("F4_CODE", os.path.join(HERE, "..", "f4"))   # directory with rootsys.py, fsym_bin.py
sys.path.insert(0, F4DIR)
from rootsys import RootSystem          # root data only
import fsym_bin                          # reader of the binary table

NB, COND = 24, 52


# ---------------------------------------------------------------- arithmetic in Q(zeta_52)
def _phi52():
    # Phi_52(x) = Phi_13(-x^2) = sum_{k=0}^{12} (-1)^k x^{2k}
    c = [0] * 25
    for k in range(13):
        c[2 * k] = (-1) ** k
    return c


PHI = _phi52()


def reduce(poly):
    a = list(poly) + [0] * max(0, NB - len(poly))
    for d in range(len(a) - 1, NB - 1, -1):
        t = a[d]
        if t:
            for j in range(25):
                a[d - 24 + j] -= t * PHI[j]
            assert a[d] == 0
    return a[:NB]


def zpow(k):
    v = [0] * (k % COND) + [1]
    return reduce(v)


ZP = [zpow(k) for k in range(COND)]


def add(a, b):
    return [x + y for x, y in zip(a, b)]


def scal(c, a):
    return [c * x for x in a]


def mul(a, b):
    out = [0] * (2 * NB - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y:
                    out[i + j] += x * y
    return reduce(out)


ZERO = [0] * NB
ONE = [1] + [0] * (NB - 1)


def zeta_minus_inv(a):
    """zeta^a - zeta^-a."""
    return add(ZP[a % COND], scal(-1, ZP[(-a) % COND]))


def real_to_zeta(c):

    out = [Fr(0)] * NB
    out[0] += c[0]
    for k in range(1, 12):
        if c[k]:
            for j in range(NB):
                out[j] += c[k] * (ZP[k][j] + ZP[COND - k][j])
    return out


def is_zero(a):
    return all(x == 0 for x in a)


# ---------------------------------------------------------------- reading R
def read_R(path):
    N, R, Ri, obj, unit = {}, {}, {}, None, None
    import lzma
    opener = lzma.open if path.endswith(".xz") else open
    for line in opener(path, "rt"):
        if line.startswith("#"):
            t = line[1:].split()
            if t[:1] == ["objects:"]:
                s = t.index(";")
                obj, unit = t[1:s], t[s + 2]
            elif t[:1] == ["N"] and len(t) == 5:
                N[(t[1], t[2], t[3])] = int(t[4])
            continue
        parts = line.split("|")
        head = parts[0].split()
        ij = tuple(int(x) for x in parts[1].split())
        cs = [Fr(x) for x in parts[2].split()]
        key = tuple(head[1:4])
        (R if head[0] == "R" else Ri).setdefault(key, {})[ij] = cs
    return obj, unit, N, R, Ri


def mat(entries, n):
    return [[entries.get((i, j), ZERO) for j in range(n)] for i in range(n)]


def matmul(A, B):
    n, m, l = len(A), len(B), len(B[0])
    out = [[ZERO] * l for _ in range(n)]
    for i in range(n):
        for j in range(l):
            acc = [Fr(0)] * NB
            for k in range(m):
                if not is_zero(A[i][k]) and not is_zero(B[k][j]):
                    acc = add(acc, mul(A[i][k], B[k][j]))
            out[i][j] = acc
    return out


def trace(A):
    acc = [Fr(0)] * NB
    for i in range(len(A)):
        acc = add(acc, A[i][i])
    return acc


def main():
    rpath, b1, b2, outpath = sys.argv[1:5]
    t0 = time.time()
    Rs = RootSystem("F4")
    k = 4
    obj, unit, N, R, Ri = read_R(rpath)
    n = lambda a, b, c: N.get((a, b, c), 0)
    wt = {o: tuple(int(c) for c in o) for o in obj}
    V = "0001"
    # ---------------- 1. basic data
    order = 52
    rho = Rs.rho
    pos = Rs.positive_roots
    def pairings(lam):
        return [int(Rs.ip(tuple(l + 1 for l in lam), a)) for a in pos]
    base = pairings((0, 0, 0, 0))
    D = ONE
    for b in base:
        D = mul(D, zeta_minus_inv(b))
    Nx = {}
    for o in obj:
        v = ONE
        for a in pairings(wt[o]):
            v = mul(v, zeta_minus_inv(a))
        Nx[o] = v                                     # d_X = Nx[X] / D
    import cmath
    def emb(x, u=1):
        return sum(complex(c) * cmath.exp(2j * math.pi * u * j / COND) for j, c in enumerate(x))
    dnum = {o: (emb(Nx[o]) / emb(D)).real for o in obj}
    # depth
    depth = {unit: 0}
    layer = [unit]
    while layer:
        nxt = []
        for a in layer:
            for c in obj:
                if n(a, V, c) and c not in depth:
                    depth[c] = depth[a] + 1
                    nxt.append(c)
        layer = nxt
    # valuation at primes above 13 by parity count: v([m]) = 1/12 iff m even, 13 does not divide m (0 < m < 26)
    def v13(lam):
        cnt = lambda ps: sum(1 for m in ps if m % 2 == 0 and m % 13 != 0)
        return Fr(cnt(pairings(lam)) - cnt(base), 12)
    # norm check: v_13(Norm_{K/Q} d_X) = 12 (v_P1 + v_P2) = 24 v
    import mpmath
    mpmath.mp.dps = 80
    def norm(lam):
        val = mpmath.mpf(1)
        for u in range(1, COND):
            if math.gcd(u, COND) != 1:
                continue
            for a, b in zip(pairings(lam), base):
                val *= mpmath.sin(mpmath.pi * u * a / 26) / mpmath.sin(mpmath.pi * u * b / 26)
        return val
    print("1. objects (Dynkin labels, Bourbaki numbering, alpha1 alpha2 long):")
    print(" X     Weyl dim   d_X              h_X     theta exp  depth  v13(d_X)  v13(Norm d_X)")
    hX = {}
    for o in obj:
        cas = int(Rs.casimir(wt[o]))
        hX[o] = cas
        nm = norm(wt[o])
        q = mpmath.identify(nm) if False else None
        # the norm is a rational number; find it with a continued fraction
        fr = Fr(str(mpmath.nstr(nm, 60))).limit_denominator(10 ** 30)
        v = 0
        x = fr
        while x.numerator % 13 == 0:
            x /= 13; v += 1
        while x.denominator % 13 == 0:
            x *= 13; v -= 1
        print(f" {o}  {Rs.weyl_dim(wt[o]):8d}   {dnum[o]:.10f}   {str(Rs.h(wt[o], k)):6s}  {cas % order:4d}       {depth[o]}      "
              f"{str(v13(wt[o])):5s}     {v}  (24 v13 = {24 * v13(wt[o])})")
    # ---------------- 2. fusion graph of V
    print("\n2. fusion graph of V = 0001:  X (x) V = sum N_{X,V}^Z Z")
    for o in obj:
        terms = [f"{n(o, V, z)}*{z}" if n(o, V, z) > 1 else z for z in obj if n(o, V, z)]
        print(f" {o} (x) V = " + " + ".join(terms))
    print(" largest fusion multiplicity:", max(N.values()), [key for key, v in N.items() if v == max(N.values())])
    print(" self-dual (N_{XX}^1 = 1 for all X):", all(n(o, o, unit) == 1 for o in obj))
    # ---------------- 3. modular data from exact R
    Rm = {key: mat(R[key], n(*key)) for key in R}
    Rim = {key: mat(Ri[key], n(*key)) for key in Ri}
    ok_inv = all(matmul(Rm[key], Rim[key]) == [[ONE if i == j else ZERO for j in range(len(Rm[key]))] for i in range(len(Rm[key]))]
                 for key in Rm)
    print("\n3. exact checks from the R-symbols")
    print(" R Rinv = 1 exactly for all", len(Rm), "blocks:", ok_inv)
    ok_twist = True
    for X in obj:
        s = [Fr(0)] * NB
        for Z in obj:
            if n(X, X, Z):
                s = add(s, mul(Nx[Z], trace(Rm[(X, X, Z)])))
        # theta_X d_X = s / D  <=>  theta_X Nx[X] = s
        if not is_zero(add(mul(ZP[hX[X] % order], Nx[X]), scal(-1, s))):
            ok_twist = False
    print(" theta_X d_X = sum_Z d_Z Tr R[(X,X,Z)] with theta_X = q^(X,X+2rho), for all 16 X:", ok_twist)
    ok_mono = True
    for (X, Y, Z) in Rm:
        M = matmul(Rm[(Y, X, Z)], Rm[(X, Y, Z)])
        lam = ZP[(hX[Z] - hX[X] - hX[Y]) % order]
        m = len(M)
        if M != [[lam if i == j else ZERO for j in range(m)] for i in range(m)]:
            ok_mono = False
    print(" R[(Y,X,Z)] R[(X,Y,Z)] = theta_Z/(theta_X theta_Y) * identity for all", len(Rm), "blocks:", ok_mono)
    # Kac-Peterson
    Q = [[Fr(Rs.Ainv[j][i]) * Rs.d[j] for i in range(4)] for j in range(4)]
    def ipw(l, m):     # (l, m) for weights in Dynkin labels
        return sum(l[i] * Q[j][i] * m[j] for i in range(4) for j in range(4))
    def orbit_signed(lam):
        seen = {lam: 1}
        fr = [lam]
        while fr:
            new = []
            for x in fr:
                for i in range(4):
                    y = tuple(x[kk] - x[i] * Rs.A[kk][i] for kk in range(4))
                    if y not in seen:
                        seen[y] = -seen[x]
                        new.append(y)
            fr = new
        return seen
    orbs = {o: orbit_signed(tuple(c + 1 for c in wt[o])) for o in obj}
    assert all(len(v) == 1152 for v in orbs.values())
    def kp(X, Y):
        acc = [0] * NB
        my = tuple(c + 1 for c in wt[Y])
        for pnt, sg in orbs[X].items():
            e = ipw(pnt, my)
            assert e.denominator == 1
            acc = add(acc, scal(sg, ZP[(-2 * int(e)) % order]))
        return acc
    K00 = kp(unit, unit)
    ok_S = True
    St = {}
    for X in obj:
        for Y in obj:
            s = [Fr(0)] * NB
            for Z in obj:
                if n(X, Y, Z):
                    s = add(s, mul(Nx[Z], trace(matmul(Rm[(Y, X, Z)], Rm[(X, Y, Z)]))))
            St[(X, Y)] = s        # = D * S~_{XY}
            # S~_{XY} = s / D  versus  kp(X,Y)/K00
            if not is_zero(add(mul(s, K00), scal(-1, mul(D, kp(X, Y))))):
                ok_S = False
    print(" S~_{XY} = sum_Z d_Z Tr(R[(Y,X,Z)] R[(X,Y,Z)]) equals the Kac-Peterson ratio S_{XY}/S_{00} for all 256 pairs:", ok_S)
    # ---------------- 4. F checks from the binary table
    readers = [fsym_bin.Reader(b1), fsym_bin.Reader(b2)]
    where = {kk: r for r in readers for kk in r.keys()}
    def entry(key, rl, cl):
        rows, cols, Dd, signs, limbs = where[key].block(key)
        i, j = rows.index(rl), cols.index(cl)
        t = i * len(cols) + j
        nl = limbs.shape[2]
        return [Fr(int(signs[t, kk]) * sum(int(limbs[t, kk, a]) << (32 * a) for a in range(nl)), Dd) for kk in range(12)]
    ok_unit_mid = all(entry((X, unit, X, unit), (X, 0, 0), (X, 0, 0)) == [1] + [0] * 11 for X in obj)
    print("\n4. exact checks from the F-table")
    print(" F[(X,1,X,1)] = 1 for all 16 X:", ok_unit_mid)
    ok_fs = True
    for X in obj:
        e = entry((X, X, X, X), (unit, 0, 0), (unit, 0, 0))
        z = real_to_zeta(e)
        if not is_zero(add(mul(z, Nx[X]), scal(-1, D))):
            ok_fs = False
            print("  FAIL", X)
    print(" F[(X,X,X,X)][(1,0,0),(1,0,0)] * d_X = 1 exactly for all 16 X:", ok_fs)
    eV = entry((V, V, V, V), (unit, 0, 0), (unit, 0, 0))
    print(" example: F[(V,V,V,V)][(1,0,0),(1,0,0)] in the basis 1, c_1, ..., c_11:", [str(x) for x in eV])
    # unit-leg identities (all blocks with a unit among the first three legs)
    bad = 0
    cnt = 0
    for key, r in where.items():
        if unit in key[:3]:
            cnt += 1
            rows, cols, Dd, signs, limbs = r.block(key)
            # identity under the natural identification: the binary stores rows/cols in label order;
            # check: exactly one nonzero entry per row, equal to 1 (coefficient vector (1,0,...,0))
            m = len(rows)
            nl = limbs.shape[2]
            for i in range(m):
                nzs = []
                for j in range(m):
                    t = i * m + j
                    v = [int(signs[t, kk]) * sum(int(limbs[t, kk, a]) << (32 * a) for a in range(nl)) for kk in range(12)]
                    if any(v):
                        nzs.append((j, v))
                a = tuple(rows[i])
                X_, Y_, W_, T_ = key
                if Y_ == unit:            # (X,1,W,T): row (X,0,j) <-> column (W,0,j)
                    want = (W_, 0, a[2])
                elif W_ == unit:          # (X,Y,1,T): row (T,i,0) <-> column (Y,0,i)
                    want = (Y_, 0, a[1])
                else:                     # (1,Y,W,T): row (Y,0,j) <-> column (T,j,0)
                    want = (T_, a[2], 0)
                if len(nzs) != 1 or nzs[0][1] != [Dd] + [0] * 11 or tuple(cols[nzs[0][0]]) != want:
                    bad += 1
    print(f" F-matrices with a unit among the first three legs: {cnt}; each is the identity under the "
          f"natural identification of labels: {bad == 0}")
    # ---------------- 5. summary file
    with open(outpath, "w") as f:
        f.write("# C(F4,4): gauge-invariant data computed exactly from F4_4_R_exact.txt (exact_checks.py)\n")
        f.write("# K = Q(zeta), zeta = exp(2 pi i/52), q = zeta; elements as 24 coefficients of 1, zeta, ..., zeta^23\n")
        f.write("# objects in the order of the data files; weights are Dynkin labels (Bourbaki), V = 0001\n")
        f.write("# lines:  OBJ X weyl_dim d_X(numerical) h_X depth theta_exponent   (theta_X = q^theta_exponent)\n")
        f.write("#         DIM X | coefficients of D d_X ; DEN | coefficients of D   (d_X = (D d_X)/D, D = prod_{alpha>0} (q^(rho,alpha) - q^-(rho,alpha)))\n")
        f.write("#         S X Y | coefficients of D S~_{XY},  S~_{XY} = S_{XY}/S_{00} = sum_Z d_Z Tr(R[(Y,X,Z)] R[(X,Y,Z)])\n")
        fmt = lambda v: " ".join(str(x) for x in v)
        for o in obj:
            f.write(f"OBJ {o} {Rs.weyl_dim(wt[o])} {dnum[o]:.12f} {Rs.h(wt[o], k)} {depth[o]} {hX[o] % order}\n")
        f.write("DEN | " + fmt(D) + "\n")
        for o in obj:
            f.write(f"DIM {o} | " + fmt(Nx[o]) + "\n")
        for X in obj:
            for Y in obj:
                f.write(f"S {X} {Y} | " + fmt(St[(X, Y)]) + "\n")
    print(f"\nwrote {outpath}; {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
