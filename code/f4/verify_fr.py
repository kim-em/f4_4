#!/usr/bin/env python3
"""Standalone exact verification of F-, R- and inverse R-symbols over K = Q(zeta_52).

Reads a text data file with F, R and Rinv (format in its header lines; for example
../metaplectic/data/B6_2_FR_exact.txt for C(so(13),2)) and checks, exactly in K:
  * the fusion ring is associative and the unit behaves,
  * R[(X,Y,Z)] Rinv[(X,Y,Z)] = 1,
  * the pentagon equation (in the layout of pentagon-reduction.md) for all five-tuples,
  * both hexagon equations for all four-tuples, in the form
      (H1)  P1 F[(Y,X,W,T)] P2 = F[(X,Y,W,T)] P3 F[(Y,W,X,T)]      with Rinv-blocks P1, P2, P3,
      (H2)  Q5 F[(X,W,Y,T)] Q4 = F[(W,X,Y,T)] Q6 F[(X,Y,W,T)]      with R-blocks Q4, Q5, Q6
    (see hexagon_tuple below for the exact index placement),
  * every F-matrix is invertible.
It uses only the Python standard library (numpy, if present, only speeds up matrix products).
It shares no code with the programs that produced the data.

Method (the height bound of the paper, Section "Verification").  Let L be the lcm of all denominators and H
the largest coefficient of L x, for x running over all entries written in the basis
1, zeta, ..., zeta^23 (F entries are converted exactly).  Each identity, multiplied by a power
of L, becomes an identity between elements of Z[zeta] whose coefficients are bounded by an
explicit B (computed below from L, H, the matrix sizes and the growth constant c of
multiplication modulo Phi_52).  At a prime p = 1 mod 52 not dividing L, Z[zeta]/p is the
product of 24 copies of F_p (zeta -> the 24 primitive 52nd roots); if the identity holds at
all 24 roots for primes p_1 ... p_s with product > 2B, it holds exactly.

usage: verify_fr.py DATA [--nprimes N] [--jobs J] [--no-pentagon] [--no-hexagon]
"""
import itertools, lzma, math, os, sys, time
from fractions import Fraction

try:
    import numpy as np
except ImportError:
    np = None

NB = 24
CONDUCTOR = 52


# ------------------------------------------------------------------ cyclotomic arithmetic
def polymul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def polydiv_exact(a, b):
    a = list(a)
    q = [0] * (len(a) - len(b) + 1)
    for i in range(len(a) - len(b), -1, -1):
        c = a[i + len(b) - 1] // b[-1]
        q[i] = c
        for j in range(len(b)):
            a[i + j] -= c * b[j]
    assert not any(a)
    return q


def cyclotomic(n):
    p = [-1] + [0] * (n - 1) + [1]
    for d in range(1, n):
        if n % d == 0:
            p = polydiv_exact(p, cyclotomic(d))
    return p


PHI = cyclotomic(CONDUCTOR)          # low degree first, monic of degree 24
assert len(PHI) == NB + 1


def reduce(t):
    t = list(t)
    for k in range(len(t) - 1, NB - 1, -1):
        c = t[k]
        if c:
            t[k] = 0
            for j in range(NB + 1):
                t[k - NB + j] -= c * PHI[j]
    return t[:NB] + [0] * (NB - len(t[:NB]))


ZPOW = [reduce([0] * k + [1]) for k in range(CONDUCTOR)]       # zeta^k in the basis, k < 52


def real_to_zeta(a):
    """a0 + sum a_k (zeta^k + zeta^-k) -> coefficients in 1, ..., zeta^23."""
    out = [Fraction(0)] * NB
    out[0] += a[0]
    for k in range(1, len(a)):
        if a[k]:
            for j in range(NB):
                out[j] += a[k] * (ZPOW[k][j] + ZPOW[CONDUCTOR - k][j])
    return out


def growth_constant():
    """least c from the triangle inequality with |ab| <= c |a| |b| (sup norms) in the basis."""
    return max(sum((min(k, 2 * NB - 2 - k) + 1) * abs(ZPOW[k][j]) for k in range(2 * NB - 1)) for j in range(NB))


# ------------------------------------------------------------------ reading
class Data:
    def __init__(self, path):
        opener = lzma.open if path.endswith(".xz") else open
        self.N, self.F, self.R, self.Ri = {}, {}, {}, {}
        with opener(path, "rt") as f:
            for line in f:
                if line.startswith("#"):
                    t = line[1:].split()
                    if t[:1] == ["objects:"]:
                        s = t.index(";")
                        self.obj, self.unit = t[1:s], t[s + 2]
                    elif t[:1] == ["N"] and len(t) == 5:
                        self.N[(t[1], t[2], t[3])] = int(t[4])
                    continue
                parts = line.split("|")
                head = parts[0].split()
                cs = [Fraction(x) for x in parts[-1].split()]
                if head[0] == "F":
                    a = parts[1].split(); b = parts[2].split()
                    assert len(cs) == 12
                    self.F[(tuple(head[1:]), (a[0], int(a[1]), int(a[2])), (b[0], int(b[1]), int(b[2])))] = real_to_zeta(cs)
                else:
                    ij = tuple(int(x) for x in parts[1].split())
                    assert len(cs) == NB
                    (self.R if head[0] == "R" else self.Ri)[(tuple(head[1:]), ij)] = cs
        n = self.n
        self.blocks = {}
        for X, Y, W, T in itertools.product(self.obj, repeat=4):
            A = [(Z, i, j) for Z in self.obj for i in range(n(X, Y, Z)) for j in range(n(Z, W, T))]
            B = [(Z, k, l) for Z in self.obj for k in range(n(Y, W, Z)) for l in range(n(X, Z, T))]
            assert len(A) == len(B), "fusion rules not associative"
            if A:
                self.blocks[(X, Y, W, T)] = (A, B)
        for (key, a, b) in self.F:
            A, B = self.blocks[key]
            assert a in A and b in B
        for (key, ij) in list(self.R) + list(self.Ri):
            assert max(ij) < n(*key)
        for X, Y in itertools.product(self.obj, repeat=2):
            assert n(self.unit, X, Y) == n(X, self.unit, Y) == (X == Y)

    def n(self, X, Y, Z):
        return self.N.get((X, Y, Z), 0)

    def entries(self):
        return list(self.F.values()) + list(self.R.values()) + list(self.Ri.values())


def height_data(data):
    L = 1
    for c in data.entries():
        for x in c:
            L = math.lcm(L, x.denominator)
    H = max(abs(x.numerator) * (L // x.denominator) for c in data.entries() for x in c)
    c = growth_constant()
    O, n = data.obj, data.n
    n5 = max(sum(n(A, B, P) * n(P, C, Q) * n(Q, D, E) for P in O for Q in O)
             for A, B, C, D, E in itertools.product(O, repeat=5))
    n3 = max(len(A) for A, B in data.blocks.values())
    m = max(n(*k) for k in data.N)
    bounds = dict(pentagon=n5 * L * c * H * H + n5 * n5 * c * c * H ** 3,
                  hexagon=2 * n3 * n3 * c * c * H ** 3,
                  inverse=m * c * H * H + L * L)
    return L, H, c, n5, n3, bounds


# ------------------------------------------------------------------ primes and roots
def is_prime(n):
    if n < 2:
        return False
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % q == 0:
            return n == q
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2; s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def primitive_roots_52(p):
    fac, m, d = [], p - 1, 2
    while d * d <= m:
        if m % d == 0:
            fac.append(d)
            while m % d == 0:
                m //= d
        d += 1
    if m > 1:
        fac.append(m)
    g = next(g for g in range(2, p) if all(pow(g, (p - 1) // f, p) != 1 for f in fac))
    z = pow(g, (p - 1) // CONDUCTOR, p)
    roots = [pow(z, w, p) for w in range(1, CONDUCTOR) if math.gcd(w, CONDUCTOR) == 1]
    for r in roots:
        assert sum(c * pow(r, k, p) for k, c in enumerate(PHI)) % p == 0
    assert len(set(roots)) == NB
    return roots


def primes_1_mod_52(count, below=2 ** 31):
    p = below - (below % CONDUCTOR) + 1
    out = []
    while len(out) < count:
        p -= CONDUCTOR
        if is_prime(p):
            out.append(p)
    return out


# ------------------------------------------------------------------ modular checks
def mm(A, B, p):
    if np is not None:
        A = np.asarray(A, dtype=np.int64); B = np.asarray(B, dtype=np.int64)
        if A.size == 0 or B.size == 0:
            return np.zeros((A.shape[0], B.shape[1]), dtype=np.int64)
        A0 = (A & 0xFFFF).astype(np.float64); A1 = (A >> 16).astype(np.float64)
        B0 = (B & 0xFFFF).astype(np.float64); B1 = (B >> 16).astype(np.float64)
        lo = (A0 @ B0).astype(np.int64) % p
        mid = ((A0 @ B1).astype(np.int64) % p + (A1 @ B0).astype(np.int64) % p) % p
        hi = (A1 @ B1).astype(np.int64) % p
        return (lo + (mid << 16) % p + (((hi << 16) % p) << 16) % p) % p
    return [[sum(a * b for a, b in zip(row, col)) % p for col in zip(*B)] for row in A]


def zeros(r, c):
    return np.zeros((r, c), dtype=np.int64) if np is not None else [[0] * c for _ in range(r)]


def eq(A, B):
    if np is not None:
        return np.array_equal(np.asarray(A) % P_CURRENT[0], np.asarray(B) % P_CURRENT[0])
    return all(x % P_CURRENT[0] == y % P_CURRENT[0] for ra, rb in zip(A, B) for x, y in zip(ra, rb))


P_CURRENT = [0]


class Mod:
    def __init__(self, data, p, r):
        self.p = p
        pw = [pow(r, k, p) for k in range(NB)]
        inv = {}

        def ev(c):
            s = 0
            for k, x in enumerate(c):
                if x:
                    d = x.denominator
                    if d not in inv:
                        assert d % p, f"p = {p} divides a denominator"
                        inv[d] = pow(d, -1, p)
                    s += x.numerator * inv[d] % p * pw[k]
            return s % p
        self.F = {}
        for key, (A, B) in data.blocks.items():
            M = zeros(len(A), len(B))
            ia = {a: i for i, a in enumerate(A)}; ib = {b: i for i, b in enumerate(B)}
            self.F[key] = (ia, ib, M)
        for (key, a, b), c in data.F.items():
            ia, ib, M = self.F[key]
            M[ia[a]][ib[b]] = ev(c)
        self.R, self.Ri = {}, {}
        for tgt, src in ((self.R, data.R), (self.Ri, data.Ri)):
            for key in data.N:
                m = data.N[key]
                tgt[key] = zeros(m, m)
            for (key, (i, j)), c in src.items():
                tgt[key][i][j] = ev(c)


def assemble(rows, cols, pieces):
    ri = {x: i for i, x in enumerate(rows)}; ci = {x: i for i, x in enumerate(cols)}
    M = zeros(len(rows), len(cols))
    for (fi, fj), v in pieces:
        M[ri[fi]][ci[fj]] = v
    return M


def block_pieces(mod, key, rmap, cmap):
    if key not in mod.F:
        return []
    ia, ib, M = mod.F[key]
    return [((rmap(a), cmap(b)), M[ia[a]][ib[b]]) for a in ia for b in ib]


def pentagon_tuple(data, mod, A, B, C, D, E):
    O, n, p = data.obj, data.n, mod.p
    I0 = [(P, i1, Q, i2, i3) for P in O for i1 in range(n(A, B, P)) for Q in O for i2 in range(n(P, C, Q)) for i3 in range(n(Q, D, E))]
    if not I0:
        return True, 0
    I1 = [(P, i1, R, k, l) for P in O for i1 in range(n(A, B, P)) for R in O for k in range(n(C, D, R)) for l in range(n(P, R, E))]
    I2 = [(R, k, S, m, nn) for R in O for k in range(n(C, D, R)) for S in O for m in range(n(B, R, S)) for nn in range(n(A, S, E))]
    J1 = [(T, a, b, Q, i3) for T in O for a in range(n(B, C, T)) for Q in O for b in range(n(A, T, Q)) for i3 in range(n(Q, D, E))]
    J2 = [(T, a, S, c, nn) for T in O for a in range(n(B, C, T)) for S in O for c in range(n(T, D, S)) for nn in range(n(A, S, E))]
    pc = lambda key, f, g: block_pieces(mod, key, f, g)
    M1 = assemble(I0, I1, [x for P in O for i in range(n(A, B, P)) for x in pc((P, C, D, E), lambda a, P=P, i=i: (P, i) + a, lambda b, P=P, i=i: (P, i) + b)])
    M2 = assemble(I1, I2, [x for R in O for k in range(n(C, D, R)) for x in pc((A, B, R, E), lambda a, R=R, k=k: (a[0], a[1], R, k, a[2]), lambda b, R=R, k=k: (R, k) + b)])
    N1 = assemble(I0, J1, [x for Q in O for i in range(n(Q, D, E)) for x in pc((A, B, C, Q), lambda a, Q=Q, i=i: (a[0], a[1], Q, a[2], i), lambda b, Q=Q, i=i: b + (Q, i))])
    N2 = assemble(J1, J2, [x for T in O for t in range(n(B, C, T)) for x in pc((A, T, D, E), lambda a, T=T, t=t: (T, t, a[1], a[0], a[2]), lambda b, T=T, t=t: (T, t) + b)])
    N3 = assemble(J2, I2, [x for S in O for nn in range(n(A, S, E)) for x in pc((B, C, D, S), lambda a, S=S, nn=nn: (a[0], a[1], S, a[2], nn), lambda b, S=S, nn=nn: (b[0], b[1], S, b[2], nn))])
    return eq(mm(M1, M2, p), mm(mm(N1, N2, p), N3, p)), len(I0)


def hexagon_tuple(data, mod, X, Y, W, T):
    """Both hexagons for (X, Y, W, T); returns (ok1, ok2)."""
    if (X, Y, W, T) not in data.blocks:
        return True, True
    p, n = mod.p, data.n
    bl = data.blocks
    A_xyw, B_xyw = bl[(X, Y, W, T)]
    A_yxw, B_yxw = bl[(Y, X, W, T)]
    A_ywx, B_ywx = bl[(Y, W, X, T)]
    A_xwy, B_xwy = bl[(X, W, Y, T)]
    A_wxy, B_wxy = bl[(W, X, Y, T)]
    Fm = lambda key: mod.F[key][2]          # indexed by the orders of data.blocks (A rows, B columns)
    Ri, R = mod.Ri, mod.R
    # H1:  P1 F^{YXW} P2 = F^{XYW} P3 F^{YWX}
    P1 = assemble(A_xyw, A_yxw, [(((Z, i, j), (Z, i2, j)), Ri[(X, Y, Z)][i][i2]) for (Z, i, j) in A_xyw for i2 in range(n(Y, X, Z))])
    P2 = assemble(B_yxw, B_ywx, [(((Z, k, l), (Z, k2, l)), Ri[(X, W, Z)][k][k2]) for (Z, k, l) in B_yxw for k2 in range(n(W, X, Z))])
    P3 = assemble(B_xyw, A_ywx, [(((Z, k, l), (Z, k, l2)), Ri[(X, Z, T)][l][l2]) for (Z, k, l) in B_xyw for l2 in range(n(Z, X, T))])
    ok1 = eq(mm(mm(P1, Fm((Y, X, W, T)), p), P2, p), mm(mm(Fm((X, Y, W, T)), P3, p), Fm((Y, W, X, T)), p))
    # H2:  Q5 F^{XWY} Q4 = F^{WXY} Q6 F^{XYW}
    Q5 = assemble(A_wxy, A_xwy, [(((Z, i2, j), (Z, i, j)), R[(X, W, Z)][i2][i]) for (Z, i2, j) in A_wxy for i in range(n(X, W, Z))])
    Q4 = assemble(B_xwy, B_xyw, [(((Z, k2, l), (Z, k, l)), R[(Y, W, Z)][k2][k]) for (Z, k2, l) in B_xwy for k in range(n(Y, W, Z))])
    Q6 = assemble(B_wxy, A_xyw, [(((Z, i, j2), (Z, i, j)), R[(Z, W, T)][j2][j]) for (Z, i, j2) in B_wxy for j in range(n(Z, W, T))])
    ok2 = eq(mm(mm(Q5, Fm((X, W, Y, T)), p), Q4, p), mm(mm(Fm((W, X, Y, T)), Q6, p), Fm((X, Y, W, T)), p))
    return ok1, ok2


def inverse_ok(mod):
    p = mod.p
    for key, M in mod.R.items():
        m = len(M)
        I = [[int(i == j) for j in range(m)] for i in range(m)]
        if not eq(mm(M, mod.Ri[key], p), I):
            return False
    return True


def det_nonzero(M, p):
    M = [[int(x) % p for x in r] for r in M]
    n = len(M)
    for c in range(n):
        piv = next((r for r in range(c, n) if M[r][c]), None)
        if piv is None:
            return False
        M[c], M[piv] = M[piv], M[c]
        inv = pow(M[c][c], -1, p)
        for r in range(c + 1, n):
            f = M[r][c] * inv % p
            if f:
                M[r] = [(x - f * y) % p for x, y in zip(M[r], M[c])]
    return True


_DATA, _FLAGS = None, None


def job(args):
    p, r = args
    P_CURRENT[0] = p
    t0 = time.time()
    mod = Mod(_DATA, p, r)
    res = dict(inverse=inverse_ok(mod))
    if _FLAGS["pentagon"]:
        res["pentagon"] = all(pentagon_tuple(_DATA, mod, *t)[0] for t in itertools.product(_DATA.obj, repeat=5))
    if _FLAGS["hexagon"]:
        h = [hexagon_tuple(_DATA, mod, *t) for t in itertools.product(_DATA.obj, repeat=4)]
        res["hexagon1"] = all(a for a, b in h); res["hexagon2"] = all(b for a, b in h)
    if _FLAGS["invertible"]:
        res["F invertible"] = all(det_nonzero(M, p) for _, _, M in mod.F.values())
    return p, r, res, time.time() - t0


def opt(argv, name, default):
    return argv[argv.index(name) + 1] if name in argv else default


def main(argv):
    global _DATA, _FLAGS
    if len(argv) < 2:
        print(__doc__); return 2
    t0 = time.time()
    data = Data(argv[1])
    _DATA = data
    print(f"read {argv[1]}: {len(data.obj)} objects, {len(data.blocks)} F-matrices, {len(data.F)} nonzero F entries, "
          f"{len(data.R)} + {len(data.Ri)} nonzero R, Rinv entries ({time.time()-t0:.1f} s)")
    L, H, c, n5, n3, bounds = height_data(data)
    need = max(2 * b for b in bounds.values()).bit_length()
    print(f"height data: L {L.bit_length()} bits, H {H.bit_length()} bits, c = {c}, n5 = {n5}, n3 = {n3}; "
          + ", ".join(f"2 x bound({k}) {(2 * b).bit_length()} bits" for k, b in bounds.items()))
    nprimes = int(opt(argv, "--nprimes", need // 30 + 2))
    jobs = int(opt(argv, "--jobs", os.cpu_count() or 1))
    _FLAGS = dict(pentagon="--no-pentagon" not in argv, hexagon="--no-hexagon" not in argv, invertible=True)
    primes = [p for p in primes_1_mod_52(nprimes) if L % p]
    tasks = [(p, r) for p in primes for r in primitive_roots_52(p)]
    print(f"{len(primes)} primes p = 1 mod 52 below 2^31 (none divides L), 24 roots each: {len(tasks)} runs")
    passed = {}
    if jobs > 1:
        import multiprocessing as mpc
        with mpc.get_context("fork").Pool(jobs) as pool:
            results = list(pool.imap_unordered(job, tasks))
    else:
        results = [job(t) for t in tasks]
    for p, r, res, dt in results:
        passed.setdefault(p, []).append(res)
    keys = list(results[0][2].keys())
    allok = {k: all(res[k] for p, rs in passed.items() for res in rs) for k in keys}
    prod = math.prod(p for p, rs in passed.items() if len(rs) == NB and all(all(res.values()) for res in rs))
    print("results over all primes and roots:", allok)
    print(f"product of primes passing at all 24 roots: {prod.bit_length()} bits")
    for k in ("pentagon", "hexagon", "inverse"):
        if k == "pentagon" and not _FLAGS["pentagon"] or k == "hexagon" and not _FLAGS["hexagon"]:
            continue
        verdict = all(v for kk, v in allok.items() if kk.startswith(k)) and prod > 2 * bounds[k]
        print(f"CONCLUSION {k}: {'holds exactly in K' if verdict else 'NOT established'}")
    print("all F-matrices invertible over K:", allok.get("F invertible"))
    return 0 if all(allok.values()) else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
