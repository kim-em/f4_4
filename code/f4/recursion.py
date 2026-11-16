"""All F-symbols from the H family G[(X, Y, V, T)] by the pentagon equation, modulo p.

Splitting convention (hfamily.py): L_a = sum_b G[a, b] R_b.  For the five bases of
Hom(E, A B C D) (labels as in pentagon-reduction.md) the pentagon reads M1 M2 = N1 N2 N3 with
  M1 = (+) G[(P,C,D,E)],  M2 = (+) G[(A,B,R,E)],  N1 = (+) G[(A,B,C,Q)],
  N2 = (+) G[(A,T,D,E)],  N3 = (+) G[(B,C,D,S)].
Step 1 (second leg V):  (A,B,C,D,E) = (X, V, W', V, E),  M2 = M1^-1 N1 N2 N3 gives G[(X,V,R,E)]
        for R in W' (x) V from G with smaller third leg and from H.
Step 2 (all):          (A,B,C,D,E) = (X, Y', V, W, E),  N2 = N1^-1 M1 M2 N3^-1 gives G[(X,T,W,E)]
        for T in Y' (x) V from G with smaller middle leg and from steps 1 and H.
Each step produces the blocks for every R (resp. T) in W' (x) V (resp. Y' (x) V) and every
multiplicity copy; the copies must agree, off-diagonal blocks must vanish, and blocks already
known (smaller depth, or H itself) must be reproduced.  These are the checks reported.
"""
import numpy as np
from lin import Fp


class Recursion:
    def __init__(self, alc, fus, G, p, unit, gen, depth):
        self.alc, self.fus, self.G, self.p, self.unit, self.gen = alc, fus, G, p, unit, gen
        self.F = Fp(p)
        self.depth = depth
        self._inv = {}
        self.checks = dict(copies=0, offdiag=0, known=0)
        self._add_units()

    def N(self, a, b, c):
        return self.fus.get((a, b, c), 0)

    def _add_units(self):
        u, alc, N = self.unit, self.alc, self.N
        for X in alc:
            for Y in alc:
                for T in alc:
                    n = N(X, Y, T)
                    if not n:
                        continue
                    # (X, 1, Y, T): rows (X,0,j), cols (Y,0,j)
                    self.G.setdefault((X, u, Y, T), ([(X, 0, j) for j in range(n)], [(Y, 0, j) for j in range(n)],
                                                     np.eye(n, dtype=np.int64)))
                    # (X, Y, 1, T): rows (T,i,0), cols (Y,0,i)
                    self.G.setdefault((X, Y, u, T), ([(T, i, 0) for i in range(n)], [(Y, 0, i) for i in range(n)],
                                                     np.eye(n, dtype=np.int64)))
                    # (1, X, Y, T): rows (X,0,j), cols (T,j,0)
                    self.G.setdefault((u, X, Y, T), ([(X, 0, j) for j in range(n)], [(T, j, 0) for j in range(n)],
                                                     np.eye(n, dtype=np.int64)))

    def ginv(self, key):
        if key not in self._inv:
            rl, cl, M = self.G[key]
            self._inv[key] = (cl, rl, self.F.inverse(M))
        return self._inv[key]

    def labels(self, A, B, C, D, E):
        N, O = self.N, self.alc
        I0 = [(P, i1, Q, i2, i3) for P in O for i1 in range(N(A, B, P)) for Q in O for i2 in range(N(P, C, Q)) for i3 in range(N(Q, D, E))]
        I1 = [(P, i1, R, k, l) for P in O for i1 in range(N(A, B, P)) for R in O for k in range(N(C, D, R)) for l in range(N(P, R, E))]
        I2 = [(R, k, S, m, n) for R in O for k in range(N(C, D, R)) for S in O for m in range(N(B, R, S)) for n in range(N(A, S, E))]
        J1 = [(T, a, b, Q, i3) for T in O for a in range(N(B, C, T)) for Q in O for b in range(N(A, T, Q)) for i3 in range(N(Q, D, E))]
        J2 = [(T, a, S, c, n) for T in O for a in range(N(B, C, T)) for S in O for c in range(N(T, D, S)) for n in range(N(A, S, E))]
        return {name: {x: i for i, x in enumerate(L)} for name, L in zip(("I0", "I1", "I2", "J1", "J2"), (I0, I1, I2, J1, J2))}

    # block lists: (row label list, col label list, matrix) in pentagon coordinates
    def blocks(self, which, A, B, C, D, E, inverse=False):
        N, O, out = self.N, self.alc, []
        get = self.ginv if inverse else (lambda key: self.G[key])
        if which == "M1":
            for P in O:
                for i in range(N(A, B, P)):
                    key = (P, C, D, E)
                    if key in self.G:
                        rl, cl, M = get(key)
                        out.append(([(P, i) + a for a in rl], [(P, i) + b for b in cl], M))
        elif which == "M2":
            for R in O:
                for k in range(N(C, D, R)):
                    key = (A, B, R, E)
                    if key in self.G:
                        rl, cl, M = get(key)
                        f = (lambda a, R=R, k=k: (a[0], a[1], R, k, a[2]))
                        g = (lambda b, R=R, k=k: (R, k) + b)
                        if inverse:
                            f, g = g, f
                        out.append(([f(a) for a in rl], [g(b) for b in cl], M))
        elif which == "N1":
            for Q in O:
                for i in range(N(Q, D, E)):
                    key = (A, B, C, Q)
                    if key in self.G:
                        rl, cl, M = get(key)
                        f = (lambda a, Q=Q, i=i: (a[0], a[1], Q, a[2], i))
                        g = (lambda b, Q=Q, i=i: b + (Q, i))
                        if inverse:
                            f, g = g, f
                        out.append(([f(a) for a in rl], [g(b) for b in cl], M))
        elif which == "N2":
            for T in O:
                for x in range(N(B, C, T)):
                    key = (A, T, D, E)
                    if key in self.G:
                        rl, cl, M = get(key)
                        f = (lambda a, T=T, x=x: (T, x, a[1], a[0], a[2]))
                        g = (lambda b, T=T, x=x: (T, x) + b)
                        if inverse:
                            f, g = g, f
                        out.append(([f(a) for a in rl], [g(b) for b in cl], M))
        elif which == "N3":
            for S in O:
                for n in range(N(A, S, E)):
                    key = (B, C, D, S)
                    if key in self.G:
                        rl, cl, M = get(key)
                        f = (lambda a, S=S, n=n: (a[0], a[1], S, a[2], n))
                        g = (lambda b, S=S, n=n: (b[0], b[1], S, b[2], n))
                        if inverse:
                            f, g = g, f
                        out.append(([f(a) for a in rl], [g(b) for b in cl], M))
        return out

    def times(self, thin, src, dst, blks):
        """thin (h x |src|) times the block matrix (src -> dst)."""
        F = self.F
        out = np.zeros((thin.shape[0], len(dst)), dtype=np.int64)
        for rl, cl, M in blks:
            ri = [src[a] for a in rl]
            ci = [dst[b] for b in cl]
            sub = thin[:, ri]
            if not sub.any():
                continue
            out[:, ci] = (out[:, ci] + F.mm(sub, M)) % self.p
        return out

    def start_rows(self, rows, src, blk):
        """Rows of a block matrix (as thin matrix) at the given row labels."""
        out = np.zeros((len(rows), len(src)), dtype=np.int64)
        pos = {r: t for t, r in enumerate(rows)}
        for rl, cl, M in blk:
            for a, lab in enumerate(rl):
                if lab in pos:
                    out[pos[lab], [src[b] for b in cl]] = M[a]
        return out

    # ---------- step 1 ----------
    def step1(self, X, Wp, E):
        """Pentagon (X, V, W', V, E): returns dict R -> G[(X,V,R,E)] for R in W' (x) V."""
        A, B, C, D = X, self.gen, Wp, self.gen
        L = self.labels(A, B, C, D, E)
        I0, I1, I2, J1, J2 = (L[k] for k in ("I0", "I1", "I2", "J1", "J2"))
        if not I1:
            return {}
        rows = list(I1)                                   # all of M2's rows
        thin = self.start_rows(rows, I0, self.blocks("M1", A, B, C, D, E, inverse=True))
        thin = self.times(thin, I0, J1, self.blocks("N1", A, B, C, D, E))
        thin = self.times(thin, J1, J2, self.blocks("N2", A, B, C, D, E))
        M2 = self.times(thin, J2, I2, self.blocks("N3", A, B, C, D, E))
        out = {}
        for R in self.alc:
            n = self.N(C, D, R)
            if not n:
                continue
            cps = []
            for k in range(n):
                ri = [i for i, a in enumerate(rows) if a[2] == R and a[3] == k]
                ci = [I2[b] for b in I2 if b[0] == R and b[1] == k]
                if not ri:
                    continue
                cps.append((ri, ci))
                other = [I2[b] for b in I2 if not (b[0] == R and b[1] == k)]
                if M2[np.ix_(ri, other)].any():
                    self.checks["offdiag"] += 1
                    raise AssertionError(("off-diagonal block", X, Wp, E, R, k))
            if not cps:
                continue
            blks = [M2[np.ix_(ri, ci)] for ri, ci in cps]
            for b in blks[1:]:
                assert np.array_equal(b, blks[0]), ("copies disagree", X, Wp, E, R)
                self.checks["copies"] += 1
            ri, ci = cps[0]
            rl = [(rows[i][0], rows[i][1], rows[i][4]) for i in ri]
            I2l = list(I2)
            cl = [(I2l[j][2], I2l[j][3], I2l[j][4]) for j in ci]
            out[R] = (rl, cl, blks[0])
        return out

    # ---------- step 2 ----------
    def step2(self, X, Yp, W, E):
        """Pentagon (X, Y', V, W, E): returns dict T -> G[(X,T,W,E)] for T in Y' (x) V."""
        A, B, C, D = X, Yp, self.gen, W
        L = self.labels(A, B, C, D, E)
        I0, I1, I2, J1, J2 = (L[k] for k in ("I0", "I1", "I2", "J1", "J2"))
        if not J1:
            return {}
        rows = list(J1)
        thin = self.start_rows(rows, I0, self.blocks("N1", A, B, C, D, E, inverse=True))
        thin = self.times(thin, I0, I1, self.blocks("M1", A, B, C, D, E))
        thin = self.times(thin, I1, I2, self.blocks("M2", A, B, C, D, E))
        N2 = self.times(thin, I2, J2, self.blocks("N3", A, B, C, D, E, inverse=True))
        out = {}
        J2l = list(J2)
        for T in self.alc:
            n = self.N(B, C, T)
            if not n:
                continue
            cps = []
            for x in range(n):
                ri = [i for i, a in enumerate(rows) if a[0] == T and a[1] == x]
                ci = [J2[b] for b in J2 if b[0] == T and b[1] == x]
                if not ri:
                    continue
                cps.append((ri, ci))
                other = [J2[b] for b in J2 if not (b[0] == T and b[1] == x)]
                if N2[np.ix_(ri, other)].any():
                    raise AssertionError(("off-diagonal block", X, Yp, W, E, T, x))
                self.checks["offdiag"] += 1
            if not cps:
                continue
            blks = [N2[np.ix_(ri, ci)] for ri, ci in cps]
            for b in blks[1:]:
                assert np.array_equal(b, blks[0]), ("copies disagree", X, Yp, W, E, T)
                self.checks["copies"] += 1
            ri, ci = cps[0]
            rl = [(rows[i][3], rows[i][2], rows[i][4]) for i in ri]      # (Q, b, i3)
            cl = [(J2l[j][2], J2l[j][3], J2l[j][4]) for j in ci]         # (S, c, n)
            out[T] = (rl, cl, blks[0])
        return out

    def compare(self, key, val):
        """Compare a recomputed block with a stored one (same label sets, maybe permuted)."""
        rl, cl, M = self.G[key]
        rl2, cl2, M2 = val
        ri = {a: i for i, a in enumerate(rl2)}
        ci = {b: j for j, b in enumerate(cl2)}
        assert set(rl) == set(rl2) and set(cl) == set(cl2), ("label sets differ", key)
        P = M2[np.ix_([ri[a] for a in rl], [ci[b] for b in cl])]
        ok = np.array_equal(P % self.p, M % self.p)
        self.checks["known"] += 1
        return ok

    def canonical(self, val):
        """Sort labels by (alcove index, multiplicity indices)."""
        rl, cl, M = val
        key = lambda a: (self.alc.index(a[0]),) + tuple(a[1:])
        ro = sorted(range(len(rl)), key=lambda i: key(rl[i]))
        co = sorted(range(len(cl)), key=lambda j: key(cl[j]))
        return [rl[i] for i in ro], [cl[j] for j in co], M[np.ix_(ro, co)]
