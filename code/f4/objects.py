"""Simple objects as idempotents in End(V^n), multiplicity spaces, and their gauge.

Each simple S gets a spectral path gamma_S = (1, V, Z_2, ..., Z_s = S) from the unit: every step
Z_{k-1} -> Z_k has multiplicity one and the twist of Z_k differs from the twists of all other
constituents of Z_{k-1} (x) V.  The idempotent is
    p_S = prod_k  prod_{t in twists(Z_{k-1} (x) V), t != theta_{Z_k}} (Theta_k - t)/(theta_{Z_k} - t)
with Theta_k = theta_V^k (sigma_1 ... sigma_{k-1})^k the ribbon twist of V^k.  At the base 1 it
is the projection on the single path gamma_S.  (The factor 1 - cup cap, the unit projection
of V (x) V removed, is used at the second step since theta_(0,0,0,2) = theta_unit for F4.)

Gauge: the basis of Hom(T, X (x) S) is the reduced echelon basis of the image of 1_X (x) p_S
on the path space P_X(T, s) (paths ordered as in pathmodel.Paths).
"""
import numpy as np
from pathmodel import Model
from lin import Fp


class Objects:
    def __init__(self, cells, h, fusion):
        """h: dict weight -> twist exponent e (theta = q^e); fusion: (X,S,T) -> N."""
        self.fusion = fusion
        self.c = cells
        self.M = Model(cells)
        self.p = cells["p"]
        self.F = Fp(self.p)
        self.alc = cells["alcove"]
        self.unit = self.alc[0]
        q, order, p = cells["q"], cells["order"], self.p
        self.theta = {w: pow(q, e % order, p) for w, e in h.items()}
        self.nb = {}
        for (A, B), m in cells["edges"].items():
            self.nb.setdefault(A, {})[B] = m
        self.gamma = self._spectral_paths()
        self.length = {S: len(g) - 1 for S, g in self.gamma.items()}
        self._basis = {}

    def N(self, A, B):
        return self.nb.get(A, {}).get(B, 0)

    def _spectral_paths(self):
        gen = self.c["gen"]
        gamma = {self.unit: (self.unit,), gen: (self.unit, gen)}
        frontier = [gen]
        while frontier:
            new = []
            for Z in frontier:
                # at the second step the unit summand of V (x) V is removed by 1 - cup cap
                tw = [self.theta[B] for B in self.nb[Z] if not (Z == gen and B == self.unit)]
                for B, m in sorted(self.nb[Z].items(), key=lambda t: self.alc.index(t[0])):
                    if B in gamma or m != 1:
                        continue
                    if tw.count(self.theta[B]) == 1:
                        gamma[B] = gamma[Z] + (B,)
                        new.append(B)
            frontier = new
        assert set(gamma) == set(self.alc), set(self.alc) - set(gamma)
        return gamma

    def ribbon(self, X, n, T, k, x):
        """Apply Theta_k (first k strands) to columns x in the T-block of P_X(n)."""
        p = self.p
        gen = self.c["gen"]
        word = [("sigma", i) for _ in range(k) for i in range(1, k)]
        y, _ = self.M.wordT(X, n, T, word, x)
        return y * pow(int(self.theta[gen]), k, p) % p

    def project(self, X, S, T, x, n=None, start=2):
        """Apply p_S (on the first s strands) to columns x in the T-block of P_X(n)."""
        p, F = self.p, self.F
        g = self.gamma[S]
        if n is None:
            n = len(g) - 1
        for k in range(max(start, 2), len(g)):
            Z0, Z1 = g[k - 1], g[k]
            target = self.theta[Z1]
            cons = [B for B in self.nb[Z0] if B != Z1]
            if k == 2 and Z1 != self.unit:
                # remove the unit summand of V (x) V: cup o cap (raw cap, cap(cup) = 1) on strands 1, 2
                y, _ = self.M.wordT(X, n, T, [("cap", 0), ("cup", 0)], x)
                x = (x - y) % p
                cons = [B for B in cons if B != self.unit]
            others = sorted({self.theta[B] for B in cons})
            assert target not in others
            for t in others:
                y = self.ribbon(X, n, T, k, x)
                x = (y - t * x) % p * F.inv((target - t) % p) % p
        return x

    def basis(self, X, S, T):
        """Rows: reduced echelon basis of Hom(T, X (x) S) inside P_X(T, s), as a pair
        (row matrix in the coordinates of the T-block, (start, stop) of that block)."""
        key = (X, S, T)
        if key not in self._basis:
            s = self.length[S]
            P = self.M.P(X, s)
            b = P.blk.get(T)
            if b is None:
                self._basis[key] = (np.zeros((0, 0), dtype=np.int64), b)
                return self._basis[key]
            # the image has dimension N(X S -> T) <= 10: project a few random vectors
            rng = np.random.default_rng(hash(key) % 2 ** 32)
            kk = self.expected(X, S, T) + 3
            x = rng.integers(0, self.p, size=(b[1] - b[0], kk), dtype=np.int64)
            y = self.project(X, S, T, x)                 # columns span the image
            R, piv = self.F.rref(y.T.copy())
            assert R.shape[0] == self.expected(X, S, T), (X, S, T, R.shape, self.expected(X, S, T))
            self._basis[key] = (R, b)
        return self._basis[key]

    def expected(self, X, S, T):
        return self.fusion.get((X, S, T), 0)

    def dim(self, X, S, T):
        return self.basis(X, S, T)[0].shape[0]
