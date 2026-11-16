"""R-symbols, quantum dimensions and the braid group action at the base 1, modulo p.

Splitting convention: c_{X,Y} o psi^{XY}_{Z,i} = sum_j Rs[(X,Y,Z)][i, j] psi^{YX}_{Z,j}.
At the base 1, (iota_X (x) iota_Y) psi^{XY}_{Z,i} is the path gamma_X followed by the i-th basis
vector of Hom(Z, X (x) Y) inside P_X(Z, y); naturality gives
  c_{V^x, V^y} (gamma_X . b_i) = sum_j Rs[i, j] gamma_Y . b'_j,
and c_{V^x, V^y} is a positive braid word in the sigma_i.
"""
import numpy as np


def bundle_braid(x, y):
    """sigma-word (in order of application) for c_{V^x, V^y} on x + y strands."""
    word = []
    # c_{V^x,V^y} = (c_{V,V^y} (x) 1^{x-1}) o (1 (x) c_{V^{x-1},V^y})
    def rec(off, x):
        if x == 0:
            return
        rec(off + 1, x - 1)                       # applied first: strands off+1 .. off+x-1+y
        # c_{V,V^y} on strands off .. off+y : sigma_off, ..., sigma_{off+y-1}
        for i in range(off, off + y):
            word.append(("sigma", i))
    rec(1, x)
    return word


class RSymbols:
    def __init__(self, O):
        self.O, self.M, self.p, self.F = O, O.M, O.p, O.F

    def path_vector(self, X, Y, Z, rows):
        """gamma_X . b for rows b of the basis of Hom(Z, X (x) Y) (in P_X(Z, y)), as columns of
        the Z-block of P_1(Z, x + y)."""
        O, M = self.O, self.M
        unit = O.unit
        x, y = O.length[X], O.length[Y]
        Pxy = M.P(unit, x + y)
        bZ = Pxy.blk[Z]
        Py = M.P(X, y)
        _, bb = O.basis(X, Y, Z)
        gX = O.gamma[X]
        out = np.zeros((bZ[1] - bZ[0], rows.shape[0]), dtype=np.int64)
        for a in range(rows.shape[1]):
            vs, es = Py.paths[bb[0] + a]
            gi = Pxy.index[(gX + vs[1:], tuple([0] * x) + es)] - bZ[0]
            out[gi, :] = rows[:, a]
        return out

    def R(self, X, Y, Z):
        O, F, p = self.O, self.F, self.p
        n = O.expected(X, Y, Z)
        if not n:
            return None
        x, y = O.length[X], O.length[Y]
        B1, _ = O.basis(X, Y, Z)
        B2, _ = O.basis(Y, X, Z)
        v = self.path_vector(X, Y, Z, B1)
        w, _ = self.M.wordT(O.unit, x + y, Z, bundle_braid(x, y), v)
        tgt = self.path_vector(Y, X, Z, B2)            # columns: gamma_Y . b'_j
        _, piv = F.rref(tgt.T.copy())
        sol = F.mm(w[piv, :].T, F.inverse(tgt[piv, :]).T)   # w = tgt sol^T  -> rows i, cols j
        assert np.array_equal(F.mm(tgt, sol.T), w % p), ("braided vector not in Y (x) X", X, Y, Z)
        return sol                                      # Rs[i, j]

    def qdim(self, X):
        """Trace of p_X at the base 1 with nested cups and caps (cap normalised by the zigzag)."""
        O, M, p = self.O, self.M, self.p
        s = O.length[X]
        unit = O.unit
        z = self.O.c["zigzag"][0]
        v = np.ones((1, 1), dtype=np.int64)
        n = 0
        # nested cups: cup at vertex n/2... build 1 -> V^{2s} as X-part (first s) and its dual (last s)
        for k in range(s):
            v, n = M.wordT(unit, n, unit, [("cup", k)], v)
        v = O.project(unit, X, unit, v, n=n)
        for k in range(s):
            v, n = M.wordT(unit, n, unit, [("cap", s - 1 - k)], v)
        zi = pow(int(z), -1, p)
        return int(v[0, 0]) * pow(zi, s, p) % p
