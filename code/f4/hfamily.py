"""F-symbols with third leg V, G[(X, Y, V, T)], in the path model, modulo p.

Splitting convention: for the bases of Hom(T, X (x) Y (x) V)
   L_(R,i,e)  = (psi^{XY}_{R,i} (x) 1) psi^{RV}_{T,e}           (concatenation of paths)
   R_(S,e',l) = (1 (x) psi^{YV}_{S,e'}) psi^{XS}_{T,l}
we compute G with L_a = sum_b G[a,b] R_b.  (The fusion-map F of export format is G^{-T}.)

psi^{YV}_{S,e'} is an edge of the path model.  To apply 1_X (x) psi^{YV}_{S,e'} to an element of
Hom(T, X (x) S) = image of 1 (x) p_S, we need the morphism
   Phi_{e'} = (iota_Y (x) 1) o e' o pi_S : V^s -> V^{y+1}
as an operator at every base.  It is a combination of compressed words
   (p_Y (x) 1) w p_S,   w a word in braidings, splits and cups,
normalised at the base 1, where Phi_{e'} gamma_S = gamma_Y . e'.
"""
import numpy as np


class HFamily:
    def __init__(self, O, seed=7):
        self.O, self.M, self.p, self.F = O, O.M, O.p, O.F
        self.alc = O.alc
        self.rng = np.random.default_rng(seed)
        self._phi = {}

    # ---------- words ----------
    def random_word(self, s, t):
        """A random word V^s -> V^t, t - s in {-1,...,2}: braid(s) . changer . braid(t)."""
        rng = self.rng
        w = []
        for _ in range(2 * s):
            if s >= 2:
                w.append(("sigma", int(rng.integers(1, s))))
        n = s
        tri = bool(self.O.c["split"])          # is there a trivalent vertex V -> V (x) V ?
        while n < t:
            if n == 0 or not tri or (t - n >= 2 and rng.random() < 0.5):
                w.append(("cup", int(rng.integers(0, n + 1)))); n += 2
            else:
                w.append(("split", int(rng.integers(1, n + 1)))); n += 1
        while n > t:
            if tri:
                w.append(("merge", int(rng.integers(1, n)))); n -= 1
            else:
                w.append(("cap", int(rng.integers(0, n - 1)))); n -= 2
        assert n == t
        for _ in range(2 * t):
            if t >= 2:
                w.append(("sigma", int(rng.integers(1, t))))
        return w

    def apply_word_projected(self, X, T, Y, word, x, s):
        """(p_Y (x) 1) word applied to columns x in the T-block of P_X(s)."""
        y, n = self.M.wordT(X, s, T, word, x)
        return self.O.project(X, Y, T, y, n=n)

    def phi(self, Y, S):
        """List over edges e' (Y -> S) of lists of (coefficient, word)."""
        key = (Y, S)
        if key in self._phi:
            return self._phi[key]
        O, p, F = self.O, self.p, self.F
        unit = O.unit
        y, s = O.length[Y], O.length[S]
        N = O.N(Y, S)
        P1 = self.M.P(unit, y + 1)
        b = P1.blk[S]
        gS = O.gamma[S]
        # gamma_S as a path vector in P_1(S, s): all edges along gamma are the single edge 0
        PS = self.M.P(unit, s)
        g_idx = PS.index[(gS, tuple([0] * s))]
        bS = PS.blk[S]
        x0 = np.zeros((bS[1] - bS[0], 1), dtype=np.int64)
        x0[g_idx - bS[0], 0] = 1
        # target paths gamma_Y . e''
        gY = O.gamma[Y]
        tgt = [P1.index[(gY + (S,), tuple([0] * y) + (e,))] - b[0] for e in range(N)]
        words, vals = [], []
        tries = 0
        while True:
            tries += 1
            assert tries < 200, ("no spanning words", Y, S)
            w = self.random_word(s, y + 1)
            v = self.apply_word_projected(unit, S, Y, w, x0, s)[:, 0]
            rest = v.copy(); rest[tgt] = 0
            assert not rest.any(), "compressed word leaves the image of p_Y (x) 1 at the base 1"
            cand = vals + [v[tgt]]
            if F.rank(np.array(cand).T) == len(cand):
                words.append(w); vals.append(v[tgt])
            if len(words) == N:
                break
        A = np.array(vals).T % p              # column k = coordinates of word k
        Ainv = F.inverse(A)                   # Phi_e = sum_k Ainv[k, e] word_k
        out = [[(int(Ainv[k, e]), words[k]) for k in range(N)] for e in range(N)]
        self._phi[key] = out
        return out

    def apply_phi(self, X, T, Y, S, e, x):
        """(1_X (x) Phi_e) applied to columns x in the T-block of P_X(s) (x in image of p_S)."""
        s = self.O.length[S]
        acc = None
        for c, w in self.phi(Y, S)[e]:
            if c == 0:
                continue
            y = self.apply_word_projected(X, T, Y, w, x, s)
            acc = y * c % self.p if acc is None else (acc + y * c) % self.p
        return acc

    # ---------- the F-matrix ----------
    def G(self, X, Y, T):
        """Splitting F-matrix G[(X, Y, V, T)] with row labels (R,i,e), column labels (S,e',l)."""
        O, M, p, F = self.O, self.M, self.p, self.F
        y = O.length[Y]
        P = M.P(X, y + 1)
        if T not in P.blk:
            return None
        bT = P.blk[T]
        dimT = bT[1] - bT[0]
        Py = M.P(X, y)
        rows, rl = [], []
        for R in self.alc:
            if not O.expected(X, Y, R) or not O.N(R, T):
                continue
            Bm, bR = O.basis(X, Y, R)
            for i in range(Bm.shape[0]):
                for e in range(O.N(R, T)):
                    v = np.zeros(dimT, dtype=np.int64)
                    for a in np.nonzero(Bm[i])[0]:
                        vs, es = Py.paths[bR[0] + a]
                        v[P.index[(vs + (T,), es + (e,))] - bT[0]] = Bm[i, a]
                    rows.append(v); rl.append((R, i, e))
        if not rows:
            return None
        cols, cl = [], []
        for S in self.alc:
            if not O.N(Y, S) or not O.expected(X, S, T):
                continue
            Bm, _ = O.basis(X, S, T)
            for e in range(O.N(Y, S)):
                img = self.apply_phi(X, T, Y, S, e, Bm.T.copy())
                for l in range(Bm.shape[0]):
                    cols.append(img[:, l]); cl.append((S, e, l))
        Lm = np.array(rows) % p
        Rm = np.array(cols) % p
        assert Lm.shape[0] == Rm.shape[0], (X, Y, T, Lm.shape, Rm.shape)
        # G: Lm = G Rm.  Use pivot columns of Rm.
        _, piv = F.rref(Rm)
        assert len(piv) == Rm.shape[0], ("right basis not independent", X, Y, T)
        Gm = F.mm(Lm[:, piv], F.inverse(Rm[:, piv]))
        assert np.array_equal(F.mm(Gm, Rm), Lm), ("left basis not in the span of the right basis", X, Y, T)
        return rl, cl, Gm
