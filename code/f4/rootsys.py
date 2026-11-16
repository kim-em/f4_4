"""Root data, weights, Weyl dimensions and level-k fusion rules for simple Lie algebras.

Weights are tuples of Dynkin labels.  The invariant form is normalised so that SHORT roots
have squared length 2 (so long roots have length 4 in types B, C, F and 6 in G2); then
(lambda, alpha_i) = lambda_i * d_i with d_i = (alpha_i, alpha_i)/2 in {1, 2, 3}, and all
pairings of weights in the root lattice are integers.  In this normalisation the quantum
group parameter that gives the twists theta_lambda = exp(2 pi i h_lambda) of the affine
category at level k is q = exp(2 pi i / (2 m (k + h^vee))), where m = max d_i is the
lacing number: h_lambda = (lambda, lambda + 2 rho) / (2 m (k + h^vee)).
"""
from fractions import Fraction as Fr
from functools import lru_cache
import itertools


class RootSystem:
    def __init__(self, cartan_type):
        self.type = cartan_type
        X, n = cartan_type[0], int(cartan_type[1:])
        self.rank = n
        if X == "F" and n == 4:
            # Bourbaki: alpha1, alpha2 long; alpha3, alpha4 short
            A = [[2, -1, 0, 0], [-1, 2, -1, 0], [0, -2, 2, -1], [0, 0, -1, 2]]
            d = [2, 2, 1, 1]
        elif X == "G" and n == 2:
            # Bourbaki: alpha1 short, alpha2 long
            A = [[2, -3], [-1, 2]]
            d = [1, 3]
        elif X == "B":
            # alpha_1..alpha_{n-1} long, alpha_n short
            A = [[0] * n for _ in range(n)]
            for i in range(n):
                A[i][i] = 2
                if i + 1 < n:
                    A[i][i + 1] = -1
                    A[i + 1][i] = -1
            A[n - 2][n - 1] = -1
            A[n - 1][n - 2] = -2
            d = [2] * (n - 1) + [1]
        elif X == "A":
            A = [[2 if i == j else (-1 if abs(i - j) == 1 else 0) for j in range(n)] for i in range(n)]
            d = [1] * n
        else:
            raise ValueError(cartan_type)
        # A[i][j] = <alpha_j, alpha_i^vee> : Dynkin label i of alpha_j
        self.A = A
        self.d = d
        # symmetrised form on simple roots: (alpha_i, alpha_j) = d_i A[i][j]
        self.B = [[d[i] * A[i][j] for j in range(n)] for i in range(n)]
        assert all(self.B[i][j] == self.B[j][i] for i in range(n) for j in range(n))
        self.m = max(d)
        self._inverse_cartan()
        self.positive_roots = self._positive_roots()
        self.rho = tuple([1] * n)
        hi = max(self.positive_roots, key=lambda r: sum(self.root_coords(r)))
        self.theta = hi                                  # highest root (long)
        # comarks: theta^vee = sum a_i^vee alpha_i^vee ;  (lambda, theta^vee) = sum a_i^vee lambda_i
        c = self.root_coords(hi)
        self.comarks = tuple(Fr(c[i] * d[i], self.m) for i in range(n))
        assert all(x.denominator == 1 for x in self.comarks)
        self.comarks = tuple(int(x) for x in self.comarks)
        self.hdual = 1 + sum(self.comarks)

    def _inverse_cartan(self):
        n = self.rank
        M = [[Fr(self.A[i][j]) for j in range(n)] + [Fr(int(i == j)) for j in range(n)] for i in range(n)]
        for c in range(n):
            p = next(r for r in range(c, n) if M[r][c] != 0)
            M[c], M[p] = M[p], M[c]
            piv = M[c][c]
            M[c] = [x / piv for x in M[c]]
            for r in range(n):
                if r != c and M[r][c] != 0:
                    f = M[r][c]
                    M[r] = [x - f * y for x, y in zip(M[r], M[c])]
        self.Ainv = [row[n:] for row in M]

    def root_coords(self, lam):
        """Coordinates of the weight lam in the basis of simple roots (rationals)."""
        n = self.rank
        # lam_i = sum_j A[i][j] c_j
        return tuple(sum(self.Ainv[j][i] * lam[i] for i in range(n)) for j in range(n))

    def alpha(self, j):
        return tuple(self.A[i][j] for i in range(self.rank))

    def ip(self, lam, mu):
        """(lam, mu) with short roots of length 2."""
        c = self.root_coords(mu)
        return sum(c[j] * lam[j] * self.d[j] for j in range(self.rank))

    def casimir(self, lam):
        return self.ip(lam, tuple(l + 2 for l in lam))   # (lam, lam + 2 rho)

    def _positive_roots(self):
        n = self.rank
        simple = [self.alpha(j) for j in range(n)]
        roots = set(simple)
        frontier = list(simple)
        while frontier:
            new = []
            for r in frontier:
                for i in range(n):
                    s = tuple(r[k] - r[i] * self.A[k][i] for k in range(n))
                    if s not in roots and all(x >= 0 for x in self.root_coords(s)) and any(x != 0 for x in s):
                        roots.add(s)
                        new.append(s)
            frontier = new
        return sorted(roots, key=lambda r: sum(self.root_coords(r)))

    def coroot_pairing(self, lam, alpha):
        """<lam, alpha^vee> = 2 (lam, alpha) / (alpha, alpha)."""
        return Fr(2) * self.ip(lam, alpha) / self.ip(alpha, alpha)

    def weyl_dim(self, lam):
        num, den = Fr(1), Fr(1)
        lr = tuple(l + 1 for l in lam)
        for a in self.positive_roots:
            num *= self.ip(lr, a)
            den *= self.ip(self.rho, a)
        r = num / den
        assert r.denominator == 1
        return int(r)

    def level(self, lam):
        return sum(c * l for c, l in zip(self.comarks, lam))

    def alcove(self, k):
        n = self.rank
        out = []
        for lam in itertools.product(range(k + 1), repeat=n):
            if self.level(lam) <= k:
                out.append(lam)
        return sorted(out, key=lambda w: (self.level(w), w))

    def h(self, lam, k):
        """Conformal weight h_lambda (exact rational)."""
        return Fr(self.casimir(lam), 2 * self.m * (k + self.hdual))

    # ---------- weights of irreducible modules (Freudenthal) ----------
    def dominant_weights(self, lam):
        """Dominant weights mu <= lam."""
        n = self.rank
        out = {tuple(lam)}
        frontier = [tuple(lam)]
        while frontier:
            new = []
            for mu in frontier:
                for a in self.positive_roots:
                    nu = tuple(x - y for x, y in zip(mu, a))
                    if all(x >= 0 for x in nu) and nu not in out:
                        out.add(nu)
                        new.append(nu)
            frontier = new
        return out

    @lru_cache(maxsize=None)
    def dominant_multiplicities(self, lam):
        lam = tuple(lam)
        dom = sorted(self.dominant_weights(lam), key=lambda mu: -sum(self.root_coords(tuple(a - b for a, b in zip(mu, lam)))))
        # order by height of lam - mu increasing
        dom = sorted(dom, key=lambda mu: sum(self.root_coords(tuple(a - b for a, b in zip(lam, mu)))))
        mult = {}
        lr = tuple(l + 1 for l in lam)
        c_lam = self.ip(lr, lr)

        def m_of(nu):
            dn = self.dominant_rep(nu)
            return mult.get(dn, 0)
        for mu in dom:
            if mu == lam:
                mult[mu] = 1
                continue
            s = Fr(0)
            for a in self.positive_roots:
                j = 1
                while True:
                    nu = tuple(x + j * y for x, y in zip(mu, a))
                    if not self._le(nu, lam):
                        break
                    mm = m_of(nu)
                    s += mm * self.ip(nu, a)
                    j += 1
            mr = tuple(x + 1 for x in mu)
            den = c_lam - self.ip(mr, mr)
            val = 2 * s / den
            assert val.denominator == 1
            mult[mu] = int(val)
        return mult

    def _le(self, nu, lam):
        c = self.root_coords(tuple(a - b for a, b in zip(lam, nu)))
        return all(x >= 0 and x.denominator == 1 for x in c)

    def dominant_rep(self, mu):
        mu = list(mu)
        while True:
            i = next((i for i, x in enumerate(mu) if x < 0), None)
            if i is None:
                return tuple(mu)
            c = mu[i]
            mu = [mu[k] - c * self.A[k][i] for k in range(self.rank)]

    def orbit(self, mu):
        mu = tuple(mu)
        out = {mu}
        frontier = [mu]
        while frontier:
            new = []
            for x in frontier:
                for i in range(self.rank):
                    if x[i] != 0:
                        y = tuple(x[k] - x[i] * self.A[k][i] for k in range(self.rank))
                        if y not in out:
                            out.add(y)
                            new.append(y)
            frontier = new
        return out

    def character(self, lam):
        """dict weight -> multiplicity for L(lam)."""
        if not hasattr(self, "_chcache"):
            self._chcache = {}
        lam = tuple(lam)
        if lam in self._chcache:
            return self._chcache[lam]
        ch = {}
        for mu, m in self.dominant_multiplicities(tuple(lam)).items():
            for nu in self.orbit(mu):
                ch[nu] = m
        assert sum(ch.values()) == self.weyl_dim(lam)
        self._chcache[lam] = ch
        return ch

    # ---------- level-k fusion by the Kac-Walton / Racah-Speiser algorithm ----------
    def reflect_shifted(self, p, k):
        """p = rho-shifted weight.  Return (sign, dominant alcove weight) or (0, None)."""
        p = list(p)
        sign = 1
        kap = k + self.hdual if k is not None else None
        th = self.theta
        for _ in range(10000):
            if 0 in p:
                return 0, None
            i = next((i for i, x in enumerate(p) if x < 0), None)
            if i is not None:
                c = p[i]
                p = [p[r] - c * self.A[r][i] for r in range(self.rank)]
                sign = -sign
                continue
            if kap is not None:
                ex = self.level(p) - kap
                if ex == 0:
                    return 0, None
                if ex > 0:
                    # reflection in the affine wall: p -> p - ex * theta  (theta long, <theta,theta^vee>=2)
                    p = [x - ex * t for x, t in zip(p, th)]
                    sign = -sign
                    continue
            return sign, tuple(x - 1 for x in p)
        raise RuntimeError("reflection loop")

    def tensor(self, lam, mu, k=None):
        """Decomposition of L(lam) (x) L(mu) (classical if k is None, else level-k fusion)."""
        out = {}
        for nu, m in self.character(mu).items():
            p = tuple(a + b + 1 for a, b in zip(lam, nu))
            s, w = self.reflect_shifted(p, k)
            if s:
                out[w] = out.get(w, 0) + s * m
        out = {w: m for w, m in out.items() if m != 0}
        assert all(m > 0 for m in out.values()), (lam, mu, out)
        return out
