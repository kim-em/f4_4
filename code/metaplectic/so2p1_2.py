"""F- and R-symbols of the metaplectic categories (fusion rules of so(2p+1)_2), exactly.

Source: E. Ardonne, P. E. Finch, M. Titsworth, "Classification of metaplectic fusion
categories", arXiv:1608.03762v2, Symmetry 13 (2021) 2102, Section 4 ("AFT" below; equation
names refer to the arXiv source).  The F- and R-symbols are transcribed as printed and need no
correction: "-1^{x}" is read as (-1)^x, the fusion rule phi_i (x) psi_pm as psi_+ + psi_-, and
the blocks AFT do not print are obtained from their rotation property (rotate), the remaining
1x1 blocks being (1) by their rule (O.W.).  Misprints in the surrounding equations (pentagon,
second hexagon, R gauge law, S-matrix) are listed in README.md.

Parameters (AFT notation): p >= 1, n = 2p + 1; r odd, 1 <= r < n, gcd(r, n) = 1; kappa = +-1;
lam = +1 (the R-symbols as printed) or -1 (their inverses).

Labels, in the AFT order:  '1', 'eps', 'phi1', ..., 'phi<p>', 'psi+', 'psi-'.
Dimensions 1, 1, 2, ..., 2, sqrt(n), sqrt(n) (with pivotal signs eps_psi = kappa).

Conventions (AFT, Section 2.2), multiplicity free:
  * F[(a,b,c,d)][(e,f)] = F_{abc}^{d;ef}, e in a (x) b, f in b (x) c: the change of basis from
    ((a b)_e c)_d to (a (b c)_f)_d.  All F-matrices are real orthogonal, so the splitting and
    the fusion-map conventions agree.  Pentagon (AFT (F3), with its summation variable renamed
    so that it is the one not occurring on the other side; README.md):
      F_{fcd}^{e;gl} F_{abl}^{e;fk} = sum_h F_{abc}^{g;fh} F_{ahd}^{e;gk} F_{bcd}^{k;hl}.
  * R[(a,b,c)] = R_{ab}^c; hexagons (AFT (R1), and (R1) for the reverse-inverse braiding):
      R_{ac}^e F_{acb}^{d;eg} R_{bc}^g = sum_f F_{cab}^{d;ef} R_{fc}^d F_{abc}^{d;fg}
      and the same with R_{xy}^z replaced by (R_{yx}^z)^-1 (AFT (R2) keeps the lower indices,
      which is equivalent here because R_{ab}^c = R_{ba}^c in AFT's gauge; README.md).
  * Gauge: under fusion-space rescalings g, F_{abc}^{d;ef} -> g_{ab}^e g_{ec}^d F / (g_{bc}^f g_{af}^d)
    and R_{ab}^c -> g_{ba}^c R_{ab}^c / g_{ab}^c (AFT print the inverse ratio; see README.md).

Everything is computed in Q(zeta_N), N = 8n, with zeta_N = exp(2 pi i/N); AFT's q is
exp(pi i/n) = zeta_N^4.
"""
import itertools
from fractions import Fraction

from cyclo import Field


def sg(k):
    """(-1)^k for an integer k."""
    return -1 if k % 2 else 1


def jacobi(a, n):
    assert n > 0 and n % 2 == 1
    a %= n
    t = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5):
                t = -t
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            t = -t
        a %= n
    return t if n == 1 else 0


class Metaplectic:
    def __init__(self, p, r=1, kappa=1, lam=1, field=None):
        n = 2 * p + 1
        assert p >= 1 and r % 2 == 1 and 1 <= r < n and __import__("math").gcd(r, n) == 1
        assert kappa in (1, -1) and lam in (1, -1)
        self.p, self.n, self.r, self.kappa, self.lam = p, n, r, kappa, lam
        self.K = K = field or Field(8 * n)
        assert K.N % (8 * n) == 0
        self.labels = ["1", "eps"] + [f"phi{i}" for i in range(1, p + 1)] + ["psi+", "psi-"]
        self.order = {x: k for k, x in enumerate(self.labels)}
        self.sqrt2 = K.sqrt2()
        self.sqrtn = K.gauss_sqrt(n)
        self.N = self._fusion()
        self.rotation_conflicts = []
        self.F = self._fmatrices()
        self.R = self._rsymbols()

    # ---------------------------------------------------------------- helpers
    def g(self, a):
        """AFT's g: Z_n -> {0, ..., p}, g(a) = |a|."""
        a %= self.n
        return min(a, self.n - a)

    def phi(self, i):
        return f"phi{i}"

    @staticmethod
    def idx(x):
        """i for 'phi<i>'; None otherwise."""
        return int(x[3:]) if x.startswith("phi") else None

    @staticmethod
    def psisign(x):
        return {"psi+": 1, "psi-": -1}.get(x)

    def q(self, t):
        """AFT's q^t, q = exp(pi i/n)."""
        return self.K.root_of_unity(t, 2 * self.n)

    def re_q(self, t):
        return (self.q(t) + self.q(-t)) * Fraction(1, 2)

    def im_q(self, t):
        return (self.q(t) - self.q(-t)) * self.K.root_of_unity(-1, 4) * Fraction(1, 2)

    # ---------------------------------------------------------------- fusion rules (AFT 4.1)
    def _fusion(self):
        p = self.p
        L = self.labels
        prod = {}
        for a in L:
            prod[("1", a)] = prod[(a, "1")] = [a]
        prod[("eps", "eps")] = ["1"]
        for i in range(1, p + 1):
            prod[("eps", self.phi(i))] = prod[(self.phi(i), "eps")] = [self.phi(i)]
            for j in range(1, p + 1):
                if i == j:
                    prod[(self.phi(i), self.phi(i))] = ["1", "eps", self.phi(self.g(2 * i))]
                else:
                    prod[(self.phi(i), self.phi(j))] = [self.phi(self.g(i - j)), self.phi(self.g(i + j))]
            for s in ("psi+", "psi-"):
                # AFT print "phi_i (x) psi_pm = psi_pm (x) psi_mp"; the intended rule is psi_+ + psi_-
                prod[(self.phi(i), s)] = prod[(s, self.phi(i))] = ["psi+", "psi-"]
        phis = [self.phi(j) for j in range(1, p + 1)]
        prod[("eps", "psi+")] = prod[("psi+", "eps")] = ["psi-"]
        prod[("eps", "psi-")] = prod[("psi-", "eps")] = ["psi+"]
        prod[("psi+", "psi+")] = prod[("psi-", "psi-")] = ["1"] + phis
        prod[("psi+", "psi-")] = prod[("psi-", "psi+")] = ["eps"] + phis
        N = {}
        for (a, b), cs in prod.items():
            for c in cs:
                N[(a, b, c)] = N.get((a, b, c), 0) + 1
        assert max(N.values()) == 1
        return N

    def n_(self, a, b, c):
        return self.N.get((a, b, c), 0)

    def rows(self, a, b, c, d):
        return [e for e in self.labels if self.n_(a, b, e) and self.n_(e, c, d)]

    def cols(self, a, b, c, d):
        return [f for f in self.labels if self.n_(b, c, f) and self.n_(a, f, d)]

    def blocks(self):
        for key in itertools.product(self.labels, repeat=4):
            if self.rows(*key):
                yield key

    # ---------------------------------------------------------------- AFT 4.2 matrices
    def A(self, s1, s2, s3, s4):
        h = self.sqrt2 * Fraction(1, 2)
        return [[h * s1, h * s2], [h * s3, h * s4]]

    def D(self, t, s1, s2, s3, s4):
        c, s = self.re_q(self.r * t), self.im_q(self.r * t)
        return [[c * s1, s * s2], [s * s3, c * s4]]

    def E(self, t, s1, s2, s3, s4):
        c, s = self.re_q(self.r * t), self.im_q(self.r * t)
        return [[s * s1, c * s2], [c * s3, s * s4]]

    def J(self, i, j, part):
        """Re or Im of AFT's J(i,j;r,kappa) = 2^zeta(i,j) kappa q^(r i j)/sqrt(n), zeta = (2 - d_i0 - d_j0)/2."""
        z2 = (i != 0) + (j != 0)                      # 2 zeta
        pref = self.sqrt2 if z2 == 1 else self.K.from_int(1 if z2 == 0 else 2)
        pref = pref * self.kappa / self.sqrtn
        return pref * (self.re_q(self.r * i * j) if part == "re" else self.im_q(self.r * i * j))

    def Gmat(self):
        p = self.p
        return [[self.J(i, j, "im") * sg((i - 1) * (j - 1)) for j in range(1, p + 1)] for i in range(1, p + 1)]

    def Hmat(self, prime=False):
        p = self.p
        H = [[self.J(i, j, "re") * sg(i * j) for j in range(p + 1)] for i in range(p + 1)]
        if prime:
            H = [[H[i][j] * sg((i == 0) + (j == 0) + 1) for j in range(p + 1)] for i in range(p + 1)]
        return H

    @staticmethod
    def neg(M, s=-1):
        return [[x * s for x in row] for row in M]

    @staticmethod
    def T(M):
        return [list(r) for r in zip(*M)]

    # ---------------------------------------------------------------- the F-symbols (AFT 4.2.2)
    def explicit(self, a, b, c, d):
        """The matrix F_{abc}^d if AFT give it explicitly, else None (rows/cols in label order)."""
        K, p = self.K, self.p
        one = K.from_int(1)
        m1 = K.from_int(-1)
        if "1" in (a, b, c, d):
            return [[one]]
        i, j, k, l = (self.idx(x) for x in (a, b, c, d))
        sa, sb, sc, sd = (self.psisign(x) for x in (a, b, c, d))
        E_ = "eps"
        # Z_2 part and the action of eps on phi
        if (a, b, c, d) == (E_, E_, E_, E_):
            return [[one]]
        if a == E_ and j and c == E_ and l == j:
            return [[one]]
        if a == E_ and b == E_ and k and l == k:
            return [[m1]]
        if a == E_ and j and k and l:
            ii, jj, kk = j, k, l                     # AFT: F_{eps phi_i phi_j}^{phi_k}
            if jj <= ii:
                return [[K.from_int(sg(jj))]]
            if kk == self.g(ii + jj):
                return [[K.from_int(sg(jj))]]
            assert kk == self.g(ii - jj)
            return [[K.from_int(sg(jj - 1))]]
        # the rest of Rep(D_{2p+1})
        if i and j and k and l:
            h = Fraction(1, 2)
            s2 = self.sqrt2 * Fraction(1, 2)
            if i == j == k == l:
                return [[one * h, m1 * h, s2], [m1 * h, one * h, s2], [s2, s2, K.zero]]   # C
            if i == k and j == l and i != j:
                return [[K.zero, one], [one, K.zero]]                                      # B
            if i == j and k == l and i != k:
                jj, ii = k, i
                return self.A(1, 1, sg(jj - ii + 1), sg(jj - ii))
            return None            # "(O.W.)" = (1), filled in after rotation
        # four psi
        if sa and sb and sc and sd:
            if sa == sb == sc == sd:
                return self.Hmat()
            if sa == sc and sb == sd and sa != sb:
                return self.neg(self.Hmat())
            if sa == sb and sc == sd and sa != sc:
                return self.Hmat(prime=True)
            if sa == sb == sc and sd == -sa:
                return self.neg(self.Gmat(), sa)
            return None
        # two psi: phi psi phi psi
        if i and sb and k and sd:
            ii, jj = i, k
            if sd == sb:
                return self.neg(self.D(ii * jj, -1, sg(ii + jj), sg(ii + jj), 1), -sb * sg(ii * jj))
            return self.neg(self.E(ii * jj, sg(ii + jj), 1, 1, sg(ii + jj + 1)), -sg(ii * jj))
        # two psi: phi phi psi psi
        if i and j and sc and sd:
            ii, jj, s = i, j, sc
            if sd == sc:
                if ii == jj:
                    return self.A(1, 1, 1, -1)
                if (ii - jj) % 2 == 0:
                    return self.A(s, s, 1, -1)
                return self.A(1, -1, s, s)
            if ii == jj:
                return self.A(-1, -1, sg(ii), sg(ii + 1))
            if ii < jj and sc == 1:
                return self.A(sg(jj + 1), sg(ii + 1), sg(jj), sg(ii + 1))
            if ii < jj and sc == -1:
                return self.A(sg(ii), sg(jj), sg(ii), sg(jj + 1))
            if ii > jj and sc == 1:
                return self.A(sg(jj), sg(ii), sg(jj), sg(ii + 1))
            if ii > jj and sc == -1:
                return self.A(sg(ii + 1), sg(jj + 1), sg(ii), sg(jj + 1))
        # two psi with eps
        if a == E_:
            if j and sc and sd:
                return [[one]]                              # F_{eps phi psi}^{psi'}
            if sb and k and sd:
                return [[one]] if sb == sd else [[m1]]      # F_{eps psi phi}^{psi'}
            if sb and sc and l:
                return [[one]]                              # F_{eps psi psi'}^{phi}
            if b == E_ and sc and sd:
                return [[m1]]                               # F_{eps eps psi}^{psi}
            if sb and c == E_ and sd:
                return [[m1]]                               # F_{eps psi eps}^{psi}
        return None

    def _fmatrices(self):
        F = {}
        src = {}
        for key in self.blocks():
            M = self.explicit(*key)
            if M is not None:
                F[key] = self._todict(key, M)
                src[key] = "explicit"
        # close under AFT (rotate): F_{bcd}^a = (F_{abc}^d)^T
        changed = True
        while changed:
            changed = False
            for (a, b, c, d) in list(F):
                rk = (b, c, d, a)
                rot = {(f, e): v for (e, f), v in F[(a, b, c, d)].items()}
                if rk in F:
                    if not all(F[rk][x] == rot[x] for x in rot) and ((a, b, c, d), rk) not in self.rotation_conflicts:
                        self.rotation_conflicts.append(((a, b, c, d), rk))
                else:
                    F[rk] = rot
                    src[rk] = "rotation"
                    changed = True
        # (O.W.) = (1) for what is left
        for key in self.blocks():
            if key not in F:
                rows, cols = self.rows(*key), self.cols(*key)
                assert len(rows) == 1, ("unspecified block larger than 1x1", key)
                F[key] = {(rows[0], cols[0]): self.K.from_int(1)}
                src[key] = "O.W."
        self.source = src
        return F

    def _todict(self, key, M):
        rows, cols = self.rows(*key), self.cols(*key)
        assert len(M) == len(rows) and all(len(r) == len(cols) for r in M), (key, rows, cols, len(M))
        return {(e, f): M[x][y] for x, e in enumerate(rows) for y, f in enumerate(cols)}

    # ---------------------------------------------------------------- R-symbols (AFT 4.3)
    def h(self, a):
        p, n, r, kappa = self.p, self.n, self.r, self.kappa
        if a == "1":
            return Fraction(0)
        if a == "eps":
            return Fraction(1)
        if a.startswith("phi"):
            i = self.idx(a)
            return Fraction(r * i * (n - i), 2 * n)
        s_p = -sg(p * (p + 1) // 2)
        assert s_p == -jacobi(2, n)
        base = p + kappa * s_p - jacobi(n, r)
        return Fraction(r * (base + (2 if a == "psi+" else 6)), 8)

    def sigma1_rule(self, a, b, c):
        p = self.p
        i, j, k = self.idx(a), self.idx(b), self.idx(c)
        if i and self.psisign(b) and self.psisign(c):
            m = i % 4
            if (m in (1, 2) and p % 2 == 1) or (m in (2, 3) and p % 2 == 0):
                return -1
            return 1
        if i and j and k and k == self.g(i + j):
            return sg(i * j)
        return None

    def sigma2_rule(self, a, b, c):
        sb, sc = self.psisign(b), self.psisign(c)
        if (a == "eps" or self.idx(a)) and sb and sc and sb == -sc:
            return sg((self.r - 1) // 2)
        return None

    def _sigma(self, rule, a, b, c):
        """apply AFT's rotation property sigma_{bc}^a = sigma_{ab}^c; unspecified = +1."""
        vals = {rule(*t) for t in ((a, b, c), (b, c, a), (c, a, b))} - {None}
        assert len(vals) <= 1, ("inconsistent sign rule", a, b, c, vals)
        return vals.pop() if vals else 1

    def _rsymbols(self):
        K = self.K
        R = {}
        for (a, b, c) in self.N:
            x = self.h(a) + self.h(b) - self.h(c)
            # (-1)^x = exp(pi i x)
            assert (x * K.N / 2).denominator == 1
            v = K.root_of_unity(int(x * K.N / 2), K.N) * (self._sigma(self.sigma1_rule, a, b, c) * self._sigma(self.sigma2_rule, a, b, c))
            R[(a, b, c)] = v if self.lam == 1 else v.inverse()
        return R

    # ---------------------------------------------------------------- derived data
    def dim(self, a):
        """quantum dimension with the pivotal structure eps_psi = kappa (positive dimensions)."""
        return self.sqrtn if a.startswith("psi") else self.K.from_int(1 if a in ("1", "eps") else 2)

    def fentry(self, a, b, c, d, e, f):
        return self.F.get((a, b, c, d), {}).get((e, f), self.K.zero)


# -------------------------------------------------------------------- the quantum group case
def quantum_group_parameters(p):
    """(r, kappa, lam) and the label dictionary for C(so(2p+1), 2), p >= 2.

    C(so(2p+1), 2) here is the category with twists theta_X = exp(2 pi i h_X), h_X the conformal
    weight (X, X + 2 rho)/(2(2p+1)) (long roots of length^2 2), i.e. the quantum group at
    q = exp(2 pi i/(4(2p+1))) in the conventions of the paper (Section 2).  Its twists and S-matrix are
    those of AFT (r, kappa, lam) = (1, (2 | 2p+1), -1) with
        1 = 0, eps = 2 w_1, phi_i = w_i (i < p), phi_p = 2 w_p,
        spinor w_p = psi+ if kappa = 1 else psi-,  w_1 + w_p = the other psi.
    For p = 6 the braided equivalence with our quantum group computation is proved by an explicit
    gauge (compare_so13.py); for p = 2, ..., 6 the data agree exactly with SageMath's solution
    (compare_sage.py --aft); for every p this is the choice singled out by the twists (README.md).
    Returns (r, kappa, lam, dict Dynkin label string -> AFT label).
    """
    assert p >= 2
    n = 2 * p + 1
    kappa = jacobi(2, n)
    def w(*pairs):
        a = [0] * p
        for i, c in pairs:
            a[i - 1] += c
        return "".join(map(str, a))
    lab = {w(): "1", w((1, 2)): "eps", w((p, 2)): f"phi{p}"}
    for i in range(1, p):
        lab[w((i, 1))] = f"phi{i}"
    spin, other = ("psi+", "psi-") if kappa == 1 else ("psi-", "psi+")
    lab[w((p, 1))] = spin
    lab[w((1, 1), (p, 1))] = other
    return 1, kappa, -1, lab


def write_text(D, path, names=None):
    """Write F and R in a text format: coefficients in the power basis of zeta_N, N = 8(2p+1)."""
    names = names or {x: x for x in D.labels}
    K = D.K
    with open(path, "w") as fh:
        fh.write(f"# metaplectic category, AFT formulas (arXiv:1608.03762), p = {D.p}, r = {D.r}, kappa = {D.kappa}, lam = {D.lam}\n")
        fh.write(f"# field Q(zeta), zeta = exp(2 pi i/{K.N}), entries sum_k b_k zeta^k, k = 0 .. {K.d - 1}\n")
        fh.write("# F a b c d | e f | b_0 .. : F_{abc}^{d;ef} (e in a b, f in b c; zero entries omitted)\n")
        fh.write("# R a b c | b_0 .. : R_{ab}^c;  conventions as in so2p1_2.py and the paper, Section 3 (export convention)\n")
        fh.write("# objects: " + " ".join(names[x] for x in D.labels) + " ; unit: " + names["1"] + "\n")
        for (a, b, c) in sorted(D.N, key=lambda t: [D.order[x] for x in t]):
            fh.write(f"# N {names[a]} {names[b]} {names[c]} 1\n")
        for key in sorted(D.F, key=lambda t: [D.order[x] for x in t]):
            for (e, f), v in sorted(D.F[key].items(), key=lambda t: (D.order[t[0][0]], D.order[t[0][1]])):
                if not v.is_zero():
                    fh.write("F " + " ".join(names[x] for x in key) + f" | {names[e]} {names[f]} | " + " ".join(map(str, v.coeffs())) + "\n")
        for key in sorted(D.R, key=lambda t: [D.order[x] for x in t]):
            fh.write("R " + " ".join(names[x] for x in key) + " | " + " ".join(map(str, D.R[key].coeffs())) + "\n")


if __name__ == "__main__":
    import sys
    a = sys.argv[1:]
    if len(a) == 2 and a[0].isdigit():
        p = int(a[0])
        r, kappa, lam, lab = quantum_group_parameters(p)
        D = Metaplectic(p, r, kappa, lam)
        write_text(D, a[1], names={v: k for k, v in lab.items()})
    elif len(a) == 5:
        write_text(Metaplectic(int(a[0]), int(a[1]), int(a[2]), int(a[3])), a[4])
    else:
        print("usage: so2p1_2.py P OUT              (C(so(2p+1),2), Dynkin labels)\n"
              "       so2p1_2.py P R KAPPA LAM OUT  (any AFT parameters, AFT labels)")
