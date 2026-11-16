"""Dense linear algebra modulo a prime p < 2^30 with numpy int64 arrays (entries in [0, p))."""
import numpy as np


class Fp:
    def __init__(self, p):
        assert p < 2 ** 30
        self.p = p

    def inv(self, x):
        return pow(int(x) % self.p, -1, self.p)

    def mm(self, A, B):
        """A @ B mod p (inner dimension < 2^17)."""
        p = self.p
        A = np.asarray(A, dtype=np.int64)
        B = np.asarray(B, dtype=np.int64)
        if A.size == 0 or B.size == 0:
            return np.zeros((A.shape[0], B.shape[1]) if B.ndim == 2 else (A.shape[0],), dtype=np.int64)
        assert A.shape[-1] < 2 ** 17
        hi = (A @ (B >> 15)) % p
        lo = (A @ (B & 0x7FFF)) % p
        return ((hi << 15) + lo) % p

    def rref(self, M):
        """Reduced row echelon form: (nonzero rows, pivot columns)."""
        p = self.p
        M = np.array(M, dtype=np.int64) % p
        rows, cols = M.shape
        piv = []
        r = 0
        for c in range(cols):
            if r == rows:
                break
            nz = np.nonzero(M[r:, c])[0]
            if len(nz) == 0:
                continue
            k = r + nz[0]
            if k != r:
                M[[r, k]] = M[[k, r]]
            M[r] = M[r] * self.inv(M[r, c]) % p
            col = M[:, c].copy()
            col[r] = 0
            idx = np.nonzero(col)[0]
            if len(idx):
                M[idx] = (M[idx] - col[idx, None] * M[r][None, :]) % p
            piv.append(c)
            r += 1
        return M[:r], piv

    def rank(self, M):
        return len(self.rref(M)[1])

    def nullspace(self, M):
        """Columns: reduced echelon basis of the right kernel of M."""
        p = self.p
        n = M.shape[1]
        if M.shape[0] == 0:
            return np.eye(n, dtype=np.int64)
        R, piv = self.rref(M)
        free = [c for c in range(n) if c not in set(piv)]
        N = np.zeros((n, len(free)), dtype=np.int64)
        for j, c in enumerate(free):
            N[c, j] = 1
            for i, pc in enumerate(piv):
                N[pc, j] = (-R[i, c]) % p
        return N

    def inverse(self, M):
        n = M.shape[0]
        R, piv = self.rref(np.concatenate([M % self.p, np.eye(n, dtype=np.int64)], axis=1))
        assert piv[:n] == list(range(n)) and len(piv) >= n, "singular matrix"
        return R[:n, n:]

    def solve(self, A, B):
        """X with A X = B, A of full column rank; asserts consistency."""
        p = self.p
        m, n = A.shape
        R, piv = self.rref(np.concatenate([A % p, B % p], axis=1))
        assert piv[:n] == list(range(n)) and (len(piv) == n), "inconsistent or rank deficient"
        return R[:n, n:]

    def solve_check(self, A, B):
        """Like solve but returns None if inconsistent; A must have full column rank."""
        p = self.p
        m, n = A.shape
        R, piv = self.rref(np.concatenate([A % p, B % p], axis=1))
        assert piv[:n] == list(range(n)), "rank deficient"
        if len(piv) > n:
            return None
        return R[:n, n:]

    def left_inverse(self, A):
        """L with L A = I for A of full column rank (rows of A selected by rref of A^T)."""
        m, n = A.shape
        R, piv = self.rref(A.T.copy())      # pivot columns of A^T = independent rows of A
        assert len(piv) == n, "not full column rank"
        sub = A[piv, :]
        Li = self.inverse(sub)
        L = np.zeros((n, m), dtype=np.int64)
        L[:, piv] = Li
        return L


def find_primes(order, count, below=2 ** 30):
    """Primes p = 1 mod order below the bound, largest first."""
    import sympy
    out = []
    p = below - (below % order) + 1
    while len(out) < count:
        p -= order
        if sympy.isprime(p):
            out.append(p)
    return out


def primitive_root_of_unity(p, order, which=1):
    """A primitive order-th root of unity in F_p (deterministic: g^((p-1)/order) for the least
    generator g), raised to the power 'which'."""
    import sympy
    g = sympy.primitive_root(p)
    z = pow(g, (p - 1) // order, p)
    return pow(z, which, p)
