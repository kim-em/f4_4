"""Exact arithmetic in the cyclotomic field Q(zeta_N).

An element is stored as num / den with num an integer vector of length d = phi(N) (coefficients
in the power basis 1, zeta, ..., zeta^(d-1)) and den a positive integer.  zeta is identified with
exp(2 pi i / N) for numerical evaluation.  Standard library plus numpy; products use int64 when
an a priori bound shows that no overflow is possible and Python integers otherwise.
"""
import math
from fractions import Fraction

import numpy as np

_BIG = 2 ** 62


def _polydiv_exact(a, b):
    a = list(a)
    q = [0] * (len(a) - len(b) + 1)
    for i in range(len(a) - len(b), -1, -1):
        c = a[i + len(b) - 1] // b[-1]
        q[i] = c
        for j in range(len(b)):
            a[i + j] -= c * b[j]
    assert not any(a)
    return q


def cyclotomic_poly(n, _cache={}):
    """Coefficients (low degree first) of the n-th cyclotomic polynomial."""
    if n not in _cache:
        p = [-1] + [0] * (n - 1) + [1]
        for e in range(1, n):
            if n % e == 0:
                p = _polydiv_exact(p, cyclotomic_poly(e))
        _cache[n] = p
    return _cache[n]


def _maxabs(v):
    return int(np.max(np.abs(v))) if len(v) else 0


class Field:
    """Q(zeta_N)."""

    def __init__(self, N):
        self.N = N
        self.phi = cyclotomic_poly(N)
        self.d = d = len(self.phi) - 1
        # RED[k] = x^(d+k) mod Phi_N, k = 0 .. N-d-1
        red, cur = [], [-c for c in self.phi[:d]]
        for k in range(N - d):
            red.append(cur)
            top = cur[-1]
            cur = [0] + cur[:-1]
            cur = [x - top * c for x, c in zip(cur, self.phi[:d])]
        self.RED = np.array(red, dtype=object)
        self.RED64 = np.array(red, dtype=np.int64)
        self.redmax = max(abs(x) for r in red for x in r) if red else 0
        self.zero = Elt(self, np.zeros(d, dtype=np.int64), 1)
        self.one = self.from_int(1)
        self._pow = {}

    def __repr__(self):
        return f"Q(zeta_{self.N})"

    def from_int(self, a):
        v = np.zeros(self.d, dtype=np.int64)
        v[0] = a
        return Elt(self, v, 1)

    def from_fraction(self, a):
        a = Fraction(a)
        return self.from_int(a.numerator) * Fraction(1, a.denominator)

    def reduce_poly(self, c):
        """integer polynomial (any length, low first) -> vector modulo Phi_N."""
        d = self.d
        c = np.asarray(c)
        if len(c) <= d:
            out = np.zeros(d, dtype=c.dtype if c.dtype != np.int64 else np.int64)
            out = out.astype(c.dtype)
            out[:len(c)] = c
            return out
        if len(c) > self.N:                 # fold with x^N = 1 first
            m = np.zeros(self.N, dtype=object)
            for s in range(0, len(c), self.N):
                blk = c[s:s + self.N].astype(object)
                m[:len(blk)] = m[:len(blk)] + blk
            c = m
        hi = c[d:]
        if c.dtype == np.int64 and _maxabs(c) * self.redmax * len(hi) + _maxabs(c) < _BIG:
            return c[:d] + hi @ self.RED64[:len(hi)]
        c = c.astype(object)
        return c[:d] + hi.astype(object) @ self.RED[:len(hi)]

    def zeta(self, k=1):
        """zeta_N^k."""
        k %= self.N
        if k not in self._pow:
            c = np.zeros(k + 1, dtype=np.int64)
            c[k] = 1
            self._pow[k] = Elt(self, self.reduce_poly(c), 1)
        return self._pow[k]

    def root_of_unity(self, a, m):
        """exp(2 pi i a / m), m must divide N."""
        assert self.N % m == 0, (a, m, self.N)
        return self.zeta(a * (self.N // m))

    def gauss_sqrt(self, n):
        """sqrt(n) for odd n > 0, as an element (requires 4n | N)."""
        g = self.zero
        for a in range(n):
            g = g + self.root_of_unity(a * a, n)
        # sum_a exp(2 pi i a^2/n) = sqrt(n) (n = 1 mod 4) or i sqrt(n) (n = 3 mod 4)
        return g if n % 4 == 1 else g * self.root_of_unity(-1, 4)

    def sqrt2(self):
        return self.root_of_unity(1, 8) + self.root_of_unity(-1, 8)

    def from_vector(self, v, den=1):
        return Elt(self, np.array([int(x) for x in v], dtype=object), den).normalize()


class Elt:
    __slots__ = ("F", "num", "den")

    def __init__(self, F, num, den):
        self.F, self.num, self.den = F, num, den

    # -- normalisation
    def normalize(self):
        num = self.num
        g = math.gcd(self.den, *[int(x) for x in num]) if len(num) else self.den
        if g < 0:
            g = -g
        if g > 1:
            num = num // g if num.dtype == np.int64 else np.array([int(x) // g for x in num], dtype=object)
            self.den //= g
        if num.dtype == object and _maxabs(num) < 2 ** 40:
            num = num.astype(np.int64)
        self.num = num
        return self

    def is_zero(self):
        return not np.any(self.num)

    def __eq__(self, other):
        if not isinstance(other, Elt):
            other = self.F.from_fraction(other)
        return (self - other).is_zero()

    def __hash__(self):
        raise TypeError("unhashable")

    # -- ring operations
    def _coerce(self, o):
        if isinstance(o, Elt):
            return o
        return self.F.from_fraction(o)

    def __add__(self, o):
        o = self._coerce(o)
        if self.den == o.den:
            return Elt(self.F, _add(self.num, o.num), self.den).normalize()
        return Elt(self.F, _add(_scale(self.num, o.den), _scale(o.num, self.den)), self.den * o.den).normalize()

    __radd__ = __add__

    def __neg__(self):
        return Elt(self.F, -self.num, self.den)

    def __sub__(self, o):
        return self + (-self._coerce(o))

    def __rsub__(self, o):
        return self._coerce(o) - self

    def __mul__(self, o):
        if isinstance(o, (int, Fraction)):
            o = Fraction(o)
            return Elt(self.F, _scale(self.num, o.numerator), self.den * o.denominator).normalize()
        a, b = self.num, o.num
        d = len(a)
        if a.dtype == np.int64 and b.dtype == np.int64 and _maxabs(a) * _maxabs(b) * d < _BIG // (self.F.redmax * d + 1):
            c = np.convolve(a, b)
        else:
            c = np.convolve(a.astype(object), b.astype(object))
        return Elt(self.F, self.F.reduce_poly(c), self.den * o.den).normalize()

    __rmul__ = __mul__

    def galois(self, k):
        """the automorphism zeta -> zeta^k, gcd(k, N) = 1."""
        F = self.F
        assert math.gcd(k, F.N) == 1
        c = np.zeros(F.N, dtype=object)
        for j, x in enumerate(self.num):
            if x:
                c[(j * k) % F.N] += int(x)
        return Elt(F, F.reduce_poly(c), self.den).normalize()

    def conj(self):
        return self.galois(-1)

    def inverse(self):
        assert not self.is_zero(), "division by zero"
        F = self.F
        prod = F.one
        for k in range(2, F.N):
            if math.gcd(k, F.N) == 1:
                prod = prod * self.galois(k)
        nrm = self * prod
        assert not np.any(nrm.num[1:]), "norm is not rational"
        return prod * Fraction(nrm.den, int(nrm.num[0]))

    def __truediv__(self, o):
        if isinstance(o, (int, Fraction)):
            return self * (1 / Fraction(o))
        return self * o.inverse()

    def __rtruediv__(self, o):
        return self._coerce(o) * self.inverse()

    def __pow__(self, e):
        if e < 0:
            return self.inverse() ** (-e)
        r, b = self.F.one, self
        while e:
            if e & 1:
                r = r * b
            b = b * b
            e >>= 1
        return r

    # -- output
    def coeffs(self):
        return [Fraction(int(x), self.den) for x in self.num]

    def __complex__(self):
        z = np.exp(2j * np.pi * np.arange(self.F.d) / self.F.N)
        return complex(sum(complex(int(x)) * zz for x, zz in zip(self.num, z)) / self.den)

    def height_bits(self):
        return max(max(abs(int(x)) for x in self.num).bit_length(), self.den.bit_length())

    def __repr__(self):
        return f"<{complex(self):.6g} in {self.F}>"


def _add(a, b):
    if a.dtype == np.int64 and b.dtype == np.int64 and max(_maxabs(a), _maxabs(b)) < 2 ** 61:
        return a + b
    return a.astype(object) + b.astype(object)


def _scale(a, k):
    k = int(k)
    if a.dtype == np.int64 and abs(k) < _BIG and _maxabs(a) * abs(k) < _BIG:
        return a * k
    return a.astype(object) * k
