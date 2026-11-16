"""Exact arithmetic in K = Q(zeta), zeta = exp(2 pi i/52), elements as tuples of 24 Fractions
(coefficients of 1, zeta, ..., zeta^23), plus conversions to the real basis of Q(zeta_52)^+ and
to the integral basis eta_j = zeta13^j + zeta13^-j (j = 1..6) of Q(zeta_13)^+, zeta13 = zeta^4.
"""
from fractions import Fraction

NB, COND = 24, 52


def _cyclotomic(n):
    def div(a, b):
        a = list(a); q = [0] * (len(a) - len(b) + 1)
        for i in range(len(a) - len(b), -1, -1):
            c = a[i + len(b) - 1] // b[-1]; q[i] = c
            for j in range(len(b)):
                a[i + j] -= c * b[j]
        assert not any(a)
        return q
    p = [-1] + [0] * (n - 1) + [1]
    for d in range(1, n):
        if n % d == 0:
            p = div(p, _cyclotomic(d))
    return p


PHI = _cyclotomic(COND)
ZERO = tuple([Fraction(0)] * NB)
ONE = tuple([Fraction(1)] + [Fraction(0)] * (NB - 1))


def reduce(t):
    t = list(t)
    for k in range(len(t) - 1, NB - 1, -1):
        c = t[k]
        if c:
            t[k] = 0
            for j in range(NB + 1):
                t[k - NB + j] -= c * PHI[j]
    t = t[:NB] + [0] * (NB - len(t[:NB]))
    return tuple(Fraction(x) for x in t)


ZPOW = [reduce([0] * k + [1]) for k in range(COND)]


def zpow(k):
    return ZPOW[k % COND]


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def scal(c, a):
    return tuple(c * x for x in a)


def mul(a, b):
    out = [Fraction(0)] * (2 * NB - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y:
                    out[i + j] += x * y
    return reduce(out)


def is_zero(a):
    return not any(a)


def inv(a):
    """inverse by solving the multiplication-matrix system exactly."""
    assert not is_zero(a)
    # columns: a * zeta^j
    cols = [mul(a, ZPOW[j]) for j in range(NB)]
    M = [[cols[j][i] for j in range(NB)] + [Fraction(int(i == 0))] for i in range(NB)]
    n = NB
    for c in range(n):
        piv = next(r for r in range(c, n) if M[r][c])
        M[c], M[piv] = M[piv], M[c]
        pv = M[c][c]
        M[c] = [x / pv for x in M[c]]
        for r in range(n):
            if r != c and M[r][c]:
                f = M[r][c]
                M[r] = [x - f * y for x, y in zip(M[r], M[c])]
    return tuple(M[i][n] for i in range(n))


def power(a, e):
    if e < 0:
        a, e = inv(a), -e
    r = ONE
    while e:
        if e & 1:
            r = mul(r, a)
        a = mul(a, a)
        e >>= 1
    return r


def sigma(a, k):
    """the automorphism zeta -> zeta^k."""
    out = [Fraction(0)] * NB
    for j, x in enumerate(a):
        if x:
            z = ZPOW[(j * k) % COND]
            for i in range(NB):
                if z[i]:
                    out[i] += x * z[i]
    return tuple(out)


def _solve(cols, target):
    """express target (a 24-vector) as a rational combination of cols; None if impossible."""
    m = len(cols)
    M = [[cols[j][i] for j in range(m)] + [target[i]] for i in range(NB)]
    piv_cols = []
    r = 0
    for c in range(m):
        piv = next((i for i in range(r, NB) if M[i][c]), None)
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        pv = M[r][c]
        M[r] = [x / pv for x in M[r]]
        for i in range(NB):
            if i != r and M[i][c]:
                f = M[i][c]
                M[i] = [x - f * y for x, y in zip(M[i], M[r])]
        piv_cols.append(c)
        r += 1
    if any(M[i][m] for i in range(r, NB)):
        return None
    sol = [Fraction(0)] * m
    for i, c in enumerate(piv_cols):
        sol[c] = M[i][m]
    return sol


REAL_BASIS = [ONE] + [add(ZPOW[k], ZPOW[COND - k]) for k in range(1, 12)]      # 1, c_1..c_11
ETA_BASIS = [add(ZPOW[4 * j], ZPOW[(COND - 4 * j) % COND]) for j in range(1, 7)]  # eta_1..eta_6


def to_real(a):
    s = _solve(REAL_BASIS, a)
    assert s is not None, "not in Q(zeta_52)^+"
    return s


def from_real(cs):
    out = ZERO
    for c, b in zip(cs, REAL_BASIS):
        if c:
            out = add(out, scal(Fraction(c), b))
    return out


def to_eta(a):
    """coefficients in eta_1..eta_6 (a Z-basis of the integers of Q(zeta_13)^+); None if a is not in that field."""
    return _solve(ETA_BASIS, a)


def from_eta(cs):
    out = ZERO
    for c, b in zip(cs, ETA_BASIS):
        if c:
            out = add(out, scal(Fraction(c), b))
    return out


def to_float(a, w=1):
    import cmath
    z = cmath.exp(2j * cmath.pi * w / COND)
    return sum(float(x) * z ** j for j, x in enumerate(a))


UNITS = [k for k in range(1, COND) if k % 2 and k % 13]
