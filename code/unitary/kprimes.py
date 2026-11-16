"""Degree-one primes of k = Q(zeta_13)^+ and quadratic characters of elements of k.

Primes p = 1 mod 52 split completely in Q(zeta_52), so their characters cannot see the classes
c in k^*/k^*2 with sqrt(c) in Q(zeta_52) (the classes of -1 and of the c0 with
Q(zeta_52)^+ = k(sqrt c0)).  Here all p = +-1 mod 13 are used (eta_1 has a root mod p), which
detects every class.  Elements of k are given by their 6 eta coefficients (Fractions).
"""
import os, sys
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "fields"))
sys.dont_write_bytecode = True
import cyc

POL = [1, 1, -5, -4, 6, 3, -1]          # minimal polynomial of eta_1, leading coefficient first


def _is_prime(n):
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


def _roots(p):
    # brute force is too slow; use the fact that eta_1 = z + 1/z for z a 13th root of unity in F_p or F_p^2
    roots = set()
    if p % 13 == 1:
        g = 2
        while True:
            z = pow(g, (p - 1) // 13, p)
            if z != 1:
                break
            g += 1
        for k in range(1, 13):
            zk = pow(z, k, p)
            roots.add((zk + pow(zk, -1, p)) % p)
    else:
        # p = -1 mod 13: roots via a root of x^2 - eta x + 1 in F_p^2; find eta by factoring POL with sympy
        import sympy
        roots = {int(r) % p for r in sympy.Poly(POL, sympy.Symbol("x"), modulus=p).ground_roots().keys()}
    return sorted(roots)


def kprimes(count, start=1 << 24):
    """yields (p, eta) with eta = [eta_1..eta_6] mod p, for degree-one primes of k (6 per p)."""
    p, n = start, 0
    while n < count:
        p += 1
        if p % 13 not in (1, 12) or not _is_prime(p):
            continue
        for r in _roots(p):
            e = [2, r]
            for j in range(2, 7):
                e.append((r * e[-1] - e[-2]) % p)
            yield p, e[1:]
            n += 1
            if n >= count:
                return


def red(cs, p, eta):
    v = 0
    for c, e in zip(cs, eta):
        if c:
            assert c.denominator % p
            v = (v + c.numerator * pow(c.denominator, -1, p) * e) % p
    return v


def char(x, p):
    assert x % p
    return 0 if pow(x, (p - 1) // 2, p) == 1 else 1


def candidates():
    """named elements of k in the eta basis: -1, 2, 3, 13, [n] (n odd), [2][n] (n even), q = zeta_52."""
    z, zi = cyc.zpow(1), cyc.zpow(-1)
    def qi(n):
        return cyc.mul(cyc.sub(cyc.zpow(n), cyc.zpow(-n)), cyc.inv(cyc.sub(z, zi)))
    out = {"-1": [Fraction(1)] * 6, "2": [Fraction(-2)] * 6, "3": [Fraction(-3)] * 6, "13": [Fraction(-13)] * 6}
    for n in range(3, 26, 2):
        out[f"[{n}]"] = cyc.to_eta(qi(n))
    for n in range(4, 26, 2):
        out[f"[2][{n}]"] = cyc.to_eta(cyc.mul(qi(2), qi(n)))
    for k, v in out.items():
        assert v is not None, k
    return out
