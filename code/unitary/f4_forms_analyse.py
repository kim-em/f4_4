#!/usr/bin/env python3
"""Positivity of the forms G of C(F4,4) at the 6 real embeddings of k, and the square classes
needed for an exact orthogonal gauge.

usage: f4_forms_analyse.py G_exact.pkl NPRIMES
  part 1 (positivity): G evaluated at 256 bits; eigenvalues of every G^{XY}_Z; the sign pattern
         eps_t of the definite spaces tested for the form s_X s_Y s_Z (torus signs).
  part 2 (square classes): for NPRIMES primes p = 1 mod 52 and the 6 roots w (6 NPRIMES degree-one
         primes of k), the quadratic characters of the LDL^T diagonals of every G^{XY}_Z; F_2 ranks of
         (a) invariant classes of one-dimensional spaces (cycle space of the torus action),
         (b) the same with the first LDL diagonal of every space, (c) plus the ratios D_i / D_1.
"""
import os, pickle, sys, time
sys.dont_write_bytecode = True
from fractions import Fraction
import numpy as np
import mpmath
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "f4"))
from lin import find_primes, primitive_root_of_unity

UNIT = "0000"
KREP = [1, 3, 5, 7, 9, 11]


def load(path):
    d = pickle.load(open(path, "rb"))
    keys, shapes = d["layout"]["keys"], d["layout"]["shapes"]
    flat = [[Fraction(a, b) for (a, b) in r] for r in d["G_eta"]]
    G, pos = {}, 0
    for t, s in zip(keys, shapes):
        n = s[0] * s[1]
        G[t] = [flat[pos + i * s[1]: pos + (i + 1) * s[1]] for i in range(s[0])]
        pos += n
    return G


def numeric(G, w, prec=256):
    mpmath.mp.prec = prec
    eta = [2 * mpmath.cos(2 * mpmath.pi * 4 * j * w / 52) for j in range(1, 7)]
    ev = lambda c: sum(mpmath.mpf(x.numerator) / x.denominator * e for x, e in zip(c, eta))
    return {t: mpmath.matrix([[ev(c) for c in row] for row in M]) for t, M in G.items()}


def f2_rank(vecs):
    """rank over F_2 of a list of python ints (bit vectors)."""
    basis = {}
    for v in vecs:
        while v:
            h = v.bit_length() - 1
            if h in basis:
                v ^= basis[h]
            else:
                basis[h] = v
                break
    return len(basis), basis


def reduce_vec(v, basis):
    while v:
        h = v.bit_length() - 1
        if h not in basis:
            return v
        v ^= basis[h]
    return 0


def solve_f2(rows, rhs, nvar):
    """rows: list of int bitmasks over nvar variables; returns a solution or None."""
    piv = {}
    for r, b in zip(rows, rhs):
        while r:
            h = r.bit_length() - 1
            if h in piv:
                r ^= piv[h][0]; b ^= piv[h][1]
            else:
                piv[h] = (r, b); break
        else:
            if b:
                return None
    x = 0
    for h in sorted(piv):
        r, b = piv[h]
        low = r & ~(1 << h)
        val = b ^ (bin(low & x).count("1") & 1)
        if val:
            x |= 1 << h
    return x


def part1(G):
    objs = sorted({t[0] for t in G} | {t[2] for t in G})
    oi = {X: i for i, X in enumerate(objs)}
    spaces = [t for t in G if UNIT not in t[:2]]
    print(f"{len(G)} spaces ({len(spaces)} without unit leg); dims {sorted(set(len(G[t]) for t in G))}")
    for w in KREP:
        t0 = time.time()
        Gn = numeric(G, w)
        unit_ok = all(mpmath.mnorm(Gn[t] - mpmath.eye(Gn[t].rows), 1) < mpmath.mpf(2) ** -200 for t in G if UNIT in t[:2])
        eps, indef, minrel, cond = {}, 0, mpmath.inf, 0
        for t in spaces:
            E = mpmath.eigsy(Gn[t], eigvals_only=True)
            E = [E[i] for i in range(len(E))]
            if all(e > 0 for e in E):
                eps[t] = 0
            elif all(e < 0 for e in E):
                eps[t] = 1
            else:
                indef += 1
                continue
            a = [abs(e) for e in E]
            minrel = min(minrel, min(a) / max(a))
        # torus signs: eps_t = s_X + s_Y + s_Z (mod 2), s_unit = 0
        rows = [(1 << oi[t[0]]) ^ (1 << oi[t[1]]) ^ (1 << oi[t[2]]) for t in eps]
        rows = [r & ~(1 << oi[UNIT]) for r in rows]
        sol = solve_f2(rows, list(eps.values()), len(objs)) if not indef else None
        npos = sum(1 for v in eps.values() if v == 0)
        verdict = "UNITARY (positive after a torus sign change)" if sol is not None else "not unitary"
        print(f"w={w:2d}: unit-leg spaces identity: {unit_ok}; definite {len(eps)} (positive {npos}), indefinite {indef}; "
              f"smallest |eigenvalue| ratio {mpmath.nstr(minrel, 3)}; torus sign pattern solvable: {sol is not None} -> {verdict} "
              f"({time.time() - t0:.0f} s)", flush=True)
        if sol is not None:
            print(f"      signs s_X = -1 on {[X for X in objs if sol >> oi[X] & 1]}")


def ldl_modp(M, p):
    n = len(M)
    A = [list(r) for r in M]
    D = []
    for i in range(n):
        d = A[i][i] % p
        if d == 0:
            return None
        D.append(d)
        inv = pow(d, -1, p)
        for r in range(i + 1, n):
            f = A[r][i] * inv % p
            for c in range(i + 1, n):
                A[r][c] = (A[r][c] - f * A[i][c]) % p
    return D


def part2(G, nprimes):
    spaces = sorted(t for t in G if UNIT not in t[:2])
    objs = sorted({t[0] for t in G} | {t[2] for t in G})
    oi = {X: i for i, X in enumerate(objs)}
    primes = find_primes(52, nprimes, below=2 ** 31)
    chars = {}          # (t, i) -> list of bits
    ndet = {}
    pos = 0
    t0 = time.time()
    for p in primes:
        z = primitive_root_of_unity(p, 52, 1)
        for w in KREP:
            zw = pow(z, w, p); zi = pow(zw, -1, p)
            eta = [(pow(zw, 4 * j, p) + pow(zi, 4 * j, p)) % p for j in range(1, 7)]
            red = lambda c: sum(x.numerator * e * pow(x.denominator, -1, p) for x, e in zip(c, eta)) % p
            res = {}
            ok = True
            for t in spaces:
                M = [[red(c) for c in row] for row in G[t]]
                D = ldl_modp(M, p)
                if D is None:
                    ok = False; break
                res[t] = [0 if pow(d, (p - 1) // 2, p) == 1 else 1 for d in D]
            if not ok:
                print(f"  skip p={p} w={w} (a pivot vanishes)"); continue
            for t in spaces:
                for i, b in enumerate(res[t]):
                    chars.setdefault((t, i), 0)
                    if b:
                        chars[(t, i)] |= 1 << pos
            pos += 1
    print(f"square classes: {pos} degree-one primes of k ({len(primes)} rational primes x 6), {time.time() - t0:.0f} s")
    # cycle space of the torus action on spaces: a with sum_t a_t (e_X + e_Y + e_Z) = 0
    def inv_span(ts, vec):
        # kernel of the incidence map on the set ts, then the span of sum a_t vec[t]
        inc = [((1 << oi[t[0]]) ^ (1 << oi[t[1]]) ^ (1 << oi[t[2]])) & ~(1 << oi[UNIT]) for t in ts]
        # gaussian elimination tracking combinations
        piv = {}
        kernel_vals = []
        for idx, (r, v) in enumerate(zip(inc, (vec[t] for t in ts))):
            comb_v = v
            while r:
                h = r.bit_length() - 1
                if h in piv:
                    r ^= piv[h][0]; comb_v ^= piv[h][1]
                else:
                    piv[h] = (r, comb_v); break
            else:
                kernel_vals.append(comb_v)
        return kernel_vals
    one = [t for t in spaces if len(G[t]) == 1]
    first = {t: chars[(t, 0)] for t in spaces}
    rk_a, _ = f2_rank(inv_span(one, first))
    rk_b, basis_b = f2_rank(inv_span(spaces, first))
    diffs = [chars[(t, i)] ^ chars[(t, 0)] for t in spaces for i in range(1, len(G[t]))]
    rk_c, _ = f2_rank(list(basis_b.values()) + diffs)
    rk_d, _ = f2_rank(diffs)
    dets = [0 for _ in []]
    for t in spaces:
        if len(G[t]) % 2 == 0:
            v = 0
            for i in range(len(G[t])):
                v ^= chars[(t, i)]
            dets.append(v)
    rk_det, _ = f2_rank(dets + inv_span(one, first))
    rk_all, _ = f2_rank([chars[k] for k in chars])
    # determinant classes: det(c G) = c^d det G, so for odd d the det class is torus covariant like a
    # one-dimensional norm, for even d it is invariant
    det = {}
    for t in spaces:
        v = 0
        for i in range(len(G[t])):
            v ^= chars[(t, i)]
        det[t] = v
    odd = [t for t in spaces if len(G[t]) % 2 == 1]
    ev = [det[t] for t in spaces if len(G[t]) % 2 == 0]
    rk_e, basis_e = f2_rank(inv_span(odd, det) + ev)
    print(f"  (e) torus-invariant determinant classes (odd-dimensional spaces covariantly, even ones invariantly): {rk_e}")
    # candidates: -1, 2, 3, 13, quantum integers [n] (n odd) and [2][n] (n even), at q = zeta_52^w
    cand = {}
    pos = 0
    for p in primes:
        z = primitive_root_of_unity(p, 52, 1)
        for w in KREP:
            zw = pow(z, w, p)
            qi = lambda n: (pow(zw, n, p) - pow(zw, -n, p)) * pow(zw - pow(zw, -1, p), -1, p) % p
            els = {"-1": p - 1, "2": 2, "3": 3, "13": 13}
            for n in range(2, 26):
                els[f"[{n}]" if n % 2 else f"[2][{n}]"] = qi(n) if n % 2 else qi(2) * qi(n) % p
            for name, x in els.items():
                cand.setdefault(name, 0)
                if pow(x, (p - 1) // 2, p) != 1:
                    cand[name] |= 1 << pos
            pos += 1
    rk_cand, basis_cand = f2_rank(list(cand.values()))
    inv_a, basis_a = f2_rank(inv_span(one, first))
    inside_a = all(reduce_vec(v, basis_cand) == 0 for v in basis_a.values())
    inside_e = all(reduce_vec(v, basis_cand) == 0 for v in basis_e.values())
    print(f"  candidates -1, 2, 3, 13, [n] (n odd), [2][n] (n even), n < 26: span rank {rk_cand}; "
          f"(a) inside their span: {inside_a}; (e) inside their span: {inside_e}")
    # which candidates generate (a)
    gen = [name for name, v in cand.items() if reduce_vec(v, basis_a) == 0 and v]
    print(f"  candidates whose class lies in (a): {gen}")
    print(f"F_2 ranks (characters at {pos} primes; a rank close to {pos} would mean too few primes):")
    print(f"  (a) torus-invariant classes from the one-dimensional spaces ({len(one)} spaces): {rk_a}")
    print(f"  (b) torus-invariant classes from the first LDL diagonal of all {len(spaces)} spaces: {rk_b}")
    print(f"  (c) (b) plus all ratios D_i/D_1 within spaces (upper bound for LDL+sqrt gauges, best torus): {rk_c}")
    print(f"      ratios D_i/D_1 alone: {rk_d}")
    print(f"  (d) (a) plus det classes of even-dimensional spaces (lower bound for the field of an orthogonal gauge): {rk_det}")
    print(f"  all LDL diagonals, no torus optimisation: {rk_all}")
    return chars


if __name__ == "__main__":
    G = load(sys.argv[1])
    if "--no-pos" not in sys.argv:
        part1(G)
    part2(G, int(sys.argv[2]))
