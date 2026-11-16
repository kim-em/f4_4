#!/usr/bin/env python3
"""Exact comparison of the quantum group data for C(so(13),2) with the AFT formulas (p = 6).

usage: compare_so13.py [OUT_GAUGE_FILE]
       compare_so13.py --check GAUGE_FILE     (only step 4, with the gauge read from the file)

1. Reads data/B6_2_FR_exact.txt (read-only) over K = Q(zeta_104).
2. Builds the AFT data with r = 1, kappa = -1, lam = -1 and the label dictionary LABELS.
3. Solves for a gauge g (one scalar per admissible triple, in K) with
       F_ours^g = F_AFT  and  R_ours^g = R_AFT,
   where  F^g_{XYW}^{T;ZZ'} = g(X,Y,Z) g(Z,W,T) F / (g(Y,W,Z') g(X,Z',T))  and
          R^g(X,Y,Z) = g(Y,X,Z) R(X,Y,Z) / g(X,Y,Z),
   by elimination with pivots of exponent +-1 only (so no roots are taken).
4. Applies g to every F- and R-symbol of ours and checks equality with AFT entrywise, exactly
   (including the zero pattern).  This is the proof of the equivalence.
5. Also checks the second dictionary (r = 1, kappa = -1, lam = +1, phi_i -> phi_{g(5i)}).
"""
import sys
import time

from cyclo import Field
from ours import Ours
from so2p1_2 import Metaplectic

# ours (Dynkin labels of B6) -> AFT labels, for (r, kappa, lam) = (1, -1, -1)
LABELS = {"000000": "1", "200000": "eps",
          "100000": "phi1", "010000": "phi2", "001000": "phi3", "000100": "phi4",
          "000010": "phi5", "000002": "phi6",
          "100001": "psi+", "000001": "psi-"}


def second_dictionary():
    """(r, kappa, lam) = (1, -1, +1): phi_i -> phi_{g(5 i)}, psi+ = spinor (5^2 = -1 mod 13)."""
    out = {}
    for k, v in LABELS.items():
        if v.startswith("phi"):
            i = 5 * int(v[3:]) % 13
            out[k] = "phi%d" % min(i, 13 - i)
        else:
            out[k] = {"psi+": "psi-", "psi-": "psi+"}.get(v, v)
    return out


def equations(O):
    """(exponent dict over triples, ('F', key, e, f) or ('R', key))."""
    eqs = []
    for (X, Y, W, T), blk in O.F.items():
        for (Z, Zp) in O.rows_cols[(X, Y, W, T)]:
            m = {}
            for t, s in (((X, Y, Z), 1), ((Z, W, T), 1), ((Y, W, Zp), -1), ((X, Zp, T), -1)):
                m[t] = m.get(t, 0) + s
            eqs.append(({t: s for t, s in m.items() if s}, ("F", (X, Y, W, T), Z, Zp)))
    for (X, Y, Z) in O.R:
        m = {}
        for t, s in (((Y, X, Z), 1), ((X, Y, Z), -1)):
            m[t] = m.get(t, 0) + s
        eqs.append(({t: s for t, s in m.items() if s}, ("R", (X, Y, Z))))
    return eqs


def solve_gauge(O, A, lab):
    """g with ours^g = AFT o lab, by +-1-pivot elimination; returns dict triple -> K element."""
    K = O.K
    def target(tag):
        if tag[0] == "F":
            _, (X, Y, W, T), Z, Zp = tag
            a = A.F[(lab[X], lab[Y], lab[W], lab[T])].get((lab[Z], lab[Zp]), K.zero)
            b = O.F[(X, Y, W, T)].get((Z, Zp), K.zero)
        else:
            _, (X, Y, Z) = tag
            a, b = A.R[(lab[X], lab[Y], lab[Z])], O.R[(X, Y, Z)]
        return a, b
    # sol[t] = (const, {free var: exponent})
    sol = {}
    inv_cache = {}
    def inv(x, key):
        if key not in inv_cache:
            inv_cache[key] = x.inverse()
        return inv_cache[key]
    pending = equations(O)
    progress = True
    while progress:
        progress = False
        rest = []
        for mono, tag in pending:
            a, b = target(tag)
            if a.is_zero() or b.is_zero():
                continue                      # zero pattern: checked in the verification
            # substitute solved variables: monomial(g) = a/b  ->  const * monomial(free) = a/b
            const, free = K.one, {}
            for t, e in mono.items():
                if t in sol:
                    c, f = sol[t]
                    const = const * (c ** e if e > 0 else inv(c, ("sol", t)) ** (-e))
                    for s, x in f.items():
                        free[s] = free.get(s, 0) + e * x
                else:
                    free[t] = free.get(t, 0) + e
            free = {s: x for s, x in free.items() if x}
            if not free:
                continue                      # a consequence; checked in the verification
            piv = next((s for s, x in free.items() if abs(x) == 1), None)
            if piv is None:
                rest.append((mono, tag))
                continue
            # g_piv^x * const * prod_{others} g^y = a/b
            x = free.pop(piv)
            rhs = a * inv(b, tag) * inv(const, ("c", tag, len(sol)))
            if x == 1:
                new = (rhs, {s: -y for s, y in free.items()})
            else:
                new = (inv(rhs, ("r", tag)), dict(free))
            # substitute piv into the existing solutions
            for t in list(sol):
                c, f = sol[t]
                if piv in f:
                    e = f.pop(piv)
                    c = c * (new[0] ** e if e > 0 else inv(new[0], ("sol", piv)) ** (-e))
                    for s, y in new[1].items():
                        f[s] = f.get(s, 0) + e * y
                    sol[t] = (c, {s: y for s, y in f.items() if y})
                    inv_cache.pop(("sol", t), None)
            sol[piv] = new
            progress = True
        pending = rest
    unresolved = [tag for mono, tag in pending]
    triples = sorted(O.N)
    g = {}
    for t in triples:
        if t in sol:
            c, f = sol[t]
            g[t] = c                          # free variables set to 1
        else:
            g[t] = K.one
    # g -> g u(a) u(b) / u(c) acts trivially on F and R; with u(1) = 1/g(1,1,1), u = 1 otherwise,
    # every triple with a unit label gets g = 1 (they all carry the same value, forced by the unit-leg F's)
    one = O.unit
    c = g[(one, one, one)]
    ci = c.inverse()
    for (a, b, e) in triples:
        k = (a == one) + (b == one) - (e == one)
        if k:
            g[(a, b, e)] = g[(a, b, e)] * (ci ** k if k > 0 else c ** (-k))
    return g, unresolved


def apply_and_compare(O, A, lab, g):
    K = O.K
    ginv = {t: v.inverse() for t, v in g.items()}
    bad = []
    nF = nR = 0
    for key, (rows, cols) in O.blocks.items():
        X, Y, W, T = key
        ak = tuple(lab[x] for x in key)
        for Z in rows:
            for Zp in cols:
                v = O.F[key].get((Z, Zp), K.zero) * g[(X, Y, Z)] * g[(Z, W, T)] * ginv[(Y, W, Zp)] * ginv[(X, Zp, T)]
                w = A.F[ak].get((lab[Z], lab[Zp]), K.zero)
                nF += 1
                if not (v - w).is_zero():
                    bad.append(("F", key, Z, Zp))
    for (X, Y, Z), r in O.R.items():
        v = r * g[(Y, X, Z)] * ginv[(X, Y, Z)]
        nR += 1
        if not (v - A.R[(lab[X], lab[Y], lab[Z])]).is_zero():
            bad.append(("R", (X, Y, Z)))
    return bad, nF, nR


def prepare(O):
    """admissible rows/cols of every block of ours (including zero entries)."""
    P = {}
    for (a, b, c) in O.N:
        P.setdefault((a, b), []).append(c)
    O.blocks, O.rows_cols = {}, {}
    for key in O.F:
        X, Y, W, T = key
        rows = [Z for Z in O.labels if (X, Y, Z) in O.N and (Z, W, T) in O.N]
        cols = [Z for Z in O.labels if (Y, W, Z) in O.N and (X, Z, T) in O.N]
        O.blocks[key] = (rows, cols)
        O.rows_cols[key] = [(Z, Zp) for Z in rows for Zp in cols]


def run(lab, r, kappa, lam, K, O, out=None):
    t0 = time.time()
    A = Metaplectic(6, r, kappa, lam, field=K)
    assert {(lab[a], lab[b], lab[c]) for (a, b, c) in O.N} == set(A.N), "fusion rules differ"
    assert set(tuple(lab[x] for x in k) for k in O.F) == set(A.F), "blocks differ"
    g, unresolved = solve_gauge(O, A, lab)
    bad, nF, nR = apply_and_compare(O, A, lab, g)
    print(f"(r, kappa, lam) = ({r}, {kappa}, {lam}), labels {lab}")
    print(f"  gauge solved by +-1 pivots; equations left with only non-unit exponents: {len(unresolved)}")
    print(f"  after the gauge: {nF} F entries (all admissible pairs, zeros included) and {nR} R-symbols compared exactly; mismatches: {len(bad)}  ({time.time() - t0:.0f} s)")
    if bad:
        print("  first mismatches:", bad[:10])
    if out:
        with open(out, "w") as fh:
            fh.write("# gauge g(X,Y,Z) in K = Q(zeta_104), zeta = exp(2 pi i/104), power basis coefficients b_0 .. b_47\n")
            fh.write(f"# takes data/B6_2_FR_exact.txt (rescaled fusion maps g(X,Y,Z) phi^(XY->Z)) to AFT (r, kappa, lam) = ({r}, {kappa}, {lam})\n")
            fh.write("# labels (ours -> AFT): " + " ".join(f"{k}->{v}" for k, v in lab.items()) + "\n")
            fh.write("# F^g = g(X,Y,Z) g(Z,W,T) F / (g(Y,W,Z') g(X,Z',T)),  R^g(X,Y,Z) = g(Y,X,Z) R(X,Y,Z) / g(X,Y,Z)\n")
            for t in sorted(g):
                fh.write("g %s %s %s | %s\n" % (t + (" ".join(str(c) for c in g[t].coeffs()),)))
    return not bad, g


def read_gauge(path, K):
    from fractions import Fraction
    g = {}
    for line in open(path):
        if line.startswith("g "):
            head, cs = line.split("|")
            v = K.zero
            for k, x in enumerate(cs.split()):
                if x != "0":
                    v = v + K.zeta(k) * Fraction(x)
            g[tuple(head.split()[1:])] = v
    return g


def main(argv):
    K = Field(104)
    O = Ours(field=K)
    prepare(O)
    if argv[1:2] == ["--check"]:
        g = read_gauge(argv[2], K)
        A = Metaplectic(6, 1, -1, -1, field=K)
        bad, nF, nR = apply_and_compare(O, A, LABELS, g)
        print(f"gauge from {argv[2]} applied to ours: {nF} F entries and {nR} R-symbols compared exactly with AFT (1, -1, -1); mismatches {len(bad)}")
        print("CONCLUSION:", "agreement" if not bad else "DISAGREEMENT")
        return
    out = argv[1] if len(argv) > 1 else None
    ok1, g = run(LABELS, 1, -1, -1, K, O, out)
    ok2, _ = run(second_dictionary(), 1, -1, 1, K, O)
    hb = max(v.height_bits() for v in g.values())
    print(f"gauge heights: at most {hb} bits")
    print("CONCLUSION:", "agreement (explicit gauge, exact)" if ok1 and ok2 else "DISAGREEMENT")


if __name__ == "__main__":
    main(sys.argv)
