"""Exact pentagon and hexagon checks for multiplicity-free F- and R-symbols.

Data interface: .labels, .N (dict (a,b,c) -> 1 for admissible triples), .F (dict (a,b,c,d) ->
{(e,f): value}), .R (dict (a,b,c) -> value), values exact (cyclo.Elt).  Conventions as in
so2p1_2.py:
  pentagon  F_{fcd}^{e;gl} F_{abl}^{e;fk} = sum_h F_{abc}^{g;fh} F_{ahd}^{e;gk} F_{bcd}^{k;hl}
  hexagons  R_{ac}^e F_{acb}^{d;eg} R_{bc}^g = sum_f F_{cab}^{d;ef} R_{fc}^d F_{abc}^{d;fg}   (H+)
            and (H+) for the reverse-inverse braiding, R'_{ab}^c = (R_{ba}^c)^-1               (H-)
(H+) is AFT (R1).  AFT print the second hexagon (R2) as (H+) with every R inverted but the
lower indices kept; that form is not gauge covariant together with (R1) and agrees with (H-)
only in gauges with R_{ab}^c = R_{ba}^c (such as AFT's); `mode="aft-literal"` checks it as printed.
These coincide with the conventions of ../f4/verify_fr.py for F[(X,Y,W,T)][(Z,0,0),(Z',0,0)]
= F_{XYW}^{T;ZZ'} and R[(X,Y,Z)] = R_{XY}^Z: the pentagon is its pentagon_tuple, (H+) is its
(H2) and (H-) is its (H1).
"""
import itertools
import sys
from multiprocessing import Pool


def _prod(D):
    """(a,b) -> list of c."""
    out = {}
    for (a, b, c) in D.N:
        out.setdefault((a, b), []).append(c)
    return out


def pentagon_failures(D, first=None, stop_after=20):
    """All pentagon equations with first label a in `first` (default: all): one scalar equation for
    every basis vector (f, g) of Hom(((ab)c)d, e) and (l, k) of Hom(a(b(cd)), e), i.e. every
    entry of the matrix identity of ../f4/verify_fr.py, including those whose left side is 0."""
    P = _prod(D)
    F = D.F
    zero = D.K.zero
    fe = lambda a, b, c, d, e, f: F[(a, b, c, d)].get((e, f), zero) if (a, b, c, d) in F else zero
    fails, count = [], 0
    L = D.labels
    for a in (first or L):
        for b, c, d in itertools.product(L, repeat=3):
            for f in P.get((a, b), []):
                for g in P.get((f, c), []):
                    for e in P.get((g, d), []):
                        for l in P.get((c, d), []):
                            for k in P.get((b, l), []):
                                if (a, k, e) not in D.N:
                                    continue
                                # if (f, l, e) is not admissible the left side is 0 and the right side must vanish
                                lhs = fe(f, c, d, e, g, l) * fe(a, b, l, e, f, k) if (f, l, e) in D.N else zero
                                rhs = zero
                                for h in P.get((b, c), []):
                                    if (a, h, g) in D.N and (h, d, k) in D.N:
                                        rhs = rhs + fe(a, b, c, g, f, h) * fe(a, h, d, e, g, k) * fe(b, c, d, k, h, l)
                                count += 1
                                if not (lhs - rhs).is_zero():
                                    fails.append((a, b, c, d, e, f, g, k, l))
                                    if len(fails) >= stop_after:
                                        return fails, count
    return fails, count


def hexagon_failures(D, mode="+", stop_after=20):
    """mode "+" (H+), "-" (H-), "aft-literal" (AFT (R2) as printed)."""
    P = _prod(D)
    F = D.F
    zero = D.K.zero
    fe = lambda a, b, c, d, e, f: F[(a, b, c, d)].get((e, f), zero) if (a, b, c, d) in F else zero
    if mode == "+":
        R = D.R
    elif mode == "-":
        R = {(a, b, c): D.R[(b, a, c)].inverse() for (a, b, c) in D.R}
    else:
        assert mode == "aft-literal"
        R = {k: v.inverse() for k, v in D.R.items()}
    fails, count = [], 0
    L = D.labels
    for a, b, c, d in itertools.product(L, repeat=4):
        if (a, c, b, d) not in F:
            continue
        for e in P.get((a, c), []):
            for g in P.get((c, b), []):
                if (e, b, d) not in D.N or (a, g, d) not in D.N:
                    continue
                lhs = R[(a, c, e)] * fe(a, c, b, d, e, g) * R[(b, c, g)]
                rhs = zero
                for f in P.get((a, b), []):
                    if (c, f, d) in D.N and (f, c, d) in D.N:
                        rhs = rhs + fe(c, a, b, d, e, f) * R[(f, c, d)] * fe(a, b, c, d, f, g)
                count += 1
                if not (lhs - rhs).is_zero():
                    fails.append((a, b, c, d, e, g))
                    if len(fails) >= stop_after:
                        return fails, count
    return fails, count


def orthogonality_failures(D):
    """blocks whose matrix M (rows, cols in label order) does not satisfy M M^T = 1."""
    bad = []
    for key, blk in D.F.items():
        rows = sorted({e for e, f in blk}, key=D.labels.index)
        cols = sorted({f for e, f in blk}, key=D.labels.index)
        for r1, r2 in itertools.product(rows, repeat=2):
            s = D.K.zero
            for f in cols:
                s = s + blk.get((r1, f), D.K.zero) * blk.get((r2, f), D.K.zero)
            if not (s - (1 if r1 == r2 else 0)).is_zero():
                bad.append(key)
                break
    return bad


def _pent_job(args):
    D, a = args
    return pentagon_failures(D, first=[a])


def pentagon_parallel(D, jobs=8):
    with Pool(jobs) as pool:
        res = pool.map(_pent_job, [(D, a) for a in D.labels])
    fails = [x for f, c in res for x in f]
    return fails, sum(c for f, c in res)
