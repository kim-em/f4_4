"""Gauge theory of multiplicity-free F-symbols as a torus action.

For multiplicity-free data the gauge group is the torus G = prod_{(X,Y,Z)} K^*, acting by
  F'[(X,Y,W,T)][Z,Z'] = g^{XY}_Z g^{ZW}_T / (g^{YW}_{Z'} g^{XZ'}_T) * F[(X,Y,W,T)][Z,Z']
(fusion-map convention of the paper, Section 3, phi'^{XY->Z} = g^{XY}_Z phi^{XY->Z}).
So log F' = log F + D log g with D an integer matrix (rows: nonzero F entries, columns: triples).

`unit_pivots` chooses a set P of entries by integer elimination with pivots +-1 only; if
rank(D_P) = rank(D), then D_P is onto Z^P and coker(D) is torsion free, and the normal form
"F'_e = 1 for e in P" is reached by a unique gauge up to the stabiliser.  Its entries are then
gauge-invariant Laurent monomials in F (a complete invariant of the gauge class).
"""
import random


def d_row(key, a, b):
    X, Y, W, T = key
    Z, Zp = a[0], b[0]
    row = {}
    for t, s in (((X, Y, Z), 1), ((Z, W, T), 1), ((Y, W, Zp), -1), ((X, Zp, T), -1)):
        row[t] = row.get(t, 0) + s
    return {t: s for t, s in row.items() if s}


def rank_mod_p(rows, cols, p=2147483647):
    idx = {c: i for i, c in enumerate(cols)}
    basis = {}          # pivot col -> row (dict) normalised to 1 at pivot
    r = 0
    for row in rows:
        v = {idx[c]: s % p for c, s in row.items() if s % p}
        while v:
            c = min(v)
            if c in basis:
                f = v[c]
                for cc, ss in basis[c].items():
                    v[cc] = (v.get(cc, 0) - f * ss) % p
                    if not v[cc]:
                        del v[cc]
            else:
                inv = pow(v[c], p - 2, p)
                basis[c] = {cc: ss * inv % p for cc, ss in v.items()}
                r += 1
                break
    return r


def unit_pivots(rows):
    """rows: list of dicts col->int.  Returns (pivots, deferred) where pivots is a list of
    (entry index, pivot col, reduced row, U) with reduced row = sum_e U[e] rows[e], pivot coefficient 1,
    and each reduced row zero at the pivot columns of the earlier pivots."""
    pivots = []
    pcol = {}
    deferred = []
    for e, row in enumerate(rows):
        v = dict(row)
        U = {e: 1}
        for (e2, c2, r2, U2) in pivots:
            f = v.get(c2, 0)
            if f:
                for cc, ss in r2.items():
                    nv = v.get(cc, 0) - f * ss
                    if nv:
                        v[cc] = nv
                    else:
                        v.pop(cc, None)
                for ee, uu in U2.items():
                    nu = U.get(ee, 0) - f * uu
                    if nu:
                        U[ee] = nu
                    else:
                        U.pop(ee, None)
        if not v:
            continue
        unit = [c for c, s in v.items() if abs(s) == 1]
        if not unit:
            deferred.append((e, v))
            continue
        c = min(unit, key=str)
        if v[c] == -1:
            v = {cc: -ss for cc, ss in v.items()}
            U = {ee: -uu for ee, uu in U.items()}
        # keep earlier reduced rows reduced?  Not needed: we only require each row to vanish
        # at the earlier pivot columns, which holds by construction.
        pivots.append((e, c, v, U))
        pcol[c] = len(pivots) - 1
    return pivots, deferred
