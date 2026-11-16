#!/usr/bin/env python3
"""Canonical (gauge-invariant) normal form of the H family of C(F4,4), modulo p at a real root.

The H family is the set of F-matrices F[(X,Y,V,T)], V = 0001, in the export convention
(read from a binary F-table, by default the archived normal-form table F4_4_F_part{1,2}.bin in
the directory $F4_ARCHIVE, default archive/ of the repository; read only).  A unit-preserving gauge g (g = 1 on
spaces with a unit leg) is fixed in three steps; every step is a rational rule applied to
F-entries whose bases are already fixed, with index choices taken from a PLAN that is computed
once (at a reference prime) and then used verbatim, so the rule is the same over K = Q(zeta_52)
and modulo every prime.

 1. one-dimensional V-edges Hom(C, A (x) V) (A != 1): torus normal form.  Entries used:
    F[(A,V,V,C)][(Z,0,0),(Z',0,0)] with Z' in {1, V} and all four spaces one dimensional.
    Unit-pivot elimination gives a pivot set P; the gauge with F'_e = 1 on P and g = 1 on the
    non-pivot variables.  Residual freedom = stabiliser torus of this subsystem; it is checked
    to have dimension 15 = dim {g^{XY}_Z = u_X u_Y / u_Z, u_1 = 1}, which acts trivially on F.
 2. the two-dimensional loops E_A = Hom(A, A (x) V), A in {0011, 0012}: two vectors
    x^C_i = F[(A,V,V,C)][(A,i,0),(V,0,0)] (C a one-dimensional neighbour, fixed in the plan)
    transform as x -> g_E x / g^{VV}_V; we set g_E = g^{VV}_V [x^C1 x^C2]^-1.
 3. second legs S of depth >= 2, by depth: with P = parent(S) (N_{PV}^S = 1), the column block
    S of F[(A,P,V,C)] after the known row gauges is M0 (rows x N_{AS}^C); we set g^{AS}_C =
    M0[piv] for a pivot row set fixed in the plan, so that the new block has identity on piv.
Then F_can = g.F on the whole H family.  Gauge action (fusion-map convention):
  F'[(Z,i,j),(Z',k',l')] = sum g^{XY}_Z[i,a] g^{ZW}_T[j,b] F[(Z,a,b),(Z',k,l)] ginv^{YW}_{Z'}[k,k'] ginv^{XZ'}_T[l,l'].

usage:
  f4_canon.py plan OUT.plan.pkl p w             compute the plan at (p, zeta -> z^w)
  f4_canon.py run PLAN OUTDIR p w [w ...]       residues of F_can and g at the given roots
"""
import os, pickle, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
F4DIR = os.path.join(HERE, "..", "f4")
sys.path.insert(0, HERE)
sys.path.insert(1, F4DIR)
import fsym_bin                      # read only
from lin import Fp, primitive_root_of_unity
import torus

UNIT, V = "0000", "0001"
# the binary F-table read by readers(): $F4_ARCHIVE/F4_4_F_part{1,2}.bin (default: archive/ of the repository)
ARCHIVE = os.environ.get("F4_ARCHIVE", os.path.join(HERE, "..", "..", "archive"))
PARTS = [os.path.join(ARCHIVE, f"F4_4_F_part{i}.bin") for i in (1, 2)]
LOOPS = ("0011", "0012")


def readers():
    rs = [fsym_bin.Reader(p) for p in PARTS]
    where = {}
    for r in rs:
        for k in r.keys():
            where[k] = r
    return where


def cvals(p, w):
    z = primitive_root_of_unity(p, 52, w)
    zi = pow(z, -1, p)
    return [1] + [(pow(z, k, p) + pow(zi, k, p)) % p for k in range(1, 12)]


def load_H(p, w, where):
    cv = cvals(p, w)
    H = {}
    for k, r in where.items():
        if k[2] == V:
            rows, cols, _, _, _ = r.block(k)
            H[k] = (rows, cols, r.residues(k, p, cv))
    return H


def structure(H):
    N = {}
    for (X, Y, W, T), (rows, cols, M) in H.items():
        for (Z, i, j) in rows:
            N[(X, Y, Z)] = max(N.get((X, Y, Z), 0), i + 1)
            N[(Z, W, T)] = max(N.get((Z, W, T), 0), j + 1)
        for (Z, k, l) in cols:
            N[(Y, W, Z)] = max(N.get((Y, W, Z), 0), k + 1)
            N[(X, Z, T)] = max(N.get((X, Z, T), 0), l + 1)
    obj = sorted({t[0] for t in N})
    n = lambda a, b, c: N.get((a, b, c), 0)
    depth = {UNIT: 0}
    frontier = [UNIT]
    while frontier:
        nxt = []
        for A in frontier:
            for B in obj:
                if n(A, V, B) and B not in depth:
                    depth[B] = depth[A] + 1
                    nxt.append(B)
        frontier = nxt
    assert len(depth) == len(obj)
    parent = {}
    for S in obj:
        if depth[S] >= 2:
            parent[S] = min(P for P in obj if depth[P] == depth[S] - 1 and n(P, V, S) == 1)
    return N, obj, depth, parent


def torus_entries(H, N):
    """entries (key, row, col) of step 1 and their D-rows."""
    n = lambda a, b, c: N.get((a, b, c), 0)
    ents = []
    for (A, Y, W, C), (rows, cols, M) in sorted(H.items()):
        if Y != V or A == UNIT:
            continue
        for ci, (Zp, k, l) in enumerate(cols):
            if Zp not in (UNIT, V):
                continue
            if Zp == V and n(A, V, C) != 1:
                continue
            for ri, (Z, i, j) in enumerate(rows):
                if n(A, V, Z) == 1 and n(Z, V, C) == 1:
                    row = {}
                    for t, s in (((A, V, Z), 1), ((Z, V, C), 1), ((V, V, Zp), -1), ((A, Zp, C), -1)):
                        if UNIT in (t[0], t[1]):
                            continue
                        row[t] = row.get(t, 0) + s
                    row = {t: s for t, s in row.items() if s}
                    ents.append(((A, Y, W, C), ri, ci, row))
    return ents


def make_plan(H, F):
    N, obj, depth, parent = structure(H)
    ents = torus_entries(H, N)
    ents = [e for e in ents if H[e[0]][2][e[1], e[2]] % F.p]          # nonzero entries only
    variables = sorted({t for e in ents for t in e[3]})
    edge_vars = sorted(t for t in N if t[1] == V and t[0] != UNIT and N[t] == 1)
    assert set(variables) <= set(edge_vars)
    rows = [e[3] for e in ents]
    pivots, deferred = torus.unit_pivots(rows)
    rk = torus.rank_mod_p(rows, edge_vars)
    print(f"[plan] step 1: {len(edge_vars)} one-dim V-edge variables, {len(ents)} entries, rank {rk}, "
          f"unit pivots {len(pivots)}, deferred {len(deferred)}, residual torus dim {len(edge_vars) - rk}")
    assert len(pivots) in (rk, rk - 1)
    plan = dict(N=N, obj=obj, depth=depth, parent=parent, edge_vars=edge_vars,
                torus=[(ents[e][:3], c, r, {ents[ee][:3]: uu for ee, uu in U.items()}) for (e, c, r, U) in pivots],
                residual=len(edge_vars) - len(pivots), rank=rk)
    return plan


def apply_gauge(F, key, rows, cols, M, g):
    """F' = Lg M Rg with the index conventions of the module docstring."""
    X, Y, W, T = key
    p = F.p
    ri = {r: a for a, r in enumerate(rows)}
    ci = {c: a for a, c in enumerate(cols)}
    L = np.zeros((len(rows), len(rows)), dtype=np.int64)
    for a, (Z, i, j) in enumerate(rows):
        g1, g2 = g[(X, Y, Z)], g[(Z, W, T)]
        for aa in range(g1.shape[1]):
            for bb in range(g2.shape[1]):
                L[a, ri[(Z, aa, bb)]] = g1[i, aa] * g2[j, bb] % p
    R = np.zeros((len(cols), len(cols)), dtype=np.int64)
    for b, (Z, k, l) in enumerate(cols):
        h1, h2 = g[("inv", Y, W, Z)], g[("inv", X, Z, T)]
        for kk in range(h1.shape[1]):
            for ll in range(h2.shape[1]):
                R[b, ci[(Z, kk, ll)]] = h1[k, kk] * h2[l, ll] % p
    return F.mm(F.mm(L, M), R)


def setg(g, F, t, mat):
    g[t] = np.asarray(mat, dtype=np.int64) % F.p
    g[("inv",) + t] = F.inverse(g[t])


def canon(H, plan, F, make=False, lam=None):
    p = F.p
    N, obj, depth, parent = plan["N"], plan["obj"], plan["depth"], plan["parent"]
    g = {}
    for t, m in N.items():
        if UNIT in (t[0], t[1]):
            setg(g, F, t, np.eye(m, dtype=np.int64))
    # step 1
    val = {t: 1 for t in plan["edge_vars"]}
    if plan.get("free") is not None:
        if lam is None:
            lam = lam_from_extra(H, plan, F)
        if not (isinstance(plan["free"], tuple) and plan["free"][0] == "loop"):
            val[plan["free"]] = lam % p
    ent = lambda e: int(H[e[0]][2][e[1], e[2]])
    for (e, c, r, U) in reversed(plan["torus"]):
        beta = 1
        for ee, uu in U.items():
            x = ent(ee)
            assert x % p, "torus pivot entry vanishes mod p"
            beta = beta * pow(x, -uu, p) % p
        for cc, ss in r.items():
            if cc != c:
                beta = beta * pow(val[cc], -ss, p) % p
        val[c] = beta
    for t in plan["edge_vars"]:
        setg(g, F, t, [[val[t]]])
    # step 2: loops
    if make:
        plan["loops"] = {}
    for A in LOOPS:
        t = (A, V, A)
        assert N[t] == 2
        cands = []
        for C in obj:
            if C != A and N.get((A, V, C)) == 1 and (A, V, V, C) in H:
                rows, cols, M = H[(A, V, V, C)]
                ci = cols.index((V, 0, 0))
                s = pow(int(g[(V, V, V)][0, 0]), -1, p)
                x = [int(M[rows.index((A, i, 0)), ci]) * s % p for i in range(2)]
                cands.append((C, x))
        for C in obj:
            if C != A and N.get((C, V, A)) == 1 and (C, V, V, A) in H:
                rows, cols, M = H[(C, V, V, A)]
                if (V, 0, 0) not in cols:
                    continue
                ci = cols.index((V, 0, 0))
                s = pow(int(g[(V, V, V)][0, 0]), -1, p)
                x = [int(M[rows.index((A, 0, j)), ci]) * s % p for j in range(2)]
                cands.append(("in:" + C, x))
        if cands and (A, V, V, A) in H:
            # y_i = sum_{j,l} T_{ijl} u_j x_l with u = B^-1 x: covariant like x (paper, Section 4.7)
            rows, cols, M = H[(A, V, V, A)]
            sV = pow(int(g[(V, V, V)][0, 0]), -1, p)
            s1 = pow(int(g[(V, V, UNIT)][0, 0]), -1, p)
            B = np.array([[int(M[rows.index((A, i, j)), cols.index((UNIT, 0, 0))]) * s1 % p for j in range(2)] for i in range(2)], dtype=np.int64)
            T = [[[int(M[rows.index((A, i, j)), cols.index((V, 0, l))]) * sV % p for l in range(2)] for j in range(2)] for i in range(2)]
            x0 = np.array(cands[0][1], dtype=np.int64)
            if F.rank(B) == 2:
                u = F.mm(F.inverse(B), x0.reshape(2, 1)).ravel()
                y = [sum(T[i][j][l] * int(u[j]) % p * int(x0[l]) for j in range(2) for l in range(2)) % p for i in range(2)]
                cands.append(("BT:" + str(cands[0][0]), y))
                # Op[i,l] = sum_j T_ijl u_j is an endomorphism of the loop space (conjugated by the
                # gauge); x is an eigenvector of Op in the examples; the other eigenvector, with a
                # free scale fixed later by an extra entry (plan['free'] = ('loop', A)), is the
                # second basis vector when no two covariant vectors are independent.
                Op = np.array([[sum(T[i][j][l] * int(u[j]) for j in range(2)) % p for l in range(2)] for i in range(2)], dtype=np.int64)
                tr = int(Op[0, 0] + Op[1, 1]) % p
                Ox = F.mm(Op, x0.reshape(2, 1)).ravel()
                k0 = 0 if x0[0] % p else 1
                mu1 = int(Ox[k0]) * pow(int(x0[k0]), -1, p) % p
                mu2 = (tr - mu1) % p
                if mu1 != mu2 and all((int(Ox[i]) - mu1 * int(x0[i])) % p == 0 for i in range(2)):
                    Kmat = (Op - mu2 * np.eye(2, dtype=np.int64)) % p
                    nsp = F.nullspace(Kmat)
                    y2 = nsp[:, 0] % p
                    sc = lam if (isinstance(plan.get("free"), tuple) and plan["free"] == ("loop", A) and lam is not None) else 1
                    cands.append(("EIG:" + str(cands[0][0]), [int(v) * sc % p for v in y2]))
        if make:
            for a in range(len(cands)):
                for b in range(a + 1, len(cands)):
                    X = np.array([cands[a][1], cands[b][1]], dtype=np.int64).T
                    if F.rank(X) == 2:
                        plan["loops"][A] = (cands[a][0], cands[b][0])
                        break
                if A in plan["loops"]:
                    break
            print(f"[plan] loop {A}: candidate neighbours {[c for c, _ in cands]}, chosen {plan['loops'][A]}")
        C1, C2 = plan["loops"][A]
        x = dict(cands)[C1]; y = dict(cands)[C2]
        X = np.array([x, y], dtype=np.int64).T % p
        setg(g, F, t, F.inverse(X))
    # step 3
    if make:
        plan["piv"] = {}
    for S in sorted(obj, key=lambda s: (depth[s], s)):
        if depth[S] < 2:
            continue
        P = parent[S]
        for A in obj:
            if A == UNIT:
                continue
            for C in obj:
                m = N.get((A, S, C), 0)
                if not m:
                    continue
                key = (A, P, V, C)
                rows, cols, M = H[key]
                Mg = apply_gauge_partial(F, key, rows, cols, M, g)
                cS = [cols.index((S, 0, l)) for l in range(m)]
                M0 = Mg[:, cS] * pow(int(g[(P, V, S)][0, 0]), -1, p) % p
                if make:
                    _, piv = F.rref(M0.T.copy())
                    assert len(piv) == m
                    plan["piv"][(A, S, C)] = piv
                piv = plan["piv"][(A, S, C)]
                setg(g, F, (A, S, C), M0[piv, :])
    assert all(t in g for t in N)
    out = {}
    for key, (rows, cols, M) in H.items():
        out[key] = (rows, cols, apply_gauge(F, key, rows, cols, M, g))
    return g, out


def lam_from_extra(H, plan, F):
    """the value of the free variable for which F_can[extra] = 1: F_can depends on it monomially,
    F_can[e](lam) = lam^s F_can[e](1) with s = +-1 (checked in the plan)."""
    p = F.p
    (key, r, c), s = plan["extra"]
    _, out = canon(H, plan, F, lam=1)
    x = int(out[key][2][r, c])
    assert x % p
    return pow(x, -s, p)


def find_extra(H, plan, F):
    """choose the free variable and an entry of F_can depending on it with exponent +-1."""
    p = F.p
    pivcols = {c for (_, c, _, _) in plan["torus"]}
    nonpiv = [t for t in plan["edge_vars"] if t not in pivcols]
    loopfree = [("loop", A) for A, (C1, C2) in plan.get("loops", {}).items() if C2.startswith("EIG:")]
    if loopfree:
        nonpiv = loopfree
    for v in nonpiv:
        plan["free"] = v
        _, o1 = canon(H, plan, F, lam=1)
        _, o2 = canon(H, plan, F, lam=2)
        pw = {pow(2, s, p): s for s in range(-8, 9)}
        found = None
        exps = {}
        for key in sorted(o1):
            A, B = o1[key][2], o2[key][2]
            nz = np.argwhere(A % p)
            for (r, c) in nz:
                ratio = int(B[r, c]) * pow(int(A[r, c]), -1, p) % p
                s = pw.get(ratio)
                assert s is not None, "non-monomial dependence"
                exps[s] = exps.get(s, 0) + 1
                if found is None and abs(s) == 1:
                    found = ((key, int(r), int(c)), s)
        print(f"[plan] free variable {v}: exponent distribution {dict(sorted(exps.items()))}")
        if found:
            plan["extra"] = found
            return
    raise RuntimeError("no entry with exponent +-1: genuine square-root ambiguity")


def apply_gauge_partial(F, key, rows, cols, M, g):
    """row gauge only (columns untouched): Lg M."""
    X, Y, W, T = key
    p = F.p
    ri = {r: a for a, r in enumerate(rows)}
    L = np.zeros((len(rows), len(rows)), dtype=np.int64)
    for a, (Z, i, j) in enumerate(rows):
        g1, g2 = g[(X, Y, Z)], g[(Z, W, T)]
        for aa in range(g1.shape[1]):
            for bb in range(g2.shape[1]):
                L[a, ri[(Z, aa, bb)]] = g1[i, aa] * g2[j, bb] % p
    return F.mm(L, M)


def check_canon(out, plan, F):
    """the normalisation conditions hold in the output (a self-check)."""
    p = F.p
    N = plan["N"]
    for (e, c, r, U) in plan["torus"]:
        assert out[e[0]][2][e[1], e[2]] % p == 1
    for (A, S, C), piv in plan["piv"].items():
        P = plan["parent"][S]
        key = (A, P, V, C)
        rows, cols, M = out[key]
        m = N[(A, S, C)]
        cS = [cols.index((S, 0, l)) for l in range(m)]
        assert (M[np.ix_(piv, cS)] % p == np.eye(m, dtype=np.int64)).all()
    for A in LOOPS:
        C1, C2 = plan["loops"][A]
        for b, C in enumerate((C1, C2)):
            if ":" in C:
                continue
            rows, cols, M = out[(A, V, V, C)]
            assert [int(M[rows.index((A, i, 0)), cols.index((V, 0, 0))]) for i in range(2)] == [int(i == b) for i in range(2)]
    return True


def main(argv):
    mode = argv[1]
    where = readers()
    if mode == "plan":
        out, p, w = argv[2], int(argv[3]), int(argv[4])
        F = Fp(p)
        t0 = time.time()
        H = load_H(p, w, where)
        print(f"[plan] loaded {len(H)} H blocks in {time.time() - t0:.1f} s")
        plan = make_plan(H, F)
        plan["free"] = None
        g, Fc = canon(H, plan, F, make=True)
        if plan["residual"] == 16:
            find_extra(H, plan, F)
            print(f"[plan] extra normalisation entry {plan['extra']}")
            g, Fc = canon(H, plan, F)
        check_canon(Fc, plan, F)
        with open(out, "wb") as f:
            pickle.dump(plan, f)
        print(f"[plan] wrote {out}; {len(plan['piv'])} echelon spaces; {time.time() - t0:.1f} s")
    elif mode == "run":
        plan = pickle.load(open(argv[2], "rb"))
        outdir, p = argv[3], int(argv[4])
        F = Fp(p)
        os.makedirs(outdir, exist_ok=True)
        for w in map(int, argv[5:]):
            t0 = time.time()
            H = load_H(p, w, where)
            g, Fc = canon(H, plan, F)
            check_canon(Fc, plan, F)
            keys = sorted(Fc)
            flatF = np.concatenate([Fc[k][2].ravel() for k in keys]).astype(np.int64)
            gkeys = sorted(t for t in g if t[0] != "inv")
            flatg = np.concatenate([g[t].ravel() for t in gkeys]).astype(np.int64)
            flatH = np.concatenate([H[k][2].ravel() for k in keys]).astype(np.int64)
            np.savez(os.path.join(outdir, f"canon_{p}_{w}.npz"), F=flatF, g=flatg, H=flatH)
            if not os.path.exists(os.path.join(outdir, "layout.pkl")):
                with open(os.path.join(outdir, "layout.pkl"), "wb") as f:
                    pickle.dump(dict(keys=keys, shapes=[Fc[k][2].shape for k in keys],
                                     rows=[H[k][0] for k in keys], cols=[H[k][1] for k in keys],
                                     gkeys=gkeys, gshapes=[g[t].shape for t in gkeys]), f)
            print(f"[run] p={p} w={w}: {len(flatF)} entries, {len(flatg)} gauge entries, {time.time() - t0:.1f} s", flush=True)


if __name__ == "__main__":
    main(sys.argv)
