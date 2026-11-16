"""Local ("cell") data of C(g,k) for the generator V, modulo p, from the quantum group.

For simple objects A, C the path space P(A -> C, 2) has basis the pairs of edges
A -> B -> C (B in the alcove), an edge B -> C being a basis vector of Hom(C, B (x) V).
This script computes the action on these path spaces of the braiding c_{V,V}, the cup, the
cap, and the trivalent split V -> V (x) V and merge V (x) V -> V, all modulo negligible
morphisms.

Method.  A highest weight vector y of weight C in A (x) V (x) V is determined by its
component tau(y) along the top vector of A, which lies in (V (x) V)_{C-A}; tau commutes with
every operator acting on the last two factors.  A (x) V is semisimple at the root (all its
constituents lie in the closed alcove), so the space of such y is the direct sum of
  paths:      (e1 (x) 1) e2  with e1: B -> A (x) V, e2: C -> B (x) V,  B in the open alcove,
  negligible: (e  (x) 1) f   with e: lam -> A (x) V, lam on the first affine wall.
An operator applied to a path vector is expanded in this basis and the wall part dropped.
The tau-image of (e1 (x) 1) x is computed with the U^- module map tau o e1 : L(B) -> V.

Gauge: an edge basis of Hom(C, B (x) V) is the reduced echelon basis of the tau-images in
V_{C-B} (dimension at most 2).

usage: cells.py TYPE LEVEL p which out.pkl [generator]   (q = (primitive root)^which)
"""
import sys, time, pickle, itertools
import numpy as np
from rootsys import RootSystem
from qg import QG, Module, TensorV, uminus_map, wadd
from lin import primitive_root_of_unity


def build(typ, k, p, which, gen=None, verbose=True):
    R = RootSystem(typ)
    order = 2 * R.m * (k + R.hdual)
    q = primitive_root_of_unity(p, order, which)
    G = QG(R, p, q, order)
    F = G.F
    r = R.rank
    if gen is None:
        gen = {"F4": (0, 0, 0, 1), "G2": (1, 0)}.get(typ) or tuple([1] + [0] * (r - 1))
    alc = R.alcove(k)
    fus = {}
    for a in alc:
        for b, m in R.tensor(a, gen, k).items():
            fus[(a, b)] = m
    t0 = time.time()
    V = Module(G, gen)
    TV = TensorV(G, V)
    n = TV.n
    depth = 2 * int(sum(R.root_coords(R.dominant_rep(gen))) ) + 2   # enough below every top
    # weights met: alcove and the classical constituents of A (x) V
    need = set(alc)
    for a in alc:
        need |= set(R.tensor(a, gen).keys())
    wall = sorted(w for w in need if w not in alc)
    mods = {}
    for lam in sorted(need, key=lambda w: (R.level(w), w)):
        mods[lam] = Module(G, lam, depth)
    if verbose:
        print(f"[cells] p={p} q={q}: built {len(mods)} modules ({len(wall)} outside the open alcove) in {time.time()-t0:.1f}s", flush=True)

    # global V basis: index off[w] + a
    off = TV.off
    vw = list(off)

    # full matrices of E_i, F_i, K_i on V
    Ef, Ff, Kf = [], [], []
    for i in range(r):
        E = np.zeros((n, n), dtype=np.int64)
        Fm = np.zeros((n, n), dtype=np.int64)
        K = np.zeros(n, dtype=np.int64)
        for w in V.dim:
            K[off[w]:off[w] + V.dim[w]] = G.K(i, w)
            up = wadd(w, G.alpha[i])
            if up in V.dim:
                E[off[up]:off[up] + V.dim[up], off[w]:off[w] + V.dim[w]] = V.Emat(i, w)
            dn = wadd(w, G.alpha[i], -1)
            if dn in V.dim:
                Fm[off[dn]:off[dn] + V.dim[dn], off[w]:off[w] + V.dim[w]] = V.Fmat(i, w)
        Ef.append(E); Ff.append(Fm); Kf.append(K)
    I = np.eye(n, dtype=np.int64)
    Kinv = [np.array([pow(int(x), -1, p) for x in K], dtype=np.int64) for K in Kf]
    # coproducts on V (x) V  (index i1*n + i2)
    import scipy.sparse as sps
    Is = sps.identity(n, dtype=np.int64, format="csr")
    dF = []
    for i in range(r):
        M = (sps.kron(sps.csr_matrix(Ff[i]), sps.csr_matrix(np.diag(Kinv[i]))) + sps.kron(Is, sps.csr_matrix(Ff[i]))).tocsr()
        M = M.astype(np.int64)
        M.data %= p
        dF.append(M)

    def spmm(M, x):
        hi = (M @ (x >> 15)) % p
        lo = (M @ (x & 0x7FFF)) % p
        return ((hi << 15) + lo) % p

    def vec_of(d, Cw):
        """dict mu -> (dim V_mu x dim V_{C-mu}) into a vector of V (x) V."""
        out = np.zeros(n * n, dtype=np.int64)
        for mu, x in d.items():
            w2 = wadd(Cw, mu, -1)
            a0, b0 = off[mu], off[w2]
            for a in range(x.shape[0]):
                out[(a0 + a) * n + b0:(a0 + a) * n + b0 + x.shape[1]] = x[a]
        return out

    # ---- decomposition of V (x) V ----
    VV = R.tensor(gen, gen)
    cols, labels = [], []
    hwvec = {}
    for W, mult in sorted(VV.items(), key=lambda t: (R.level(t[0]), t[0])):
        vecs, tau = TV.hw(V, W)
        assert len(vecs) == mult == 1, (W, len(vecs), mult)
        x = vec_of(vecs[0], W)
        hwvec[W] = x
        MW = Module(G, W)
        img = uminus_map(MW, x, lambda j, src, c: spmm(dF[j], c), list(MW.dim))
        for mu in sorted(MW.dim, key=lambda m: (MW.depth[m], m)):
            for c in range(img[mu].shape[1]):
                cols.append(img[mu][:, c]); labels.append((W, mu, c))
    Bfull = np.array(cols).T % p
    assert Bfull.shape == (n * n, n * n)
    # Bfull is block diagonal for the weight grading: invert weight space by weight space
    wt_of = {}
    for w1 in V.dim:
        for a in range(V.dim[w1]):
            for w2 in V.dim:
                for b in range(V.dim[w2]):
                    wt_of[(off[w1] + a) * n + off[w2] + b] = wadd(w1, w2)
    Binv = np.zeros_like(Bfull)
    for nu in set(wt_of.values()):
        ri = [t for t in range(n * n) if wt_of[t] == nu]
        ci = [t for t, l in enumerate(labels) if l[1] == nu]
        assert len(ri) == len(ci)
        Binv[np.ix_(ci, ri)] = F.inverse(Bfull[np.ix_(ri, ci)])
    xr = np.random.default_rng(0).integers(0, p, (n * n, 2))
    assert np.array_equal(F.mm(Bfull, F.mm(Binv, xr)), xr % p)
    sign = classical_signs(R, gen, VV)
    cV = R.casimir(gen)
    Rhat = np.zeros((n * n, n * n), dtype=np.int64)
    eig, ecol = {}, np.zeros(n * n, dtype=np.int64)
    for W in VV:
        idx = [t for t, l in enumerate(labels) if l[0] == W]
        ex = (R.casimir(W) - 2 * cV)
        assert ex % 2 == 0
        e = sign[W] * G.qp(ex // 2) % p
        eig[W] = (sign[W], ex // 2)
        ecol[idx] = e
    # Rhat = Bfull diag(e) Binv, weight space by weight space
    for nu in set(wt_of.values()):
        ri = [t for t in range(n * n) if wt_of[t] == nu]
        ci = [t for t, l in enumerate(labels) if l[1] == nu]
        Rhat[np.ix_(ri, ri)] = F.mm(Bfull[np.ix_(ri, ci)] * ecol[ci][None, :] % p, Binv[np.ix_(ci, ri)])
    zero = tuple([0] * r)
    cup = hwvec[zero]
    i0 = [t for t, l in enumerate(labels) if l[0] == zero]
    ev = Binv[i0[0], :]                      # ev(cup) = 1
    iv = [t for t, l in enumerate(labels) if l[0] == gen]
    tmat = Bfull[:, iv]                      # V -> V (x) V, columns in V-basis order of labels
    mmat = Binv[iv, :]                       # V (x) V -> V, m t = 1
    # order of labels for gen block vs V basis order: build permutation
    perm = [off[l[1]] + l[2] for t, l in enumerate(labels) if l[0] == gen]
    tm = np.zeros((n * n, n), dtype=np.int64); tm[:, perm] = tmat
    mm_ = np.zeros((n, n * n), dtype=np.int64); mm_[perm, :] = mmat
    tmat, mmat = tm, mm_
    # zigzags
    z1 = zigzag(F, ev, cup, n, left=True)
    z2 = zigzag(F, ev, cup, n, left=False)
    if verbose:
        print("[cells] V(x)V:", {str(W): eig[W] for W in eig}, " zigzags (scalars):", z1, z2, flush=True)
    # Yang-Baxter check on a random vector of V^3
    rng = np.random.default_rng(1)
    x = rng.integers(0, p, n ** 3)
    R12 = lambda y: F.mm(Rhat, y.reshape(n * n, n)).reshape(-1)
    R23 = lambda y: F.mm(y.reshape(n, n * n), Rhat.T).reshape(-1)
    ybe = np.array_equal(R12(R23(R12(x))), R23(R12(R23(x))))
    if verbose:
        print("[cells] Yang-Baxter on V^3:", ybe, flush=True)
    assert ybe

    # ---- edges ----
    edges = {}          # (A, B) -> (list of full hw vectors in L(A)(x)V, tau matrix n x k)
    for A in alc:
        for B, mult in R.tensor(A, gen).items():
            vecs, tau = TV.hw(mods[A], B)
            assert len(vecs) == mult, (A, B, len(vecs), mult)
            if B in alc:
                assert mult == fus.get((A, B), 0), (A, B, mult, fus.get((A, B)))
            edges[(A, B)] = (vecs, tau)
    if verbose:
        print(f"[cells] edges done {time.time()-t0:.1f}s", flush=True)

    def tau_map(A, B, kk, weights):
        """tau o e_kk : L(B) -> V as dict mu -> (dim V_{mu-A} x dim L(B)_mu)."""
        tau = edges[(A, B)][1][:, kk]
        w0 = wadd(B, A, -1)
        u = tau[off[w0]:off[w0] + V.dim[w0]]
        return uminus_map(mods[B], u, lambda j, src, c: F.mm(V.Fmat(j, wadd(src, A, -1)), c), weights)

    def tau_compose(A, B, kk, x, Cw):
        """tau((e_kk (x) 1) x) for x a hw vector of weight Cw in L(B) (x) V (dict)."""
        tm_ = tau_map(A, B, kk, list(x.keys()))
        out = np.zeros(n * n, dtype=np.int64)
        for mu, blk in x.items():
            y = F.mm(tm_[mu], blk)            # (dim V_{mu-A}) x (dim V_{C-mu})
            w1, w2 = wadd(mu, A, -1), wadd(Cw, mu, -1)
            for a in range(y.shape[0]):
                s = (off[w1] + a) * n + off[w2]
                out[s:s + y.shape[1]] = (out[s:s + y.shape[1]] + y[a]) % p
        return out

    paths, sigma, cupc, capc, splitc, mergec, basisdim = {}, {}, {}, {}, {}, {}, {}
    for A in alc:
        Bs = [B for B in R.tensor(A, gen)]
        targets = set()
        for B in Bs:
            targets |= set(R.tensor(B, gen))
        for C in sorted(targets, key=lambda w: (R.level(w), w)):
            if C not in alc:
                continue
            P, cols_p = [], []
            for B in Bs:
                if B not in alc or C not in R.tensor(B, gen):
                    continue
                vecs2 = edges[(B, C)][0]
                for i1 in range(len(edges[(A, B)][0])):
                    for i2, x in enumerate(vecs2):
                        P.append((B, i1, i2))
                        cols_p.append(tau_compose(A, B, i1, x, C))
            cols_n = []
            for B in Bs:
                if B in alc or C not in R.tensor(B, gen):
                    continue
                for i1 in range(len(edges[(A, B)][0])):
                    for x in edges[(B, C)][0] if (B, C) in edges else hw_cache(TV, mods, B, C):
                        cols_n.append(tau_compose(A, B, i1, x, C))
            if not P:
                continue
            Bm = np.array(cols_p + cols_n).T % p
            npth = len(P)
            classical = sum(m1 * R.tensor(B, gen).get(C, 0) for B, m1 in R.tensor(A, gen).items())
            rk = F.rank(Bm)
            assert rk == Bm.shape[1] == classical, (A, C, rk, Bm.shape, classical)
            paths[(A, C)] = P
            basisdim[(A, C)] = (npth, len(cols_n))
            coords = lambda Y: F.solve(Bm, Y % p)[:npth]
            sigma[(A, C)] = coords(F.mm(Rhat, Bm[:, :npth]))
            if C == A:
                cupc[A] = coords(cup.reshape(-1, 1))[:, 0]
                capc[A] = F.mm(ev.reshape(1, -1), Bm[:, :npth])[0]
            # merge: tau(m y) in V_{C-A}, expressed in the edge basis of Hom(C, A (x) V)
            if (A, C) in edges and C in alc:
                tauE = edges[(A, C)][1]
                mergec[(A, C)] = F.solve(tauE, F.mm(mmat, Bm[:, :npth]))
                splitc[(A, C)] = coords(F.mm(tmat, tauE))
        if verbose:
            print(f"[cells] base {A} done {time.time()-t0:.1f}s", flush=True)
    return dict(type=typ, k=k, p=p, q=q, which=which, order=order, gen=gen, alcove=alc, fus=fus,
                paths=paths, sigma=sigma, cup=cupc, cap=capc, split=splitc, merge=mergec,
                basisdim=basisdim, eig=eig, zigzag=(z1, z2),
                edges={key: len(v[0]) for key, v in edges.items() if key[1] in alc})


_hwc = {}


def hw_cache(TV, mods, B, C):
    if (B, C) not in _hwc:
        _hwc[(B, C)] = TV.hw(mods[B], C)[0]
    return _hwc[(B, C)]


def zigzag(F, ev, cup, n, left=True):
    """(ev (x) 1)(1 (x) cup) [left] or (1 (x) ev)(cup (x) 1) as a scalar (asserted scalar)."""
    p = F.p
    E = ev.reshape(n, n)
    Cc = cup.reshape(n, n)
    if left:
        # x_a -> sum_{b,c} ev[a,b] cup[b,c] x_c
        M = F.mm(E, Cc).T
    else:
        M = F.mm(Cc, E)
    s = M[0, 0]
    assert np.array_equal(M % p, (s * np.eye(n, dtype=np.int64)) % p), "zigzag not scalar"
    return int(s)


def classical_signs(R, gen, VV):
    """+1 on constituents of S^2(gen), -1 on Lambda^2(gen) (classical characters)."""
    ch = R.character(gen)
    # character of S^2 and Lambda^2 : weights mu+nu (mu<=nu ordered) and 2mu
    wl = []
    for mu, m in ch.items():
        wl += [mu] * m
    sym, alt = {}, {}
    for a in range(len(wl)):
        for b in range(a, len(wl)):
            w = wadd(wl[a], wl[b])
            sym[w] = sym.get(w, 0) + 1
            if b > a:
                alt[w] = alt.get(w, 0) + 1
    out = {}
    rem = dict(alt)
    for W in sorted(VV, key=lambda w: -sum(R.root_coords(w))):
        chW = R.character(W)
        if all(rem.get(w, 0) >= m for w, m in chW.items()):
            out[W] = -1
            for w, m in chW.items():
                rem[w] -= m
        else:
            out[W] = 1
    assert all(v == 0 for v in rem.values())
    return out


if __name__ == "__main__":
    typ, k, p, which, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
    gen = tuple(int(x) for x in sys.argv[6].split(",")) if len(sys.argv) > 6 else None
    res = build(typ, k, p, which, gen)
    pickle.dump(res, open(out, "wb"))
    print("[cells] wrote", out)
