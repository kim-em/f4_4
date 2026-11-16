"""Compute the whole H family G[(X, Y, V, T)] modulo p, in parallel over X.

usage: run_h.py cells.pkl out.pkl [nproc]
"""
import sys, pickle, time, os
import multiprocessing as mp
from rootsys import RootSystem
from objects import Objects
from hfamily import HFamily

cells_path, out = sys.argv[1], sys.argv[2]
nproc = int(sys.argv[3]) if len(sys.argv) > 3 else 16
c = pickle.load(open(cells_path, "rb"))
R = RootSystem(c["type"])
alc = c["alcove"]
h = {w: int(R.casimir(w)) for w in alc}
fus = {(X, S, T): m for X in alc for S in alc for T, m in R.tensor(X, S, c["k"]).items()}


def work(X):
    t = time.time()
    O = Objects(c, h, fus)
    H = HFamily(O, seed=alc.index(X) + 1)
    res = {}
    for Y in alc:
        for T in alc:
            g = H.G(X, Y, T)
            if g is not None:
                res[(X, Y, c["gen"], T)] = g
    bases = {key: v[0] for key, v in O._basis.items()}
    print(f"[run_h] base {X}: {len(res)} matrices, {time.time()-t:.0f}s", flush=True)
    return res, bases


if __name__ == "__main__":
    t = time.time()
    with mp.get_context("fork").Pool(nproc) as pool:
        parts = pool.map(work, sorted(alc, key=lambda X: -fus.get((X, X, X), 0)))
    G, bases = {}, {}
    for r, b in parts:
        G.update(r)
        bases.update(b)
    pickle.dump(dict(G=G, bases=bases, p=c["p"], q=c["q"], which=c["which"]), open(out, "wb"))
    print(f"[run_h] wrote {out}: {len(G)} matrices, {sum(v[2].size for v in G.values())} entries, {time.time()-t:.0f}s")
