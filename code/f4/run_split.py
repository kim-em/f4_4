"""The pentagon recursion of recursion.py (class Recursion), split by the first leg X.

Step 2 for G[(X, T, W, E)] only uses H, the second-leg-V family G1 (step 1, all first legs) and
blocks G[(X, ., ., .)] with the same X, so the table can be computed one X at a time.

usage:
  run_split.py g1 cells.pkl H.pkl out.G1.pkl
      step 1 (all X); writes H, unit legs and G1.
  run_split.py x cells.pkl G1.pkl X1,X2,... out.npy [--structure s.pkl]
      step 2 for the listed first legs (Dynkin labels like 0011); writes the residues of all
      G[(X, Y, W, T)] for those X as int32, blocks in sorted key order, each block row-major,
      labels in Recursion.canonical order; --structure also writes the block list.
"""
import sys, pickle, time
import numpy as np
from rootsys import RootSystem
from objects import Objects
from recursion import Recursion


def setup(cells_path, Gdict):
    c = pickle.load(open(cells_path, "rb"))
    R = RootSystem(c["type"])
    alc = c["alcove"]
    fus = {(X, S, T): m for X in alc for S in alc for T, m in R.tensor(X, S, c["k"]).items()}
    h = {w: int(R.casimir(w)) for w in alc}
    O = Objects(c, h, fus)
    depth = O.length
    parent = {S: O.gamma[S][-2] for S in alc if S != alc[0]}
    REC = Recursion(alc, fus, Gdict, c["p"], alc[0], c["gen"], depth)
    return c, alc, fus, depth, parent, REC


def keys_for(alc, fus, X):
    return sorted((X, b, cc, d) for b in alc for cc in alc for d in alc
                  if sum(fus.get((X, b, e), 0) * fus.get((e, cc, d), 0) for e in alc))


if __name__ == "__main__":
    mode = sys.argv[1]
    t0 = time.time()
    if mode == "g1":
        H = pickle.load(open(sys.argv[3], "rb"))
        c, alc, fus, depth, parent, REC = setup(sys.argv[2], dict(H["G"]))
        gen = c["gen"]
        nk = 0
        for d in range(2, max(depth.values()) + 1):
            new = {}
            for X in alc:
                for Rr in [w for w in alc if depth[w] == d]:
                    for E in alc:
                        for R2, val in REC.step1(X, parent[Rr], E).items():
                            key = (X, gen, R2, E)
                            if R2 == Rr:
                                new[key] = REC.canonical(val)
                            elif key in REC.G:
                                assert REC.compare(key, val), ("step 1 disagrees", key)
                                nk += 1
            REC.G.update(new)
        pickle.dump(dict(G=REC.G, p=c["p"], q=c["q"], which=c["which"]), open(sys.argv[4], "wb"))
        print(f"[split g1] {len(REC.G)} blocks (H, unit legs, G1), {nk} known blocks reproduced, {time.time()-t0:.0f}s")
    elif mode == "x":
        G1 = pickle.load(open(sys.argv[3], "rb"))
        c, alc, fus, depth, parent, REC = setup(sys.argv[2], dict(G1["G"]))
        byname = {"".join(map(str, w)): w for w in alc}
        Xs = [byname[x] for x in sys.argv[4].split(",")]
        out = sys.argv[5]
        nk = 0
        for d in range(2, max(depth.values()) + 1):
            parents = sorted({parent[T] for T in alc if depth[T] == d})
            for X in Xs:
                new = {}
                for Yp in parents:
                    for W in alc:
                        for E in alc:
                            for T2, val in REC.step2(X, Yp, W, E).items():
                                key = (X, T2, W, E)
                                if depth[T2] == d and parent[T2] == Yp:
                                    new[key] = REC.canonical(val)
                                elif key in REC.G:
                                    assert REC.compare(key, val), ("step 2 disagrees", key)
                                    nk += 1
                REC.G.update(new)
        keys = [k for X in Xs for k in keys_for(alc, fus, X)]
        missing = [k for k in keys if k not in REC.G]
        assert not missing, missing[:3]
        flat = np.concatenate([REC.G[k][2].reshape(-1) for k in keys]).astype(np.int32)
        np.save(out, flat)
        if "--structure" in sys.argv:
            st = [(k, REC.G[k][0], REC.G[k][1]) for k in keys]
            pickle.dump(st, open(sys.argv[sys.argv.index("--structure") + 1], "wb"))
        print(f"[split x] {','.join(sys.argv[4].split(','))}: {len(keys)} blocks, {flat.size} entries, "
              f"{nk} recomputed known blocks agree, {time.time()-t0:.0f}s")
