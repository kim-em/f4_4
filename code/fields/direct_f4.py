#!/usr/bin/env python3
"""Direct computation of C(F4,4) in the normal-form gauge, modulo p at one root, using the
construction scripts of ../f4 (read only; called as subprocesses with outputs in WORK).

For (p, w), zeta -> z^w with z = g^((p-1)/52) (g the least primitive root):
  1. ../f4/cells.py F4 4 p w                      quantum group cells (construction gauge)
  2. ../f4/run_h.py                               H family G[(X,Y,V,T)], splitting convention
  3. normal form (f4_canon.canon) on F = (G^-1)^T of the H family; the index choices are
     recomputed at this (p, w) by the rule of the paper (Section 4.7, "The normal form") and
     asserted equal to the stored plan data/F4_4_canon.plan.pkl; then G_nf = (F_nf^-1)^T replaces H
  4. ../f4/run_split.py g1 / x                    the pentagon recursion, run once, on G_nf:
     it only solves pentagons, so its output is the whole table in the normal-form gauge
  5. export convention F = (G^-1)^T, written in the key order of the binary tables
     (sorted keys, rows and columns sorted by (object index, i, j)) as WORK/F_<p>_<w>.npy
  6. ../f4/run_r.py; R in the export convention (R = Rs^T), regauged by the normal-form gauge g
     of step 3 (R_nf = g^{YX} R (g^{XY})^-1) and its inverse, as WORK/R_<p>_<w>.npz; the same at
     the conjugate root 52 - w (cells and R only: F and g are real), WORK/R_<p>_<52-w>.npz
No exact data and no exact gauge transformation is used; step 3 and 6 are gauge changes modulo p.

usage: direct_f4.py p w WORK [--nproc n] [--keep]
  WORK may be relative; it is made absolute (the ../f4 scripts run with working directory ../f4).
"""
import os, pickle, subprocess, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
F4 = os.path.join(HERE, "..", "f4")
sys.path.insert(0, HERE)
sys.path.insert(1, F4)
import f4_canon as fc
from lin import Fp

OBJ = ["0000", "0001", "0002", "0010", "1000", "0003", "0011", "0100", "1001", "0004", "0012",
       "0020", "0101", "1002", "1010", "2000"]          # the alcove order of ../f4 (labels are sorted by it)
PLAN = os.path.join(HERE, "data", "F4_4_canon.plan.pkl")


def nm(w):
    return "".join(map(str, w))


def run(cmd, log):
    t = time.time()
    with open(log, "a") as f:
        subprocess.run(cmd, cwd=F4, check=True, stdout=f, stderr=subprocess.STDOUT)
    return time.time() - t


def order_key(a):
    return (OBJ.index(a[0]),) + tuple(a[1:])


def main(argv):
    p, w, work = int(argv[1]), int(argv[2]), os.path.abspath(argv[3])
    nproc = int(argv[argv.index("--nproc") + 1]) if "--nproc" in argv else 4
    os.makedirs(work, exist_ok=True)
    S = os.path.join(work, f"F4_{p}_{w}")
    log = S + ".log"
    T = {}
    py = sys.executable
    T["cells"] = run([py, "cells.py", "F4", "4", str(p), str(w), S + ".cells.pkl", "0,0,0,1"], log)
    T["H"] = run([py, "run_h.py", S + ".cells.pkl", S + ".H.pkl", str(nproc)], log)
    # ---- normal form of the H family
    t0 = time.time()
    F = Fp(p)
    Hd = pickle.load(open(S + ".H.pkl", "rb"))
    H = {}
    back = {}
    for key, (rl, cl, G) in Hd["G"].items():
        k2 = tuple(nm(x) for x in key)
        rl2 = [(nm(a[0]),) + tuple(a[1:]) for a in rl]
        cl2 = [(nm(b[0]),) + tuple(b[1:]) for b in cl]
        ro = sorted(range(len(rl2)), key=lambda i: order_key(rl2[i]))
        co = sorted(range(len(cl2)), key=lambda j: order_key(cl2[j]))
        Fm = F.inverse(np.asarray(G, dtype=np.int64) % p).T
        H[k2] = ([rl2[i] for i in ro], [cl2[j] for j in co], Fm[np.ix_(ro, co)] % p)
        back[k2] = (key, [rl[i] for i in ro], [cl[j] for j in co])
    # unit legs of the H family are not produced by run_h; add identities (strict unit)
    plan = pickle.load(open(PLAN, "rb"))
    fresh = fc.make_plan(H, F)
    fresh["free"] = None
    fc.canon(H, fresh, F, make=True)
    same = all(fresh[k] == plan[k] for k in ("torus", "loops", "piv", "parent", "edge_vars"))
    assert same, "the plan rule gives different choices at this (p, w): bad prime"
    g, Hn = fc.canon(H, plan, F)
    G2 = {}
    for k2, (rows, cols, M) in Hn.items():
        key, rl, cl = back[k2]
        G2[key] = (rl, cl, F.inverse(M % p).T % p)
    Hd2 = dict(Hd, G=G2)
    pickle.dump(Hd2, open(S + ".Hnf.pkl", "wb"))
    T["normal form"] = time.time() - t0
    # ---- recursion in the normal-form gauge
    T["g1"] = run([py, "run_split.py", "g1", S + ".cells.pkl", S + ".Hnf.pkl", S + ".G1.pkl"], log)
    groups = [OBJ[i::nproc] for i in range(nproc)]
    t0 = time.time()
    procs = []
    for i, grp in enumerate(groups):
        if not grp:
            continue
        cmd = [py, "run_split.py", "x", S + ".cells.pkl", S + ".G1.pkl", ",".join(grp), f"{S}.x{i}.npy", "--structure", f"{S}.x{i}.st.pkl"]
        procs.append(subprocess.Popen(cmd, cwd=F4, stdout=open(log, "a"), stderr=subprocess.STDOUT))
    for pr in procs:
        assert pr.wait() == 0
    T["recursion x"] = time.time() - t0
    # ---- export convention, binary key order
    t0 = time.time()
    blocks = {}
    for i, grp in enumerate(groups):
        if not grp:
            continue
        a = np.load(f"{S}.x{i}.npy").astype(np.int64)
        st = pickle.load(open(f"{S}.x{i}.st.pkl", "rb"))
        pos = 0
        for key, rl, cl in st:
            n = len(rl)
            Gm = a[pos:pos + n * n].reshape(n, n); pos += n * n
            Fm = F.inverse(Gm % p).T
            k2 = tuple(nm(x) for x in key)
            rl2 = [(nm(x[0]),) + tuple(x[1:]) for x in rl]
            cl2 = [(nm(x[0]),) + tuple(x[1:]) for x in cl]
            ro = sorted(range(n), key=lambda i: order_key(rl2[i]))
            co = sorted(range(n), key=lambda j: order_key(cl2[j]))
            blocks[k2] = Fm[np.ix_(ro, co)] % p
    keys = sorted(blocks)
    np.save(f"{work}/F_{p}_{w}.npy", np.concatenate([blocks[k].ravel() for k in keys]))
    pickle.dump(keys, open(f"{work}/F_{p}_{w}.npy.keys.pkl", "wb"))
    T["export"] = time.time() - t0
    # ---- R in the normal-form gauge
    T["R"] = run([py, "run_r.py", S + ".cells.pkl", S + ".R.pkl"], log)
    Rd = pickle.load(open(S + ".R.pkl", "rb"))["G"]
    Rn, Rin = {}, {}
    for key, (rl, cl, M) in Rd.items():
        X, Y, Z = (nm(x) for x in key)
        Ru = np.asarray(M, dtype=np.int64).T % p                 # export R[j, i]
        Rc = F.mm(F.mm(g[(Y, X, Z)], Ru), g[("inv", X, Y, Z)])
        Rn[(X, Y, Z)] = Rc
        Rin[(X, Y, Z)] = F.inverse(Rc)
    tk = sorted(Rn)
    np.savez(f"{work}/R_{p}_{w}.npz", R=np.concatenate([Rn[t].ravel() for t in tk]), Rinv=np.concatenate([Rin[t].ravel() for t in tk]))
    pickle.dump(tk, open(f"{work}/R_{p}_{w}.keys.pkl", "wb"))
    # R at the conjugate root 52 - w (F and the gauge g are real, so only cells and R are needed)
    wc = 52 - w
    Sc = os.path.join(work, f"F4_{p}_{wc}")
    T["R conj"] = run([py, "cells.py", "F4", "4", str(p), str(wc), Sc + ".cells.pkl", "0,0,0,1"], log)
    T["R conj"] += run([py, "run_r.py", Sc + ".cells.pkl", Sc + ".R.pkl"], log)
    Rd = pickle.load(open(Sc + ".R.pkl", "rb"))["G"]
    Rn, Rin = {}, {}
    for key, (rl, cl, M) in Rd.items():
        X, Y, Z = (nm(x) for x in key)
        Rc = F.mm(F.mm(g[(Y, X, Z)], np.asarray(M, dtype=np.int64).T % p), g[("inv", X, Y, Z)])
        Rn[(X, Y, Z)] = Rc
        Rin[(X, Y, Z)] = F.inverse(Rc)
    np.savez(f"{work}/R_{p}_{wc}.npz", R=np.concatenate([Rn[t].ravel() for t in tk]), Rinv=np.concatenate([Rin[t].ravel() for t in tk]))
    pickle.dump(tk, open(f"{work}/R_{p}_{wc}.keys.pkl", "wb"))
    for suf in (".cells.pkl", ".R.pkl"):
        if os.path.exists(Sc + suf) and "--keep" not in argv:
            os.remove(Sc + suf)
    if "--keep" not in argv:
        for suf in [".cells.pkl", ".H.pkl", ".Hnf.pkl", ".G1.pkl", ".R.pkl"] + [f".x{i}.npy" for i in range(nproc)] + [f".x{i}.st.pkl" for i in range(nproc)]:
            if os.path.exists(S + suf):
                os.remove(S + suf)
    print(f"[direct] p={p} w={w}: {len(keys)} blocks, {sum(b.size for b in blocks.values())} entries; plan rule reproduced; "
          + ", ".join(f"{k} {v:.0f}s" for k, v in T.items()), flush=True)


if __name__ == "__main__":
    main(sys.argv)
