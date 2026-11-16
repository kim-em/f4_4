"""Consistency checks of the cell data in the path model, at every base object.

usage: check_cells.py cells.pkl
"""
import sys, pickle
import numpy as np
from pathmodel import Model
from lin import Fp

c = pickle.load(open(sys.argv[1], "rb"))
p = c["p"]
F = Fp(p)
M = Model(c)
alc = c["alcove"]
q, order = c["q"], c["order"]
qp = lambda e: pow(q, e % order, p)
eigs = {W: (s * qp(e)) % p for W, (s, e) in c["eig"].items()}
zero = tuple([0] * len(alc[0]))
gen = c["gen"]


def dense(X, n, word):
    P = M.P(X, n)
    I = np.eye(len(P), dtype=np.int64)
    return M.word(X, n, word, I)[0]


def scalar(A):
    s = A[0, 0] if A.size else 0
    return int(s) if np.array_equal(A % p, (s * np.eye(A.shape[0], dtype=np.int64)) % p) else None


fails = 0
results = {}
for X in alc:
    # braid relation
    a = dense(X, 3, [("sigma", 1), ("sigma", 2), ("sigma", 1)])
    b = dense(X, 3, [("sigma", 2), ("sigma", 1), ("sigma", 2)])
    br = np.array_equal(a, b)
    # minimal polynomial of sigma on P_X(2)
    S = dense(X, 2, [("sigma", 1)])
    Z = np.eye(S.shape[0], dtype=np.int64)
    for lam in eigs.values():
        Z = F.mm(Z, (S - lam * np.eye(S.shape[0], dtype=np.int64)) % p)
    mp = not Z.any()
    # zigzags, loop, curl
    z1 = scalar(dense(X, 1, [("cup", 1), ("cap", 0)]))
    z2 = scalar(dense(X, 1, [("cup", 0), ("cap", 1)]))
    loop = scalar(dense(X, 0, [("cup", 0), ("cap", 0)]))
    curl = scalar(dense(X, 1, [("cup", 1), ("sigma", 1), ("cap", 1)]))
    curl2 = scalar(dense(X, 1, [("cup", 0), ("sigma", 2), ("cap", 0)]))
    # naturality of merge and cup with sigma
    trivalent = gen in eigs           # V is a summand of V (x) V
    if trivalent:
        m1 = dense(X, 2, [("sigma", 1), ("merge", 1)])
        m0 = dense(X, 2, [("merge", 1)])
        nat_m = np.array_equal(m1, (eigs[gen] * m0) % p)
    else:
        nat_m = True
    c1 = dense(X, 0, [("cup", 0), ("sigma", 1)])
    c0 = dense(X, 0, [("cup", 0)])
    nat_c = np.array_equal(c1, (eigs[zero] * c0) % p)
    # split followed by merge is a scalar (the theta graph)
    tg = scalar(dense(X, 1, [("split", 1), ("merge", 1)])) if trivalent else 1
    # sliding a cap through a split (pivotality of the trivalent vertex)
    results[X] = dict(braid=br, minpoly=mp, z1=z1, z2=z2, loop=loop, curl=curl, curl2=curl2, nat_m=nat_m,
                      nat_c=nat_c, theta=tg)
    ok = br and mp and nat_m and nat_c and None not in (z1, z2, loop, curl, curl2, tg)
    fails += not ok
    print(X, "paths(3) =", len(M.P(X, 3)), results[X], "OK" if ok else "FAIL")
vals = {k: {r[k] for r in results.values()} for k in ("z1", "z2", "loop", "curl", "curl2", "theta")}
print("trivalent vertex V (x) V -> V present:", gen in eigs)
print("scalars independent of the base:", {k: len(v) == 1 for k, v in vals.items()})
curl = vals["curl"].pop()
for e in range(order):
    if qp(e) == curl:
        print("curl = q^%d" % e)
    if (p - qp(e)) % p == curl:
        print("curl = -q^%d" % e)
print("failures:", fails)
