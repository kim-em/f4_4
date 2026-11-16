"""R-symbols only (splitting convention Rs, see rmatrix.py), modulo p.

usage: run_r.py cells.pkl out.pkl
"""
import sys, pickle, time
from rootsys import RootSystem
from objects import Objects
from rmatrix import RSymbols

c = pickle.load(open(sys.argv[1], "rb"))
R = RootSystem(c["type"])
alc = c["alcove"]
h = {w: int(R.casimir(w)) for w in alc}
fus = {(X, S, T): m for X in alc for S in alc for T, m in R.tensor(X, S, c["k"]).items()}
O = Objects(c, h, fus)
RS = RSymbols(O)
t = time.time()
G = {}
for X in alc:
    for Y in alc:
        for Z in alc:
            if fus.get((X, Y, Z), 0):
                M = RS.R(X, Y, Z)
                G[(X, Y, Z)] = ([(Z, i) for i in range(M.shape[0])], [(Z, j) for j in range(M.shape[1])], M)
pickle.dump(dict(G=G, p=c["p"], q=c["q"], which=c["which"]), open(sys.argv[2], "wb"))
print(f"[run_r] {len(G)} R-matrices in {time.time()-t:.0f}s")
