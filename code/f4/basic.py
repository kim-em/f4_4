"""Basic exact data for C(g, k): alcove weights, dimensions, twists, fusion rules.

usage: basic.py [TYPE LEVEL GENERATOR [OUT.pkl]]      default: F4 4 0,0,0,1
Writes OUT.pkl (default: work/basic_<TYPE>_<LEVEL>.pkl in the repository, or in $F4_WORK if set)
and prints a summary.
"""
import sys, pickle, math, cmath, itertools, os
from fractions import Fraction as Fr
from rootsys import RootSystem

typ = sys.argv[1] if len(sys.argv) > 1 else "F4"
k = int(sys.argv[2]) if len(sys.argv) > 2 else 4
gen = tuple(int(x) for x in sys.argv[3].split(",")) if len(sys.argv) > 3 else (0, 0, 0, 1)

R = RootSystem(typ)
kap = k + R.hdual
order = 2 * R.m * kap            # q = exp(2 pi i / order) in the short-root-2 normalisation
W = R.alcove(k)
n = len(W)
print(f"{typ} level {k}: h^vee={R.hdual}, lacing m={R.m}, comarks={R.comarks}, k+h^vee={kap}, q of order {order}")
print(f"rank of C = {n}")


def qdim(lam):
    val = 1.0
    for a in R.positive_roots:
        x, y = R.ip(tuple(l + 1 for l in lam), a), R.ip(R.rho, a)
        val *= math.sin(math.pi * x / (R.m * kap)) / math.sin(math.pi * y / (R.m * kap))
    return val


dims = [qdim(w) for w in W]
hs = [R.h(w, k) for w in W]
D2 = sum(d * d for d in dims)
c = Fr(k * len(R.positive_roots) * 2 + k * R.rank, 1) / kap   # k dim g / (k + h^vee)
print(f"central charge c = {c},  D^2 = {D2:.10f}, D = {math.sqrt(D2):.10f}")
fus = {}
for a in W:
    for b in W:
        for cc, m in R.tensor(a, b, k).items():
            fus[(a, b, cc)] = m
N = lambda a, b, cc: fus.get((a, b, cc), 0)
# duality
dual = {}
for a in W:
    ds = [b for b in W if N(a, b, W[0]) == 1]
    assert len(ds) == 1
    dual[a] = ds[0]
print("all objects self-dual:", all(dual[a] == a for a in W))
# associativity of the fusion ring
import numpy as np
Nm = {a: np.array([[N(a, b, cc) for cc in W] for b in W]) for a in W}
for a, b in itertools.product(W, W):
    lhs = Nm[a] @ Nm[b]
    rhs = sum(N(a, b, cc) * Nm[cc] for cc in W)
    assert (lhs == rhs).all()
# dimensions are a character
for a, b in itertools.product(W, W):
    assert abs(dims[W.index(a)] * dims[W.index(b)] - sum(N(a, b, cc) * dims[W.index(cc)] for cc in W)) < 1e-8
print("fusion ring associative, commutative, dimensions multiplicative: ok")
maxm = max(fus.values())
print("largest fusion multiplicity:", maxm, [key for key, v in fus.items() if v == maxm][:6])
print(" #  weight        dim(Weyl)  qdim            h        level")
for i, w in enumerate(W):
    print(f"{i:2d} {str(w):14s} {R.weyl_dim(w):9d}  {dims[i]:14.10f}  {str(hs[i]):8s} {R.level(w)}")

# generator graph
g = gen
print(f"generator {g}: dim {R.weyl_dim(g)}, qdim {qdim(g):.10f}, h = {R.h(g, k)}")
depth = {W[0]: 0}
layer = [W[0]]
while layer:
    nxt = []
    for a in layer:
        for cc in W:
            if N(a, g, cc) and cc not in depth:
                depth[cc] = depth[a] + 1
                nxt.append(cc)
    layer = nxt
print("depth (first tensor power of the generator containing the object):",
      {str(w): depth.get(w) for w in W})
loops = [(a, cc, N(a, g, cc)) for a in W for cc in W if N(a, g, cc) > 1]
print("generator edges with multiplicity > 1:", loops)

# sizes of the F-symbol tables
def msize(a, b, cc, d):
    return sum(N(a, b, e) * N(e, cc, d) for e in W)
tot_mats = tot_entries = biggest = 0
for a, b, cc, d in itertools.product(W, repeat=4):
    m = msize(a, b, cc, d)
    if m:
        tot_mats += 1
        tot_entries += m * m
        biggest = max(biggest, m)
print(f"full F table: {tot_mats} matrices, {tot_entries} entries, largest {biggest}")
_work = os.environ.get("F4_WORK", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "work"))
_out = sys.argv[4] if len(sys.argv) > 4 else os.path.join(_work, f"basic_{typ}_{k}.pkl")
os.makedirs(os.path.dirname(os.path.abspath(_out)), exist_ok=True)
pickle.dump(dict(type=typ, k=k, weights=W, dims=dims, h=hs, fus=fus, gen=g, order=order, depth=depth,
                 c=c), open(_out, "wb"))

# twists through the fusion ring: Verlinde check against the Kac-Peterson S matrix
theta = [cmath.exp(2j * math.pi * float(h)) for h in hs]
S = np.array([[sum(N(dual[a], b, cc) * theta[W.index(cc)] / (theta[W.index(a)] * theta[W.index(b)]) * dims[W.index(cc)]
                    for cc in W) for b in W] for a in W]) / math.sqrt(D2)
U = S @ S.conj().T
print("S from (fusion, twists, dims) is unitary:", np.allclose(U, np.eye(n), atol=1e-9))
T = np.diag(theta)
ST = S @ T
lhs = np.linalg.matrix_power(ST, 3)
ratio = lhs[0, 0] / (S @ S)[0, 0]
print("(ST)^3 = p S^2 with p =", ratio, " expected exp(2 pi i c/8) =", cmath.exp(2j * math.pi * float(c) / 8))
print("(ST)^3 proportional to S^2:", np.allclose(lhs, ratio * (S @ S), atol=1e-9))
