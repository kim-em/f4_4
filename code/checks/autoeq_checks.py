#!/usr/bin/env python3
"""Exact checks for the paper's Section "Autoequivalences and braidings", from the R-file only.

  1. fusion ring automorphisms (backtracking over permutations preserving N), their orders;
  2. no product X (x) Y of non-unit simples is simple;
  3. the Mueger centre is trivial: for each X != 1 some monodromy scalar
     theta_Z / (theta_X theta_Y) (Z in X (x) Y) differs from 1; the monodromy matrices
     R[Y,X,Z] R[X,Y,Z] are these scalars times the identity (exact_checks.py);
  4. for each automorphism pi and each a in (Z/52)^*, whether the gauge-invariant numbers
     R[X,X,Z] on one-dimensional spaces satisfy  R[pi X, pi X, pi Z] = sigma_a(R[X,X,Z])
     ("braided": pi transports c to c^{sigma_a}) or = sigma_a(R[X,X,Z])^{-1} ("reversed").
     Both values are roots of unity q^e (checked), so this is a comparison of exponents mod 52.

usage: autoeq_checks.py R.txt[.xz]
"""
import sys, os, math, itertools
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from exact_checks import read_R, ZP, NB, COND

obj, unit, N, R, Ri = read_R(sys.argv[1])
n = lambda a, b, c: N.get((a, b, c), 0)
UNITS = [a for a in range(1, COND) if math.gcd(a, COND) == 1]

# 1. automorphisms: permutations pi of obj with pi(unit) = unit and N(pi a, pi b, pi c) = N(a,b,c)
row = {a: tuple(sorted(n(a, b, c) for b in obj for c in obj)) for a in obj}
auts = []
def extend(m):
    if len(m) == len(obj):
        auts.append(dict(m)); return
    a = obj[len(m)]
    for b in obj:
        if b in m.values() or row[b] != row[a]:
            continue
        m[a] = b
        if all(n(x, y, z) == n(m[x], m[y], m[z]) for x in m for y in m for z in m):
            extend(m)
        del m[a]
extend({})
def order(p):
    k, q = 1, dict(p)
    while any(q[x] != x for x in obj):
        q = {x: p[q[x]] for x in obj}; k += 1
    return k
def cycles(p):
    seen, out = set(), []
    for x in obj:
        if x in seen or p[x] == x:
            continue
        c, y = [], x
        while y not in seen:
            seen.add(y); c.append(y); y = p[y]
        out.append("(" + " ".join(c) + ")")
    return " ".join(out) or "id"
print(f"1. fusion ring automorphisms: {len(auts)}, orders {sorted(order(p) for p in auts)}")
for p in auts:
    print("   ", cycles(p))
# 2. simple products
simp = [(x, y) for x in obj for y in obj if x != unit and y != unit and sum(n(x, y, z) for z in obj) == 1]
print("2. pairs of non-unit simples with simple tensor product:", simp)

# exponents of one-dimensional R[X,X,Z] = zeta^e
def exponent(v):
    for e in range(COND):
        if [Fr(x) for x in ZP[e]] == [Fr(x) for x in v]:
            return e
    raise ValueError("not a root of unity")
E = {(X, Z): exponent(R[(X, X, Z)][(0, 0)]) for X in obj for Z in obj if n(X, X, Z) == 1}
print(f"   one-dimensional R[X,X,Z]: {len(E)}, all of the form zeta^e")

# 3. Mueger centre via twists theta = zeta^{h52}, h52 = (X, X + 2 rho) mod 52, read from R: theta_X^{-1} = R[X,X,1] for nu_X = 1
th = {X: (-E[(X, unit)]) % COND for X in obj}          # R[X,X,1] = theta_X^{-1} (self-dual, indicator +1)
transparent = [X for X in obj if X != unit and all((th[Z] - th[X] - th[Y]) % COND == 0
               for Y in obj for Z in obj if n(X, Y, Z))]
print("3. objects X != 1 with all monodromy scalars theta_Z/(theta_X theta_Y) = 1 (transparent):", transparent)

# 4. transport of the braiding
print("4. pi, and the a in (Z/52)^* with R[pi X,pi X,pi Z] = sigma_a(R[X,X,Z]) (braided) or its inverse (reversed):")
for p in auts:
    br = [a for a in UNITS if all(E[(p[X], p[Z])] == (a * e) % COND for (X, Z), e in E.items())]
    rv = [a for a in UNITS if all(E[(p[X], p[Z])] == (-a * e) % COND for (X, Z), e in E.items())]
    print(f"   {cycles(p)}: braided for a in {br}; reversed for a in {rv}")
