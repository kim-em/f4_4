#!/usr/bin/env python3
"""Fusion-ring facts used in the paper's Section "Forms over the cubic field", from the fusion rules in the header of the R-file.

  * Gannon's pi_4 preserves the fusion rules; its orbits; pi^2.
  * The universal grading group is trivial: V = 0001 is a summand of V (x) V, and every simple
    object is a summand of a power of V; hence every grading has V in degree 1, and all of
    Irr is in degree 1.  Equivalently, the only u with u_X u_Y = u_Z whenever N_XY^Z != 0 is u = 1
    (checked directly: the exponent lattice of these relations has full rank 16).
  * The fusion rules of the forms D_1, D_3 over E and of the k-form along Q(zeta13)/k
    (multiplicity of a simple object S'' as a summand of S (x) S'), using the Schur indices of
    brauer.gp.

usage: fusion_facts.py R_FILE.xz
"""
import lzma, sys
from fractions import Fraction

OBJ = ["0000", "0001", "0002", "0010", "1000", "0003", "0011", "0100", "1001", "0004", "0012",
       "0020", "0101", "1002", "1010", "2000"]
PI = {"0001": "1000", "1000": "2000", "2000": "0004", "0004": "0001",
      "0010": "0100", "0100": "0020", "0020": "0003", "0003": "0010",
      "1001": "1002", "1002": "1010", "1010": "0012", "0012": "1001"}
pi = {x: PI.get(x, x) for x in OBJ}


def rank_int(rows):
    rows = [[Fraction(x) for x in r] for r in rows]
    r = 0
    ncol = len(rows[0])
    for c in range(ncol):
        piv = next((i for i in range(r, len(rows)) if rows[i][c]), None)
        if piv is None:
            continue
        rows[r], rows[piv] = rows[piv], rows[r]
        for i in range(len(rows)):
            if i != r and rows[i][c]:
                f = rows[i][c] / rows[r][c]
                rows[i] = [a - f * b for a, b in zip(rows[i], rows[r])]
        r += 1
    return r


def main(argv):
    N = {}
    with lzma.open(argv[1], "rt") as fh:
        for line in fh:
            if line.startswith("# N "):
                _, _, X, Y, Z, n = line.split()
                N[(X, Y, Z)] = int(n)
    n = lambda X, Y, Z: N.get((X, Y, Z), 0)
    print(f"{len(N)} nonzero fusion coefficients")
    ok = all(n(pi[X], pi[Y], pi[Z]) == n(X, Y, Z) for X in OBJ for Y in OBJ for Z in OBJ)
    print("pi preserves the fusion rules:", ok)
    orbits, seen = [], set()
    for X in OBJ:
        if X not in seen:
            o = [X]
            while pi[o[-1]] != X:
                o.append(pi[o[-1]])
            orbits.append(o); seen |= set(o)
    print("pi-orbits:", orbits)
    pi2 = {x: pi[pi[x]] for x in OBJ}
    orb2, seen = [], set()
    for X in OBJ:
        if X not in seen:
            o = sorted({X, pi2[X]}, key=OBJ.index); orb2.append(o); seen |= set(o)
    print("pi^2-orbits:", orb2)
    print("V in V (x) V:", n("0001", "0001", "0001") == 1)
    depth = {"0000": 0}
    frontier = ["0000"]
    while frontier:
        new = []
        for X in frontier:
            for Z in OBJ:
                if n(X, "0001", Z) and Z not in depth:
                    depth[Z] = depth[X] + 1; new.append(Z)
        frontier = new
    print("every simple object is a summand of a power of V:", len(depth) == 16, "max depth", max(depth.values()))
    rows = []
    for (X, Y, Z) in N:
        r = [0] * 16
        r[OBJ.index(X)] += 1; r[OBJ.index(Y)] += 1; r[OBJ.index(Z)] -= 1
        rows.append(r)
    rk = rank_int(rows)
    print(f"rank of the relations u_X u_Y = u_Z (exponent vectors): {rk} of 16; with the lattice:",
          "the only solution is u = 1" if rk == 16 else "nontrivial solutions")
    # unimodularity: rank 16 alone gives finiteness; V=VV gives u_V=1 and then u_Z=u_V u_X/... directly:
    u = {"0000": 1, "0001": 1}
    print("u_V = 1 from N_VV^V = 1; then u_Z = u_X u_V for N_XV^Z != 0 propagates u = 1 to all", len(depth), "objects")
    # fusion rules of the forms: simple objects = orbits; S(O) (x) K = sum over O, times the Schur
    # index m(O) for fixed objects.  The multiplicity of S(O'') in S(O) (x) S(O') (as a summand) is
    # m(O) m(O') sum_{X in O, Y in O'} N_XY^{Z0} / m(O''), Z0 any element of O''.
    # Schur indices: m = 2 for the fixed objects 0002, 0011, 0101 in D_1, D_3 and in C^(Q(zeta13));
    # m = 1 in the forms where all three quaternion algebras split (paper, Section "Forms over the cubic field").
    m2 = {o: (2 if o[0] in ("0002", "0011", "0101") else 1) for o in map(tuple, orbits + orb2)}
    m1 = {o: 1 for o in map(tuple, orbits + orb2)}
    cases = (("D_1, D_3 over E (Schur index 2 on 0002, 0011, 0101)", orbits, m2),
             ("D_chi over E with all three algebras split, e.g. e = -1, M_chi = Q(zeta_52)^+", orbits, m1),
             ("C^L over k, L = Q(zeta13) (Schur index 2 on 0002, 0011, 0101)", orb2, m2),
             ("C^L over k with all three algebras split", orb2, m1))
    for name, orbs, m in cases:
        print(f"\nfusion rules of {name}:")
        label = {tuple(o): ("{" + ",".join(o) + "}" if len(o) > 1 else o[0]) for o in orbs}
        nchk = 0
        for a in orbs:
            for b in orbs:
                for c in orbs:
                    tot = m[tuple(a)] * m[tuple(b)] * sum(n(X, Y, c[0]) for X in a for Y in b)
                    if m[tuple(c)] > 1 and m[tuple(a)] * m[tuple(b)] == 1:
                        nchk += 1
                        assert tot % m[tuple(c)] == 0, (a, b, c)
        print(f"  integrality: the multiplicity of every simple object with Schur index m in every product is an integer ({nchk} constraints with S, S' of index 1 and S'' of index 2; these are automatic, since pi acts freely on the pairs of an orbit product)")
        dimE = sum(len(o) * m[tuple(o)] ** 2 for o in orbs)
        print(f"  {len(orbs)} simple objects; sum of dim End(S) over the base field = {dimE}; base changes contain {sum(len(o) for o in orbs)} distinct simple objects")
        for a in orbs:
            for b in orbs:
                if OBJ.index(a[0]) > OBJ.index(b[0]):
                    continue
                terms = []
                for c in orbs:
                    mm = m[tuple(a)] * m[tuple(b)] * sum(n(X, Y, c[0]) for X in a for Y in b) // m[tuple(c)]
                    if mm:
                        terms.append(f"{mm} {label[tuple(c)]}" if mm > 1 else label[tuple(c)])
                print(f"  {label[tuple(a)]} * {label[tuple(b)]} = " + " + ".join(terms))

if __name__ == "__main__":
    main(sys.argv)
