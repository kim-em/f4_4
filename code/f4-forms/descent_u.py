#!/usr/bin/env python3
"""The descent data of C(F4,4) along Q(zeta_13)/E and the scalars u_X (exact).

Input: the exact gauges h of relabel_exact_nf.py (pickles), proved there to satisfy
      h . (sigma_5 . F) = pi^j . F        on every block with third leg V = 0001,
hence on the whole table (paper, Lemma "H determines").  Here F is the normal-form table of the
paper (entries in k = Q(zeta_13)^+), (pi.F)[X,Y,W,T] = F[pi^-1 X, ...], and pi is Gannon's
pi_4 (0001 -> 1000 -> 2000 -> 0004).

Conventions (paper, Section "Forms over the cubic field").  For a field automorphism g, a fusion ring automorphism
p and a gauge c indexed by target spaces, Psi = (g, p, c) is the g-semilinear monoidal functor
with Psi(X) = pX on simple objects and
      Psi(phi^{XY->Z}_i) = sum_a c^{pX pY}_{pZ}[i, a] phi^{pX pY -> pZ}_a ;
it is compatible with the associator iff c . F = p . (g . F).  Composition:
      (g2, p2, c2) o (g1, p1, c1) = (g2 g1, p2 p1, c),  c_s = g2(c1_{p2^-1 s}) c2_s  (matrices).
From h . (sigma_5 F) = pi^j F we get c . F = pi^j . (sigma_5^-1 F) with c = sigma_5^-1(h); on k,
sigma_5^-1 = sigma_21 and sigma_5 agree (sigma_25 fixes k).  So Psi_j = (sigma_5, pi^j, c) is a
sigma_5-semilinear monoidal autoequivalence of C over Q(zeta_13) for j = 1, 3.  Psi_j^4 is
linear, fixes the objects, and has gauge c4 on every space; a monoidal isomorphism
theta: Psi_j^4 => id with theta_X = u_X exists iff c4^{XY}_Z = u_X u_Y / u_Z.

This script, exactly in Q(zeta_52) (fields/cyc.py):
  * converts h to c, checks that every entry of c lies in k;
  * computes c4 on every one-dimensional space, the unique u with u_1 = 1 from the spaces
    (A, V, C), and checks c4 = u_X u_Y / u_Z on every one-dimensional space;
  * checks the cocycle identity u_{pi^j X} = sigma_5(u_X) and that u_X lies in E for pi-fixed X;
  * modulo a prime, at the 12 real roots, checks that c4 is the scalar u_X u_Y / u_Z also on
    the spaces of dimension >= 2, and checks the braiding test R[pX,pX,pZ] = sigma_5(R[X,X,Z])
    (or its inverse) on the 100 one-dimensional spaces Hom(X X, Z);
  * writes u (eta basis and as a polynomial in the generator of E) and a PARI script.

usage: descent_u.py H_PI1.pkl H_PI3.pkl R_FILE.xz OUTDIR
"""
import lzma, os, pickle, sys, time
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "fields"))
sys.path.insert(1, os.path.join(HERE, "..", "f4"))
import cyc
from lin import find_primes, primitive_root_of_unity

UNIT, V = "0000", "0001"
PI = {"0001": "1000", "1000": "2000", "2000": "0004", "0004": "0001",
      "0010": "0100", "0100": "0020", "0020": "0003", "0003": "0010",
      "1001": "1002", "1002": "1010", "1010": "0012", "0012": "1001"}
OBJ = ["0000", "0001", "0002", "0010", "1000", "0003", "0011", "0100", "1001", "0004", "0012",
       "0020", "0101", "1002", "1010", "2000"]


def perm_power(j):
    p = {x: x for x in OBJ}
    for _ in range(j % 4):
        p = {x: PI.get(p[x], p[x]) for x in OBJ}
    return p


def load_h(path):
    d = pickle.load(open(path, "rb"))
    h = {}
    for t, (o, sh) in d["gpos"].items():
        ents = []
        for i in range(sh[0] * sh[1]):
            ents.append(cyc.from_real([Fraction(x, y) for (x, y) in d["hr"][o + i]]))
        h[t] = [ents[r * sh[1]:(r + 1) * sh[1]] for r in range(sh[0])]
    return d, h


def matmul(A, B):
    n, m, l = len(A), len(B), len(B[0])
    return [[ (lambda acc: acc)(sum_list([cyc.mul(A[i][k], B[k][j]) for k in range(m)])) for j in range(l)] for i in range(n)]


def sum_list(xs):
    out = cyc.ZERO
    for x in xs:
        out = cyc.add(out, x)
    return out


def sig(M, a):
    return [[cyc.sigma(x, a) for x in r] for r in M]


def eta_str(a):
    e = cyc.to_eta(a)
    return None if e is None else [str(x) for x in e]


def main(argv):
    out_dir = argv[4]
    os.makedirs(out_dir, exist_ok=True)
    t0 = time.time()
    results = {}
    gp_lines = []
    for j, path in ((1, argv[1]), (3, argv[2])):
        d, h = load_h(path)
        p = perm_power(j)
        assert d["pi"] == p, "the pickle is not the expected relabelling"
        assert d["a"] == 5
        pinv = {y: x for x, y in p.items()}
        # c = sigma_5^-1(h) = sigma_21(h)
        c = {t: sig(M, 21) for t, M in h.items()}
        notk = sum(1 for M in c.values() for r in M for x in r if cyc.to_eta(x) is None)
        print(f"[pi^{j}] {len(c)} spaces; entries of c outside k: {notk}  ({time.time() - t0:.0f} s)", flush=True)
        assert notk == 0
        one = [t for t, M in c.items() if len(M) == 1 and len(M[0]) == 1]
        # c4 on one-dimensional spaces: c_{n+1, s} = sigma_5(c_{n, p^-1 s}) c_s
        cn = {t: c[t][0][0] for t in one}
        for n in range(3):
            cn = {s: cyc.mul(cyc.sigma(cn[tuple(pinv[x] for x in s)], 5), c[s][0][0]) for s in one}
        c4 = cn
        print(f"[pi^{j}] c4 on {len(one)} one-dimensional spaces ({time.time() - t0:.0f} s)", flush=True)
        u = {UNIT: cyc.ONE, V: c4[(V, V, V)]}
        changed = True
        while changed:
            changed = False
            for (A, B, C) in one:
                if B == V and A in u and C not in u:
                    u[C] = cyc.mul(cyc.mul(u[A], u[V]), cyc.inv(c4[(A, B, C)]))
                    changed = True
        assert len(u) == len(OBJ)
        bad = 0
        uinv = {X: cyc.inv(x) for X, x in u.items()}
        for (A, B, C) in one:
            if cyc.sub(c4[(A, B, C)], cyc.mul(cyc.mul(u[A], u[B]), uinv[C])) != cyc.ZERO:
                bad += 1
        print(f"[pi^{j}] c4 = u_X u_Y / u_Z on one-dimensional spaces: {len(one) - bad} hold, {bad} fail", flush=True)
        assert bad == 0
        coc = all(cyc.sigma(u[X], 5) == u[p[X]] for X in OBJ)
        print(f"[pi^{j}] cocycle identity u_(pX) = sigma_5(u_X) for all X: {coc}")
        fixed = [X for X in OBJ if p[X] == X]
        for X in OBJ:
            e = cyc.to_eta(u[X])
            inE = cyc.sigma(u[X], 5) == u[X]
            print(f"[pi^{j}] u_{X} = {' + '.join(f'({x})*eta{i + 1}' for i, x in enumerate(e) if x)}   in E: {inE}   ~ {cyc.to_float(u[X]).real:.6g}")
        results[j] = {"u": {X: cyc.to_eta(u[X]) for X in OBJ}, "c": c, "fixed": fixed}
    # modular checks: full scalar property of c4 (all spaces) and the braiding test
    q = find_primes(52, 60)[-1]
    z = primitive_root_of_unity(q, 52, 1)

    def ev(a, w):
        zz = pow(z, w, q)
        acc = 0
        for k, x in enumerate(a):
            if x:
                acc = (acc + x.numerator % q * pow(x.denominator, -1, q) % q * pow(zz, k, q)) % q
        return acc

    def mmq(A, B):
        return [[sum(A[i][k] * B[k][jj] for k in range(len(B))) % q for jj in range(len(B[0]))] for i in range(len(A))]

    REALW = [1, 3, 5, 7, 9, 11, 15, 17, 19, 21, 23, 25]
    cls = lambda w: (w % 52) if (w % 52) < 26 else 52 - (w % 52)
    for j in (1, 3):
        p = perm_power(j); pinv = {y: x for x, y in p.items()}
        c = results[j]["c"]
        cw = {w: {t: [[ev(x, w) for x in r] for r in M] for t, M in c.items()} for w in REALW}
        cur = cw
        for n in range(3):
            cur = {w: {s: mmq(cur[cls(5 * w)][tuple(pinv[x] for x in s)], cw[w][s]) for s in c} for w in REALW}
        uw = {w: {X: ev(cyc.from_eta(results[j]["u"][X]), w) for X in OBJ} for w in REALW}
        nbad = 0
        for w in REALW:
            for (X, Y, Z), M in cur[w].items():
                sc = uw[w][X] * uw[w][Y] % q * pow(uw[w][Z], -1, q) % q
                n_ = len(M)
                if any((M[a][b] - (sc if a == b else 0)) % q for a in range(n_) for b in range(n_)):
                    nbad += 1
        print(f"[pi^{j}] mod {q} at 12 real roots: c4 = (u_X u_Y / u_Z) Id on all {len(c)} spaces: {nbad == 0} ({nbad} failures)", flush=True)
    # braiding test with the R-symbols of the paper (entries in Q(zeta_13), power basis of zeta_52)
    R = {}
    with lzma.open(argv[3], "rt") as fh:
        for line in fh:
            if line.startswith("R "):
                head, idx, coef = line.split("|")
                _, X, Y, Z = head.split()
                if X == Y and idx.split() == ["0", "0"]:
                    R[(X, Z)] = [Fraction(s) for s in coef.split()]
    N = {}
    with lzma.open(argv[3], "rt") as fh:
        for line in fh:
            if line.startswith("# N "):
                _, _, X, Y, Z, n = line.split()
                N[(X, Y, Z)] = int(n)
    keys = [(X, Z) for (X, Z) in R if N.get((X, X, Z)) == 1]
    for j in (1, 3):
        p = perm_power(j)
        same = all(ev(R[(p[X], p[Z])], w) == ev(R[(X, Z)], 5 * w) for (X, Z) in keys for w in range(1, 52, 2) if w != 13 and w != 39)
        rev = all(ev(R[(p[X], p[Z])], w) * ev(R[(X, Z)], 5 * w) % q == 1 for (X, Z) in keys for w in range(1, 52, 2) if w != 13 and w != 39)
        print(f"[pi^{j}] {len(keys)} one-dimensional R[X,X,Z], mod {q}, all 24 roots: R[pX,pX,pZ] = sigma_5(R[X,X,Z]): {same};  = sigma_5(R[X,X,Z])^-1: {rev}")
        results[j]["braided"] = same
    with open(os.path.join(out_dir, "u_values.txt"), "w") as fh:
        fh.write("# u_X for Psi_j = (sigma_5, pi^j, c_j): Psi_j^4 => id is theta_X = u_X.  eta basis of k, eta_i = zeta13^i + zeta13^-i\n")
        for j in (1, 3):
            for X in OBJ:
                fh.write(f"j={j} {X} " + " ".join(str(x) for x in results[j]["u"][X]) + "\n")
    pickle.dump({j: {"u": results[j]["u"], "fixed": results[j]["fixed"], "braided": results[j]["braided"]} for j in (1, 3)},
                open(os.path.join(out_dir, "u_values.pkl"), "wb"))
    print(f"done ({time.time() - t0:.0f} s)")


if __name__ == "__main__":
    main(sys.argv)
