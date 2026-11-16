#!/usr/bin/env python3
"""Third check: SageMath's F-matrix solver for C(B6,2) against the quantum group data.

usage: compare_sage.py SAGE_OUTPUT              Sage B6 level 2 against data/B6_2_FR_exact.txt
       compare_sage.py SAGE_OUTPUT --aft P      Sage B<P> level 2 against the AFT formulas with
                                                so2p1_2.quantum_group_parameters(P)
(SAGE_OUTPUT written by sage_fmatrix.sage B<P> 2 SAGE_OUTPUT)

Sage (FusionRing("B<P>", 2).get_fmatrix().find_orthogonal_solution(), SageMath 10.9) solves the
pentagon from the fusion rules alone; its F-symbols are exact elements of Q(zeta_N),
N = 8(2P+1).  Sage's generator zetaN is embedded in C as exp(2 pi i k/N) for some k, found from
the numerical values Sage also writes (k = N/2 - 1 in all runs here; -k fits equally, the
F-symbols being real), so its F-symbols are mapped by zetaN -> zeta^k into our Q(zeta),
zeta = exp(2 pi i/N).  Its R-symbols come from a closed formula, are written numerically to 60
digits, and each is recognised here exactly as an element of mu_N.  Sage F_{abc}^d(x, y),
x in a(x)b, y in b(x)c, is compared with F_{abc}^{d;xy} of verify.py.  Missing Sage F-symbols
are those with a unit label, which are 1.

The check: solve for a gauge g in Q(zeta_N) with (ours or AFT)^g = Sage (F and R), as in
compare_so13.py, apply it to every entry and compare exactly.
"""
import cmath
import lzma
import math
import re
import sys

from cyclo import Field
from ours import Ours
from compare_so13 import prepare, solve_gauge, apply_and_compare


def dynkin(s):
    m = re.fullmatch(r"B\d2\(([\d,]+)\)", s)
    return "".join(m.group(1).split(","))


def parse_cyclo(K, s):
    """'a/b*zetaN^k + ...' (Sage repr) -> element of K = Q(zeta_N)."""
    from fractions import Fraction
    s = s.replace(" ", "")
    if not s.startswith("-"):
        s = "+" + s
    v = K.zero
    for sign, body in re.findall(r"([+-])([^+-]+)", s):
        if "zeta" in body:
            if "*" in body:
                c, z = body.split("*")
            else:
                c, z = "1", body
            k = int(z.split("^")[1]) if "^" in z else 1
        else:
            c, k = body, 0
        term = K.zeta(k) * Fraction(c)
        v = v + (term if sign == "+" else -term)
    return v


class SageData:
    def __init__(self, path, K):
        self.K = K
        self.F, self.R = {}, {}
        N = K.N
        raw = []
        with (lzma.open(path, "rt") if path.endswith(".xz") else open(path)) as fh:
            for line in fh:
                if line.startswith("F"):
                    head, num, ex = line.split("|")
                    raw.append((parse_cyclo(K, ex.strip()), complex(*map(float, num.split()))))
        # Sage's generator zetaN is embedded in C (get_qqbar_embedding) as exp(2 pi i k/N) for some k;
        # find k from the numerical values (k and -k both fit: the F-symbols are real)
        ks = [k for k in range(1, N) if math.gcd(k, N) == 1 and all(abs(complex(v.galois(k)) - z) < 1e-9 for v, z in raw)]
        assert ks, "no embedding fits"
        self.embedding = k = min(ks)
        with (lzma.open(path, "rt") if path.endswith(".xz") else open(path)) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                head, *rest = line.split("|")
                t = head.split()
                if t[0] == "F":
                    a, b, c, d, x, y = (dynkin(u) for u in t[1:7])
                    val = parse_cyclo(K, rest[1].strip()).galois(k)
                    re_, im_ = (float(u) for u in rest[0].split())
                    assert abs(complex(val) - complex(re_, im_)) < 1e-9
                    self.F.setdefault((a, b, c, d), {})[(x, y)] = val
                else:
                    a, b, c = (dynkin(u) for u in t[1:4])
                    re_, im_ = (float(u) for u in rest[0].split())
                    z = complex(re_, im_)
                    j = round(cmath.phase(z) / (2 * math.pi) * N) % N
                    assert abs(z - cmath.exp(2j * math.pi * j / N)) < 1e-12, (a, b, c, z)
                    self.R[(a, b, c)] = K.zeta(j)


class AFTDynkin:
    """the AFT data for C(so(2p+1),2) relabelled by Dynkin labels (interface of ours.Ours)."""
    def __init__(self, p):
        from so2p1_2 import Metaplectic, quantum_group_parameters
        r, kappa, lam, lab = quantum_group_parameters(p)
        D = Metaplectic(p, r, kappa, lam)
        inv = {v: k for k, v in lab.items()}
        self.K, self.unit = D.K, inv["1"]
        self.labels = [inv[x] for x in D.labels]
        self.N = {tuple(inv[x] for x in t): 1 for t in D.N}
        self.F = {tuple(inv[x] for x in k): {(inv[e], inv[f]): v for (e, f), v in b.items()} for k, b in D.F.items()}
        self.R = {tuple(inv[x] for x in k): v for k, v in D.R.items()}
        self.params = (p, r, kappa, lam)


def main(argv):
    if "--aft" in argv:
        O = AFTDynkin(int(argv[argv.index("--aft") + 1]))
        K = O.K
        print(f"AFT (p, r, kappa, lam) = {O.params} in Dynkin labels, against Sage")
    else:
        K = Field(104)
        O = Ours(field=K)
    prepare(O)
    S = SageData(argv[1], K)
    print(f"Sage generator embedded as exp(2 pi i {S.embedding}/{K.N})")
    one = K.one
    unit = O.unit
    for key, (rows, cols) in O.blocks.items():
        if unit in key:
            S.F.setdefault(key, {})
            for Z in rows:
                for Zp in cols:
                    S.F[key].setdefault((Z, Zp), one)
    extra = [k for k in S.F if k not in O.blocks]
    assert not extra, extra[:5]
    nS = sum(len(b) for b in S.F.values())
    print(f"Sage: {len(S.F)} blocks, {nS} F-symbols (unit-leg ones filled in as 1), {len(S.R)} R-symbols")
    lab = {x: x for x in O.labels}
    for name, R in (("Sage R", dict(S.R)), ("Sage R inverted", {k: v.inverse() for k, v in S.R.items()})):
        S.R = R
        g, unresolved = solve_gauge(O, S, lab)
        bad, nF, nR = apply_and_compare(O, S, lab, g)
        print(f"{name}: gauge with +-1 pivots, {len(unresolved)} equations left over; after the gauge {nF} F and {nR} R entries compared exactly, mismatches {len(bad)}")
        if bad:
            nf = sum(1 for b in bad if b[0] == "F")
            print(f"   ({nf} F mismatches, {len(bad) - nf} R mismatches; first: {bad[:4]})")


if __name__ == "__main__":
    main(sys.argv)
