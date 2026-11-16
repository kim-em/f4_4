#!/usr/bin/env python3
"""Exact verification of the AFT formulas (so2p1_2.py) for given p, all parameters.

usage: check_formulas.py P [P ...] [--jobs J] [--r R] [--kappa K]

For every r (odd, 1 <= r < 2p+1, coprime to 2p+1) and kappa = +-1, exactly in Q(zeta_{8(2p+1)}):
  * the rotation property used to complete the tables is consistent, every F-matrix is orthogonal,
  * the pentagon for all instances,
  * the hexagons (H+) and (H-) of verify.py for lam = +1 and lam = -1, and AFT's (R2) as printed,
  * the twists theta_a = d_a^-1 sum_c d_c R_{aa}^c (d_psi = +sqrt(2p+1)) against AFT (tw1),
  * S_{ab} = theta_a^-1 theta_b^-1 sum_c N_{ab}^c d_c theta_c against AFT (sm1)-(sm4), against
    (sm1)-(sm4) with the sign of the psi-psi block flipped, and against D Shat D from AFT (smat)
    with its second F-symbol replaced by the inverse F-matrix entry.
"""
import math
import sys
import time

from so2p1_2 import Metaplectic, jacobi
from verify import hexagon_failures, orthogonality_failures, pentagon_parallel


def twists(D):
    K = D.K
    th = {}
    for a in D.labels:
        s = K.zero
        for c in D.labels:
            if (a, a, c) in D.N:
                s = s + D.dim(c) * D.R[(a, a, c)]
        th[a] = s / D.dim(a)
    return th


def printed_twists(D):
    """AFT (tw1), with (-1)^x read as exp(pi i x)."""
    K, p, n, r, kappa, lam = D.K, D.p, D.n, D.r, D.kappa, D.lam
    e = lambda x: K.root_of_unity(int(x * K.N / 2), K.N)        # exp(pi i x), x rational
    from fractions import Fraction as Fr
    th = {"1": K.one, "eps": K.one}
    for i in range(1, p + 1):
        th[f"phi{i}"] = e(Fr(i)) * e(Fr(lam * r * i * i, n))
    for s, lab in ((1, "psi+"), (-1, "psi-")):
        th[lab] = e(Fr(kappa * lam * r, 2)) * e(Fr(lam * r, 4) * (jacobi(n, r) + kappa * jacobi(2, n) - p)) * (-s)
    return th


def smatrix(D, th):
    K = D.K
    S = {}
    for a in D.labels:
        for b in D.labels:
            s = K.zero
            for c in D.labels:
                if (a, b, c) in D.N:
                    s = s + D.dim(c) * th[c]
            S[(a, b)] = s / (th[a] * th[b])
    return S


def printed_smatrix(D, corrected=False):
    """AFT (sm1)-(sm4); corrected=True flips the sign of the four psi-psi entries (see README)."""
    K, p, n, r, kappa = D.K, D.p, D.n, D.r, D.kappa
    q = D.sqrtn
    two = K.from_int(2)
    S = {}
    for a in D.labels:
        for b in D.labels:
            ia, ib = D.idx(a), D.idx(b)
            pa, pb = D.psisign(a), D.psisign(b)
            if a == "1" or b == "1":
                S[(a, b)] = D.dim(b if a == "1" else a)
            elif a == b == "eps":
                S[(a, b)] = K.one
            elif (a == "eps" and ib) or (ia and b == "eps"):
                S[(a, b)] = two
            elif ia and ib:
                S[(a, b)] = (K.root_of_unity(ia * ib * r, n) + K.root_of_unity(-ia * ib * r, n)) * 2
            elif (ia and pb) or (pa and ib):
                S[(a, b)] = K.zero
            elif (a == "eps" and pb) or (pa and b == "eps"):
                S[(a, b)] = -q
            elif pa and pb:
                S[(a, b)] = q * (-kappa * jacobi(2, n) * (1 if pa == pb else -1) * (-1 if corrected else 1))
    return S


def shat(D):
    """D Shat D with AFT (smat) corrected: Shat_ab = sum_c F_{abb}^{a;c1} R_{ab}^c R_{ba}^c (F^-1)_{abb}^{a;1c}
    (all objects are self-dual and F is orthogonal, so (F^-1)^{1c} = F^{c1})."""
    K = D.K
    out = {}
    for a in D.labels:
        for b in D.labels:
            s = K.zero
            for c in D.labels:
                if (a, b, c) in D.N:
                    x = D.fentry(a, b, b, a, c, "1")
                    s = s + x * D.R[(a, b, c)] * D.R[(b, a, c)] * x
            out[(a, b)] = D.dim(a) * s * D.dim(b)
    return out


def main(argv):
    jobs = int(argv[argv.index("--jobs") + 1]) if "--jobs" in argv else 8
    only_r = int(argv[argv.index("--r") + 1]) if "--r" in argv else None
    only_k = int(argv[argv.index("--kappa") + 1]) if "--kappa" in argv else None
    ps = [int(x) for x in argv[1:] if x.isdigit() and argv[argv.index(x) - 1] not in ("--jobs", "--r", "--kappa")]
    allok = True
    for p in ps:
        n = 2 * p + 1
        for r in range(1, n, 2):
            if math.gcd(r, n) != 1 or (only_r and r != only_r):
                continue
            for kappa in (1, -1):
                if only_k and kappa != only_k:
                    continue
                t0 = time.time()
                D = Metaplectic(p, r, kappa, 1)
                Dm = Metaplectic(p, r, kappa, -1, field=D.K)
                orth = orthogonality_failures(D)
                pf, pc = pentagon_parallel(D, jobs)
                hx = {}
                for name, X in (("lam=+1", D), ("lam=-1", Dm)):
                    for mode in ("+", "-", "aft-literal"):
                        f, c = hexagon_failures(X, mode)
                        hx[(name, mode)] = len(f)
                tw = {}
                for name, X in (("lam=+1", D), ("lam=-1", Dm)):
                    th, pt = twists(X), printed_twists(X)
                    tw[name] = [a for a in X.labels if not (th[a] - pt[a]).is_zero()]
                    S, PS, CS = smatrix(X, th), printed_smatrix(X), printed_smatrix(X, corrected=True)
                    Sh = shat(X)
                    tw[name + " S entries differing from (sm1)-(sm4)"] = sorted(k for k in S if not (S[k] - PS[k]).is_zero())
                    tw[name + " ... from corrected (sm3),(sm4)"] = len([k for k in S if not (S[k] - CS[k]).is_zero()])
                    tw[name + " ... D Shat D (corrected (smat)) vs S"] = len([k for k in S if not (Sh[k] - S[k]).is_zero()])
                ok = not D.rotation_conflicts and not orth and not pf and not any(hx[(nm, m)] for nm in ("lam=+1", "lam=-1") for m in ("+", "-"))
                allok &= ok
                print(f"p={p} r={r} kappa={kappa:+d}: blocks {len(D.F)} ({sum(v == 'explicit' for v in D.source.values())} explicit, "
                      f"{sum(v == 'rotation' for v in D.source.values())} by rotation, {sum(v == 'O.W.' for v in D.source.values())} O.W.); "
                      f"rotation conflicts {len(D.rotation_conflicts)}; non-orthogonal {len(orth)}; "
                      f"pentagon {pc} equations, failures {len(pf)}; hexagon failures {hx}; "
                      f"twists differing from (tw1) {tw}  [{time.time() - t0:.0f} s]", flush=True)
    print("ALL OK" if allok else "FAILURES")


if __name__ == "__main__":
    main(sys.argv)
