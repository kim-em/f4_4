#!/usr/bin/env python3
"""Evidence for misprint 3 of README.md: which R gauge law and which second hexagon go with
AFT's (fgauge) and (R1).

usage: gauge_covariance.py [P]      (default P = 2)

Applies a random gauge g (values in Q(zeta_N), g = 1 on triples with a unit label) to AFT's data
(r, kappa, lam) = (1, 1, 1): F by (fgauge), R once by (rgauge) as printed,
R -> g_{ab}^c R / g_{ba}^c, and once by the inverse ratio, R -> g_{ba}^c R / g_{ab}^c; then
counts failures (up to 20) of the pentagon, of (R1) = (H+), of the corrected second hexagon
(H-) and of (R2) as printed ("aft-literal").  Also reports on the quantum group data for p = 6.
"""
import random
import sys

from so2p1_2 import Metaplectic
from verify import hexagon_failures, pentagon_failures


def main(argv):
    p = int(argv[1]) if len(argv) > 1 else 2
    random.seed(1)
    D = Metaplectic(p, 1, 1, 1)
    K = D.K
    g = {t: K.one if "1" in t[:2] else K.from_int(random.randint(2, 9)) + K.zeta(random.randint(0, K.N - 1)) for t in D.N}
    gi = {t: v.inverse() for t, v in g.items()}
    for (a, b, c, d), blk in D.F.items():
        for (e, f) in blk:
            blk[(e, f)] = blk[(e, f)] * g[(a, b, e)] * g[(e, c, d)] * gi[(b, c, f)] * gi[(a, f, d)]
    R0 = dict(D.R)
    print(f"p = {p}, random gauge applied to F by (fgauge); pentagon failures: {len(pentagon_failures(D)[0])}")
    for name, law in (("(rgauge) as printed, g_ab R / g_ba", lambda a, b, c: g[(a, b, c)] * gi[(b, a, c)]),
                      ("inverse ratio,       g_ba R / g_ab", lambda a, b, c: g[(b, a, c)] * gi[(a, b, c)])):
        D.R = {(a, b, c): R0[(a, b, c)] * law(a, b, c) for (a, b, c) in R0}
        res = {m: len(hexagon_failures(D, m)[0]) for m in ("+", "-", "aft-literal")}
        print(f"  R by {name}: failures (R1) {res['+']}, corrected (R2) {res['-']}, printed (R2) {res['aft-literal']}")
    if p == 2:
        from ours import Ours
        O = Ours()
        res = {m: len(hexagon_failures(O, m)[0]) for m in ("+", "-", "aft-literal")}
        print(f"quantum group data so(13)_2 (data/B6_2_FR_exact.txt): failures (R1) {res['+']}, corrected (R2) {res['-']}, printed (R2) {res['aft-literal']}")


if __name__ == "__main__":
    main(sys.argv)
