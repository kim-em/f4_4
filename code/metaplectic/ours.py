"""Reader for the quantum group data data/B6_2_FR_exact.txt (C(so(13),2)), read-only.

The file stores F in the basis 1, c_k = zeta^k + zeta^-k and R, Rinv in the power basis of
zeta = exp(2 pi i/52), all in the fusion-map ("export") convention of the paper (Section 3).  The data
are returned in the interface of verify.py over a field Q(zeta_M) with 52 | M, mapping
zeta_52 -> zeta_M^(M/52).  For this multiplicity-free category the export convention
F[(X,Y,W,T)][(Z,0,0),(Z',0,0)] is F_{XYW}^{T;ZZ'} of so2p1_2.py and R[(X,Y,Z)][0,0] is
R_{XY}^Z, with the same pentagon and hexagons (see verify.py).
"""
import os
from fractions import Fraction

from cyclo import Field

DEFAULT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "B6_2_FR_exact.txt")


class Ours:
    def __init__(self, path=DEFAULT, field=None):
        self.K = K = field or Field(104)
        assert K.N % 52 == 0
        m = K.N // 52
        z = [K.zeta(m * k) for k in range(52)]
        self.N, self.F, self.R, self.Rinv = {}, {}, {}, {}
        with open(path) as fh:
            for line in fh:
                if line.startswith("#"):
                    t = line[1:].split()
                    if t[:1] == ["objects:"]:
                        s = t.index(";")
                        self.labels, self.unit = t[1:s], t[s + 2]
                    elif t[:1] == ["N"] and len(t) == 5:
                        assert int(t[4]) == 1
                        self.N[(t[1], t[2], t[3])] = 1
                    continue
                parts = line.split("|")
                head = parts[0].split()
                cs = [Fraction(x) for x in parts[-1].split()]
                if head[0] == "F":
                    a, b = parts[1].split(), parts[2].split()
                    assert a[1:] == ["0", "0"] and b[1:] == ["0", "0"]
                    v = K.from_fraction(cs[0])
                    for k in range(1, 12):
                        if cs[k]:
                            v = v + (z[k] + z[52 - k]) * cs[k]
                    self.F.setdefault(tuple(head[1:]), {})[(a[0], b[0])] = v
                else:
                    assert parts[1].split() == ["0", "0"]
                    v = K.zero
                    for k in range(24):
                        if cs[k]:
                            v = v + z[k] * cs[k]
                    (self.R if head[0] == "R" else self.Rinv)[tuple(head[1:])] = v
