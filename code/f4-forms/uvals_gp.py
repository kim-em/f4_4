#!/usr/bin/env python3
"""Write uvals.gp: the scalars u^(j)_X (j = 1, 3; X = 0002, 0011, 0101) of data/u_values.pkl as
elements of Q(zeta13) (z = Mod(y, polcyclo(13, y))), for endo.gp.  usage: uvals_gp.py u_values.pkl OUT"""
import pickle, sys
from fractions import Fraction
d = pickle.load(open(sys.argv[1], "rb"))
out = []
for j in (1, 3):
    for X in ("0002", "0011", "0101"):
        e = d[j]["u"][X]
        s = " + ".join(f"({Fraction(x)})*(z^{i + 1} + z^{12 - i})" for i, x in enumerate(e) if Fraction(x))
        out.append(f"uz{j}_{X} = {s};")
open(sys.argv[2], "w").write("\n".join(out) + "\n")
