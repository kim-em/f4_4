#!/usr/bin/env python3
"""Write the PARI/GP script brauer.gp: endomorphism algebras of the simple objects of the forms.

Reads u_values.pkl (descent_u.py).  For each pi-fixed simple object X and j = 1, 3:
  * u_X lies in E = Q(theta), theta = zeta13 + zeta13^5 + zeta13^8 + zeta13^12 (written there in
    the eta basis; here as a + b theta + c theta^2, found by exact linear algebra and checked);
  * the signs of u_X at the three real places of E;
  * whether u_X is a norm from Q(zeta_13) (rnfisnorm; when it is, the preimage is checked by an
    exact norm computation);
  * the local Hasse invariants of the cyclic algebra A_X = (Q(zeta_13)/E, s, u_X), s: zeta13 ->
    zeta13^5, and its index (PARI alginit / alghassei / alghassef / algindex);
  * over k = Q(zeta_13)^+: the quaternion algebra (eta^2 - 4, u_X)_k with eta = zeta13 + zeta13^-1
    (so k(sqrt(eta^2 - 4)) = Q(zeta_13)), its Hasse invariants, and the classes of u_X^(1)/u_X^(3)
    in k^*/k^*2 (expected trivial).
Elements b of E are made integral by multiplying by d^4 (a norm) before alginit.

usage: brauer.py u_values.pkl OUT.gp
"""
import pickle, sys
from fractions import Fraction
from math import lcm

ORDER = ["0000", "0001", "0002", "0010", "1000", "0003", "0011", "0100", "1001", "0004", "0012",
         "0020", "0101", "1002", "1010", "2000"]


def gp_eta(e):
    """element of Q(zeta13) (variable z = Mod(y, polcyclo(13, y))) from eta coefficients."""
    return "(" + " + ".join(f"({Fraction(x)})*(z^{i + 1} + z^{13 - i - 1})" for i, x in enumerate(e) if Fraction(x)) + ")"


def main(argv):
    d = pickle.load(open(argv[1], "rb"))
    L = []
    w = L.append
    w("default(parisize, \"2G\");")
    w("\\\\ fields: z = zeta13, theta generates E, eta = z + 1/z generates k")
    w("z = Mod(y, polcyclo(13, y));")
    w("th = z + z^5 + z^8 + z^12;")
    w("Epol = subst(minpoly(th), 'x, 't);")
    w("print(\"E: minimal polynomial of theta: \", Epol, \"  disc \", factor(poldisc(Epol)));")
    w("nfE = nfinit(Epol);")
    w("facs = nffactor(nfE, polcyclo(13, 'x))[,1];")
    w("f = 0; for(i = 1, #facs, if(lift(lift(Mod(Mod(1, nfE.pol) * ('x + 'x^5 + 'x^8 + 'x^12), facs[i]))) == 't, f = facs[i]));")
    w("if(f == 0, error(\"no factor matching theta\"));")
    w("print(\"Q(zeta13) = E[x]/f, f = \", f);")
    w("rnf = rnfinit(nfE, f);")
    w("s = Mod('x^5, f);   \\\\ s: zeta13 -> zeta13^5, a generator of Gal(Q(zeta13)/E), s^2 = complex conjugation")
    w("Tn = rnfisnorminit(Epol, f, 1);")
    w("vc(v) = my(p = lift(v)); Col(vector(12, i, polcoef(p, i - 1, 'y)));")
    w("M3 = Mat([vc(th^0), vc(th), vc(th^2)]);")
    w("toE(v) = my(c = matinverseimage(M3, vc(v))); if(#c == 0 || Mod(1, polcyclo(13, y)) * (c[1] + c[2] * th + c[3] * th^2) != v, error(\"not in E\")); c[1] + c[2] * 't + c[3] * 't^2;")
    w("rE = polroots(Epol);")
    w("\\\\ k = Q(eta), eta = z + 1/z")
    w("et = z + z^12; kpol = subst(minpoly(et), 'x, 'v); nfk = nfinit(kpol);")
    w("M6 = matconcat(vector(6, i, vc(et^(i - 1))));")
    w("tok(v) = my(c = matinverseimage(M6, vc(v))); if(#c == 0 || Mod(1, polcyclo(13, y)) * sum(i = 1, 6, c[i] * et^(i - 1)) != v, error(\"not in k\")); sum(i = 1, 6, c[i] * 'v^(i - 1));")
    w("dk = 'v^2 - 4;  \\\\ k(sqrt(dk)) = Q(zeta13)")
    w("print(\"k: \", kpol, \";  d = eta^2 - 4, totally negative: \", vecmax(apply(r -> real(r)^2 - 4, polroots(kpol))) < 0);")
    for j in (1, 3):
        u = d[j]["u"]
        fixed = [X for X in ORDER if X in d[j]["fixed"]]
        w(f"print(\"\\n== j = {j}: Psi_{j} = (s, pi^{j}, c_{j}); braided: {d[j]['braided']}\");")
        for X in fixed:
            if X == "0000":
                continue
            w(f"uz = {gp_eta(u[X])};")
            w("b = toE(uz); den = denominator(content(b)); bi = b * den^4;")
            w(f"print(\"u_{X} = \", b, \"  in E;  (integral multiple by d^4, d = \", den, \")\");")
            w(f"print(\"  signs at the real places of E: \", vector(3, i, sign(real(subst(b, 't, rE[i])))));")
            w("fa = idealfactor(nfE, bi); print(\"  factorisation of the integral multiple: \", vector(#fa~, i, [fa[i,1].p, fa[i,1].f, fa[i,1].e, fa[i,2], if(fa[i,1].p == 13, \"ramified\", znorder(Mod(fa[i,1].p, 13)) / fa[i,1].f)]));")
            w("print(\"  (columns: p, f, e of the prime of E, exponent, local degree of Q(zeta13)/E at the prime (p != 13: unramified, degree = ord_13(p)/f))\");")
            w("nr = rnfisnorm(Tn, Mod(bi, Epol));")
            w("isnorm = (nr[2] == 1);")
            w("if(isnorm, al = Mod(lift(nr[1]), f); chk = Mod(1, f) * al * subst(lift(al), 'x, 'x^5) * subst(lift(al), 'x, 'x^25) * subst(lift(al), 'x, 'x^125); print(\"  u is a norm from Q(zeta13): exact check N(a) = u: \", chk == Mod(bi, f)), print(\"  u is not a norm from Q(zeta13) (rnfisnorm, flag 1: guaranteed for Galois extensions given the bnf)\"));")
            w("A = alginit(rnf, [s, bi]);")
            w(f"print(\"  cyclic algebra (Q(zeta13)/E, s, u_{X}): Hasse invariants at the real places \", alghassei(A), \", at finite places \", alghassef(A), \", index \", algindex(A));")
            w("pr = alghassef(A)[1]; if(#pr, print(\"    finite places: \", vector(#pr, i, [pr[i].p, pr[i].f, pr[i].e])));")
            w("bk = tok(uz);")
            w("print(\"  real embeddings of E agree with the ordering of nfE.roots: \", vector(3,i,real(nfE.roots[i])), \" vs \", vector(3,i,real(rE[i])));")
            w("Q = alginit(nfk, [dk, bk * denominator(content(bk))^2]);")
            w(f"print(\"  quaternion algebra (eta^2-4, u_{X})_k: Hasse invariants real \", alghassei(Q), \", finite \", alghassef(Q)[2], \", split: \", algissplit(Q));")
            w(f"bk{j}_{X} = bk;")
        w("")
    for X in ("0002", "0011", "0101"):
        w(f"rq = lift(Mod(bk1_{X}, kpol) / Mod(bk3_{X}, kpol)); print(\"u_{X}^(1)/u_{X}^(3) is a square in k: \", #nfroots(nfk, 'x^2 - rq) > 0);")
    w("quit;")
    open(argv[2], "w").write("\n".join(L) + "\n")


if __name__ == "__main__":
    main(sys.argv)
