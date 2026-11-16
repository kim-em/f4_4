\\ Endomorphism algebras of the forms of C(F4,4) over E and over k (paper, Section "Forms over the cubic field").
\\ Reads uvals.gp (u^(j)_X for j = 1, 3 and X = 0002, 0011, 0101; written by uvals_gp.py).
\\ Run: nix-shell -p pari --run 'gp -q endo.gp'
default(parisize, "2G");
z = Mod(y, polcyclo(13, y));
read("uvals.gp");
XS = ["0002", "0011", "0101"];
U1 = [uz1_0002, uz1_0011, uz1_0101];
U3 = [uz3_0002, uz3_0011, uz3_0101];
th = z + z^5 + z^8 + z^12;
Epol = subst(minpoly(th), 'x, 't);
bnfE = bnfinit(Epol, 1); nfE = bnfE.nf;
vc(v) = my(p = lift(v)); Col(vector(12, i, polcoef(p, i - 1, 'y)));
M3 = Mat([vc(th^0), vc(th), vc(th^2)]);
toE(v) = my(c = matinverseimage(M3, vc(v))); if(#c == 0 || Mod(1, polcyclo(13, y)) * (c[1] + c[2] * th + c[3] * th^2) != v, error("not in E")); Mod(c[1] + c[2] * 't + c[3] * 't^2, Epol);
et = z + z^12; kpol = subst(minpoly(et), 'x, 'v); bnfk = bnfinit(kpol, 1); nfk = bnfk.nf;
M6 = matconcat(vector(6, i, vc(et^(i - 1))));
tok(v) = my(c = matinverseimage(M6, vc(v))); if(#c == 0 || Mod(1, polcyclo(13, y)) * sum(i = 1, 6, c[i] * et^(i - 1)) != v, error("not in k")); Mod(sum(i = 1, 6, c[i] * 'v^(i - 1)), kpol);
uE = vector(3, i, toE(U1[i])); uE3 = vector(3, i, toE(U3[i]));
uk = vector(3, i, tok(U1[i]));
rE = vector(3, i, real(nfE.roots[i]));
sgnE(b) = vector(3, i, sign(real(subst(lift(b), 't, rE[i]))));
P2 = idealprimedec(nfE, 2); P13 = idealprimedec(nfE, 13);
print("E = Q(t), ", Epol, "; real places t = ", rE, "; primes above 2: ", #P2, " (f = ", P2[1].f, "), above 13: ", #P13, " (e = ", P13[1].e, ")");
print("class number of E: ", bnfE.no, ", of k: ", bnfk.no, "; bnfcertify(E): ", bnfcertify(bnfE), ", bnfcertify(k): ", bnfcertify(bnfk));
\\ invariant vector of the quaternion algebra (a, b)_E: [real places 1..3, P2, P13, other finite places with symbol -1]
qinv(a, b) = {
  my(r = vector(3, i, (sign(real(subst(lift(a), 't, rE[i]))) < 0 && sign(real(subst(lift(b), 't, rE[i]))) < 0)),
     f2 = (nfhilbert(nfE, lift(a), lift(b), P2[1]) < 0), f13 = (nfhilbert(nfE, lift(a), lift(b), P13[1]) < 0),
     other = List(), fa = idealfactor(nfE, lift(a) * lift(b)));
  for(i = 1, #fa~, my(pr = fa[i, 1]); if(pr.p != 2 && pr.p != 13 && nfhilbert(nfE, lift(a), lift(b), pr) < 0, listput(other, [pr.p, pr.f])));
  [r, f2, f13, Vec(other)];
}
print("\n== 1. consistency of the two base points");
for(i = 1, 3, print(XS[i], ": u1/u3 a square in E: ", #nfroots(nfE, 'x^2 - lift(uE[i] / uE3[i])) > 0, ";  (13, u1)_E split: ", nfhilbert(nfE, 13, lift(uE[i])) == 1));
print("\n== 2. B_X = (-1, u_X)_E (same local invariants as the cyclic algebra of brauer.gp)");
for(i = 1, 3, print(XS[i], ": signs of u at real places ", sgnE(uE[i]), ";  invariants of (-1,u)_E [real; P2; P13; others]: ", qinv(Mod(-1, Epol), uE[i])));
print("\n== 3. the pairing (g, u_X)_E for generators g of E^*/E^*2 restricted to {2, 13, infinity}-units");
fu = bnfE.fu; pi13 = nfbasistoalg(nfE, bnfisprincipal(bnfE, P13[1])[2]);
print("fundamental units: ", fu, ";  generator of P13: ", pi13, " (norm ", norm(pi13), ")");
GEN = [Mod(-1, Epol), Mod(lift(fu[1]), Epol), Mod(lift(fu[2]), Epol), Mod(2, Epol), Mod(lift(pi13), Epol)];
GN = ["-1", "eps1", "eps2", "2", "p13"];
for(g = 1, #GEN, print("  g = ", GN[g], " signs ", sgnE(GEN[g])); for(i = 1, 3, print("    (g, u_", XS[i], ")_E invariants [real; P2; P13; others]: ", qinv(GEN[g], uE[i]))));
print("\n== 4. examples e: End(S_X) in D_(chi_1 + 2 psi_e) is (-e, u_X)_E");
EX = [[1, "1 (D_1)"], [13, "13 (D_3)"], [-1, "-1"], [-13, "-13"], [2, "2"], [-2, "-2"], [26, "26"], [-26, "-26"]];
for(j = 1, #EX, my(e = Mod(EX[j][1], Epol)); print("  e = ", EX[j][2], ": ", vector(3, i, my(q = qinv(-e, uE[i])); if(vecsum(q[1]) + q[2] + q[3] + #q[4] == 0, "split", Str("division, ram ", q)))));
print("\n== 5. quadratic characters psi_e of E unramified outside 13 and infinity (e in <-1, eps1, eps2, p13>, nonsquare)");
{
for(m = 1, 15,
  my(e = prod(i = 1, 4, if(bittest(m, i - 1), [GEN[1], GEN[2], GEN[3], GEN[5]][i], 1)));
  my(D = rnfdisc(nfE, 'x^2 - lift(e))[1]);
  if(idealval(nfE, D, P2[1]) == 0,
    my(sg = sgnE(e), nsplit = sum(i = 1, 3, qinv(-e, uE[i])[1] == [0,0,0] && qinv(-e, uE[i])[2] == 0 && qinv(-e, uE[i])[3] == 0));
    print("  e = ", lift(e), "  [", m, "]  signs ", sg, "; ramified at 13: ", idealval(nfE, D, P13[1]) > 0,
          "; M_chi totally real (e totally negative): ", vecmax(sg) < 0, "; number of X with (-e,u_X) split: ", nsplit)));
}
print("\n== 6. parity check: N_{0002,0011}^{0101} = 3 is odd, so u_0002 u_0011 / u_0101 must be a square in k");
for(j = 1, 2, my(w = if(j == 1, U1, U3)); print("  j = ", 2 * j - 1, ": square in k: ", #nfroots(nfk, 'x^2 - lift(tok(w[1] * w[2] / w[3]))) > 0));
print("\n== 7. real places of k: embedding a (zeta13 -> exp(2 pi i a/13)), E-place, signs of u_X");
{
for(a = 1, 6,
  my(ev = 2 * cos(2 * Pi * a / 13), thv = 2 * cos(2 * Pi * a / 13) + 2 * cos(2 * Pi * 5 * a / 13),
     ip = vecsearch(vecsort(rE), thv) , iE = 0);
  for(i = 1, 3, if(abs(rE[i] - thv) < 1e-20, iE = i));
  print("  a = ", a, ": eta = ", ev, ", theta = ", thv, ", E-place ", iE, ", signs (u_0002, u_0011, u_0101) = ", vector(3, i, sign(real(subst(lift(uk[i]), 'v, ev))))));
}
print("\n== 8. the 2-torsion order of E^*/E^*2 at the generators (sanity: GEN independent mod squares)");
print("  ", matrank(Mat(vector(5, g, Col(concat([vector(3, i, sgnE(GEN[g])[i] < 0), [idealval(nfE, lift(GEN[g]), P2[1]) % 2, idealval(nfE, lift(GEN[g]), P13[1]) % 2]])))) * Mod(1, 2)), " (rank over F_2 of signs and valuations)");

print("\n== 9. all classes e in <-1, eps1, eps2, p13> mod squares (2 pairs trivially with every u_X): split pattern of (-e, u_X)_E, X = 0002, 0011, 0101");
{
for(m = 0, 15,
  my(e = Mod(1, Epol) * prod(i = 1, 4, if(bittest(m, i - 1), [GEN[1], GEN[2], GEN[3], GEN[5]][i], 1)));
  print("  e = ", lift(e), "  (-1,eps1,eps2,p13 exponents ", vector(4, i, bittest(m, i - 1)), "), signs ", sgnE(e), ": ",
        vector(3, i, my(q = qinv(-e, uE[i])); if(vecsum(q[1]) + q[2] + q[3] + #q[4] == 0, "split", Str("ram ", q[1], " P2:", q[2])))));
}
print("\n== 10. 13 and -5t+7 (the class unramified at 2) agree mod squares: ", #nfroots(nfE, 'x^2 - lift(Mod(13, Epol) / Mod(-5 * 't + 7, Epol))) > 0);
print("     u_0002 u_0011 / u_0101 is a square in E: ", #nfroots(nfE, 'x^2 - lift(uE[1] * uE[2] / uE[3])) > 0, ", 13 times a square: ", #nfroots(nfE, 'x^2 - 13 * lift(uE[1] * uE[2] / uE[3])) > 0);
print("     (2, u_X)_E split for all X: ", vector(3, i, nfhilbert(nfE, 2, lift(uE[i])) == 1));

print("\n== 11. local squares: u_X in E_P13^*2 (so (x, u_X) is trivial at 13 for every x); also at P2");
for(i = 1, 3, print("  ", XS[i], ": square at P13: ", nfislocalpower(nfE, P13[1], lift(uE[i]), 2), ", square at P2: ", nfislocalpower(nfE, P2[1], lift(uE[i]), 2)));
quit;
