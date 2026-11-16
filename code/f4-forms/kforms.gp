\\ Twisted forms C^L of C over k, L = k(sqrt d): End(S_X) = (d, u_X)_k for X = 0002, 0011, 0101.
\\ Every real place of k has some X with u_X < 0 (logs/f4-forms/endo.log, section 7), so all three split only if L
\\ is totally real.  Enumerate all totally real quadratic L/k by the norm of the conductor
\\ (= norm of the relative discriminant), up to BOUND, via ray class groups with finite modulus.
\\ Run: nix-shell -p pari --run 'gp -q kforms.gp'
default(parisize, "2G");
z = Mod(y, polcyclo(13, y));
read("uvals.gp");
U1 = [uz1_0002, uz1_0011, uz1_0101];
vc(v) = my(p = lift(v)); Col(vector(12, i, polcoef(p, i - 1, 'y)));
et = z + z^12; kpol = subst(minpoly(et), 'x, 'v); bnfk = bnfinit(kpol, 1); nfk = bnfk.nf;
M6 = matconcat(vector(6, i, vc(et^(i - 1))));
tok(v) = my(c = matinverseimage(M6, vc(v))); sum(i = 1, 6, c[i] * 'v^(i - 1));
uk = vector(3, i, tok(U1[i]));
print("k = Q(v), ", kpol, ", disc ", nfk.disc, " = ", factor(nfk.disc), "; class number ", bnfk.no, "; certified: ", bnfcertify(bnfk));
print("narrow class group of k: ", bnrinit(bnfk, [1, [1, 1, 1, 1, 1, 1]]).cyc);
BOUND = 300000;
L = ideallist(nfk, BOUND);
cnt = 0; found = List(); patt = vector(4);
{
for(n = 1, BOUND, for(j = 1, #L[n],
  my(m = L[n][j], bnr = bnrinit(bnfk, m), subs);
  if(bnr.cyc == [] || #select(c -> c % 2 == 0, bnr.cyc) == 0, next);
  subs = subgrouplist(bnr, [2], 1);
  for(h = 1, #subs,
    if(bnrconductor(bnr, subs[h])[1] != m, next);
    my(P = bnrclassfield(bnr, subs[h])[1], d = lift(poldisc(P)), s = vector(3, i, nfhilbert(nfk, d, uk[i]) == 1), ns = vecsum(s));
    cnt++; patt[ns + 1]++;
    if(n < 7000, print("  L with N(conductor) = ", n, ": number of X with (d,u_X)_k split: ", ns, ", pattern ", s));
    if(ns == 3, listput(found, [n, d]); print("all three split: N(conductor) = ", n, " = ", factor(n), ", d = ", d, ", |disc L| = disc(k)^2 * ", n)))));
}
print("totally real quadratic L/k with N(conductor) <= ", BOUND, ": ", cnt, "; number with 0, 1, 2, 3 of the (d,u_X)_k split: ", patt);
print("d = 2: ", vector(3, i, nfhilbert(nfk, 2, uk[i]) == 1), ", conductor norm ", idealnorm(nfk, rnfconductor(bnfk, 'x^2 - 2)[1][1]));
quit;
