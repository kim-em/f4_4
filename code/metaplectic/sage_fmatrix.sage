import sys, time
ct, lev, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
t0 = time.time()
R = FusionRing(ct, lev)
f = R.get_fmatrix()
f.find_orthogonal_solution(verbose=False)
emb = f.get_qqbar_embedding()
fv = f.get_fvars()
with open(out, "w") as fh:
    fh.write("# %s level %d, field %s, time %.1f s\n" % (ct, lev, f.field(), time.time() - t0))
    fh.write("# field defining polynomial %s\n" % (f.field().defining_polynomial() if hasattr(f.field(), "defining_polynomial") else ""))
    for b in R.basis():
        fh.write("# object %s weight %s twist %s qdim %s\n" % (b, b.weight(), b.twist(), b.q_dimension()))
    for key, val in fv.items():
        a, b, c, d, x, y = key
        z = CC(emb(val)) if val in f.field() else CC(val)
        z = ComplexField(200)(emb(val))
        fh.write("F %s %s %s %s %s %s | %s %s | %s\n" % (a, b, c, d, x, y, z.real(), z.imag(), val))
    for a in R.basis():
        for b in R.basis():
            for c in R.basis():
                if R.Nk_ij(a, b, c):
                    z = ComplexField(200)(R.r_matrix(a, b, c))
                    fh.write("R %s %s %s | %s %s\n" % (a, b, c, z.real(), z.imag()))
print("done", time.time() - t0)
