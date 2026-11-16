#!/usr/bin/env python3
"""Square classes for an exact orthogonal gauge of C(F4,4), from the exact forms G over k.

(1) Lower bound.  In any gauge with G = 1, an entry F'_ab whose four multiplicity spaces are one
    dimensional satisfies F'_ab^2 = F'_ab (F'^-1)_ba = I_ab (gauge invariant), and
    I_ab = F_ab^2 n_R(b) / n_L(a) in terms of the normal-form table and the forms.  So the classes
    [n_R(b) n_L(a)] must become squares in the field of every orthogonal table.
(2) Determinant classes.  An orthogonal gauge u with u^T u = t G (t in the torus) forces
    det(t G) to be a square: for even-dimensional spaces [det G] is torus invariant, for odd ones
    it moves like a one-dimensional norm.  The torus-invariant part of these classes.
(3) LDL^T diagonals (the naive upper bound: Gram-Schmidt without roots plus one square root per
    diagonal entry), with the best torus.
Characters at degree-one primes of k, all p = +-1 mod 13 (kprimes.py).

usage: f4_classes.py G_exact.pkl NCHARS
"""
import os, sys, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from f4_forms_analyse import load, f2_rank, reduce_vec, ldl_modp, UNIT
from kprimes import kprimes, red, char, candidates
import f4_forms
import f4_canon

G = load(sys.argv[1])
nch = int(sys.argv[2])
dim = {t: len(M) for t, M in G.items()}
spaces = sorted(t for t in G if UNIT not in t[:2])
objs = sorted({t[0] for t in G} | {t[2] for t in G})
oi = {X: i for i, X in enumerate(objs)}
H = f4_canon.load_H(1073741189, 1, f4_canon.readers())
ents = []
for (X, Y, W, T), (rows, cols, M) in H.items():
    for a, (Z, i, j) in enumerate(rows):
        for b, (Zp, k, l) in enumerate(cols):
            sp = [(X, Y, Z), (Z, W, T), (Y, W, Zp), (X, Zp, T)]
            if M[a, b] % 1073741189 and all(dim[s] == 1 for s in sp):
                ents.append([s for s in sp if UNIT not in s[:2]])
print(f"H family: {len(ents)} nonzero entries with four one-dimensional spaces")
cand = candidates()
t0 = time.time()
diag = {}          # (t, i) -> bit vector
cv = {k: 0 for k in cand}
pos = 0
for p, eta in kprimes(nch):
    res = {}
    for t in spaces:
        D = ldl_modp([[red(c, p, eta) for c in row] for row in G[t]], p)
        if D is None:
            break
        res[t] = [char(d, p) for d in D]
    else:
        for t in spaces:
            for i, b in enumerate(res[t]):
                diag[(t, i)] = diag.get((t, i), 0) | (b << pos)
        for k, v in cand.items():
            cv[k] |= char(red(v, p, eta), p) << pos
        pos += 1
print(f"{pos} degree-one primes of k (p = +-1 mod 13), {time.time() - t0:.0f} s")


def inv_span(ts, vec):
    """torus-invariant combinations: kernel of a -> sum_t a_t (e_X + e_Y + e_Z), mapped by vec."""
    piv, out = {}, []
    for t in ts:
        r = ((1 << oi[t[0]]) ^ (1 << oi[t[1]]) ^ (1 << oi[t[2]])) & ~(1 << oi[UNIT])
        v = vec[t]
        while r:
            h = r.bit_length() - 1
            if h in piv:
                r ^= piv[h][0]; v ^= piv[h][1]
            else:
                piv[h] = (r, v); break
        else:
            out.append(v)
    return out


one = {t: diag[(t, 0)] for t in spaces if dim[t] == 1}
cls = []
for sp in ents:
    v = 0
    for s in sp:
        v ^= one[s]
    cls.append(v)
rk1, b1 = f2_rank(cls)
rk1c, _ = f2_rank(inv_span([t for t in spaces if dim[t] == 1], one))
det = {}
for t in spaces:
    v = 0
    for i in range(dim[t]):
        v ^= diag[(t, i)]
    det[t] = v
rk2, b2 = f2_rank(inv_span([t for t in spaces if dim[t] % 2], det) + [det[t] for t in spaces if dim[t] % 2 == 0])
rk12, _ = f2_rank(list(b1.values()) + list(b2.values()))
first = {t: diag[(t, 0)] for t in spaces}
rk3, b3 = f2_rank(inv_span(spaces, first) + [diag[(t, i)] ^ diag[(t, 0)] for t in spaces for i in range(1, dim[t])])
ndiag = sum(dim[t] for t in spaces if dim[t] > 1)
print(f"(1) entry invariants I_ab: F_2 rank {rk1} (torus-invariant classes of all one-dimensional norms: {rk1c})")
print(f"(2) torus-invariant determinant classes: rank {rk2}; together with (1): {rk12}")
print(f"(3) LDL^T diagonals with the best torus: rank {rk3} of at most {pos} detectable "
      f"({ndiag} diagonal entries in the {sum(1 for t in spaces if dim[t] > 1)} spaces of dimension > 1)")
for name, basis in (("(1)", b1), ("(2)", b2)):
    chosen, cb = [], {}
    for k, v in cv.items():
        if reduce_vec(v, basis) == 0 and reduce_vec(v, cb) != 0:
            chosen.append(k)
            _, cb = f2_rank(list(cb.values()) + [v])
    print(f"{name}: basis among -1, 2, 3, 13, [n], [2][n]: {chosen}; spans: {all(reduce_vec(v, cb) == 0 for v in basis.values())}")
