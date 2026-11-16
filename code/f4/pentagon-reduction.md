# Reducing the pentagon to pentagons whose fourth object is a generator

## Statement

**Lemma.** Let `I` be a finite set of simple objects with fusion multiplicities `N_{ab}^c`
(associative: `sum_z N_{ab}^z N_{zc}^t = sum_z N_{bc}^z N_{az}^t`, so every F-matrix is square),
and let `F[(a,b,c,t)]` be invertible matrices indexed as in the paper (Section 3, the fusion-map convention). Let `V in I` and let
`S` be the set of `d in I` that are summands of `V^{(x) n}` for some `n >= 1` (equivalently: reachable
in `n >= 1` steps from the unit in the fusion graph of right multiplication by `V`). If the
indexed pentagon equation holds for all five-tuples `(a, b, c, V, e)` (all `a, b, c, e in I`),
then it holds for all five-tuples `(a, b, c, d, e)` with `d in S`.

In particular, if `S = I`, the pentagon holds for all five-tuples. `S = I` holds when every
simple object is a summand of a power of `V` *and the unit is a summand of some `V^{(x) n}`
with `n >= 1`* (for instance `V` self-dual: `1 < V (x) V`). If instead the unit is reachable only
with `n = 0`, the remaining case `d = 1` holds as soon as the F-matrices with unit third leg are
identities (`F[(x,y,1,t)][(t,i,0),(y,0,i)] = 1`, the normalisation of our data): the pentagon
`(a,b,c,1,e)` then reads `F[(a,b,c,e)] = F[(a,b,c,e)]`.

So the lemma as proposed holds, slightly more generally than stated: **no unit normalisation is
needed** when the unit is a summand of a positive power of `V`, and the hypothesis must include
the pentagons `(a, b, c, V, e)` for **all totals `e`** and all `a, b, c`. The hypothesis "every
simple object is reachable" cannot simply be dropped: the conclusion is only for `d in S`.

Nothing else is used: neither the triangle axiom, nor any pentagon used in constructing the
data, nor any property of the F-matrices beyond invertibility.

## The indexed pentagon (convention)

In all checkers here (`pent_reduced.py`, `pent_exact_subset.py`, `verify_fr.py`, `recursion.py`) the five-tuple `(A, B, C, D, E)` means the
equation `M1 M2 = N1 N2 N3` of square matrices on the label sets
`I0 = {(P,i1,Q,i2,i3)}` (bracketing `((AB)C)D -> E`), `I1 = {(P,i1,R,k,l)}` (`(AB)(CD)`),
`I2 = {(R,k,S,m,n)}` (`A(B(CD))`), `J1 = {(T,a,b,Q,i3)}` (`(A(BC))D`), `J2 = {(T,a,S,c,n)}`
(`A((BC)D)`), with
`M1 = (+)_{P,i1} F[(P,C,D,E)]`, `M2 = (+)_{R,k} F[(A,B,R,E)]`, `N1 = (+)_{Q,i3} F[(A,B,C,Q)]`,
`N2 = (+)_{T,x} F[(A,T,D,E)]`, `N3 = (+)_{S,n} F[(B,C,D,S)]`.
`D` is the **fourth tensor factor** and `E` the total.

## The model category

Let `C` be the category whose objects are `I`-graded finite-dimensional vector spaces
`X = (X_a)_{a in I}` (think `X = sum_a X_a (x) a`) and whose morphisms are graded linear maps.
For `a, b, c in I` let `H^{ab}_c` be the space with basis the multiplicity indices `0..N_{ab}^c - 1`.
Define `(X (x) Y)_c = sum_{a,b} X_a (x) Y_b (x) H^{ab}_c`, and on morphisms `f (x) g` componentwise;
this is a bifunctor. For simple `a, b, c` and total `t` the label sets of `((ab)c)_t` and
`(a(bc))_t` are `{(z,i,j)}` and `{(z',k,l)}`; define the associator
`alpha_{a,b,c}` on the `t`-component as the linear map with matrix `F[(a,b,c,t)]`, acting on
coefficient row vectors (`x -> x F`), so that composites are matrix products in the order
written. For arbitrary objects define
`alpha_{X,Y,W} = sum_{a,b,c} id_{X_a (x) Y_b (x) W_c} (x) alpha_{a,b,c}`. Then `alpha` is a
natural isomorphism `(X (x) Y) (x) W -> X (x) (Y (x) W)` (natural because it acts as the identity on
the argument spaces `X_a, Y_b, W_c`; an isomorphism because every `F` is invertible), whether
or not any pentagon holds.

For objects `X, Y, W, U` write `P(X,Y,W,U)` for the statement
`(alpha_{X,Y,W} (x) 1_U) . alpha_{X,Y (x) W,U} . (1_X (x) alpha_{Y,W,U}) = alpha_{X (x) Y,W,U} . alpha_{X,Y,W (x) U}`
(maps `((XY)W)U -> X(Y(WU))`, composed left to right as matrices act on row vectors).

**(i) Dictionary.** For simple `a, b, c, d`: on the `e`-component, `alpha_{a (x) b, c, d}` has matrix
`M1`, `alpha_{a, b, c (x) d}` has `M2`, `alpha_{a,b,c} (x) 1_d` has `N1`, `alpha_{a, b (x) c, d}` has `N2`
and `1_a (x) alpha_{b,c,d}` has `N3` (the direct sums over `P, R, Q, T, S` are exactly the additive
extension over the simple components of `a (x) b`, `c (x) d`, `b (x) c`, ..., with the identity on their
multiplicity spaces, and the label maps are those of the code). Hence `P(a,b,c,d)` holds if and
only if the indexed pentagon `(a, b, c, d, e)` holds for **every** total `e`.

**(ii) Additivity.** For arbitrary objects, all five maps in `P(X,Y,W,U)` are, by definition,
`sum_{a,b,c,d} id_{X_a (x) Y_b (x) W_c (x) U_d} (x) (the corresponding map for a, b, c, d)`. Hence
`P(X,Y,W,U)` holds iff `P(a,b,c,d)` holds for all simple `a, b, c, d` with `X_a, Y_b, W_c, U_d`
all nonzero. Multiplicities play no role beyond being nonzero. Also `P(a,b,c,d)` implies
`P(a,b,c,d) (x) 1_U` and `1_U (x) P(a,b,c,d)` (functoriality of `(x)`).

## The associahedron argument

For objects `X1, ..., X5` the boundary of the associahedron `K_5` is a polyhedral 2-sphere
whose 14 vertices are the bracketings of `X1 X2 X3 X4 X5`, whose 21 edges are single applications
of `alpha` (tensored with identities), and whose 9 faces are
* six pentagons: `P(X1,X2,X3,X4) (x) 1_{X5}`, `1_{X1} (x) P(X2,X3,X4,X5)`, and
  `P(X1 X2, X3, X4, X5)`, `P(X1, X2 X3, X4, X5)`, `P(X1, X2, X3 X4, X5)`, `P(X1, X2, X3, X4 X5)`;
* three squares, each an instance of naturality of `alpha` with respect to an `alpha` in one
  argument, or of bifunctoriality of `(x)` (two moves on disjoint parts of the bracketing).
  (Mac Lane, *Categories for the working mathematician*, VII.1-2.)

The squares commute in `C` for any invertible `F`. All edges are isomorphisms. A sphere minus
one open face is a polygonal disk; the boundary of a polytope is shellable (Bruggesser-Mani), so
the disk can be dismantled one face at a time along boundary arcs, and in a groupoid a polygonal
disk all of whose faces commute has commuting boundary (induction on the number of faces: the
boundary loop of the disk equals, up to conjugation by invertible edges, the boundary loop of
the disk with one face removed composed with that face's loop). Hence: **if eight of the nine
faces commute, so does the ninth.**

## Proof of the lemma

For simple `d` let `n(d)` be the least `n >= 1` with `d` a summand of `V^{(x) n}` (defined for
`d in S`). We show by induction on `n` that `P(a,b,c,d)` holds for all simple `a,b,c` and all
`d in S` with `n(d) = n`; by (i) this is the indexed pentagon for all totals.

`n = 1`: then `d = V`; this is the hypothesis (with (i)).

`n -> n + 1`: let `n(d') = n + 1`. Since `V^{(x)(n+1)} = V^{(x) n} (x) V`, `d'` is a summand of
`d (x) V` for some simple summand `d` of `V^{(x) n}`, so `n(d) <= n`. Apply the associahedron to
`(a, b, c, d, V)`. The faces other than `P(a, b, c, d (x) V)`:
* `P(a,b,c,d) (x) 1_V`: induction hypothesis and (ii);
* `1_a (x) P(b,c,d,V)`: fourth object `V`, hypothesis and (ii);
* `P(a (x) b, c, d, V)`, `P(a, b (x) c, d, V)`, `P(a, b, c (x) d, V)`: fourth object `V`; by (ii) they
  reduce to `P(x, y, z, V)` for simple `x, y, z`, the hypothesis;
* the three squares: commute in `C`.
So `P(a, b, c, d (x) V)` holds, and by (ii) `P(a, b, c, d')` holds for every simple summand `d'` of
`d (x) V`, in particular for our `d'`.

This proves the lemma. The answers to the four questions: (i) `D` is the fourth tensor factor,
and the hypothesis is needed for all totals `E` (and all `A, B, C`); (ii) yes, with multiplicities
(only nonvanishing matters); (iii) the unit case is not needed when `1` is a summand of a
positive power of `V`, and otherwise follows from identities `F[(x,y,1,t)]`; (iv) no pentagon from
the construction is used, only the hypothesis and invertibility of the F-matrices.

## Use for C(F4,4) and C(so(13),2)

`V = 0001` (the 26) is self-dual, `1 < V (x) V`, and every simple object is a summand of `V^{(x) n}`,
`n <= 4`; so `S` is everything. Invertibility of all 58,740 F-matrices of the exported F4 data is
checked modulo one prime (`pent_reduced.py`), which implies invertibility over the field. For
so(13)_2, `V = 000001` (the spinor) likewise generates and `1 < V (x) V`.

## Result for C(F4,4)

`pent_reduced.py` proves the pentagon exactly for all 63,841 five-tuples `(A, B, C, 0001, E)` of
the normal-form table (largest dimension 3380, height bounds up to 1176 bits, 1 to 38 primes per
five-tuple at the 12 real roots, 39,892 s with 2 workers, `deposit/logs/pent_reduced_descended.log`),
and checks that all 58,740 F-matrices are invertible. By the lemma the pentagon holds for all 1,030,051
five-tuples.
