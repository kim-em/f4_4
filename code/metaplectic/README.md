# so(2p+1)_2: F- and R-symbols from the published formulas

The F- and R-symbols of `C(so(2p+1), 2)`, and of every category with the same fusion rules
(the metaplectic categories), are taken from the closed formulas of

> E. Ardonne, P. E. Finch, M. Titsworth, "Classification of metaplectic fusion categories",
> arXiv:1608.03762 (v2), Symmetry 13 (2021) 2102 ("AFT" below; equation names are those of the
> arXiv source).

`so2p1_2.py` implements them exactly, for any `p` and all parameters `(r, kappa, lam)`, over
`Q(zeta_N)`, `N = 8(2p+1)`. For `p = 6` the result was compared with our own quantum group
computation (`data/B6_2_FR_exact.txt`) and with SageMath's pentagon solver: all three are
the same braided category, by explicit gauge transformations checked exactly entry by entry.

## Result

* **The formulas are correct as printed.** For `p = 1, ..., 6` and every `(r, kappa)`, exactly in
  `Q(zeta_N)`: the completion of the printed tables by AFT's rotation rule is consistent, every
  F-matrix is orthogonal, the pentagon holds (all equations, e.g. 294,526 for `p = 6`) and both
  hexagons hold for `lam = +1` and `lam = -1` (`check_formulas.py`, `../../logs/metaplectic/check_formulas_*.log`).
  The F- and R-symbols needed no correction. The misprints are in the surrounding equations (below).
* **C(so(13),2) agrees with AFT.** With the dictionary below, an explicit gauge `g` (one element
  of `Q(zeta_104)^+` per admissible triple, heights at most 20 bits, in
  `data/gauge_so13_ours_to_AFT.txt`) takes all 6,454 F-symbols (zeros included) and all 190
  R-symbols of `data/B6_2_FR_exact.txt` exactly to AFT's `(r, kappa, lam) = (1, -1, -1)`
  (`compare_so13.py`, `../../logs/metaplectic/compare_so13.log`, `../../logs/metaplectic/compare_so13_check_gauge.log`). So the two
  are equivalent as braided fusion categories (same labels, same F, same R up to gauge).
* **Third check, SageMath.** `FusionRing("B6", 2).get_fmatrix().find_orthogonal_solution()`
  (SageMath 10.9, 6.5 minutes) gives F in `Q(zeta_104)`; with Sage's own R-symbols, a gauge again
  maps our data exactly onto Sage's, F and R (`compare_sage.py`, `../../logs/metaplectic/compare_sage_B6_ours.log`),
  and maps AFT onto Sage (`../../logs/metaplectic/compare_sage_B6_aft.log`). For `p = 2, 3, 4, 5` the AFT data with
  the dictionary of `quantum_group_parameters(p)` agree exactly with Sage's solutions as well
  (`../../logs/metaplectic/compare_sage_B{2,3,4,5}_aft.log`; Sage run times 10 s, 30 s, 57 s, 3 min).

## Dictionary: AFT and the quantum group category

**Conventions are the same as ours.** AFT write `F_{abc}^{d;ef}` (`e` in `a b`, `f` in `b c`) and
`R_{ab}^c`; with `F[(X,Y,W,T)][(Z,0,0),(Z',0,0)] = F_{XYW}^{T;ZZ'}` and `R[(X,Y,Z)] = R_{XY}^Z`
(export convention of the paper, Section 3) the pentagon of `../f4/verify_fr.py` is AFT's (F3), our
hexagon (H2) is AFT's (R1), and our (H1) is (R1) for the reverse-inverse braiding. AFT's
F-matrices are real orthogonal, so splitting and fusion-map conventions give the same numbers.
`verify.py` checks exactly these equations; on our so(13) data it reproduces the 294,526
pentagon equations of `../f4/verify_fr.py` (`../../logs/f4/verify_fr_B6_2.log`) and both hexagons hold.

**Labels** (AFT order `1, eps, phi_1 .. phi_p, psi_+, psi_-`; Dynkin labels of `B_p`):

| AFT | so(2p+1) | so(13) | `h` (conformal weight) |
|---|---|---|---|
| `1` | `0` | `000000` | 0 |
| `eps` | `2 w_1` | `200000` | 1 |
| `phi_i`, `i < p` | `w_i` (exterior powers of the vector) | `100000, 010000, 001000, 000100, 000010` | `i(2p+1-i)/(2(2p+1))` |
| `phi_p` | `2 w_p` | `000002` | `p(p+1)/(2(2p+1))` |
| spinor | `w_p` | `000001` | `p/8` |
| other `psi` | `w_1 + w_p` | `100001` | `p/8 + 1/2` |

with the spinor `= psi_+` if `kappa = 1` and `= psi_-` if `kappa = -1`.

**Parameters.** `C(so(2p+1),2)` here means the category with twists `theta_X = exp(2 pi i h_X)`
(the quantum group at `q = exp(2 pi i/(4(2p+1)))` in the conventions of the paper, Section 2; for
p = 6 this is our `q = exp(2 pi i/52)`). It is AFT's

    (r, kappa, lam) = (1, (2 | 2p+1), -1),

i.e. `r = 1`, `kappa` the Jacobi symbol (`= -1` for p = 6, matching the Frobenius-Schur
indicator `-1` of the spinor: in our data `F_{VVV}^{V;11} = -1/sqrt 13`, AFT's `H_00 = kappa/sqrt n`),
and **the inverses** of the printed R-symbols (AFT's `R_{ab}^c = +- (-1)^{h_a+h_b-h_c}` is the
inverse of the usual `R^{ab}_c = +- exp(pi i (h_c - h_a - h_b))`). How it was found: for p = 6
all 288 choices of `(r, kappa, lam)` and of the 12 relabellings in `Aut(N) = (Z/13)^*/{+-1} x Z/2`
were screened numerically on gauge invariants (`R_{XX}^Z`, monodromies `R_{XY}^Z R_{YX}^Z`,
`F_{XXX}^{X;11}`, AFT's `X_p` invariants); F-invariants single out `kappa = -1`,
`r in {1, 3, 9}` (the squares mod 13, i.e. `r = z^2` up to sign), and the R-invariants
leave six matches, all related by relabellings. Two of them were proved exactly
(`compare_so13.py`): the one above with `phi_i = w_i`, and `(1, -1, +1)` with
`phi_i -> phi_{g(5i)}` and `psi_+ = spinor` (`5^2 = -1 mod 13`, which is why `lam = +1` and `-1`
give equivalent braidings for p = 6, as AFT note for `2p+1` admitting `sqrt(-1)`). For general
`p` the choice `(1, (2|2p+1), -1)` is the one whose twists are `exp(2 pi i h)`: AFT's (tw1) with
`lam = -1, r = 1` gives `theta_{phi_i} = exp(2 pi i h_{w_i})`, and `theta_{psi_pm}` equals
`exp(2 pi i p/8)` on one of `psi_pm` exactly when `kappa = (2|2p+1)`. For p = 6 the
equivalence with the quantum group category is proved (explicit gauge to our quantum group
computation). For p = 2, 3, 4, 5 the dictionary agrees exactly (explicit gauge, F and R) with
Sage's solution of the pentagon and hexagons, whose R-symbols Sage takes from the quantum group
formula; that Sage's solution is the quantum group category is Sage's claim, confirmed here for
p = 6 only.

## Misprints in AFT

The F- and R-symbol tables (Sections 4.2, 4.3) are correct. Notational points needed to read them:
`-1^{x}` means `(-1)^x` (e.g. (l1)-(l3), the `A(...)` arguments, (twotwonext)); the blocks not
printed (e.g. `F_{phi_i phi_j phi_j}^{phi_i}`, `F_{psi_+ psi_+ psi_-}^{psi_+}`, `F_{psi phi psi}^{phi}`)
are determined by (rotate) `F_{bcd}^a = (F_{abc}^d)^T`; "(O.W.) = (1)" then applies only to 1x1
blocks; `(-1)^{kappa lam r/2}` in (tw1) means `exp(pi i kappa lam r/2)`. Actual misprints:

1. **Fusion rules (Section 4.1):** `phi_i (x) psi_pm = psi_pm (x) psi_mp` should be
   `psi_+ (+) psi_-`. (The printed right side `eps + sum_j phi_j` has dimension `1 + 2p`, not
   `2 sqrt(2p+1)`.)
2. **Pentagon (F3):** the sum is printed over `h`, but `h` occurs on the right side; the summed
   variable is `g` (the one absent on the right). With that reading (F3) is the standard
   pentagon and AFT's F satisfy it.
3. **R gauge law (rgauge) and second hexagon (R2) are inconsistent with (fgauge) and (R1).**
   (rgauge) prints `R_{ab}^c -> g_{ab}^c R_{ab}^c / g_{ba}^c`; with (fgauge) the law that
   preserves (R1) is `g_{ba}^c R_{ab}^c / g_{ab}^c`. And (R2) is printed as (R1) with every R
   inverted but the lower indices kept; the correct second hexagon is (R1) for
   `R'_{ab}^c = (R_{ba}^c)^-1`. Printed (R2) and (rgauge) are consistent with each other, so the
   two hexagons as printed are covariant under opposite gauge laws. Evidence
   (`gauge_covariance.py`, `../../logs/metaplectic/gauge_covariance.log`): after a random gauge on AFT's p = 2 data,
   with R transformed by (rgauge) as printed, (R1) and the corrected (R2) fail and printed (R2)
   holds; with the inverse ratio, (R1) and the corrected (R2) hold and printed (R2) fails. AFT's
   own data satisfy both printed forms only because `R_{ab}^c = R_{ba}^c` in their gauge; our
   quantum group data (a non-symmetric gauge) satisfy (R1) and the corrected (R2) but not
   printed (R2).
4. **S-matrix formula (smat):** `Shat_ab = sum_c F_{ab*b}^{a;c1} R R F_{ab*b}^{a;1c}`: the last
   factor has upper index `1`, which is not in `a (x) b*` unless `a = b`, so as printed `Shat` is
   diagonal. It should be the inverse F-matrix entry `(F^-1)_{ab*b}^{a;1c}` (`= F_{ab*b}^{a;c1}`
   for AFT's orthogonal F). With that, `D Shat D` equals the S-matrix from the balancing formula
   `S_ab = theta_a^-1 theta_b^-1 sum_c N_ab^c d_c theta_c` for all `p <= 6` and all parameters.
5. **S-matrix (sm3), (sm4):** the sign of the `psi`-`psi` block is reversed. From AFT's own F
   and R (positive dimensions, `eps_psi = kappa`), for all `p <= 6` and all `(r, kappa, lam)`:
   `S_{psi_pm psi_pm} = +kappa (2|2p+1) q_psi` and `S_{psi_pm psi_mp} = -kappa (2|2p+1) q_psi`
   (printed: `-kappa (2|2p+1) q_psi` and `+kappa (2|2p+1) q_psi`). Independently, for p = 6 our
   data (whose S equals the Kac-Peterson S) give `S_{w6, w6} = +sqrt 13` with
   `kappa = -1`, `(2|13) = -1`. All other entries of (sm1)-(sm4) and the twists (tw1) are correct.
6. Trivial: in the arXiv source `\theta_\epsilon = 1` in (tw1) is typeset as `heta_\epsilon`.

## Field and heights of AFT's gauge (p = 6)

AFT's gauge is unitary (real orthogonal F). For p = 6, `(r, kappa, lam) = (1, -1, -1)`:

* F-symbols: fixed exactly by `{zeta -> zeta^{+-1}}` in `Gal(Q(zeta_104)/Q)`, so they generate
  `Q(zeta_104)^+ = Q(cos(pi/52))`, degree 24 (`= Q(zeta_52)^+(sqrt 2)`). In the power basis of
  `zeta_104` every coefficient is in `(1/26) Z` with numerator at most 4 in absolute value
  (denominators 1, 2, 13). `sqrt 2` is needed in 650 of the 2,880 blocks (576 with both `phi`
  and `psi` legs, 66 of the `Rep(D_13)` part (AFT's `A`, `C`), 8 with only `psi` legs); it comes
  from the unitary normalisation (`1/sqrt 2`, `2^{zeta(i,j)}`), not from the category: the
  gauge invariants of F lie in `Q(zeta_13)^+`. Sage's orthogonal solution has
  the same field and the same heights.
* R-symbols: fixed by `{1, 53}`, i.e. in `Q(zeta_52)` (degree 24); every one is `+- zeta_52^k`.
* The gauge from our data to AFT's lies in `Q(zeta_104)^+` (heights at most 20 bits).

So AFT's gauge needs degree 24 for F (ours: `Q(zeta_52)^+`, degree 12, 10 bits); neither is the
smallest field of definition.

## Generating any p

    python3 so2p1_2.py P OUT                # C(so(2P+1),2), Dynkin labels, P >= 2
    python3 so2p1_2.py P R KAPPA LAM OUT    # any AFT parameters, AFT labels

Output: fusion rules in the header, then `F a b c d | e f | b_0 ... b_{d-1}` and
`R a b c | b_0 ...`, coefficients in the power basis of `zeta = exp(2 pi i/(8(2P+1)))`,
conventions as above (export convention of the paper, Section 3). From Python:
`Metaplectic(p, r, kappa, lam)` has `.F[(a,b,c,d)][(e,f)]` and `.R[(a,b,c)]` (exact `cyclo.Elt`,
`complex(x)` for numbers), `.N`, `.labels`, `.dim(a)`; `quantum_group_parameters(p)` returns
`(r, kappa, lam, Dynkin -> AFT label)`.

## Reproducing

    python3 check_formulas.py 1 2 3 4 5 6 --jobs 10   # exact pentagon, hexagons, twists, S: all p, r, kappa
    python3 compare_so13.py data/gauge_so13_ours_to_AFT.txt   # solve and verify the gauge (25 s)
    python3 compare_so13.py --check data/gauge_so13_ours_to_AFT.txt
    sage sage_fmatrix.sage B<P> 2 sage_B<P>_2.txt && xz sage_B<P>_2.txt   # SageMath 10.9, P = 2..6
    python3 compare_sage.py sage_B6_2.txt.xz           # Sage against our data
    python3 compare_sage.py sage_B<P>_2.txt.xz --aft P

The SageMath outputs `sage_B<P>_2.txt.xz` are not included in the repository; `sage_fmatrix.sage`
regenerates them (Sage run times about 10 s, 30 s, 57 s, 3 min and 6.5 min for P = 2, ..., 6).

## Files

| file | content |
|---|---|
| `so2p1_2.py` | the AFT formulas (F, R, fusion rules, dimensions), quantum group dictionary, text export |
| `cyclo.py` | exact arithmetic in `Q(zeta_N)` (numpy integer vectors, exact inverses via norms) |
| `verify.py` | exact pentagon, hexagons (H+), (H-) and AFT's literal (R2), orthogonality |
| `check_formulas.py` | runs `verify.py` and the twist and S checks over all parameters |
| `gauge_covariance.py` | the evidence for misprint 3 |
| `ours.py` | reader for `data/B6_2_FR_exact.txt` over `Q(zeta_104)` |
| `compare_so13.py` | gauge solve (pivots of exponent +-1 only, no roots) and exact entrywise comparison |
| `compare_sage.py`, `sage_fmatrix.sage` | the Sage run and its exact comparison with ours and with AFT |
| `data/gauge_so13_ours_to_AFT.txt` | the gauge `g(X,Y,Z)` (power basis of `zeta_104`) |
| `data/B6_2_FR_exact.txt` | our C(so(13),2) data: F, R, Rinv over `Q(zeta_52)`, computed with the construction of `../f4` (generator the spinor `000001`) at all 24 roots of four primes and reconstructed exactly; the drivers of that multi-root run are not part of this repository; `python3 ../f4/verify_fr.py data/B6_2_FR_exact.txt` proves its pentagons and hexagons |
| `../../logs/metaplectic/` | outputs of the runs above |
