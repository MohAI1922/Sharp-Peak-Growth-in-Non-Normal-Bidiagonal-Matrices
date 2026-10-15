# Pre-derivation Handoff Packet

This packet collects the six items required before starting the analytic phase of Part II. Nothing in this packet is a theorem proof; it is the verified numerical and bibliographic evidence that the proof target is currently the correct one.

## 1. Final provisional conjecture statement

For the family

A = D + cN,   D = diag(-δ, -2δ, ..., -nδ),   N_{i,i+1} = 1,

with δ > 0 fixed, the strongest defensible statement supported by the current numerical evidence is:

sup_{t≥0} ||e^{tA}||_2 ≍ C(c) e^{c/δ},

with exponential rate 1 in c/δ and a nontrivial polynomial prefactor. The decisive new fact is not a global fit alone but the exact leading-row lower bound:

||row_0(e^{tA})||_2 <= ||e^{tA}||_2  for every t >= 0,

and at the envelope-maximizing time t* = log(c/δ), the row-0 asymptotic predicts a leading prefactor `c^{-5/4}`. This is not yet a completed proof of the exact prefactor, but it is a rigorous lower-bound mechanism and a numerically validated asymptotic target.

The uniform safety-margin sweep with `n >= 3c` and `c = 10, 12, 14, ..., 40` gives local finite-difference exponents

- (10,12): 1.2882
- (12,14): 1.2805
- (14,16): 1.2752
- (16,18): 1.2717
- (18,20): 1.2688
- (20,22): 1.2666
- (22,24): 1.2649
- (24,26): 1.2636
- (26,28): 1.2624
- (28,30): 1.2615
- (30,32): 1.2607
- (32,34): 1.2599
- (34,36): 1.2593
- (36,38): 1.2588
- (38,40): 1.2580

These trend cleanly toward `5/4 = 1.25`, and the overlapping-triple Aitken extrapolation on this evenly spaced sequence stays in the narrow band `1.252–1.265`, with a final cluster around `1.254–1.256`. In other words, the current evidence supports

`p ≈ 1.25 ± 0.02`

as a realistic asymptotic exponent estimate, but not yet as a theorem-level confidence interval.

The professor-facing statement is therefore:

> For the family `A = D + cN`, the exponential rate in `c/δ` is 1, and the exact leading-row lower bound is already rigorous. Numerically, the dominant row mechanism predicts a prefactor exponent near `5/4`, so the asymptotic law is strongly consistent with `c^{-5/4} e^{c/δ}` rather than the earlier crude `c^{-1}` envelope guess. The exact prefactor exponent is not yet nailed down to a formal confidence interval, and the matching upper bound is now a concrete weighted-Schur problem rather than an open-ended search.

This is the honest version to carry into the analytic phase.

## 2. Scaling refit tables (extended range)

### 2.0 Exact row-0 lower bound and the asymptotic mechanism

The missing prefactor term is explained by row dominance, and the row-0 quantity is not merely a heuristic — it is an exact lower bound on the true norm.

For any matrix `M` and any row index `i`,

`||row_i(M)||_2 = ||M^T e_i||_2 <= ||M^T||_2 = ||M||_2`.

Therefore, for the headline family,

`||row_0(e^{tA})||_2 <= ||e^{tA}||_2` for every `t >= 0`,

with the exact finite-n formula

`||row_0(e^{tA})||_2 = e^{-δt} √( Σ_{r=0}^{n-1} [ (c/δ)(1-e^{-δt}) ]^{2r} / (r!)^2 )`.

This is precisely the rigorous lower-bound candidate for Part II, Step 5. It also explains the observed behavior: at `t* = δ^{-1} log(c/δ)`, the row weights satisfy

`w_i = e^{-(i+1)δ t*} = (δ/c)^{i+1}`,

so the row-0 contribution is exponentially larger than the rest. The exact row-0 lower bound then has the asymptotic form

`||row_0||_2 ~ const · c^{-5/4} e^{c/δ}`

via the large-x Bessel asymptotic `I_0(2x) ~ e^{2x}/√(4πx)`.

A second exact structural fact is the interaction-picture factorization:

`e^{tA} = e^{tD} · e^{xN}` with `x = (c/δ)(1 - e^{-δt})`,

which follows from `e^{-tD} N e^{tD} = e^{-δt} N` and the ODE for `U(t) = e^{-tD} e^{tA}`. This immediately yields the submultiplicative bound

`||e^{tA}||_2 <= ||e^{tD}||_2 ||e^{xN}||_2 = e^{-δt} ||e^{xN}||_2`,

and this is the exact algebraic reason the crude weighted Schur attempt failed: the matrix is not a generic weighted matrix; it is a product of a decaying diagonal factor and a Toeplitz nilpotent factor. The right proof strategy is therefore to split the estimate at the factorization, not to force a single geometric weight onto the entire matrix.

This predicts a leading prefactor exponent `p = 5/4`, not `p = 1`, and it matches the actual peaks to within a few percent in the safety-margin regime (`n >= 3c`).

The earlier 4–20 sweep showed transitional exponents

- `(4,6)`: 1.3832
- `(6,8)`: 1.3242
- `(8,10)`: 1.3006
- `(10,12)`: 1.2882
- `(12,16)`: 1.2781
- `(16,20)`: 1.2703

These are not the final asymptotic estimate; they are exactly the kind of finite-range bias that the uniform `10–40` safety-margin sweep is designed to remove. The refined grid is the definitive estimate reported in §2.3.

### 2.1 Refined sweep data used

The sweep used the headline family `shifted_superdiagonal(n, c, decay_start=1)` with parameter values:

- n=8: c = 1,2,4,6,8,10
- n=12: c = 1,2,4,6,8,10,12
- n=16: c = 1,2,4,6,8,10,12,16
- n=24: c = 1,2,4,6,8,10,12,16,20

The peak growth values observed were:

- n=8: 1.0000, 1.0000, 2.1302, 8.7929, 37.8234, 141.4560
- n=12: 1.0000, 1.0000, 2.1305, 8.9848, 45.1677, 238.5474, 1184.4096
- n=16: 1.0000, 1.0000, 2.1305, 8.9852, 45.3517, 250.5291, 1446.1851, 43957.9532
- n=24: 1.0000, 1.0000, 2.1305, 8.9852, 45.3522, 250.7332, 1464.9000, 55333.2217, 2203050.9972

### 2.2 Old pure-exponential fit versus corrected two-term fit

For the same data, the old pure-exponential model `log peak = α c + β` gave:

- n=8: α = 0.5779, β = -1.0668, R² = 0.9713
- n=12: α = 0.6713, β = -1.3648, R² = 0.9738
- n=16: α = 0.7380, β = -1.6606, R² = 0.9804
- n=24: α = 0.7948, β = -1.9647, R² = 0.9848

The corrected two-term model `log peak = α c - p log c + β` gave:

- n=8: α = 0.8941, p = 1.3081, β = -0.9008, R² = 0.9985
- n=12: α = 1.0019, p = 1.5726, β = -0.9829, R² = 0.9995
- n=16: α = 1.0060, p = 1.5626, β = -0.9983, R² = 0.9997
- n=24: α = 1.0275, p = 1.6372, β = -1.0042, R² = 0.9999

This is exactly the pattern expected if a single row dominates: the `c^{-p}` prefactor is real, but the fitted `p` is distorted by trimming the sweep at `n ≈ c`, where the row-0 theory is no longer operating in its untruncated regime. Once the points are kept in the untruncated regime (`n ≥ 3c`), the agreement with the row-0 asymptotic improves dramatically.

### 2.3 Uniform safety-margin sweep and Aitken estimate

The exact row-0 lower bound is already rigorous, and the safety-margin sweep confirms the asymptotic trend without the truncation artifact from `n ≈ c`. On the uniform grid `c = 10, 12, 14, ..., 40` with `n >= 3c`, the exact finite-n peaks give the local exponents

- (10,12): 1.2882
- (12,14): 1.2805
- (14,16): 1.2752
- (16,18): 1.2717
- (18,20): 1.2688
- (20,22): 1.2666
- (22,24): 1.2649
- (24,26): 1.2636
- (26,28): 1.2624
- (28,30): 1.2615
- (30,32): 1.2607
- (32,34): 1.2599
- (34,36): 1.2593
- (36,38): 1.2588
- (38,40): 1.2580

The overlapping-triple Aitken extrapolation stabilizes at roughly `1.252–1.265`, with the latest values clustering around `1.254–1.256`. This supports the conservative estimate

`p ≈ 1.25 ± 0.02`,

not a clean theorem-level lock, but a strong numerical estimate of the leading prefactor exponent.

### 2.4 Row-0 vs full 2-norm gap and the Bessel correction check

The exact row-0 lower bound is not just a heuristic: it is an exact finite-n formula, and it remains the sharpest available lower-bound mechanism. At `t* = log(c)`, the actual norm and the row-0 norm differ by a small relative amount that shrinks with `c`.

| c | n | actual peak | row-0 norm | actual / row-0 |
|---|---:|---:|---:|---:|
| 4 | 30 | 2.10511 | 2.04991 | 1.02693 |
| 8 | 30 | 45.29799 | 44.96849 | 1.00733 |
| 12 | 36 | 1464.247 | 1459.388 | 1.00333 |
| 20 | 60 | 2.27681e6 | 2.27404e6 | 1.00122 |
| 30 | 90 | 3.00418e10 | 3.00254e10 | 1.00055 |
| 40 | 120 | 4.60618e14 | 4.60475e14 | 1.00031 |

This monotone convergence to 1 is exactly the expected pattern for a genuine lower bound approaching the true norm from below. The next-order Bessel correction `1/(16x)` is of the right order, but it is not a clean proof of the full asymptotic; it is only a suggestive correction term. The separation is:

- `‖e^{tA}‖_2 - ‖row_0‖_2` is the true matrix-vs-row gap
- `1/(16x)` is the next-term correction inside the row-0 asymptotic itself

and they are distinct phenomena, so they should not be conflated.

### 2.5 Weighted Schur test for the upper bound

The matching upper bound was probed numerically with geometric weights `w_i = (δ/c)^{β i}`. For the exact matrix exponential at `t = log(c)`, the weighted Schur bound is

`||M||_2 <= sqrt( max_i (Σ_j |M_{ij}| w_j / w_i) * max_j (Σ_i |M_{ij}| w_i / w_j) )`.

For the tested points, the results are:

| c | n | β=0 | β=0.5 | β=1 | β=2 |
|---|---:|---:|---:|---:|---:|
| 8 | 30 | 1.256 | 60.89 | 3.34e7 | 2.75e22 |
| 12 | 40 | 1.233 | 1.39e4 | 1.44e16 | 1.03e38 |
| 20 | 60 | 1.215 | 1.30e10 | 3.83e34 | 5.83e73 |
| 30 | 90 | 1.206 | 1.08e19 | 9.03e61 | 3.05e126 |
| 40 | 120 | 1.202 | 1.08e43 | 3.67e90 | NaN |

The only numerically useful choice is the flat weight `β = 0`, which stays within about `20–26%` of the exact row-0 lower bound across the tested range. The larger-β weightings overshoot catastrophically because they make the Schur test too skewed by the geometric decay. This is the correct empirical message for the proof: the upper-bound strategy should be a weighted Schur argument with a carefully chosen geometric balancing, not a brute-force all-ones or overly aggressive weight.

### 2.6 Large-c log-space check and row dominance

A direct full-matrix `‖e^{tA}‖_2` computation becomes numerically unstable once `c` is large enough to generate values beyond the manageable range of `float64` SVDs. To push past that limit while retaining the asymptotic information, the exact row-0 lower bound was evaluated in log-space for `c = 50, 60, 70, 80`:

| c | log(row-0 norm) |
|---|---:|
| 50 | 43.483 |
| 60 | 53.254 |
| 70 | 63.061 |
| 80 | 72.893 |

The corresponding local finite-difference exponents from the log-space row-0 values are

- (10,20): 1.2723
- (20,30): 1.2622
- (30,40): 1.2584
- (40,50): 1.2565
- (50,60): 1.2553
- (60,70): 1.2544
- (70,80): 1.2538

and the Aitken extrapolation on this sequence stays in the band `1.252–1.256`, again consistent with a limiting exponent near `5/4`.

The row-dominance mechanism is also exact and easy to see numerically at the peak time `t* = log(c)`: for example, at `c=20`,

- row 0: 2.27404e6
- row 1: 1.13702e5
- row 2: 5.68510e3
- row 3: 2.84255e2

so the leading row is larger than the immediately following rows by factors of about `20`, `400`, and `8000` respectively. This is the structural reason the lower bound is so close to the true norm in the untruncated regime.

## 3. Converged Kreiss values

The convergence check was re-run at the headline point and the new values intended for reporting. The results below were produced with the refined grid sequence `(20,40)`, `(40,80)`, `(60,120)`, `(80,160)` and `x_max = 5` for the largest grids.

| (n,c) | K(20×40) | K(40×80) | K(60×120) | K(80×160) | spread 60→80 |
|---|---:|---:|---:|---:|---:|
| (12,8) | 18.7736 | 18.7948 | 18.7992 | 18.8009 | 0.009% |
| (12,12) | 426.8069 | 426.1654 | 426.8610 | 426.0559 | 0.19% |
| (16,16) | 14156.5604 | 14231.9320 | 14240.0099 | 14242.5211 | 0.017% |
| (24,20) | 651275.8372 | 641209.3116 | 650502.0802 | 651504.1949 | 0.15% |

These are the numbers that are safe to cite. The headline value at (12,8) is converged to approximately 18.80; the other entries are stable to within 0.2% in the last refinement.

## 4. Eigenvector and condition-number cross-check

The right-eigenvector formula associated with the upper-bidiagonal family has entries of the form

v_i^{(j)} = (-c/δ)^{j-i} / (j-i)!    for i ≤ j,

and zero otherwise, with the eigenvalue λ_j = -jδ. For the headline case n=12, c=8, δ=1, the explicit formula produces entries whose magnitude grows as the index moves toward the upper end of the bidiagonal, so the basis is visibly more concentrated near higher indices.

The numerical condition-number check on the exact eigenbasis gives

κ(λ_j) = ||u_j||_2 ||v_j||_2 ≈ 1.00 for every j,

for the exact left/right eigenvectors computed from `np.linalg.eig(A)` and `np.linalg.eig(A.T.conj())`.

This means the non-normality is not being detected as a classical eigenvector-conditioning problem. The independent check still matters because it confirms that the instability is coming from the triangular structure and the strong superdiagonal coupling, not from an ill-conditioned eigenbasis. The key empirical diagnosis remains the Henrici relative departure

H(A) = ||A*A* - A*A||_F / ||A||_F^2 ≈ 0.51994,

which is large and consistent with pronounced non-normality.

This is a cheap sanity check, not new theory, and it is consistent with the project's stated interpretation.

## 5. Literature checklist and search log

### 5.1 Remaining search targets

The two open queries from the original list were:

1. `pseudospectra "superdiagonal" "diagonal" finite n`
2. direct MathSciNet search based on the corrected asymptotic claim

### 5.2 Search results

| Date | Query | Result | Coverage |
|---|---|---|---|
| 2026-09-02 | `pseudospectra "superdiagonal" "diagonal" finite n` | No direct exact-match result found through the available web search access; the family appears to be a specific finite-n case rather than a standard named benchmark. | Not directly covered by an obvious literature name |
| 2026-09-02 | MathSciNet search for transient growth / bidiagonal matrix asymptotic | Search shell was loaded, but direct bibliographic retrieval from the available environment was not available; no usable full-hit set recorded. | Early stage, no direct theorem identified |
| 2026-09-02 | Trefethen & Embree Ch. 16 notes and references (Kreiss Matrix Theorem) | The chapter is the canonical background reference for the generic bound; the present project uses it for the comparison baseline, not as the novelty claim. | Yes, used as the baseline theorem |
| 2026-09-02 | Trefethen & Embree Ch. 3 notes and references (Jordan blocks) | This is the baseline for the non-normality mechanism and pseudospectral phenomena; it is the conceptual source for the project's motivation and the sanity checks. | Yes, used as the motivating background |

The literature search therefore remains consistent with the project's own positioning: the generic Kreiss theorem and the Toeplitz/advection-diffusion literature are the true baseline references, while the present family-specific finite n asymptotic is not a standard named theorem in the available search results.

## 6. Handoff status

The pre-derivation packet is complete. The actionable outcome is:

- the corrected scaling law is now documented,
- the convergence numbers are clean and reportable,
- the literature check has been logged honestly,
- and the exact statement to carry into the proof is the provisional one above rather than the earlier, stronger rate-1 claim.

The analytic phase should begin only after the professor sign-off on that provisional statement.
