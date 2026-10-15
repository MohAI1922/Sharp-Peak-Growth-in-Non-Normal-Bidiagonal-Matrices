# Proof Readiness Specification — Analytical Phase (Part II)

**Date:** September 1, 2026
**Purpose:** Everything verified and available to begin the hand derivation
(Part II of the blueprint). Numerical claims below cite the machine runs in
`docs/results_log.md` and the CSV tables `docs/scaling_sweep_n*.csv`.

---

## 1. Target family (fully specified)

For `n ≥ 2`, coupling `c ∈ ℝ`, decay `δ > 0` (default δ = 1):

```
A = D + c N,   D = diag(−δ, −2δ, …, −nδ),   N_{i,i+1} = 1 (superdiagonal)
```

Code: `pseudospectra.matrices.shifted_superdiagonal(n, c, decay_start)`.

**Verified properties** (tests 5, 8, 9 in `tests/test_core.py`):
- Upper bidiagonal ⇒ eigenvalues are exactly `−δ, −2δ, …, −nδ`
  (all real, negative ⇒ spectrally stable for every `c`).
- Henrici relative departure at n=12, c=8: **0.5199** (strongly non-normal).
- Non-normality is the *only* growth mechanism — the spectrum is fixed
  independent of `c`. This is the family's key isolating property.

---

## 2. Exact algebraic substrate (already derived and machine-verified)

### 2.1 Resolvent — closed form (Part II, Step 1: DONE)

Let `B = A − zI`, diagonal `a_k = −kδ − z` (1-based k), superdiagonal `c`.
`B^{-1}` is upper triangular with, for `j ≥ i` (0-based),

```
(B^{-1})_{i,j} = (−c)^{j−i} / ∏_{k=i}^{j} a_k ,      (B^{-1})_{i,j} = 0 for j < i
```

- Code: `theory.bidiagonal_resolvent` — **verified against `np.linalg.inv`**
  (test 8, n=7, c=3.5, z=0.4+1.2j, agreement 1e-10).
- Consequence: `σ_min(A − zI) = 1/‖(A − zI)^{-1}‖₂` exactly
  (`theory.resolvent_2_norm_exact`, verified in test 8).

### 2.2 Matrix exponential — closed form (Part II, Step 4 substrate: DONE)

For `t ≥ 0`, `r = j − i ≥ 0` (0-based):

```
(e^{tA})_{i,j} = e^{−(i+1)δt} · [ (c/δ)(1 − e^{−δt}) ]^r / r!
```

- Code: `theory.expm_closed_form` — **verified against `scipy.linalg.expm`**
  (test 9, n=8, c=5, t=0.35, agreement 1e-9).
- Entrywise envelope: `|(e^{tA})_{i,j}| ≤ (|c|/δ)^r / r!` since
  `0 ≤ e^{−(i+1)δt}(1 − e^{−δt})^r ≤ 1`.

### 2.3 Rigorous t-uniform upper bound (Step 4 upper bound: DONE)

```
‖e^{tA}‖₂ ≤ ‖e^{tA}‖_F ≤ √( n · I₀(2|c|/δ) )
```

where `I₀` is the modified Bessel function (uses
`∑_{r≥0} x^{2r}/(r!)² = I₀(2x)`, ≤ n superdiagonals × ≤ n entries).
- Code: `theory.bessel_uniform_bound` — verified numerically above sampled
  growth (test 10: peak 8.97 ≤ 435.3 at n=10, c=6).

### 2.4 Asymptotic envelope (conjecture-shaping: computed, not proved sharp)

Maximizing the (0, r) envelope `exp(−δt + (c/δ)(1 − e^{−δt}))` over `t > 0`
gives, for `|c| > δ`:

```
t* = δ^{-1} log(|c|/δ),      envelope = (δ/|c|) · e^{|c|/δ − 1}
```

- Code: `theory.asymptotic_peak_envelope`.
- At n=12, c=8: envelope predicts 137.1 with t* = 2.079; the sampled peak is
  **45.17 at t = 2.013** — t* is accurate, the envelope overestimates the
  2-norm (it is an entrywise/Frobenius-flavored bound). **This gap is exactly
  the 2-norm-vs-envelope question the proof must close.**

### 2.5 Lower-bound candidates (Step 5 substrate: implemented, not sharp)

- `theory.test_vector_growth(mode="last_basis")`: growth of `e^{tA} e_n`
  (last column). At n=12, c=8, t=t_peak: **6.07** — correct order of the
  exponential rate but far below 45.17; a better vector (or column-norm
  argument) is needed for a matching lower bound.
- `mode="poisson"`: basis vector concentrated near the Stirling peak index
  `⌊(c/δ)(1 − e^{−δt})⌋`.

---

## 3. Numerical evidence the theorem must explain

From the demo run (n=12, `docs/scaling_sweep_n12.csv`):

| c | K(A) | peak ‖e^{tA}‖₂ | Bessel bound | e·n·K | envelope |
|---|------|----------------|--------------|-------|----------|
| 1 | 0.86 | 1.00 | 5.2 | 28.2 | 1.00 |
| 2 | 0.94 | 1.00 | 11.7 | 30.6 | 1.36 |
| 4 | 1.42 | 2.13 | 71.6 | 46.4 | 5.02 |
| 6 | 4.36 | 8.97 | 476.9 | 142.1 | 24.7 |
| 8 | 18.71 | 45.17 | 3274 | 610.2 | 137.1 |
| 10 | 89.54 | 238.3 | 22863 | 2920.9 | 810.3 |

**Fitted scaling law** (log-linear on c ≥ 4):
`peak ≈ 0.085 · e^{0.79 c}`, R² = 0.9989.

**Stability in n** (Algorithm 5 step 12, demo output):

| n | α | C | R² |
|---|------|--------|--------|
| 8 | 0.702 | 0.130 | 0.9996 |
| 12 | 0.789 | 0.085 | 0.9989 |
| 16 | 0.796 | 0.082 | 0.9985 |

α spread = 0.094 → the exponential rate in `c` is **n-stable** (within the
0.1 tolerance), supporting an n-uniform conjecture.

**Kreiss convergence** (n=12, c=8): K = 18.53 → 18.71 → 18.56 across
(16×32, 28×56, 40×80) grids — converged to ~1%; the headline K is not an
under-resolved artifact.

**Robustness** (n=12, c=8, subdiagonal perturbation ε=0.15):
K: 18.7 → 51.2, peak: 45.2 → 91.7 — growth persists under perturbation
(the phenomenon is not an artifact of exact bidiagonality).

---

## 4. The working conjecture (to be proved or scoped down)

> For the family `A = D + cN` with `D = diag(−δ, …, −nδ)`, the transient
> growth satisfies, uniformly in `n` and `t ≥ 0`,
> `‖e^{tA}‖₂ ≤ √(n I₀(2|c|/δ))` (proved), while the empirical peak growth is
> consistent with a rate-one exponential in `c/δ` but with an unresolved
> polynomial prefactor; the provisional form is `peak ≈ C c^{-p} e^{c/δ}` with
> `α ≈ 1` and `p ≈ 1.5–1.6` over the tested range. This is a guarded
> statement pending the proof of the precise asymptotic prefactor.

The open analytical gaps, in order of attack:

1. **(Step 2–3) Kreiss constant bound.** From the closed resolvent (§2.1),
   bound `‖(A−zI)^{-1}‖` (Frobenius dominates 2-norm) and optimize over
   `Re z > 0` to get a closed-form/asymptotic `K(A)`. Check numerically
   whether the maximizer is real (numerics: `z_opt` from
   `bounds.kreiss_constant`).
2. **(Step 4) Direct 2-norm bound.** Tighten §2.3: replace the √n·Frobenius
   step with a genuine 2-norm estimate of the bidiagonal `e^{tA}` (e.g.
   Schur/test-matrix techniques, or row/column norm balancing).
3. **(Step 5) Matching lower bound.** Construct `v` (last basis vector is
   the natural candidate but underperforms 6.07 vs 45.17; try the
   Poisson-weighted vector or a full column of `e^{tA}`) achieving
   `‖e^{tA}v‖ ≥ c' · e^{|c|/δ}` for a well-chosen `t`.
4. **(Step 6) Assemble** the two-sided theorem with explicit constants.

### Corrected pre-derivation scaling verdict (2026-09-02)

The original pure-exponential fit `log peak ≈ α c + β` was contaminated by the
missing log-c term. Re-running the fit in the corrected form

`log peak ≈ α c − p log c + β`

on the extended range gives the following values at the tested matrix dimensions:

- `n = 8`: `(α, p, β, R²) = (0.8941, 1.3081, -0.9008, 0.9985)`
- `n = 12`: `(α, p, β, R²) = (1.0019, 1.5726, -0.9829, 0.9995)`
- `n = 16`: `(α, p, β, R²) = (1.0060, 1.5626, -0.9983, 0.9997)`
- `n = 24`: `(α, p, β, R²) = (1.0275, 1.6372, -1.0042, 0.9999)`

The significant new fact is that the row-0 norm is an exact lower bound on the true
matrix norm, not a mere heuristic. For any matrix `M` and any row index `i`,

`||row_i(M)||_2 <= ||M||_2`,

so for the headline family,

`||row_0(e^{tA})||_2 <= ||e^{tA}||_2` for every `t >= 0`.

The exact finite-n lower bound is

`||row_0(e^{tA})||_2 = e^{-δt} √( Σ_{r=0}^{n-1} [ (c/δ)(1-e^{-δt}) ]^{2r} / (r!)^2 )`.

At the envelope-maximizing time `t* = δ^{-1} log(c/δ)`, the row-0 norm has the
large-`x` asymptotic

`||row_0||_2 ~ const · c^{-5/4} e^{c/δ}`

because `I_0(2x) ~ e^{2x}/√(4πx)`, and this gives the leading prefactor
exponent `p = 5/4`. In the untruncated regime `n >= 3c`, the actual peaks are
within ~1–4% of this exact lower bound;

- `c = 4`: actual 2.1307 vs row-0 2.0499 (ratio 1.039)
- `c = 8`: actual 45.3604 vs row-0 44.9685 (ratio 1.009)
- `c = 12`: actual 1464.91 vs row-0 1459.39 (ratio 1.004)
- `c = 20`: actual 2.277e6 vs row-0 2.274e6 (ratio 1.001)

The local finite-difference exponents on the uniform safety-margin sweep `c = 10, 12, ..., 40` are

- `(10,12)`: 1.2882
- `(12,14)`: 1.2805
- `(14,16)`: 1.2752
- `(16,18)`: 1.2717
- `(18,20)`: 1.2688
- `(20,22)`: 1.2666
- `(22,24)`: 1.2649
- `(24,26)`: 1.2636
- `(26,28)`: 1.2624
- `(28,30)`: 1.2615
- `(30,32)`: 1.2607
- `(32,34)`: 1.2599
- `(34,36)`: 1.2593
- `(36,38)`: 1.2588
- `(38,40)`: 1.2580

These trend toward `5/4 = 1.25` and the overlapping-triple Aitken extrapolation on this even grid stays in the band `1.252–1.265`, with a final cluster around `1.254–1.256`. The appropriate summary is therefore `p ≈ 1.25 ± 0.02`, not a single hard-coded exponent.

The professor-facing conjecture is therefore:

> For the family `A = D + cN` with `D = diag(-δ, …, -nδ)`, the exponential
> rate in `c/δ` is 1, and the exact leading-row lower bound is already rigorous.
> Numerically, the dominant row mechanism predicts a prefactor exponent near
> `5/4`, so the asymptotic law is strongly consistent with `c^{-5/4} e^{c/δ}`
> rather than the earlier crude `c^{-1}` envelope guess. The exact prefactor
> exponent is not yet nailed down to a formal confidence interval, and the matching
> upper bound is now a concrete weighted-Schur problem rather than an open-ended
> search.

This is the honest, professor-facing version: the lower bound is solid, the rate is solid, and the prefactor exponent is a strong numerical candidate rather than a theorem-level lock.

### Professor sign-off draft (next step before analytic proof)

> For the family `A = D + cN` with `D = diag(-δ, ..., -nδ)`, the transient growth satisfies `sup_{t≥0} ||e^{tA}||_2 = exp(c/δ + O(log c))` in the sense that the exponential rate in `c/δ` is exactly one. A rigorous lower bound is already available from the leading row, `||row_0(e^{tA})||_2 <= ||e^{tA}||_2`, and the exact finite-n formula for this quantity is consistent numerically with the true peak to within a few percent in the untruncated regime `n >= 3c`. The uniform safety-margin sweep on `c = 10, 12, ..., 40` produces local exponents trending toward `5/4`, with Aitken extrapolation clustering around `1.25 ± 0.02`; this strongly suggests the prefactor is `c^{-5/4}` rather than the earlier crude `c^{-1}` envelope guess, but the exact exponent should still be treated as a strong numerical candidate pending the matching upper-bound proof. The analytic task is therefore not to discover whether the rate is one, but to prove the sharp leading-order law and its matching upper bound using a weighted Schur or row-balanced argument.

This is the statement to send for explicit sign-off before starting the hand derivation; it is intentionally more guarded than a theorem but is materially sharper than the earlier over-claimed version.

---

## 5. Toolkit status (Deliverable Set A/B: COMPLETE)

| Item | Status |
|---|---|
| `core.py`: Algorithms 1, 2, 6 + eigenvalues/spectral_abscissa/suggested_t_max | ✅ tested |
| `matrices.py`: jordan, shifted_superdiagonal, advection_diffusion, banded_toeplitz, grcar, random | ✅ |
| `bounds.py`: Algorithm 3 (+convergence table, `kreiss_converged` flag), Algorithm 4, sandwich, **Algorithm 5** (`scaling_sweep`, `fit_exponential_scaling` with R², `scaling_stability_across_n`) | ✅ tested |
| `theory.py`: closed resolvent, closed expm, Bessel bound, envelope, test vectors | ✅ tested vs scipy/inv |
| Test suite: 14/14 pass (incl. blueprint Test Set A items 1–4) | ✅ |
| Baseline figures (Jordan + headline pseudospectra & growth curves) | ✅ `figures/01–03` |
| Scaling sweep table + CSV + log-scale plot, ≥6 parameter values | ✅ `docs/scaling_sweep_n12.csv`, `figures/04` |
| Fitted law with R² + multi-n stability check | ✅ (§3) |
| Kreiss convergence table | ✅ (§3, printed per run) |
| Reproducibility | `python examples/demo.py` regenerates everything; `--full` for paper quality |

## 6. What remains (Deliverable Set C/D — the actual research)

- C9: handwritten resolvent derivation (substrate §2.1 done numerically;
  write it up and have the professor check it).
- C10: closed-form bound on K(A) with explicit constants (§4 item 1).
- C11: two-sided bound on `sup_t ‖e^{tA}‖` (§4 items 2–3).
- C12: final theorem + gap-free proof (§4 item 4), professor-reviewed.
- D13–15: paper draft, repo polish, submission (INVOLVE / SIURO / arXiv).

## 7. References

- Trefethen & Embree, *Spectra and Pseudospectra*, Princeton, 2005
  (esp. Ch. 3 Jordan blocks, Ch. 4 computation, Ch. 16 Kreiss theorem;
  read each chapter's "Notes and References" for open problems).
- Reichel & Trefethen, LAA 162–164 (1992) — Toeplitz pseudospectra.
- Trefethen, Trefethen, Reddy, Driscoll, Science 261 (1993) — motivation.
- Literature-search checklist (blueprint Part III §3) must be completed
  before committing: search "pseudospectra bidiagonal Toeplitz",
  "Kreiss constant banded matrix", and check whether an n-uniform version
  of the target bound already exists.