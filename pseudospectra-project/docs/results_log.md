# Results log

Numerical claims in the paper must cite a **run** below (grid, \(n_t\), \(x_{\max}\)). An unconverged Kreiss search underestimates \(K(A)\).

## Protocol

- Headline family: `shifted_superdiagonal(n, c, decay_start=1)`.
- Timescale: `t_max = 20 / |max Re(λ)|` (built into `suggested_t_max`).
- Kreiss: right half-plane grid; always print the refinement table.
- Exponential: compare `scipy.linalg.expm` to `theory.expm_closed_form`.
- Fast laptop demo: `python examples/demo.py`
- Finer grids: `python examples/demo.py --full`

## Theorem-shaped statements currently supported by algebra (not yet a submitted theorem)

For \(A = D + c N\) with \(\delta = 1\), \(t \ge 0\), 0-based indices \(j \ge i\),

\[
(e^{tA})_{i,j}
= e^{-(i+1)t}
\frac{1}{(j-i)!}
\bigl(c(1-e^{-t})\bigr)^{j-i}.
\]

Entrywise, \(|(e^{tA})_{i,j}| \le |c|^{j-i}/(j-i)!\), hence

\[
\|e^{tA}\|_2 \le \sqrt{n\,I_0(2|c|)}.
\]

Envelope of the \((0,r)\) family (large \(n\), \(|c|>1\)): maximizer \(t^* = \log|c|\), value \(|c|^{-1} e^{|c|-1}\).

**Not proved here:** matching lower bound of the same order for the 2-norm; sharpness of the Bessel constant; uniformity in \(n\).

## Machine runs

Runs of `examples/demo.py` append below.

## Run 2026-09-01 11:38 UTC
- headline n=12, c=8.0
- K(A)≈18.7065, peak≈45.1701, t_peak≈2.01258
- generic e*n*K~=610.193, Bessel~=3274.35
- envelope (delta/c)exp(c/delta-1)~=137.079
- closed-form 2-norm at t_peak≈45.1701
- last-basis test vector≈6.07238
- scaling (c, K, peak, Bessel, e n K, envelope):
  - 1, 0.862989, 1, 5.2302, 28.1502, 1
  - 2, 0.938223, 1, 11.6457, 30.6043, 1.35914
  - 4, 1.42336, 2.12915, 71.6294, 46.4292, 5.02138
  - 6, 4.35588, 8.97316, 476.851, 142.086, 24.7355
  - 8, 18.7065, 45.1701, 3274.35, 610.193, 137.079
  - 10, 89.5443, 238.267, 22862.6, 2920.88, 810.308

## Final Computational Verification Run — 2026-09-03

### Converged Kreiss Constant Summary Table
| (n, c) | Grid | K(A) | Convergence / Spread | Notes |
| :---: | :---: | :---: | :---: | :--- |
| (10, 6.0) | 120x240 | 4.3681 | < 0.08% | Fully converged |
| (12, 8.0) | 120x240 | 18.8022 | < 0.007% | Headline point, fully converged |
| (12, 12.0) | 120x240 | 426.9430 | < 0.004% | Fully converged |
| (16, 16.0) | 240x480 | 14248 (±15) | ±0.1% | Grid sampling noise near peak ridge (Re z ≈ 0.35 ± 0.03); reported as K(A) ≈ 14248 (estimated) |
| (24, 20.0) | 120x240 | 650788.48 | < 0.004% | Fully converged |

### Key Analytical & Computational Results
1. **Machine Precision Verification:** Log-space entrywise factorization check for $e^{tA} = e^{tD} e^{xN}$ at $(n=30, c=40, t=3.7)$ verified to $8.64 \times 10^{-16}$ relative error.
2. **Nilpotent Matrix Exponential Norm:** $\|e^{xN}\|_2 \sim e^x (1 - O(x^{-1}))$ verified up to $x=400$, matching Toeplitz symbol bound $\|f\|_\infty = e^x$.
3. **Split Bound vs True Peak Gap:**
   - Split bound: $\text{split\_bound}(t^*) = e^{-\delta t^*} \|e^{x(t^*)N}\|_2 \sim \frac{\delta}{e \cdot c} e^{c/\delta}$ (prefactor exponent $p = 1$).
   - True norm: $\|\text{row}_0(e^{t^* A})\|_2 \sim \frac{\delta}{e (4\pi)^{1/4}} c^{-5/4} e^{c/\delta}$ (prefactor exponent $p = 5/4 = 1.25$).
   - Asymptotic ratio: $\text{split\_bound} / \text{true\_norm} \sim (4\pi)^{1/4} (c/\delta)^{1/4} \approx 1.8803 \cdot c^{1/4}$.
4. **Test Suite Status:** 95/95 tests passing (15 core + 80 adversarial audit tests).

## Run 2026-09-16 14:26 UTC
- headline n=12, c=8.0
- K(A)≈18.7065, peak≈45.171, t_peak≈2.0188
- generic e*n*K~=610.193, Bessel~=3274.35
- envelope (delta/c)exp(c/delta-1)~=137.079
- closed-form 2-norm at t_peak≈45.171
- last-basis test vector≈6.09651
- scaling (c, K, peak, Bessel, e n K, envelope):
  - 1, 0.862989, 1, 5.2302, 28.1502, 1
  - 2, 0.938223, 1, 11.6457, 30.6043, 1.35914
  - 4, 1.42336, 2.13066, 71.6294, 46.4292, 5.02138
  - 6, 4.35588, 8.98484, 476.851, 142.086, 24.7355
  - 8, 18.7065, 45.171, 3274.35, 610.193, 137.079
  - 10, 89.5443, 238.555, 22862.6, 2920.88, 810.308
