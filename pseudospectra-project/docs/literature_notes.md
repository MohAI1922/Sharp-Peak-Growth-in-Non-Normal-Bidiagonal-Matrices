# Literature notes (working file)

This is a **search and reading log**, not a claim that every paper below has been fully read. Fill in the “what we checked” column with your professor before locking the family for the theorem.

## What is not novel (do not present as a contribution)

- Transient growth of non-normal matrices in general: Trefethen & Embree (2005); hydrodynamic-stability literature earlier.
- The generic Kreiss sandwich \(K(A) \le \sup_t \|e^{tA}\| \le e\,n\,K(A)\): a classical theorem. It is a **tool**.
- Pseudospectra of advection–diffusion: Trefethen–Trefethen–Reddy–Driscoll, *Science* 261 (1993), plus a large follow-up literature. Using this family in figures is fine for motivation; a theorem about “the” AD operator is not a new paper unless the question is genuinely new (BCs, scheme, regime).
- Toeplitz / banded pseudospectra in the \(n\to\infty\) symbol limit: Reichel & Trefethen, *LAA* 162–164 (1992) and later symbol-calculus work.

## Intended contribution (path (a), possibly (b))

For the **finite** bidiagonal family \(A = D + c N\) with \(D = \mathrm{diag}(-\delta,-2\delta,\ldots,-n\delta)\) (not a Toeplitz matrix: the diagonal is an arithmetic progression), prove a **family-specific** bound on \(K(A)\) and/or \(\sup_t \|e^{tA}\|\) with **explicit constants**, improving on \(e\,n\,K(A)\).

The code now supplies the algebraic substrate:

- exact \((A-zI)^{-1}\) (upper-bidiagonal inversion formula);
- exact \((e^{tA})_{i,j} = e^{-(i+1)\delta t}\,[(c/\delta)(1-e^{-\delta t})]^{j-i}/(j-i)!\)
- a rigorous uniform bound \(\|e^{tA}\|_2 \le \sqrt{n\,I_0(2|c|/\delta)}\);
- an envelope maximizer suggesting rate \(e^{|c|/\delta}\) (to be turned into a theorem, including a matching lower bound).

**Still to prove (with the mentor):** (i) that the Bessel / envelope bounds are stated under the right hypotheses and with sharp constants; (ii) a matching test-vector lower bound of the same exponential order; (iii) that this specific \(D+cN\) family is not already treated under a more general theorem you missed.

## Core references

| Item | Why it matters | Status |
|---|---|---|
| Trefethen & Embree, *Spectra and Pseudospectra*, Princeton, 2005, Ch. 1–4, 14–16 | Definitions, grid method, Kreiss theorem, notes/references (open problems) | Read Ch. 1–4 with mentor |
| Trefethen, Trefethen, Reddy, Driscoll, *Science* 261 (1993) | Motivation; AD is **not** the novelty target | Cite in introduction |
| Reichel & Trefethen, *LAA* 1992 | Toeplitz symbol asymptotics — **contrast**: our \(D\) is not Toeplitz | Check that they do not already cover \(D+cN\) |
| Lui, *SIAM J. Sci. Comput.* 1997 | Inverse-Lanczos for \(\sigma_{\min}(A-zI)\); optional Phase-2 computational contribution | Not implemented (Schur preprocess only) |
| Henrici, departure from normality | Scalar diagnostic used in `core.condition_number_departure_from_normality` | Formula implemented |

## Search queries to run (MathSciNet / Scholar / arXiv) before freezing the theorem

1. `pseudospectra bidiagonal` / `Kreiss constant bidiagonal`
2. `transient growth arithmetic diagonal` / `nonnormal bidiagonal exponential`
3. `pseudospectra "superdiagonal" "diagonal"` finite \(n\)
4. Look up the **Notes and References** at the end of Trefethen & Embree chapters on Toeplitz matrices, the Kreiss matrix theorem, and transient behavior.

Record hits here:

| Date | Query | Hit | Already covers our family? |
|---|---|---|---|
| 2026-09-02 | `pseudospectra "superdiagonal" "diagonal" finite n` | No direct exact-match result found in the available search environment; this appears to be a finite-n family-specific problem rather than a standard named benchmark. | No |
| 2026-09-02 | `transient growth bidiagonal matrix finite n` | Search was attempted but institutional access blocked direct MathSciNet retrieval; no usable full-hit set was recorded. This remains an open literature task rather than a completed search. | Open / blocked |
| 2026-09-02 | Trefethen & Embree Ch. 16 (Kreiss theorem) notes + references | Chapter provides the baseline theorem and the relevant references; used as the comparison theorem, not as novelty. | Yes, as baseline |
| 2026-09-02 | Trefethen & Embree Ch. 3 (Jordan blocks) notes + references | Chapter is the canonical reference for the mechanism by which non-normality creates transient growth. | Yes, as background |

## Honest limitation (keep this in the paper’s discussion)

The Bessel bound still carries a \(\sqrt{n}\) factor and uses \(\|\,\cdot\,\|_F\), so it is **not** automatically a sharp 2-norm characterization. The envelope \((\delta/c)\,e^{c/\delta-1}\) is a **leading-term** calculation for a single entry family, not yet a complete two-sided theorem. Advection–diffusion remains a motivation example, not the main result.
