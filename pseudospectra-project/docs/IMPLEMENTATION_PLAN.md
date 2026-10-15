# Implementation Plan: Pseudospectra and Transient Growth of Non-Normal Matrices

**Author:** [your name] | **Mentor:** [your professor] | **Target:** GitHub portfolio + student-research journal paper (INVOLVE / SIURO)

## 1. The core idea, in one paragraph

A matrix `A` can have every eigenvalue safely in the left half-plane (so the ODE
`x' = Ax` is "spectrally stable" and `x(t) -> 0`) and yet `||exp(tA)||` can blow
up by orders of magnitude before it decays, if `A` is far from normal
(`A A* != A* A`). Eigenvalues alone cannot see this; **pseudospectra** can. The
project has two layers:
- **Computational**: build a correct, tested, reusable pseudospectra + transient-growth
  toolkit (already scaffolded for you below — see Section 3).
- **Theoretical**: for ONE specific structured matrix family, discover (numerically)
  and then PROVE a scaling law for how bad the transient growth gets as a
  structural parameter (coupling strength, Peclet number, matrix size) grows.
  That proof is your paper's theorem.

## 2. Required background (you already have most of this)

| Topic | Where it's used | You have it? |
|---|---|---|
| Linear algebra (SVD, eigenvalues, norms) | everything | Yes |
| ODEs (matrix exponential, stability) | transient growth | Yes |
| Complex analysis (contour integrals, residues) | proving sharp bounds | Ask your professor — 2-3 weeks of self-study if not yet covered |
| Basic numerical analysis (finite differences) | advection-diffusion discretization | Light lift given ODE background |

**Immediately read** (in this order), with your professor:
1. Trefethen & Embree, *Spectra and Pseudospectra*, Chapters 1–4 (the grid
   method, definitions, Kreiss Matrix Theorem). This is your bible — nearly
   every idea in this repo traces back to specific chapters here.
2. Trefethen, Trefethen, Reddy, Driscoll, "Hydrodynamic Stability Without
   Eigenvalues," *Science* 261 (1993) — the paper that made this field famous;
   short and very readable, gives you the "why this matters" framing for your
   introduction.
3. Reichel & Trefethen, "Eigenvalues and pseudo-eigenvalues of Toeplitz
   matrices," *Linear Algebra Appl.* 162–164 (1992) — directly relevant if you
   pursue the Toeplitz/banded family analytically (Section 5).

## 3. What's already built for you (this repository)

```
pseudospectra-project/
├── pseudospectra/
│   ├── core.py        # grid-based sigma_min(A - zI) computation (the algorithm)
│   ├── matrices.py     # structured non-normal matrix families
│   ├── bounds.py       # Kreiss constant + transient growth (the theory hooks)
│   └── plotting.py     # pseudospectra contour plots, growth curves
├── examples/demo.py     # end-to-end reproducible study, generates figures/
├── tests/test_core.py   # validated against KNOWN analytical results
├── figures/              # output plots (already generated once, re-run anytime)
└── docs/                 # this plan + a literature/results log template
```

Run it now:
```bash
pip install -r requirements.txt
python3 examples/demo.py         # generates all figures + prints numbers
python3 -m pytest tests/ -v      # confirm correctness
```

You should already see (from the included run):
- `04_scaling_with_c.png`: a clean log-linear trend, i.e. **peak growth ≈ C·e^(αc)**
  for the headline `shifted_superdiagonal` family, with the fitted values printed
  by the demo (α ≈ 0.86, C ≈ 0.054 at n=20 — **your job is to verify this is
  robust across n and then prove the exponential rate analytically**).

## 4. Milestones and timeline (aim for ~5-6 months)

### Phase 0 (Week 1-2): Orientation
- Read Trefethen & Embree Ch. 1-4 with your professor.
- Run the existing code, reproduce every number/figure yourself, understand
  every line of `core.py` and `bounds.py` (don't skip this — you need to be
  able to defend the numerics in a review).
- Deliverable: a one-page written summary of the Kreiss Matrix Theorem in
  your own words, with the proof sketch of the *lower* bound (it's short and
  a good warm-up).

### Phase 1 (Week 3-6): Validate and extend the toolkit
- Add 2-3 more structured matrix families beyond what's here (try
  `grcar_matrix`, and a Toeplitz family via `banded_toeplitz` with your own
  coefficients).
- Implement the Lui (1997) inverse-Lanczos speedup for larger matrices
  (optional but a nice "systems/algorithms" contribution — see the docstring
  in `core.py` for where it plugs in).
- Deliverable: expanded test suite; a short internal report comparing
  pseudospectra across your families, including departure-from-normality
  numbers for each (function already provided: `condition_number_departure_from_normality`).

### Phase 2 (Week 7-12): Numerical exploration -> conjecture
- Pick ONE family as your main object of study. **Recommendation: start with
  `shifted_superdiagonal`** (simple enough to hand-analyze) then, if time
  allows, extend results to `advection_diffusion_matrix` (physically
  meaningful, more literature to compare against).
- Systematically sweep the structural parameter (c, or Pe, or n) and fit
  scaling laws for: (a) the Kreiss constant K(A), (b) the true peak transient
  growth, (c) the time-to-peak.
- Write your conjecture as a precise mathematical statement, e.g.:
  > "For the family `A_n(c) = D_n + c*N` where `D_n = diag(-1,...,-n)` and `N`
  > is the nilpotent shift, `sup_t ||exp(tA_n(c))|| = Theta(exp(alpha * c))`
  > for some `alpha` depending on [...], uniformly in `n` for `n >= c`."
- Deliverable: a precise, falsifiable conjecture, plus numerical evidence
  (convergence checks: does the fitted exponent change as you increase grid
  resolution nx/ny, increase n_real/n_imag in `kreiss_constant`, or increase n?).

### Phase 3 (Week 13-20): Prove it
- This is where your professor's guidance matters most. Standard technique
  for bidiagonal/banded structured matrices: write down `(A - zI)^{-1}`
  explicitly (bidiagonal matrices have simple closed-form inverses — this is
  very tractable for `shifted_superdiagonal`), bound its norm directly using
  the explicit formula, rather than going through the generic Kreiss theorem.
- Try to prove a MATCHING lower bound too (construct an explicit vector that
  achieves near the growth you observe) — a two-sided (upper and lower) bound
  that matches up to a constant is a much stronger, more publishable result
  than an upper bound alone.
- Deliverable: a written proof (LaTeX), reviewed by your professor, of your
  main theorem — even if it's a special case (e.g., proven for the exact
  `shifted_superdiagonal` family, or for `n` above some threshold) that's fine.

### Phase 4 (Week 21-24): Write the paper + polish the repo
- Structure (standard for this genre):
  1. Introduction (why pseudospectra matter — cite Trefethen-Embree, the 1993
     Science paper)
  2. Preliminaries (definitions, Kreiss Matrix Theorem statement)
  3. The matrix family and numerical exploration (figures from your repo)
  4. Main theorem + proof
  5. Numerical validation of the theorem
  6. Discussion / future work (e.g., extending to advection-diffusion)
- Polish the GitHub repo: clean README with figures embedded, clear
  installation instructions, all tests passing, a `CITATION.cff` file.
- Submit to INVOLVE (they explicitly want faculty-mentored student papers —
  ideal fit) and post a preprint to arXiv (math.SP or math.NA) if your
  professor can endorse you.

## 5. Concrete "next experiment" ideas once you've mastered the basics

- Does the exponential rate `alpha` in the conjecture depend on `n`? Fix `c`,
  vary `n` from 5 to 100, refit, plot `alpha(n)` — does it converge?
- Compare `shifted_superdiagonal` against `banded_toeplitz(n, [-1,0,c], [-1,0,1])`
  (adds a subdiagonal too) — does adding symmetric coupling change the
  exponential rate or just the constant?
- For `advection_diffusion_matrix`: the reason growth was small at Pe=50 in
  the demo is the diffusion timescale dominates. Try much larger Pe (weaker
  relative diffusion) or Neumann boundary conditions, and see if you can find
  the parameter regime where transient growth becomes dramatic — this itself,
  done carefully, is worth reporting ("we characterize the (n, Pe) region
  where transient growth exceeds a factor of K").

## 6. Honest notes on scope

- Do not try to cover all of Trefethen-Embree or all matrix families. Depth
  on one family with a real theorem beats breadth with only numerics.
- Getting a Kreiss constant of `K(A)` accurate to 2 significant figures
  numerically is easy; PROVING a sharp scaling law is the actual hard,
  valuable part — budget most of your time there, not on more plots.
- If Phase 3 proves harder than expected and you only get a partial result
  (e.g. proven for a special case, or only an improved — not sharp — bound),
  that is still a legitimate paper. Journals like INVOLVE explicitly welcome
  well-executed partial results with honest discussion of what remains open.
