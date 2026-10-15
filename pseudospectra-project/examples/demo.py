"""
End-to-end reproducible study of pseudospectra and transient growth.

Default (``python examples/demo.py``) uses laptop-scale grids. Paper-quality
resolution: ``python examples/demo.py --full``.

Figures land in ``figures/``. A numerical summary is appended to
``docs/results_log.md``.
"""

from __future__ import annotations

import argparse
import os
import sys
from datetime import datetime, timezone

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Windows consoles often use a legacy codepage (e.g. cp1256) that cannot
# encode the mathematical symbols used in the output; force UTF-8.
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")

import numpy as np
import matplotlib.pyplot as plt

from pseudospectra.core import (
    compute_pseudospectra_grid,
    eigenvalues,
    condition_number_departure_from_normality,
    suggested_t_max,
)
from pseudospectra.matrices import (
    jordan_block,
    advection_diffusion_matrix,
    shifted_superdiagonal,
    grcar_matrix,
)
from pseudospectra.bounds import (
    kreiss_constant,
    kreiss_constant_convergence,
    max_transient_growth,
    kreiss_sandwich_bounds,
    family_sandwich_bounds,
    scaling_sweep,
    fit_exponential_scaling,
    scaling_stability_across_n,
)
from pseudospectra.plotting import plot_pseudospectra, plot_transient_growth
from pseudospectra.theory import (
    bessel_uniform_bound,
    asymptotic_peak_envelope,
    expm_closed_form,
    test_vector_growth,
)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIGDIR = os.path.join(ROOT, "figures")
LOGPATH = os.path.join(ROOT, "docs", "results_log.md")
os.makedirs(FIGDIR, exist_ok=True)


def _cfg(full: bool) -> dict:
    if full:
        return dict(n=20, nx=90, ny=90, n_t=300, n_real=50, n_imag=90, n_adv=40)
    return dict(n=12, nx=50, ny=50, n_t=160, n_real=28, n_imag=56, n_adv=24)


def sanity_check_jordan(nx: int, ny: int) -> None:
    print("\n" + "=" * 70)
    print("SANITY CHECK: Jordan block pseudospectra")
    print("=" * 70)

    n = 6
    A = jordan_block(n, lam=0.0)
    X, Y, S = compute_pseudospectra_grid(A, (-1.5, 1.5), (-1.5, 1.5), nx=nx, ny=ny)
    eigs = eigenvalues(A)
    plot_pseudospectra(
        X, Y, S, eigs=eigs,
        title=rf"Jordan block (n={n}): $\varepsilon$-pseudospectra $\sim$ disk of radius $\varepsilon^{{1/{n}}}$",
    )
    plt.savefig(os.path.join(FIGDIR, "01_jordan_sanity_check.png"), dpi=150)
    plt.close()
    print("[OK] Saved figures/01_jordan_sanity_check.png")


def main_study_headline(n: int, c: float, cfg: dict) -> dict:
    print("\n" + "=" * 70)
    print(f"HEADLINE FAMILY: shifted superdiagonal (n={n}, c={c})")
    print("=" * 70)

    A = shifted_superdiagonal(n, c=c)
    eigs = eigenvalues(A)
    dep_abs, dep_rel = condition_number_departure_from_normality(A)
    print(f"  max Re(lambda) = {eigs.real.max():.4f}  (spectral stability iff < 0)")
    print(f"  Henrici relative departure = {dep_rel:.4f}")

    xr = (float(eigs.real.min()) - 2.0, 3.0)
    yr = (-8.0, 8.0)
    X, Y, S = compute_pseudospectra_grid(A, xr, yr, nx=cfg["nx"], ny=cfg["ny"])
    plot_pseudospectra(
        X, Y, S, eigs=eigs,
        title=f"Pseudospectra of A = D + cN  (n={n}, c={c})",
    )
    plt.savefig(os.path.join(FIGDIR, "02_toeplitz_pseudospectra.png"), dpi=150)
    plt.close()
    print("  [OK] figures/02_toeplitz_pseudospectra.png")

    t_max = suggested_t_max(A, multiple=20.0)
    ts, norms, t_peak, peak = max_transient_growth(A, t_max=t_max, n_t=cfg["n_t"])
    sandwich = family_sandwich_bounds(
        A, c=c, x_max=5.0, n_real=cfg["n_real"], n_imag=cfg["n_imag"]
    )
    K = sandwich["K"]
    print(f"  peak ||exp(tA)||_2 = {peak:.4f} at t = {t_peak:.4f}")
    print(f"  K(A) ~= {K:.4f}")
    print(f"  generic upper e*n*K ~= {sandwich['upper_generic']:.4f}")
    print(f"  family Bessel upper sqrt(n I0(2c)) ~= {sandwich['upper_bessel']:.4f}")
    t_star, env = asymptotic_peak_envelope(c)
    print(f"  envelope maximizer t* ~= {t_star:.4f},  (delta/c) exp(c/delta-1) ~= {env:.4f}")

    # Closed-form vs scipy at the sampled peak (consistency of the theorem substrate)
    E_closed = expm_closed_form(n, t_peak, c)
    closed_peak = float(np.linalg.norm(E_closed, ord=2))
    print(f"  closed-form ||exp(tA)||_2 at t_peak = {closed_peak:.4f}")
    tv = test_vector_growth(n, t_peak, c, mode="last_basis")
    print(f"  last-basis test-vector growth = {tv:.4f}  (lower bound)")

    plot_transient_growth(
        ts, norms, K=K, upper_bound=sandwich["upper_generic"],
        family_upper=sandwich["upper_bessel"],
        title=f"Transient growth vs Kreiss / Bessel sandwiches (n={n}, c={c})",
    )
    plt.savefig(os.path.join(FIGDIR, "03_toeplitz_transient_growth.png"), dpi=150)
    plt.close()
    print("  [OK] figures/03_toeplitz_transient_growth.png")

    print("  Kreiss grid-convergence check:")
    table = kreiss_constant_convergence(
        A,
        resolutions=[(16, 32), (28, 56), (40, 80)],
        x_maxes=[3.0, 5.0, 8.0],
    )
    for row in table:
        print(
            f"    n_real={row['n_real']:3d}  n_imag={row['n_imag']:3d}  "
            f"x_max={row['x_max']:.1f}  K≈{row['K']:.4f}"
        )

    return {
        "n": n, "c": c, "K": K, "peak": peak, "t_peak": t_peak,
        "generic_upper": sandwich["upper_generic"],
        "bessel_upper": sandwich["upper_bessel"],
        "envelope": env, "closed_peak": closed_peak, "test_vector": tv,
        "departure_rel": dep_rel, "kreiss_table": table,
    }


def scaling_study_headline(n: int, cfg: dict, c_values=None) -> list:
    print("\n" + "=" * 70)
    print(f"SCALING vs coupling c (n={n})")
    print("=" * 70)
    if c_values is None:
        c_values = (1, 2, 4, 6, 8, 10)

    sweep = scaling_sweep(
        lambda n_, c_: shifted_superdiagonal(n_, c=c_),
        c_values,
        n,
        n_real=cfg["n_real"],
        n_imag=cfg["n_imag"],
        n_t=cfg["n_t"],
        timescale_multiple=20.0,
    )
    results = []
    print("   c        K(A)     peak     Bessel     e n K     envelope")
    for row in sweep:
        c = row["parameter"]
        K, peak = row["K"], row["peak"]
        bessel = bessel_uniform_bound(n, c)
        generic = float(np.e * n * K)
        _, env = asymptotic_peak_envelope(c)
        results.append((c, K, peak, bessel, generic, env))
        print(f"  {c:5.1f}  {K:8.2f}  {peak:8.2f}  {bessel:8.2f}  {generic:8.2f}  {env:8.2f}")

    # Save the sweep table as CSV (Deliverable B6)
    csv_path = os.path.join(ROOT, "docs", f"scaling_sweep_n{n}.csv")
    with open(csv_path, "w", encoding="utf-8") as f:
        f.write("c,K,peak,bessel,generic_enK,envelope\n")
        for row in results:
            f.write(",".join(f"{v:.10g}" for v in row) + "\n")
    print(f"  [OK] saved {csv_path}")

    cs, Ks, peaks, bessels, generics, envs = zip(*results)
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.semilogy(cs, Ks, "o-", label="Kreiss constant K(A)")
    ax.semilogy(cs, peaks, "s-", label=r"sampled peak $\|e^{tA}\|_2$")
    ax.semilogy(cs, envs, "--", label=r"envelope $(\delta/c)\,e^{c/\delta-1}$")
    ax.semilogy(cs, bessels, "-.", label=r"Bessel bound $\sqrt{n I_0(2c/\delta)}$")
    ax.set_xlabel("coupling strength c")
    ax.set_ylabel("growth measure (log scale)")
    ax.set_title(f"Headline family: scaling with c (n={n})")
    ax.legend()
    ax.grid(True, which="both", alpha=0.3)
    plt.savefig(os.path.join(FIGDIR, "04_scaling_with_c.png"), dpi=150)
    plt.close()
    print("  [OK] figures/04_scaling_with_c.png")

    # Fit log(peak) ~ a + b c on c >= 4 if enough points (Deliverable B7)
    cs_arr = np.array(cs, dtype=float)
    pk_arr = np.array(peaks, dtype=float)
    mask = cs_arr >= 4
    fit = None
    if np.count_nonzero(mask) >= 2:
        sub = [{"parameter": c, "peak": pk} for c, pk in zip(cs_arr[mask], pk_arr[mask])]
        fit = fit_exponential_scaling(sub)
        print(
            f"  log-linear fit on c>=4: peak ~= {fit['C']:.4f} * exp({fit['alpha']:.4f} c)"
            f"  (R² = {fit['r_squared']:.4f})"
        )
        print("  compare to envelope rate 1/delta = 1.0 (leading exp(c/delta)).")
    return results, fit


def robustness_subdiagonal(n: int, c: float, cfg: dict, eps: float = 0.15) -> dict:
    """Does a small subdiagonal term destroy the headline scaling?"""
    print("\n" + "=" * 70)
    print(f"ROBUSTNESS: A = D + c N + eps P  (eps={eps})")
    print("=" * 70)
    A0 = shifted_superdiagonal(n, c=c)
    A1 = A0.copy()
    A1 += eps * np.diag(np.ones(n - 1), k=-1)

    out = {}
    for name, A in (("headline", A0), (f"eps={eps}", A1)):
        t_max = suggested_t_max(A, multiple=20.0)
        _, _, _, peak = max_transient_growth(A, t_max=t_max, n_t=cfg["n_t"])
        K, _ = kreiss_constant(A, x_max=5.0, n_real=cfg["n_real"], n_imag=cfg["n_imag"])
        out[name] = {"K": K, "peak": peak}
        print(f"  {name:12s}  K~={K:.3f}  peak~={peak:.3f}")
    return out


def grcar_illustration(n: int, nx: int, ny: int) -> None:
    print("\n" + "=" * 70)
    print(f"BENCHMARK: Grcar matrix (n={n})")
    print("=" * 70)
    A = grcar_matrix(n)
    eigs = eigenvalues(A)
    X, Y, S = compute_pseudospectra_grid(A, (-4, 6), (-4, 4), nx=nx, ny=ny)
    plot_pseudospectra(X, Y, S, eigs=eigs, title=f"Grcar matrix pseudospectra (n={n})")
    plt.savefig(os.path.join(FIGDIR, "08_grcar_pseudospectra.png"), dpi=150)
    plt.close()
    print("  [OK] figures/08_grcar_pseudospectra.png")


def main_study_advdiff(n: int, Pe: float, cfg: dict) -> dict:
    print("\n" + "=" * 70)
    print(f"SECONDARY FAMILY: advection-diffusion (n={n}, Pe={Pe})")
    print("=" * 70)
    print("  Not the paper's theorem target (TTRD 1993 already covers the PDE).")
    print("  Included to keep the hydrodynamic-stability motivation honest.")

    A = advection_diffusion_matrix(n, Pe=Pe)
    eigs = eigenvalues(A)
    print(f"  max Re(lambda) = {eigs.real.max():.4e}  (timescale ~ 1/|alpha| = {1/abs(eigs.real.max()):.4e})")
    _, dep_rel = condition_number_departure_from_normality(A)
    print(f"  Henrici relative departure = {dep_rel:.4f}")

    xr = (float(eigs.real.min()) - 5.0, 2.0)
    yr = (float(eigs.imag.min()) - 5.0, float(eigs.imag.max()) + 5.0)
    X, Y, S = compute_pseudospectra_grid(A, xr, yr, nx=cfg["nx"], ny=cfg["ny"])
    plot_pseudospectra(
        X, Y, S, eigs=eigs,
        title=f"Advection-diffusion pseudospectra (n={n}, Pe={Pe})",
    )
    plt.savefig(os.path.join(FIGDIR, "05_advdiff_pseudospectra.png"), dpi=150)
    plt.close()

    t_max = suggested_t_max(A, multiple=8.0)
    ts, norms, t_peak, peak = max_transient_growth(A, t_max=t_max, n_t=cfg["n_t"])
    K, lower, upper = kreiss_sandwich_bounds(
        A, x_max=max(20.0, Pe), n_real=cfg["n_real"], n_imag=cfg["n_imag"]
    )
    print(f"  peak = {peak:.4f} at t = {t_peak:.4e};  K~={K:.4f};  generic upper~={upper:.4f}")
    plot_transient_growth(
        ts, norms, K=K, upper_bound=upper,
        title=f"Transient growth: advection-diffusion (n={n}, Pe={Pe})",
    )
    plt.savefig(os.path.join(FIGDIR, "06_advdiff_transient_growth.png"), dpi=150)
    plt.close()
    print("  [OK] figures/05_advdiff_pseudospectra.png, 06_advdiff_transient_growth.png")
    return {"n": n, "Pe": Pe, "K": K, "peak": peak, "t_peak": t_peak}


def scaling_study_advdiff(n: int, cfg: dict, Pe_values=None) -> list:
    print("\n" + "=" * 70)
    print(f"SCALING vs Peclet number (n={n}) - exploratory only")
    print("=" * 70)
    if Pe_values is None:
        Pe_values = (10, 20, 40, 80)

    results = []
    for Pe in Pe_values:
        A = advection_diffusion_matrix(n, Pe=Pe)
        t_max = suggested_t_max(A, multiple=8.0)
        K, _ = kreiss_constant(A, x_max=max(20.0, Pe), n_real=cfg["n_real"], n_imag=cfg["n_imag"])
        _, _, _, peak = max_transient_growth(A, t_max=t_max, n_t=cfg["n_t"])
        results.append((Pe, K, peak))
        print(f"  Pe={Pe:6.1f}  K~={K:8.2f}  peak~={peak:8.2f}")

    Pes, Ks, peaks = zip(*results)
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.loglog(Pes, Ks, "o-", label="Kreiss constant K(A)")
    ax.loglog(Pes, peaks, "s-", label=r"sampled peak $\|e^{tA}\|_2$")
    ax.set_xlabel("Péclet number Pe")
    ax.set_ylabel("growth measure")
    ax.set_title(f"Advection-diffusion: scaling with Pe (n={n})")
    ax.legend()
    ax.grid(True, which="both", alpha=0.3)
    plt.savefig(os.path.join(FIGDIR, "07_scaling_with_Pe.png"), dpi=150)
    plt.close()
    print("  [OK] figures/07_scaling_with_Pe.png")
    return results


def _append_log(text: str) -> None:
    os.makedirs(os.path.dirname(LOGPATH), exist_ok=True)
    with open(LOGPATH, "a", encoding="utf-8") as f:
        f.write(text)


def write_run_log(headline: dict, scaling: list) -> None:
    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    lines = [
        f"\n## Run {stamp}\n",
        f"- headline n={headline['n']}, c={headline['c']}\n",
        f"- K(A)≈{headline['K']:.6g}, peak≈{headline['peak']:.6g}, t_peak≈{headline['t_peak']:.6g}\n",
        f"- generic e*n*K~={headline['generic_upper']:.6g}, Bessel~={headline['bessel_upper']:.6g}\n",
        f"- envelope (delta/c)exp(c/delta-1)~={headline['envelope']:.6g}\n",
        f"- closed-form 2-norm at t_peak≈{headline['closed_peak']:.6g}\n",
        f"- last-basis test vector≈{headline['test_vector']:.6g}\n",
        "- scaling (c, K, peak, Bessel, e n K, envelope):\n",
    ]
    for row in scaling:
        lines.append("  - " + ", ".join(f"{v:.6g}" for v in row) + "\n")
    _append_log("".join(lines))
    print(f"  [OK] appended {LOGPATH}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Pseudospectra demo")
    parser.add_argument("--full", action="store_true", help="larger n and finer grids")
    parser.add_argument("--skip-advdiff", action="store_true")
    args = parser.parse_args()
    cfg = _cfg(args.full)

    print("=" * 70)
    print("PSEUDOSPECTRA & TRANSIENT GROWTH")
    print("mode:", "full" if args.full else "fast (laptop)")
    print("=" * 70)

    sanity_check_jordan(cfg["nx"], cfg["ny"])
    headline = main_study_headline(n=cfg["n"], c=8.0, cfg=cfg)
    scaling, fit = scaling_study_headline(n=cfg["n"], cfg=cfg)
    robustness_subdiagonal(n=cfg["n"], c=8.0, cfg=cfg)
    grcar_illustration(n=min(cfg["n"], 16), nx=cfg["nx"], ny=cfg["ny"])

    # Deliverable B7: check whether the fitted rate alpha is stable in n
    print("\n" + "=" * 70)
    print("SCALING STABILITY ACROSS n (Algorithm 5, step 12)")
    print("=" * 70)
    n_values = sorted({max(8, cfg["n"] - 4), cfg["n"], cfg["n"] + 4})
    stability = scaling_stability_across_n(
        lambda n_, c_: shifted_superdiagonal(n_, c=c_),
        (4, 6, 8, 10),
        n_values,
        n_real=cfg["n_real"],
        n_imag=cfg["n_imag"],
        n_t=cfg["n_t"],
        timescale_multiple=20.0,
        verbose=False,
    )
    for n_ in n_values:
        f_ = stability[n_]["fit"]
        print(
            f"  n={n_:3d}: alpha={f_['alpha']:.4f}  C={f_['C']:.4f}  "
            f"R²={f_['r_squared']:.4f}"
        )
    if fit is not None:
        alphas = [stability[n_]["fit"]["alpha"] for n_ in n_values]
        spread = max(alphas) - min(alphas)
        print(
            f"  alpha spread across n = {spread:.4f} "
            f"({'stable' if spread < 0.1 else 'n-dependent — report explicitly'})"
        )

    if not args.skip_advdiff:
        main_study_advdiff(n=cfg["n_adv"], Pe=50.0, cfg=cfg)
        scaling_study_advdiff(n=cfg["n_adv"], cfg=cfg)

    write_run_log(headline, scaling)
    print("\n[OK] Demo complete. See figures/ and docs/results_log.md.")


if __name__ == "__main__":
    main()
