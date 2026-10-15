"""
Kreiss constant, transient growth, and sandwich bounds.

Generic tool (not a contribution)
=================================
    K(A) ≤ sup_{t≥0} ‖e^{tA}‖_2 ≤ e · n · K(A)

with K(A) = sup_{Re(z)>0} Re(z) / σ_min(A − z I). This is the Kreiss
Matrix Theorem. Numerics here *use* it; they do not reprove it.

Family-specific bound (headline family A = D + c N)
====================================================
See ``pseudospectra.theory``: ‖e^{tA}‖_2 ≤ √(n I_0(2 |c|/δ)) uniformly
in t. Compare that against e·n·K(A) in the demo scaling study.
"""

from __future__ import annotations

from typing import List, Optional, Tuple

import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize_scalar

from .core import min_singular_value, spectral_abscissa, suggested_t_max, eigenvalues
from .theory import bessel_uniform_bound


def scaling_sweep(
    family_generator,
    parameter_values,
    n: int,
    n_real: int = 50,
    n_imag: int = 80,
    n_t: int = 300,
    timescale_multiple: float = 15.0,
    x_max: float = 5.0,
    verbose: bool = True,
) -> List[dict]:
    """
    Algorithm 5 (blueprint): sweep a structural parameter of a matrix family.

    For each parameter value p: build A = family_generator(n, p), verify
    spectral stability, set t_max = 15/|α(A)|, estimate K(A) and the peak
    transient growth, and record (p, K, peak).

    Returns a list of dicts with keys: parameter, K, peak, t_peak,
    spectral_abscissa, t_max.
    """
    results: List[dict] = []
    for p in parameter_values:
        A = family_generator(n, p)
        alpha = spectral_abscissa(A)
        if alpha >= 0:
            raise ValueError(
                f"parameter {p}: matrix is not spectrally stable (α={alpha})"
            )
        t_max = suggested_t_max(A, multiple=timescale_multiple)
        K, _ = kreiss_constant(A, x_max=x_max, n_real=n_real, n_imag=n_imag)
        _, _, t_peak, peak = max_transient_growth(A, t_max=t_max, n_t=n_t)
        results.append(
            {
                "parameter": float(p),
                "K": float(K),
                "peak": float(peak),
                "t_peak": float(t_peak),
                "spectral_abscissa": float(alpha),
                "t_max": float(t_max),
            }
        )
        if verbose:
            print(
                f"  p={p:8.3g}  K≈{K:10.4g}  peak≈{peak:10.4g}  "
                f"t_peak≈{t_peak:.4g}  (t_max={t_max:.4g})"
            )
    return results


def fit_exponential_scaling(results: List[dict]) -> dict:
    """
    Least-squares fit log(peak) ≈ α·p + β (blueprint Algorithm 5, step 11).

    Returns dict with alpha, beta, C = exp(beta), R² of the log-linear
    regression, and the number of points used.
    """
    p = np.array([r["parameter"] for r in results], dtype=float)
    log_peak = np.log(np.array([r["peak"] for r in results], dtype=float))
    alpha, beta = np.polyfit(p, log_peak, 1)
    fitted = alpha * p + beta
    ss_res = float(np.sum((log_peak - fitted) ** 2))
    ss_tot = float(np.sum((log_peak - log_peak.mean()) ** 2))
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return {
        "alpha": float(alpha),
        "beta": float(beta),
        "C": float(np.exp(beta)),
        "r_squared": r_squared,
        "n_points": int(len(p)),
    }


def scaling_stability_across_n(
    family_generator,
    parameter_values,
    n_values,
    **sweep_kwargs,
) -> dict:
    """
    Blueprint Algorithm 5, step 12: repeat the sweep at several n and
    report whether the fitted rate α depends on n.

    Returns {n: {"results": [...], "fit": {...}}}.
    """
    out = {}
    for n in n_values:
        results = scaling_sweep(family_generator, parameter_values, n, **sweep_kwargs)
        out[n] = {"results": results, "fit": fit_exponential_scaling(results)}
    return out


def kreiss_converged(
    A: np.ndarray,
    rtol: float = 0.01,
    resolutions: Optional[List[Tuple[int, int]]] = None,
    x_maxes: Optional[List[float]] = None,
    y_range: Optional[Tuple[float, float]] = None,
) -> Tuple[float, bool, np.ndarray]:
    """
    Run the Kreiss convergence check and declare convergence.

    Returns (K_finest, converged_flag, table) where converged_flag is True
    when the two finest-grid estimates agree within ``rtol`` (relative).
    """
    table = kreiss_constant_convergence(A, resolutions, x_maxes, y_range)
    K_finest = float(table["K"][-1])
    if len(table) < 2:
        return K_finest, False, table
    prev = float(table["K"][-2])
    converged = bool(abs(K_finest - prev) <= rtol * max(prev, 1e-300))
    return K_finest, converged, table


def kreiss_constant(
    A: np.ndarray,
    x_max: float = 5.0,
    n_real: int = 60,
    n_imag: int = 120,
    y_range: Optional[Tuple[float, float]] = None,
) -> Tuple[float, complex]:
    """
    Numerically estimate the Kreiss constant by a right half-plane grid search.

    Returns
    -------
    K : float
        Estimated K(A). An under-resolved grid *underestimates* K(A).
    z_opt : complex
        Grid point attaining the recorded maximum (nan+nan j if none).
    """
    if y_range is None:
        eigs = np.linalg.eigvals(A)
        y_extent = max(1.0, float(np.max(np.abs(eigs.imag))) * 1.5)
        y_range = (-y_extent, y_extent)

    xs = np.linspace(1e-3, x_max, n_real)
    ys = np.linspace(*y_range, n_imag)

    best = 0.0
    best_z: complex = complex("nan")
    for x in xs:
        for y in ys:
            z = x + 1j * y
            smin = min_singular_value(A, z)
            if smin < 1e-14:
                continue
            val = x / smin
            if val > best:
                best = val
                best_z = z

    return float(best), best_z


def kreiss_constant_convergence(
    A: np.ndarray,
    resolutions: Optional[List[Tuple[int, int]]] = None,
    x_maxes: Optional[List[float]] = None,
    y_range: Optional[Tuple[float, float]] = None,
) -> np.ndarray:
    """
    Recompute K(A) on a sequence of refined grids.

    Returns a structured array with fields n_real, n_imag, x_max, K, Re(z), Im(z).
    Report this table in the paper: if K is still climbing at the finest
    grid, the published number is an underestimate.
    """
    if resolutions is None:
        resolutions = [(20, 40), (40, 80), (80, 160)]
    if x_maxes is None:
        x_maxes = [2.0, 5.0, 10.0]

    rows = []
    for (nr, ni), xmax in zip(resolutions, x_maxes):
        K, z = kreiss_constant(A, x_max=xmax, n_real=nr, n_imag=ni, y_range=y_range)
        rows.append((nr, ni, xmax, K, float(np.real(z)), float(np.imag(z))))
    dtype = np.dtype(
        [
            ("n_real", int),
            ("n_imag", int),
            ("x_max", float),
            ("K", float),
            ("z_real", float),
            ("z_imag", float),
        ]
    )
    return np.array(rows, dtype=dtype)


def find_exact_peak_norm(
    A: np.ndarray,
    t_star: float,
    t_window: float = 0.5,
) -> Tuple[float, float]:
    """
    Finds the exact continuous maximum of ‖exp(tA)‖_2 near t_star
    using bounded scalar minimization, preventing grid-sampling errors
    that could violate theoretical lower bounds.
    """
    def neg_norm(t: float) -> float:
        if t < 0:
            return 0.0
        E = expm(t * A)
        return float(-np.linalg.norm(E, ord=2))

    t_min = max(0.0, float(t_star - t_window))
    t_max_bracket = float(t_star + t_window)
    res = minimize_scalar(neg_norm, bounds=(t_min, t_max_bracket), method='bounded')
    t_peak = float(res.x)
    peak_val = float(-res.fun)
    return t_peak, peak_val


def max_transient_growth(
    A: np.ndarray,
    t_max: Optional[float] = None,
    n_t: int = 200,
    timescale_multiple: float = 20.0,
) -> Tuple[np.ndarray, np.ndarray, float, float]:
    """
    Sample ‖e^{tA}‖_2 on [0, t_max] and refine the peak using continuous scalar optimization.

    If ``t_max`` is None, it is set to ``timescale_multiple / |α(A)|`` so
    advection-diffusion (α = O(n²)) is not sampled on a useless window.
    """
    if t_max is None:
        t_max = suggested_t_max(A, multiple=timescale_multiple)

    ts = np.linspace(0.0, t_max, n_t)
    norms = np.zeros(n_t)
    for i, t in enumerate(ts):
        E = expm(t * A)
        norms[i] = np.linalg.norm(E, ord=2)

    peak_idx = int(np.argmax(norms))
    grid_t_peak = float(ts[peak_idx])

    dt = float(ts[1] - ts[0]) if n_t > 1 else 0.1
    t_peak, peak_val = find_exact_peak_norm(A, t_star=grid_t_peak, t_window=2.0 * dt)

    if peak_val > norms[peak_idx]:
        return ts, norms, t_peak, peak_val
    return ts, norms, grid_t_peak, float(norms[peak_idx])


def kreiss_sandwich_bounds(
    A: np.ndarray, K: Optional[float] = None, **kreiss_kwargs
) -> Tuple[float, float, float]:
    """
    Generic Kreiss sandwich: K(A) ≤ sup ‖e^{tA}‖ ≤ e·n·K(A).

    This is the baseline to beat, not the paper's theorem.
    """
    n = A.shape[0]
    if K is None:
        K, _ = kreiss_constant(A, **kreiss_kwargs)
    lower_bound = K
    upper_bound_generic = float(np.e * n * K)
    return float(K), float(lower_bound), upper_bound_generic


def family_sandwich_bounds(
    A: np.ndarray,
    c: float,
    decay_start: float = 1.0,
    K: Optional[float] = None,
    **kreiss_kwargs,
) -> dict:
    """
    Place the generic e·n·K bound next to the Bessel bound for this A = D+cN.

    ``A`` must be the headline family with the given (c, decay_start);
    the Bessel formula does not apply to other structures.
    """
    n = A.shape[0]
    if K is None:
        K, z_opt = kreiss_constant(A, **kreiss_kwargs)
    else:
        z_opt = None
    generic = float(np.e * n * K)
    bessel = bessel_uniform_bound(n, c, decay_start)
    return {
        "K": float(K),
        "z_opt": z_opt,
        "lower": float(K),
        "upper_generic": generic,
        "upper_bessel": bessel,
        "spectral_abscissa": spectral_abscissa(A),
        "n_eigs_stable": bool(np.all(np.real(eigenvalues(A)) < 0)),
    }
