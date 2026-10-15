"""
Visualization helpers for pseudospectra and transient growth.

Good figures are essential for a convincing numerical linear algebra paper.
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import Optional, List


def plot_pseudospectra(
    X: np.ndarray,
    Y: np.ndarray,
    SIGMA_MIN: np.ndarray,
    eigs: Optional[np.ndarray] = None,
    levels: Optional[List[float]] = None,
    title: str = "Pseudospectra",
    ax: Optional[plt.Axes] = None,
    log_levels: bool = True,
) -> plt.Axes:
    """
    Contour plot of σ_min(A - zI) — the pseudospectra.

    Parameters
    ----------
    X : np.ndarray
        Real-axis meshgrid (ny×nx).
    Y : np.ndarray
        Imaginary-axis meshgrid (ny×nx).
    SIGMA_MIN : np.ndarray
        Smallest singular values at each grid point (ny×nx).
    eigs : np.ndarray, optional
        Eigenvalues to overlay as red × markers.
    levels : list, optional
        Epsilon levels for contours. Default: 10^k for k = -6 to 0.
    title : str, default="Pseudospectra"
        Plot title.
    ax : matplotlib.axes.Axes, optional
        Axes to plot on. If None, create new figure.
    log_levels : bool, default=True
        If True, plot log₁₀(σ_min) with logarithmic contour levels
        (natural since σ_min often spans many orders of magnitude).

    Returns
    -------
    ax : matplotlib.axes.Axes
        The axes object with the plot.

    Notes
    -----
    Eigenvalues lie at the center (σ_min smallest there).
    ε-pseudospectrum is the region where σ_min ≤ ε.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 7))

    if levels is None:
        levels = [10.0 ** k for k in range(-6, 1)]

    Z = SIGMA_MIN
    if log_levels:
        cs = ax.contour(
            X, Y, np.log10(Z + 1e-300),
            levels=np.log10(levels),
            colors="steelblue"
        )
        ax.clabel(cs, inline=True, fontsize=8, fmt=lambda v: f"1e{v:.0f}")
    else:
        cs = ax.contour(X, Y, Z, levels=levels, colors="steelblue")
        ax.clabel(cs, inline=True, fontsize=8)

    if eigs is not None:
        ax.plot(eigs.real, eigs.imag, "rx", markersize=8, label="eigenvalues")
        ax.legend(loc="best")

    ax.set_xlabel("Re(z)")
    ax.set_ylabel("Im(z)")
    ax.set_title(title)
    ax.axhline(0, color="gray", lw=0.5, alpha=0.5)
    ax.axvline(0, color="gray", lw=0.5, alpha=0.5)
    ax.set_aspect("equal")
    return ax


def plot_transient_growth(
    ts: np.ndarray,
    norms: np.ndarray,
    ax: Optional[plt.Axes] = None,
    K: Optional[float] = None,
    upper_bound: Optional[float] = None,
    family_upper: Optional[float] = None,
    title: str = "Transient growth",
) -> plt.Axes:
    """
    Plot ‖exp(tA)‖₂ vs t with Kreiss sandwich bounds overlaid.

    Parameters
    ----------
    ts : np.ndarray
        Array of time points.
    norms : np.ndarray
        ‖exp(tA)‖₂ at each time point.
    ax : matplotlib.axes.Axes, optional
        Axes to plot on. If None, create new figure.
    K : float, optional
        Kreiss constant (plotted as green dashed line, lower bound).
    upper_bound : float, optional
        Generic upper bound e·n·K(A) (plotted as red dashed line).
    family_upper : float, optional
        Family-specific uniform bound (purple dash-dot), e.g. the Bessel bound.
    title : str, default="Transient growth"
        Plot title.

    Returns
    -------
    ax : matplotlib.axes.Axes
        The axes object with the plot.

    Notes
    -----
    **CRITICAL**: verify visually that t_max captured the peak and
    that the sampling was fine enough. Do not trust the peak value blindly.
    The sandwich (green and red dashed lines) shows how loose the
    generic Kreiss bound is — your theorem should tighten this.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))

    ax.plot(ts, norms, "o-", lw=2, markersize=3, label=r"$\|e^{tA}\|_2$")
    ax.axhline(1.0, color="gray", ls=":", lw=1, label="baseline = 1")

    if K is not None:
        ax.axhline(K, color="green", ls="--", lw=1.5, label=f"Kreiss lower K≈{K:.2f}")
    if upper_bound is not None:
        ax.axhline(
            upper_bound, color="red", ls="--", lw=1.5,
            label=f"generic upper e n K ≈{upper_bound:.2f}"
        )
    if family_upper is not None:
        ax.axhline(
            family_upper, color="purple", ls="-.", lw=1.5,
            label=fr"family Bessel bound ≈{family_upper:.2f}"
        )

    ax.set_xlabel("t")
    ax.set_ylabel(r"$\|e^{tA}\|_2$")
    ax.set_title(title)
    ax.legend(loc="best")
    ax.grid(alpha=0.3)
    return ax

