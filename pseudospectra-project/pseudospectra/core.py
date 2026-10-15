"""
core.py
-------
Core routines for computing pseudospectra of a matrix A.

Mathematical background
========================
The epsilon-pseudospectrum of a matrix A (n x n) is defined as:

    Lambda_eps(A) = { z in C : sigma_min(A - zI) <= eps }

equivalently (and equivalently to a resolvent-norm definition):

    Lambda_eps(A) = { z in C : || (A - zI)^{-1} || >= 1/eps }

where sigma_min(M) is the smallest singular value of M.

For a NORMAL matrix (A A* = A* A, e.g. symmetric/Hermitian), the
pseudospectrum is just the union of eps-disks around the eigenvalues.
The interesting case -- and the whole point of this project -- is
NON-NORMAL matrices, where Lambda_eps(A) can be much larger than the
union of eigenvalue disks, revealing sensitivity to perturbation and
the potential for large transient growth in exp(tA) even when A is
spectrally stable (all eigenvalues have negative real part).

Reference: L.N. Trefethen & M. Embree, "Spectra and Pseudospectra:
The Behavior of Nonnormal Matrices and Operators", Princeton Univ.
Press, 2005. (THE reference for this entire project.)

Algorithm (grid method, Trefethen Ch. 4)
=========================================
1. Choose a rectangular grid of points z = x + iy in the complex plane.
2. For each grid point z, compute sigma_min(A - zI) via SVD.
3. The pseudospectrum boundary for level eps is the contour
   sigma_min(z) = eps.

This is O(grid_points * n^3) with the naive approach (full SVD at each
point) -- fine for n up to a few hundred and grids up to ~200x200.
For bigger problems, replace with the Lui (1997) inverse-Lanczos
speedup (see bounds.py docstring for where that would slot in --
that upgrade is itself a nice "phase 2" contribution to the project).
"""

import numpy as np
from scipy.linalg import svd, eig, schur
from typing import Tuple


def min_singular_value(A: np.ndarray, z: complex) -> float:
    """
    Compute the smallest singular value of (A - zI).

    Parameters
    ----------
    A : np.ndarray
        n×n complex or real matrix.
    z : complex
        A complex scalar (or real value).

    Returns
    -------
    float
        The smallest singular value of (A - zI).

    Notes
    -----
    Uses full SVD computation: O(n³) per call.
    """
    n = A.shape[0]
    # Promote to complex so a complex z never silently drops the imaginary part
    # against a real identity.
    M = np.asarray(A, dtype=complex) - complex(z) * np.eye(n, dtype=complex)
    s = svd(M, compute_uv=False)
    return float(np.min(s.real))


def compute_pseudospectra_grid(
    A: np.ndarray,
    x_range: Tuple[float, float],
    y_range: Tuple[float, float],
    nx: int = 150,
    ny: int = 150,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute pseudospectra on a grid using the textbook grid method.

    Parameters
    ----------
    A : np.ndarray
        n×n complex or real matrix.
    x_range : tuple
        (xmin, xmax) for real-axis range.
    y_range : tuple
        (ymin, ymax) for imaginary-axis range.
    nx : int, default=150
        Grid resolution in real direction.
    ny : int, default=150
        Grid resolution in imaginary direction.

    Returns
    -------
    X : np.ndarray
        Real-axis meshgrid (ny×nx).
    Y : np.ndarray
        Imaginary-axis meshgrid (ny×nx).
    SIGMA_MIN : np.ndarray
        Matrix of smallest singular values (ny×nx).
        SIGMA_MIN[j, i] = σ_min(A - z_ij) where z_ij = X[j,i] + i*Y[j,i].

    Notes
    -----
    Complexity: O(nx · ny · n³), so be mindful of grid size.
    For n≤100, nx,ny≤200 is reasonable on a laptop; beyond that use
    ``compute_pseudospectra_grid_fast`` and/or a coarser grid.
    """
    x = np.linspace(x_range[0], x_range[1], nx)
    y = np.linspace(y_range[0], y_range[1], ny)
    X, Y = np.meshgrid(x, y)

    SIGMA_MIN = np.zeros((ny, nx))
    for j in range(ny):
        for i in range(nx):
            z = X[j, i] + 1j * Y[j, i]
            SIGMA_MIN[j, i] = min_singular_value(A, z)

    return X, Y, SIGMA_MIN


def compute_pseudospectra_grid_fast(
    A: np.ndarray,
    x_range: Tuple[float, float],
    y_range: Tuple[float, float],
    nx: int = 150,
    ny: int = 150,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute pseudospectra on a grid using Schur form preprocessing.

    Reduces A to upper triangular Schur form T = Q^H A Q before
    computing the grid. Preserves singular values of (A - zI) via
    unitary similarity, but provides a foundation for future
    inverse-Lanczos speedup (Lui 1997).

    Parameters
    ----------
    A : np.ndarray
        n×n complex or real matrix.
    x_range : tuple
        (xmin, xmax) for real-axis range.
    y_range : tuple
        (ymin, ymax) for imaginary-axis range.
    nx : int, default=150
        Grid resolution in real direction.
    ny : int, default=150
        Grid resolution in imaginary direction.

    Returns
    -------
    X : np.ndarray
        Real-axis meshgrid (ny×nx).
    Y : np.ndarray
        Imaginary-axis meshgrid (ny×nx).
    SIGMA_MIN : np.ndarray
        Matrix of smallest singular values (ny×nx).

    Notes
    -----
    Computational complexity per point is not reduced asymptotically,
    but the Schur form is the standard preprocessing for production tools
    (EigTool). This is the natural place to plug in the Lui (1997)
    inverse-Lanczos speedup in the future.
    """
    # Compute Schur decomposition: A = Q T Q^H
    T, Q = schur(A.astype(complex), output="complex")

    # Run grid computation on triangular T instead of full A
    x = np.linspace(x_range[0], x_range[1], nx)
    y = np.linspace(y_range[0], y_range[1], ny)
    X, Y = np.meshgrid(x, y)

    SIGMA_MIN = np.zeros((ny, nx))
    for j in range(ny):
        for i in range(nx):
            z = X[j, i] + 1j * Y[j, i]
            SIGMA_MIN[j, i] = min_singular_value(T, z)

    return X, Y, SIGMA_MIN


def eigenvalues(A: np.ndarray) -> np.ndarray:
    """
    Compute eigenvalues of a matrix.

    Thin wrapper around scipy.linalg.eig; used to overlay eigenvalues
    on pseudospectra plots and to check spectral stability.

    Parameters
    ----------
    A : np.ndarray
        n×n complex or real matrix.

    Returns
    -------
    np.ndarray
        Array of n eigenvalues.

    Notes
    -----
    Spectral stability is guaranteed if max(Re(eigenvalues)) < 0.
    """
    eigs, _ = eig(A)
    return eigs


def spectral_abscissa(A: np.ndarray) -> float:
    """Return max Re(λ). Negative means spectral stability of x' = Ax."""
    return float(np.max(np.real(eigenvalues(A))))


def suggested_t_max(A: np.ndarray, multiple: float = 20.0) -> float:
    """
    Time window scaled to the spectral decay rate.

    For spectrally stable A, the natural time unit is 1/|α(A)| with
    α(A) = max Re(λ). Advection-diffusion discretizations have
    α(A) = O(n²), so a fixed t_max silently misses the transient peak.
    """
    alpha = spectral_abscissa(A)
    if abs(alpha) < 1e-14:
        return float(multiple)
    return float(multiple / abs(alpha))


def condition_number_departure_from_normality(A: np.ndarray) -> Tuple[float, float]:
    """
    Compute Henrici's departure from normality.

    Parameters
    ----------
    A : np.ndarray
        n×n complex or real matrix.

    Returns
    -------
    departure_abs : float
        Absolute departure: ‖A‖_F² - Σ|λᵢ|²
    departure_rel : float
        Relative (Frobenius-normalized) departure:
        departure_abs / ‖A‖_F²

    Notes
    -----
    For a normal matrix (e.g., Hermitian), departure is exactly zero.
    A scalar measure of "how non-normal" a matrix is.
    Used to justify why certain matrix families show dramatic transient growth.
    """
    eigs = eigenvalues(A)
    frobenius_norm_sq = np.sum(np.abs(A) ** 2)
    eigenvalue_norm_sq = np.sum(np.abs(eigs) ** 2)
    departure_abs = frobenius_norm_sq - eigenvalue_norm_sq

    if frobenius_norm_sq > 0:
        departure_rel = departure_abs / frobenius_norm_sq
    else:
        departure_rel = 0.0

    return float(departure_abs), float(departure_rel)


# Blueprint-canonical name (Part I, Algorithm 6).
departure_from_normality = condition_number_departure_from_normality
