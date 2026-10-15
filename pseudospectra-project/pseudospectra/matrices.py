"""
Structured matrix family generators for pseudospectra studies.

Each function returns an n×n real or complex ndarray with a well-defined
mathematical structure, suitable for theoretical analysis.
"""

import numpy as np
from typing import List


def jordan_block(n: int, lam: complex = 0.0) -> np.ndarray:
    """
    Create an n×n Jordan block with eigenvalue lam.

    Parameters
    ----------
    n : int
        Block size.
    lam : complex, default=0.0
        The single eigenvalue of the block.

    Returns
    -------
    np.ndarray
        n×n matrix A = λI + N where N is the nilpotent superdiagonal matrix
        (1's on the superdiagonal, 0's elsewhere).

    Notes
    -----
    For a Jordan block with eigenvalue λ, the ε-pseudospectrum is
    approximately a disk of radius ε^(1/n) centered at λ.
    This is a fundamental correctness check for the pseudospectra code.
    """
    A = lam * np.eye(n, dtype=complex)
    A += np.diag(np.ones(n - 1), k=1)
    return A


def shifted_superdiagonal(
    n: int, c: float = 8.0, decay_start: float = 1.0
) -> np.ndarray:
    """
    The headline matrix family: A = D + c·N.

    Parameters
    ----------
    n : int
        Matrix size.
    c : float, default=8.0
        Coupling strength between the diagonal and the nilpotent part.
    decay_start : float, default=1.0
        Starting value for the diagonal decay.

    Returns
    -------
    np.ndarray
        n×n matrix with:
        - Diagonal: D[i,i] = -decay_start * (i+1) for i=0,...,n-1
        - Superdiagonal: c (constant, 1's multiplied by c)
        - All other entries: 0

    Notes
    -----
    Eigenvalues are exactly the diagonal entries of D (since it's bidiagonal).
    All eigenvalues have negative real part, so spectral stability is
    guaranteed for any c > 0. This isolates the non-normality effect:
    any transient growth comes purely from the non-normal structure,
    not spectral instability.

    The bidiagonal structure makes (A - zI) analytically tractable:
    its inverse and singular values have near-closed-form expressions,
    making this family the natural target for theoretical analysis.
    """
    A = np.diag(-decay_start * np.arange(1, n + 1, dtype=float))
    A += c * np.diag(np.ones(n - 1), k=1)
    return A


def advection_diffusion_matrix(
    n: int, Pe: float = 50.0, L: float = 1.0
) -> np.ndarray:
    """
    Central-difference discretization of advection-diffusion: u_t = -c·u_x + ν·u_xx.

    Parameters
    ----------
    n : int
        Number of interior grid points.
    Pe : float, default=50.0
        Péclet number: Pe = c·L / ν (ratio of advection to diffusion).
    L : float, default=1.0
        Domain length [0, L].

    Returns
    -------
    np.ndarray
        n×n tridiagonal matrix from central-difference discretization
        with homogeneous Dirichlet boundary conditions (u=0 at x=0,L).

    Notes
    -----
    Grid spacing: h = L / (n+1).
    Main diagonal: -2ν/h²
    Sub/super-diagonals: ν/h² ± c/(2h)

    Eigenvalues scale like O(n²), so the natural timescale shrinks as n grows.
    Any transient-growth time sweep MUST use a time range relative to
    1/|max(Re(eigenvalues))|, not a fixed constant.

    This is the physical benchmark from Trefethen-Trefethen-Reddy-Driscoll
    (1993, Science paper).
    """
    h = L / (n + 1)
    nu = 1.0
    c = Pe

    main = -2.0 * nu / h**2 * np.ones(n)
    off_diffusion = nu / h**2 * np.ones(n - 1)
    off_advection = c / (2 * h) * np.ones(n - 1)

    A = (
        np.diag(main)
        + np.diag(off_diffusion + off_advection, -1)
        + np.diag(off_diffusion - off_advection, 1)
    )
    return A


def banded_toeplitz(n: int, coeffs: List[float], offsets: List[int]) -> np.ndarray:
    """
    Construct a general banded Toeplitz matrix.

    Parameters
    ----------
    n : int
        Matrix size.
    coeffs : list of float
        Coefficients to place on each diagonal.
    offsets : list of int
        Diagonal offsets for each coefficient:
        - 0: main diagonal
        - k > 0: k-th superdiagonal
        - k < 0: |k|-th subdiagonal

    Returns
    -------
    np.ndarray
        n×n banded Toeplitz matrix.

    Notes
    -----
    Allows flexible construction of band-structured matrices for robustness
    checks: e.g., testing whether a proved scaling law holds when you
    add a small symmetric sub-diagonal term to the headline family.
    """
    A = np.zeros((n, n))
    for coeff, offset in zip(coeffs, offsets):
        A += coeff * np.diag(np.ones(n - abs(offset)), k=offset)
    return A


def grcar_matrix(n: int) -> np.ndarray:
    """
    Construct the Grcar matrix (well-known benchmark for non-convex pseudospectra).

    Parameters
    ----------
    n : int
        Matrix size.

    Returns
    -------
    np.ndarray
        n×n banded non-symmetric Toeplitz matrix.

    Notes
    -----
    Structure: A[i, j] = 1 if j ∈ {i, i+1, i+2, i+3},
                        = -1 if i = j+1,
                        = 0 otherwise (clipped at matrix boundaries).

    The Grcar matrix is famous for its non-convex pseudospectrum,
    a canonical example of pseudospectral phenomena beyond simple
    eigenvalue disks.
    """
    A = np.zeros((n, n))
    for i in range(n):
        for k in [-1, 0, 1, 2, 3]:
            j = i + k
            if 0 <= j < n:
                A[i, j] = 1.0 if k >= 0 else -1.0
    return A


def random_non_normal(n: int, seed: int = 42, skew: float = 1.0) -> np.ndarray:
    """
    Generate a random non-normal matrix: A = D + skew·S.

    Parameters
    ----------
    n : int
        Matrix size.
    seed : int, default=42
        Random seed for reproducibility.
    skew : float, default=1.0
        Scaling factor for the random strictly upper-triangular part.

    Returns
    -------
    np.ndarray
        n×n matrix: D is diagonal with fixed spectrum,
        S is random strictly upper-triangular matrix.

    Notes
    -----
    Useful for exploratory studies and sanity checks. Sweeping skew
    parameter builds intuition about the relationship between
    non-normality strength and transient growth, but the matrix
    family is not analytically tractable, so results should not
    appear in the paper's main theorems (only for motivation/exploration).
    """
    rng = np.random.RandomState(seed)
    eigenvalues = -np.arange(1, n + 1, dtype=float)
    D = np.diag(eigenvalues)
    S = np.triu(rng.randn(n, n), k=1)
    return D + skew * S
