"""
Exact formulas and family-specific bounds for A = D + c N.

Headline family
===============
Let δ = decay_start > 0 and

    A = D + c N ∈ R^{n×n},

where D = diag(-δ, -2δ, …, -nδ) and N is the superdiagonal shift
(N_{i,i+1} = 1). Then A is upper bidiagonal, spectrally stable for every
c, and both the resolvent and the matrix exponential have closed forms.

Resolvent
---------
Let B = A - z I, with diagonal a_k = -k δ - z (1-based k) and constant
superdiagonal c. The inverse is upper triangular with

    (B^{-1})_{i,j} = 0                                 if i > j,
    (B^{-1})_{i,j} = (-c)^{j-i} / ∏_{k=i}^{j} a_k     if i ≤ j

(0-based indices in code). The sign (-c)^{j-i} is required by B R = I.

    σ_min(A - z I) = 1 / ‖(A - z I)^{-1}‖_2.

Matrix exponential
------------------
With 0-based indices, r = j - i ≥ 0, and t ≥ 0,

    (e^{tA})_{i,j}
        = e^{-(i+1) δ t} · [ (c/δ) (1 - e^{-δ t}) ]^r / r! .

A uniform bound follows from 0 ≤ e^{-(i+1)δ t} (1 - e^{-δ t})^r ≤ 1:

    |(e^{tA})_{i,j}| ≤ (|c|/δ)^r / r!,

hence

    ‖e^{tA}‖_2 ≤ ‖e^{tA}‖_F ≤ √( n I_0(2 |c|/δ) ),

where I_0 is the modified Bessel function of the first kind. This replaces
the generic Kreiss factor e·n (times K(A)) by an explicit function of c
and a √n factor that saturates much more slowly.

Asymptotic peak (n large enough that r ~ |c|/δ is attainable)
-------------------------------------------------------------
Maximizing the (0, r) envelope exp(-δ t + (c/δ)(1 - e^{-δ t})) over t > 0
gives, for |c| > δ,

    t* = δ^{-1} log(|c|/δ),
    leading exponential  e^{|c|/δ - 1} · (δ/|c|).

This is a conjecture-shaping asymptotic, not a complete theorem: a matching
lower bound still needs an explicit test vector and control of the 2-norm
versus the Frobenius (or entrywise) envelope. That is the intended paper
contribution; the code below is the exact algebraic substrate for the proof.
"""

from __future__ import annotations

import math
from typing import Optional, Tuple

import numpy as np
from scipy.linalg import expm
from scipy.special import iv

from .matrices import shifted_superdiagonal


def diagonal_entries(n: int, decay_start: float = 1.0) -> np.ndarray:
    """Return D = (-δ, -2δ, …, -nδ)."""
    return -decay_start * np.arange(1, n + 1, dtype=float)


def bidiagonal_resolvent(
    n: int,
    c: float,
    z: complex,
    decay_start: float = 1.0,
) -> np.ndarray:
    """
    Exact inverse of A - z I for the headline family.

    Parameters
    ----------
    n, c, decay_start
        Family parameters (same meaning as ``shifted_superdiagonal``).
    z
        Complex shift.

    Returns
    -------
    R : (n, n) complex ndarray
        (A - z I)^{-1}. Raises ZeroDivisionError if z equals an eigenvalue.
    """
    a = diagonal_entries(n, decay_start) - z
    if np.any(np.abs(a) < 1e-15):
        raise ZeroDivisionError("z is (numerically) an eigenvalue; resolvent undefined")

    R = np.zeros((n, n), dtype=complex)
    # R[i, j] = (-c)^{j-i} / prod_{k=i}^{j} a[k]  for j >= i
    for i in range(n):
        prod = 1.0 + 0.0j
        for j in range(i, n):
            prod *= a[j]
            R[i, j] = ((-c) ** (j - i)) / prod
    return R


def resolvent_frobenius_norm(
    n: int,
    c: float,
    z: complex,
    decay_start: float = 1.0,
) -> float:
    """‖(A - z I)^{-1}‖_F from the closed-form inverse."""
    R = bidiagonal_resolvent(n, c, z, decay_start)
    return float(np.linalg.norm(R, ord="fro"))


def resolvent_2_norm_exact(
    n: int,
    c: float,
    z: complex,
    decay_start: float = 1.0,
) -> float:
    """‖(A - z I)^{-1}‖_2 via the closed-form inverse (then SVD)."""
    R = bidiagonal_resolvent(n, c, z, decay_start)
    return float(np.linalg.norm(R, ord=2))


def expm_entry(i: int, j: int, t: float, c: float, decay_start: float = 1.0) -> float:
    """Closed-form (e^{tA})_{i,j} (0-based). Zero if j < i."""
    if j < i:
        return 0.0
    if t < 0:
        raise ValueError("formula stated for t >= 0")
    r = j - i
    delta = decay_start
    rho = (c / delta) * (1.0 - np.exp(-delta * t))
    # r! via gamma is fine for r up to a few hundred; we stay in n ≤ 200
    fact = float(math.factorial(r))
    return float(np.exp(-(i + 1) * delta * t) * (rho ** r) / fact)


def expm_closed_form(
    n: int,
    t: float,
    c: float,
    decay_start: float = 1.0,
) -> np.ndarray:
    """Full matrix e^{tA} from the closed-form entries."""
    E = np.zeros((n, n), dtype=float)
    delta = decay_start
    decay = np.exp(-delta * t)
    rho = (c / delta) * (1.0 - decay)
    # factorial table
    fact = np.ones(n, dtype=float)
    for r in range(1, n):
        fact[r] = fact[r - 1] * r
    # row i has prefactor exp(-(i+1) δ t) = decay^{i+1}
    prefactor = decay
    for i in range(n):
        pref = prefactor  # decay^{i+1}
        rho_pow = 1.0
        for r in range(0, n - i):
            E[i, i + r] = pref * rho_pow / fact[r]
            rho_pow *= rho
        prefactor *= decay
    return E


def interaction_picture_factorization(
    n: int,
    t: float,
    c: float,
    decay_start: float = 1.0,
) -> Tuple[np.ndarray, float, np.ndarray, np.ndarray]:
    """Exact factorization e^{tA} = e^{tD}·e^{xN} for A = D + cN.

    Let D = diag(-δ, -2δ, ..., -nδ) and N be the superdiagonal shift. Then
    e^{-tD} N e^{tD} = e^{-δt} N, so with U(t) = e^{-tD} e^{tA} one gets

        U'(t) = c e^{-δt} N U(t),

    hence U(t) = exp( x N ) with x = (c/δ)(1 - e^{-δt}); therefore

        e^{tA} = e^{tD} exp(x N).

    This is an exact algebraic identity, not a numerical approximation.
    """
    delta = decay_start
    x = (c / delta) * (1.0 - np.exp(-delta * t))
    D = np.diag(np.exp(-delta * t * np.arange(1, n + 1, dtype=float)))
    N = np.diag(np.ones(n - 1), k=1)
    factor = expm(x * N)
    return D, x, N, D @ factor


def expm_frobenius_norm(
    n: int,
    t: float,
    c: float,
    decay_start: float = 1.0,
) -> float:
    """Exact ‖e^{tA}‖_F from the closed form (no matrix exponential)."""
    E = expm_closed_form(n, t, c, decay_start)
    return float(np.linalg.norm(E, ord="fro"))


def row_norm_lower_bound(
    n: int,
    t: float,
    c: float,
    row_index: int = 0,
    decay_start: float = 1.0,
) -> float:
    """Exact 2-norm of a given row of e^{tA} for the headline family.

    This is a rigorous lower bound on the true matrix norm because for any matrix
    M and any row index i,

        ||row_i(M)||_2 = ||M^T e_i||_2 <= ||M^T||_2 = ||M||_2.

    For the headline family A = D + cN, the row-i norm is the finite sum

        ||row_i||_2 = e^{-(i+1)δt} * sqrt(Σ_{r=0}^{n-1-i} [(c/δ)(1-e^{-δt})]^{2r}/(r!)^2).
    """
    if not 0 <= row_index < n:
        raise ValueError("row_index must satisfy 0 <= row_index < n")
    delta = decay_start
    rho = (c / delta) * (1.0 - np.exp(-delta * t))
    pref = np.exp(-(row_index + 1) * delta * t)
    terms = []
    for r in range(0, n - row_index):
        terms.append((rho ** r / math.factorial(r)) ** 2)
    return float(pref * np.sqrt(sum(terms)))


def row0_norm(
    n: int,
    t: float,
    c: float,
    decay_start: float = 1.0,
) -> float:
    """Exact finite-n lower bound from the leading row of e^{tA}."""
    return row_norm_lower_bound(n, t, c, row_index=0, decay_start=decay_start)


def bessel_uniform_bound(n: int, c: float, decay_start: float = 1.0) -> float:
    """
    Rigorous t-uniform upper bound

        ‖e^{tA}‖_2 ≤ √( n I_0(2 |c|/δ) )    for all t ≥ 0.

    Derived from |(e^{tA})_{i,j}| ≤ (|c|/δ)^{j-i} / (j-i)! and
    ∑_{r≥0} x^{2r} / (r!)^2 = I_0(2x), counting at most n superdiagonals
    with at most n entries each.
    """
    x = abs(c) / decay_start
    return float(np.sqrt(n * iv(0, 2.0 * x)))


def verify_bessel_tail_bound(X: float, n: int) -> Tuple[float, bool]:
    """
    Rigorously verifies that the tail of the power series for I_0(2X)
    beyond term n is strictly less than 10^-12 for n >= 4X.
    """
    full_val = float(iv(0, 2.0 * X))
    if abs(X) < 1e-14:
        return 0.0, True
    log_X = math.log(abs(X))
    terms = [math.exp(2.0 * r * log_X - 2.0 * math.lgamma(r + 1)) for r in range(n)]
    partial_sum = float(sum(terms))
    tail_error = float(abs(full_val - partial_sum))
    rel_error = tail_error / full_val if full_val > 0 else tail_error
    return tail_error, bool(rel_error < 1e-11)


def asymptotic_peak_envelope(c: float, decay_start: float = 1.0) -> Tuple[float, float]:
    """
    Leading-order envelope max_t exp(-δ t + (|c|/δ)(1 - e^{-δ t})).

    Returns (t_star, envelope_value). For |c| ≤ δ the maximum on t ≥ 0 is
    at t = 0 and the envelope is 1.
    """
    delta = decay_start
    ac = abs(c)
    if ac <= delta:
        return 0.0, 1.0
    t_star = float(np.log(ac / delta) / delta)
    value = float((delta / ac) * np.exp(ac / delta - 1.0))
    return t_star, value


def test_vector_growth(
    n: int,
    t: float,
    c: float,
    decay_start: float = 1.0,
    mode: str = "last_basis",
) -> float:
    """
    Explicit lower bound on ‖e^{tA}‖_2 via a test vector.

    Parameters
    ----------
    mode : {'last_basis', 'poisson'}
        ``last_basis`` uses e_n (last standard basis vector) — the last
        column of e^{tA}. ``poisson`` uses a Poisson-like vector on the
        coordinates that the closed form suggests are active at time t
        (concentrated near index ~ (c/δ)(1-e^{-δ t})).
    """
    E = expm_closed_form(n, t, c, decay_start)
    if mode == "last_basis":
        v = np.zeros(n)
        v[-1] = 1.0
    elif mode == "poisson":
        delta = decay_start
        rho = max((c / delta) * (1.0 - np.exp(-delta * t)), 0.0)
        # Target the column that carries the Stirling peak, clipped to [0, n-1]
        k = int(np.clip(round(rho), 0, n - 1))
        v = np.zeros(n)
        v[k] = 1.0
    else:
        raise ValueError(f"unknown mode {mode!r}")
    return float(np.linalg.norm(E @ v, ord=2))


def family_kreiss_sandwich(
    n: int,
    c: float,
    decay_start: float = 1.0,
    K: Optional[float] = None,
) -> dict:
    """
    Compare the generic Kreiss upper bound with the family-specific Bessel bound.

    Returns a dict with K (optional), generic_upper = e n K (if K given),
    bessel_upper, asymptotic_envelope, and the ratio generic/bessel.
    """
    bessel = bessel_uniform_bound(n, c, decay_start)
    t_star, envelope = asymptotic_peak_envelope(c, decay_start)
    out = {
        "n": n,
        "c": c,
        "decay_start": decay_start,
        "bessel_upper": bessel,
        "asymptotic_t_star": t_star,
        "asymptotic_envelope": envelope,
        "K": K,
        "generic_upper": None if K is None else float(np.e * n * K),
    }
    if K is not None:
        out["generic_over_bessel"] = out["generic_upper"] / bessel
    return out


def headline_matrix(n: int, c: float, decay_start: float = 1.0) -> np.ndarray:
    """Alias used by tests and the demo."""
    return shifted_superdiagonal(n, c=c, decay_start=decay_start)
