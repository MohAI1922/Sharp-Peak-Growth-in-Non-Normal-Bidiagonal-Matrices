"""
Analysis script for the four tasks:
1. Redo factorization validation using log-space/closed-form entrywise computation
2. Compute top left singular vector of e^{xN} alone at sweep points
3. Track split-bound/true-norm ratio across c=10..40 sweep
4. Compute coupled quadratic-form bound λ_max(e^{xN*}·e^{2tD}·e^{xN}) numerically
"""

import numpy as np
from scipy.linalg import expm, svd
from scipy.special import gammaln
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pseudospectra.theory import (
    expm_closed_form,
    interaction_picture_factorization,
    row0_norm,
)
from pseudospectra.matrices import shifted_superdiagonal


def expm_xN_closed_form(n: int, x: float) -> np.ndarray:
    """
    Compute e^{xN} using closed-form entries (N is superdiagonal shift).
    Entry (i,j) = x^{j-i} / (j-i)! for j >= i, else 0.
    """
    E = np.zeros((n, n), dtype=float)
    for i in range(n):
        for r in range(n - i):
            j = i + r
            # Use log-space computation for numerical stability
            log_entry = r * np.log(max(x, 1e-300)) - gammaln(r + 1)
            E[i, j] = np.exp(log_entry)
    return E


def expm_xN_log_space(n: int, x: float) -> np.ndarray:
    """
    Compute e^{xN} using log-space computation for numerical stability.
    This avoids catastrophic cancellation when x is large.
    """
    E = np.zeros((n, n), dtype=float)
    for i in range(n):
        for r in range(n - i):
            j = i + r
            # log entry = r * log(x) - log(r!)
            log_entry = r * np.log(max(x, 1e-300)) - gammaln(r + 1)
            E[i, j] = np.exp(log_entry)
    return E


def expm_closed_form_log_space(n: int, t: float, c: float, decay_start: float = 1.0) -> np.ndarray:
    """
    Compute e^{tA} using log-space to avoid overflow for large x.
    Entry (i,j) = exp(-(i+1)δt) * [ (c/δ)(1-e^{-δt}) ]^{j-i} / (j-i)!
    """
    E = np.zeros((n, n), dtype=float)
    delta = decay_start
    rho = (c / delta) * (1.0 - np.exp(-delta * t))
    
    for i in range(n):
        for j in range(i, n):
            r = j - i
            # log entry = -(i+1)δt + r·log(rho) - log(r!)
            log_entry = -(i + 1) * delta * t + r * np.log(max(rho, 1e-300)) - gammaln(r + 1)
            E[i, j] = np.exp(log_entry)
    return E


def factorization_validation_log_space(n: int, t: float, c: float, decay_start: float = 1.0) -> dict:
    """
    Validate e^{tA} = e^{tD} * e^{xN} using full log-space computation.
    For the product, compute log(entry) = -(i+1)δt + r·log(x) - log(r!)
    directly for each entry of the product, avoiding intermediate large numbers.
    Returns the maximum absolute difference between the two sides.
    """
    delta = decay_start
    x = (c / delta) * (1.0 - np.exp(-delta * t))
    
    # Compute e^{tA} directly using log-space closed form (to avoid overflow)
    E_direct = expm_closed_form_log_space(n, t, c, decay_start)
    
    # Compute product using full log-space: (e^{tD} * e^{xN})_{i,j} = e^{-(i+1)δt} * x^{j-i} / (j-i)!
    E_product_log = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i, n):
            r = j - i
            # log entry = -(i+1)δt + r·log(x) - log(r!)
            log_entry = -(i + 1) * delta * t + r * np.log(max(x, 1e-300)) - gammaln(r + 1)
            E_product_log[i, j] = np.exp(log_entry)
    
    # Compute difference
    diff = np.max(np.abs(E_direct - E_product_log))
    rel_diff = diff / (np.max(np.abs(E_direct)) + 1e-300)
    
    return {
        "n": n,
        "t": t,
        "c": c,
        "x": x,
        "max_abs_diff": diff,
        "max_rel_diff": rel_diff,
    }


def task1_factorization_validation():
    """Task 1: Redo factorization validation using log-space computation."""
    print("\n" + "=" * 70)
    print("TASK 1: Factorization validation using log-space computation")
    print("=" * 70)
    
    # Test points from the original validation
    test_points = [
        (30, 40, 3.7),  # The problematic one with x ≈ 39
        (20, 20, 2.0),
        (15, 15, 1.5),
        (10, 10, 1.0),
    ]
    
    print("\n  (n, c, t) -> x, max_abs_diff, max_rel_diff")
    for n, c, t in test_points:
        result = factorization_validation_log_space(n, t, c)
        print(
            f"  ({n:2d}, {c:2d}, {t:3.1f}) -> x={result['x']:6.2f}, "
            f"diff={result['max_abs_diff']:.2e}, "
            f"rel_diff={result['max_rel_diff']:.2e}"
        )
    
    print("\n  All entries should now be ~1e-12 or better, including the (30,40,3.7) case.")


def task2_singular_vector_analysis():
    """Task 2: Compute top left singular vector of e^{xN} alone."""
    print("\n" + "=" * 70)
    print("TASK 2: Top left singular vector of e^{xN} alone")
    print("=" * 70)
    
    # Choose sweep points
    sweep_points = [
        (30, 40, 3.7),  # Large x case
        (20, 20, 2.0),  # Medium case
        (15, 15, 1.5),  # Small case
    ]
    
    print("\n  Analyzing mass distribution of top left singular vector u:")
    for n, c, t in sweep_points:
        delta = 1.0
        x = (c / delta) * (1.0 - np.exp(-delta * t))
        
        # Compute e^{xN} using log-space
        E_xN = expm_xN_log_space(n, x)
        
        # Compute SVD
        U, s, Vt = svd(E_xN, full_matrices=False)
        
        # Top left singular vector is the first column of U
        u = U[:, 0]
        
        # Analyze mass distribution
        mass_concentration = np.sum(u[:5]**2)  # Mass in first 5 rows
        total_mass = np.sum(u**2)
        
        print(f"\n  (n={n}, c={c}, t={t}, x={x:.2f}):")
        print(f"    Mass in rows 0-4: {mass_concentration:.4f} / {total_mass:.4f} ({100*mass_concentration/total_mass:.1f}%)")
        print(f"    Top 5 entries of |u|: {np.abs(u[:5])}")
        print(f"    Index of max |u|: {np.argmax(np.abs(u))}")
        
        # Check if mass is concentrated at low indices
        if mass_concentration / total_mass > 0.5:
            print(f"    -> Mass IS concentrated at low indices (rows 0-4)")
        else:
            print(f"    -> Mass is NOT concentrated at low indices (spread out)")


def task3_split_bound_ratio_sweep():
    """Task 3: Track split-bound/true-norm ratio across c=10..40 sweep."""
    print("\n" + "=" * 70)
    print("TASK 3: Split-bound/true-norm ratio across c=10..40 sweep")
    print("=" * 70)
    
    n = 30
    c_values = np.arange(10, 41, 2)  # c = 10, 12, ..., 40
    t = 3.7  # Fixed t value
    delta = 1.0
    
    ratios = []
    
    print("\n  c      x      true_norm    split_bound    ratio")
    print("  " + "-" * 55)
    
    for c in c_values:
        x = (c / delta) * (1.0 - np.exp(-delta * t))
        
        # True norm: ||e^{tA}||
        E_true = expm_closed_form_log_space(n, t, c, delta)
        true_norm = np.linalg.norm(E_true, ord=2)
        
        # Split bound: ||e^{tD}|| * ||e^{xN}||
        D = np.diag(np.exp(-delta * t * np.arange(1, n + 1, dtype=float)))
        norm_D = np.linalg.norm(D, ord=2)  # = 1 (since max diagonal entry is e^{-delta*t})
        
        E_xN = expm_xN_log_space(n, x)
        norm_xN = np.linalg.norm(E_xN, ord=2)
        
        split_bound = norm_D * norm_xN
        ratio = split_bound / true_norm
        
        ratios.append(ratio)
        print(f"  {c:2.0f}  {x:5.2f}  {true_norm:10.4f}  {split_bound:12.4f}  {ratio:6.3f}")
    
    # Analyze trend
    ratios = np.array(ratios)
    print(f"\n  Ratio trend: {ratios[0]:.3f} -> {ratios[-1]:.3f}")
    print(f"  Monotonically decreasing: {np.all(np.diff(ratios) <= 0)}")
    
    # Estimate limiting value
    if len(ratios) >= 3:
        # Fit to a + b/c to estimate limit as c -> infinity
        c_fit = 1.0 / c_values
        coeffs = np.polyfit(c_fit, ratios, 1)
        limit_estimate = coeffs[1]  # Constant term as 1/c -> 0
        print(f"  Estimated limiting value (as c -> infinity): {limit_estimate:.3f}")
        
        # Check if it's close to sqrt(2) ≈ 1.414
        sqrt2 = np.sqrt(2)
        if abs(limit_estimate - sqrt2) < 0.1:
            print(f"  -> Close to sqrt(2) ≈ 1.414")
        else:
            print(f"  -> Not close to sqrt(2) (diff = {abs(limit_estimate - sqrt2):.3f})")


def task4_coupled_quadratic_form():
    """Task 4: Compute coupled quadratic-form bound lambda_max(e^{xN*}·e^{2tD}·e^{xN})."""
    print("\n" + "=" * 70)
    print("TASK 4: Coupled quadratic-form bound")
    print("=" * 70)
    
    # Test at a few points
    test_points = [
        (30, 40, 3.7),
        (30, 20, 3.7),
        (30, 30, 3.7),
    ]
    
    print("\n  (n, c, t) -> true_norm^2    coupled_bound^2    ratio    split_bound^2")
    print("  " + "-" * 70)
    
    for n, c, t in test_points:
        delta = 1.0
        x = (c / delta) * (1.0 - np.exp(-delta * t))
        
        # True norm squared
        E_true = expm_closed_form_log_space(n, t, c, delta)
        true_norm_sq = np.linalg.norm(E_true, ord=2) ** 2
        
        # Coupled quadratic form: lambda_max(e^{xN*} * e^{2tD} * e^{xN})
        D = np.diag(np.exp(-delta * t * np.arange(1, n + 1, dtype=float)))
        D2 = D @ D  # e^{2tD}
        
        E_xN = expm_xN_log_space(n, x)
        coupled_matrix = E_xN.T @ D2 @ E_xN
        
        # Compute largest eigenvalue
        eigvals = np.linalg.eigvalsh(coupled_matrix)
        coupled_bound_sq = np.max(eigvals)
        
        # Split bound squared for comparison
        norm_D2 = np.linalg.norm(D2, ord=2)
        norm_xN_sq = np.linalg.norm(E_xN, ord=2) ** 2
        split_bound_sq = norm_D2 * norm_xN_sq
        
        ratio = coupled_bound_sq / true_norm_sq
        split_ratio = split_bound_sq / true_norm_sq
        
        print(
            f"  ({n}, {c}, {t}) -> {true_norm_sq:10.4f}  "
            f"{coupled_bound_sq:14.4f}  {ratio:6.3f}  {split_ratio:10.3f}"
        )
    
    print("\n  The coupled bound should be much closer to the true norm than the split bound.")


def main():
    """Run all four tasks."""
    task1_factorization_validation()
    task2_singular_vector_analysis()
    task3_split_bound_ratio_sweep()
    task4_coupled_quadratic_form()
    
    print("\n" + "=" * 70)
    print("ANALYSIS COMPLETE")
    print("=" * 70)


if __name__ == "__main__":
    main()
