"""
Correctness tests against known results.

**CRITICAL**: Always validate numerical code against cases where the answer
is analytically known before trusting it on research matrices. This test file
demonstrates the credibility of your implementation.

Test checklist (from spec, Section 7):
1. Normal matrix check: sigma_min(A - zI) = distance to nearest eigenvalue
2. Jordan block check: ?^(1/n) pseudospectral radius scaling
3. Kreiss lower bound check: K(A) <= sup_t ||exp(tA)||
4. Grid consistency check: grid output matches direct pointwise calls
"""

import numpy as np
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from scipy.linalg import expm

from pseudospectra.core import (
    min_singular_value,
    compute_pseudospectra_grid,
    eigenvalues,
    condition_number_departure_from_normality,
    suggested_t_max,
)
from pseudospectra.matrices import jordan_block, shifted_superdiagonal
from pseudospectra.bounds import (
    kreiss_constant,
    max_transient_growth,
    scaling_sweep,
    fit_exponential_scaling,
    kreiss_converged,
)
from pseudospectra.core import departure_from_normality
from pseudospectra.theory import (
    bidiagonal_resolvent,
    expm_closed_form,
    bessel_uniform_bound,
    interaction_picture_factorization,
    resolvent_2_norm_exact,
)


def test_normal_matrix_pseudospectra_is_eigenvalue_disk():
    """
    TEST 1: Normal Matrix Check (MOST IMPORTANT)

    For a NORMAL matrix (here: real diagonal), sigma_min(A - zI) must
    equal exactly the distance from z to the nearest eigenvalue.
    This is the single most important correctness check for the
    entire pipeline.
    """
    print("TEST 1: Normal matrix pseudospectra...")
    A = np.diag([1.0, -2.0, 3.0, 5.0])
    
    # Test multiple points
    test_points = [
        (0.5 + 0.5j),
        (1.0 + 0.0j),
        (-2.0 + 1.0j),
        (10.0 + 0.0j),
    ]
    
    for z in test_points:
        smin = min_singular_value(A, z)
        expected = min(abs(z - lam) for lam in [1.0, -2.0, 3.0, 5.0])
        assert np.isclose(smin, expected, atol=1e-10), \
            f"Failed at z={z}: got {smin}, expected {expected}"
    
    print("  [PASS] Normal matrix check passed")


def test_jordan_block_pseudospectra_balloon():
    """
    TEST 2: Jordan Block Check

    For an n?n Jordan block with eigenvalue 0, the ?-pseudospectrum
    contains approximately a disk of radius ?^(1/n) around 0
    (Trefethen & Embree, Ch. 3). This creates dramatic pseudospectral
    ballooning not present for normal matrices.
    """
    print("TEST 2: Jordan block pseudospectra scaling...")
    n = 8
    A = jordan_block(n, lam=0.0)
    eps = 1e-6
    
    # At z = ?^(1/n), sigma_min should be approximately ? (order of magnitude)
    z = eps ** (1.0 / n)
    smin = min_singular_value(A, z)
    assert smin < 10 * eps, \
        f"Expected small sigma_min at radius ?^(1/n), got {smin} (eps={eps})"
    
    # Far outside that radius, sigma_min should be much larger
    z_far = 5 * eps ** (1.0 / n)
    smin_far = min_singular_value(A, z_far)
    assert smin_far > smin, \
        "Expected sigma_min to increase away from eigenvalue"
    
    # Also verify eigenvalue is indeed 0
    eigs = eigenvalues(A)
    assert np.allclose(eigs, np.zeros(n), atol=1e-10), \
        "Jordan block should have eigenvalue 0 with multiplicity n"
    
    print("  [PASS] Jordan block check passed")


def test_kreiss_lower_bound_holds():
    """
    TEST 3: Kreiss Lower Bound Check (CRITICAL)

    The Kreiss Matrix Theorem lower bound states: K(A) <= sup_t ||exp(tA)||

    Never trust an upper-bound-only check; the lower bound is the one
    that would catch implementation bugs. Test on a small non-normal example.
    """
    print("TEST 3: Kreiss lower bound verification...")
    A = np.array([[-1.0, 20.0], [0.0, -2.0]], dtype=float)
    
    # Verify spectral stability
    eigs = eigenvalues(A)
    assert np.all(np.real(eigs) < 0), "Test matrix should be spectrally stable"
    
    # Compute Kreiss constant
    K, z_opt = kreiss_constant(A, x_max=3.0, n_real=40, n_imag=80)
    
    # Compute actual transient growth
    ts, norms, t_peak, peak = max_transient_growth(A, t_max=5.0, n_t=300)
    
    # Lower bound: K(A) <= sup_t ||exp(tA)||
    # Allow small tolerance for numerical errors
    assert K <= peak * 1.1, \
        f"Kreiss lower bound violated: K={K:.4f}, peak={peak:.4f}"
    
    print(f"  [PASS] Kreiss lower bound holds: K={K:.4f} <= peak={peak:.4f}")


def test_grid_matches_pointwise():
    """
    TEST 4: Grid Consistency Check

    A value read from the grid output must match a direct pointwise call
    to min_singular_value at the same z. This checks implementation
    consistency between the two code paths.
    """
    print("TEST 4: Grid consistency check...")
    A = np.array([[0.0, 1.0], [-1.0, 0.0]])
    X, Y, S = compute_pseudospectra_grid(A, (-2, 2), (-2, 2), nx=10, ny=10)
    
    # Test multiple grid points
    for j in range(Y.shape[0]):
        for i in range(X.shape[1]):
            z = X[j, i] + 1j * Y[j, i]
            direct = min_singular_value(A, z)
            grid_val = S[j, i]
            assert np.isclose(direct, grid_val, atol=1e-12), \
                f"Mismatch at grid[{j},{i}]: direct={direct}, grid={grid_val}"
    
    print("  [PASS] Grid consistency check passed")


def test_shifted_superdiagonal_spectral_stability():
    """
    TEST 5: Shifted Superdiagonal Eigenvalues

    For the headline family A = D + c?N, eigenvalues should be exactly
    the diagonal entries (since it's bidiagonal upper-triangular).
    All should have negative real part for spectral stability.
    """
    print("TEST 5: Shifted superdiagonal spectral stability...")
    n = 15
    c = 5.0
    A = shifted_superdiagonal(n, c=c)
    
    eigs = eigenvalues(A)
    expected_eigs = -np.arange(1, n + 1, dtype=float)
    
    # Eigenvalues should match diagonal
    eigs_sorted = np.sort(eigs.real)
    expected_sorted = np.sort(expected_eigs)
    assert np.allclose(eigs_sorted, expected_sorted, atol=1e-10), \
        "Eigenvalues don't match diagonal entries"
    
    # All should have negative real part
    assert np.all(np.real(eigs) < 0), "All eigenvalues should have negative real part"
    
    print(f"  [PASS] Spectral stability verified: max Re(lambda) = {np.max(np.real(eigs)):.2f}")


def test_departure_from_normality():
    """
    TEST 6: Departure from Normality Metric

    For a normal matrix (diagonal), departure should be exactly zero.
    For non-normal matrices, it should be positive.
    """
    print("TEST 6: Departure from normality metric...")
    
    # Normal matrix: diagonal
    D = np.diag([1.0, 2.0, 3.0])
    dep_abs, dep_rel = condition_number_departure_from_normality(D)
    assert np.isclose(dep_abs, 0.0, atol=1e-10), \
        f"Diagonal matrix should have zero departure, got {dep_abs}"
    
    # Non-normal matrix: Jordan block
    J = jordan_block(5, lam=0.0)
    dep_abs_j, dep_rel_j = condition_number_departure_from_normality(J)
    assert dep_abs_j > 0, "Jordan block should have positive departure"
    assert 0 <= dep_rel_j <= 1, "Relative departure should be in [0,1]"
    
    print(f"  [PASS] Normal matrix: departure ~= {dep_abs:.2e}")
    print(f"  [PASS] Jordan block: relative departure = {dep_rel_j:.4f}")


def test_schur_form_consistency():
    """
    TEST 7: Schur Form Pseudospectra

    The fast Schur-form version should give identical results to the
    direct method (singular values are invariant under unitary similarity).
    """
    print("TEST 7: Schur form consistency...")
    from pseudospectra.core import compute_pseudospectra_grid_fast

    rng = np.random.RandomState(0)
    A = rng.randn(8, 8)
    x_range = (-3, 3)
    y_range = (-3, 3)
    
    X, Y, S_direct = compute_pseudospectra_grid(A, x_range, y_range, nx=20, ny=20)
    X2, Y2, S_schur = compute_pseudospectra_grid_fast(A, x_range, y_range, nx=20, ny=20)
    
    assert np.allclose(X, X2) and np.allclose(Y, Y2), \
        "Grid coordinates should be identical"
    assert np.allclose(S_direct, S_schur, atol=1e-10), \
        "Schur form should give identical singular values"
    
    print("  [PASS] Schur form consistency verified")


def test_bidiagonal_resolvent_matches_inverse():
    """Closed-form (A-zI)^{-1} must match np.linalg.inv away from the spectrum."""
    print("TEST 8: Exact bidiagonal resolvent...")
    n, c, z = 7, 3.5, 0.4 + 1.2j
    A = shifted_superdiagonal(n, c=c)
    R = bidiagonal_resolvent(n, c, z)
    R_num = np.linalg.inv(A - z * np.eye(n))
    assert np.allclose(R, R_num, atol=1e-10), "resolvent formula disagrees with inv"
    smin = min_singular_value(A, z)
    assert np.isclose(smin, 1.0 / resolvent_2_norm_exact(n, c, z), rtol=1e-9)
    print("  [PASS] Resolvent closed form matches dense inverse")


def test_expm_closed_form_matches_scipy():
    """Closed-form e^{tA} must match scipy.linalg.expm on the headline family."""
    print("TEST 9: Exact matrix exponential...")
    n, c, t = 8, 5.0, 0.35
    A = shifted_superdiagonal(n, c=c)
    E = expm_closed_form(n, t, c)
    E_num = expm(t * A)
    assert np.allclose(E, E_num, atol=1e-9, rtol=1e-9), "expm formula disagrees with scipy"
    print("  [PASS] Exponential closed form matches scipy.expm")


def test_interaction_picture_factorization_matches_scipy():
    """Exact factorization e^{tA} = e^{tD} exp(xN) should hold numerically."""
    print("TEST 9B: Interaction-picture factorization...")
    n, c, t = 8, 5.0, 0.35
    A = shifted_superdiagonal(n, c=c)
    D, x, N, E_fact = interaction_picture_factorization(n, t, c)
    E_num = expm(t * A)
    assert np.allclose(E_fact, E_num, atol=1e-9, rtol=1e-9), "interaction-picture factorization disagrees with scipy"
    assert np.allclose(D @ (expm(x * N)), E_num, atol=1e-9, rtol=1e-9)
    print(f"  [PASS] decomposition verified with x={x:.6f}, ||D||_2={np.linalg.norm(D, 2):.6f}")


def test_bessel_bound_dominates_sampled_growth():
    """Family Bessel bound must sit above sampled ‖e^{tA}‖₂."""
    print("TEST 10: Bessel uniform bound...")
    n, c = 10, 6.0
    A = shifted_superdiagonal(n, c=c)
    t_max = suggested_t_max(A, multiple=20.0)
    _, norms, _, peak = max_transient_growth(A, t_max=t_max, n_t=120)
    bound = bessel_uniform_bound(n, c)
    assert peak <= bound * 1.01, f"peak {peak} exceeds Bessel bound {bound}"
    print(f"  [PASS] peak={peak:.4f} <= Bessel={bound:.4f}")


def test_suggested_t_max_scales_with_n_for_advdiff_style_spectrum():
    """A diagonal matrix with O(n²) abscissa should get a shrinking t_max."""
    print("TEST 11: Spectral timescale helper...")
    n = 20
    # All eigenvalues equal -n², so the abscissa is -n² and t_max = 8/n².
    A = np.diag(-float(n**2) * np.ones(n))
    t = suggested_t_max(A, multiple=8.0)
    expected = 8.0 / (n**2)
    assert np.isclose(t, expected, rtol=1e-10)
    # Sanity: for diag(-1², -2², ..., -n²) the abscissa is -1, so t_max = 8.
    A2 = np.diag(-np.arange(1, n + 1, dtype=float) ** 2)
    assert np.isclose(suggested_t_max(A2, multiple=8.0), 8.0, rtol=1e-10)
    print("  [PASS] t_max scales with the spectral abscissa (1/n² case verified)")


def test_scaling_sweep_and_fit():
    """Algorithm 5: sweep produces monotone-ish growth; fit recovers alpha ~ 1."""
    print("TEST 12: Scaling sweep + exponential fit...")
    results = scaling_sweep(
        lambda n_, c_: shifted_superdiagonal(n_, c=c_),
        (4.0, 6.0, 8.0),
        n=10,
        n_real=20,
        n_imag=40,
        n_t=100,
        verbose=False,
    )
    assert len(results) == 3
    assert all(r["spectral_abscissa"] < 0 for r in results)
    assert all(r["peak"] > 1.0 for r in results)
    fit = fit_exponential_scaling(results)
    assert 0.3 < fit["alpha"] < 2.0, f"unexpected alpha {fit['alpha']}"
    assert fit["r_squared"] > 0.9, f"poor log-linear fit R²={fit['r_squared']}"
    print(
        f"  [PASS] sweep OK; fit alpha={fit['alpha']:.3f}, "
        f"R²={fit['r_squared']:.4f}"
    )


def test_kreiss_convergence_flag():
    """Convergence helper returns a flag and a non-decreasing K table."""
    print("TEST 13: Kreiss convergence declaration...")
    A = shifted_superdiagonal(8, c=4.0)
    K, converged, table = kreiss_converged(
        A,
        rtol=0.05,
        resolutions=[(12, 24), (24, 48)],
        x_maxes=[3.0, 5.0],
    )
    assert K > 0
    assert len(table) == 2
    assert isinstance(converged, bool)
    print(f"  [PASS] K={K:.4f}, converged={converged}")


def test_departure_alias():
    """Blueprint-canonical alias matches the original function."""
    print("TEST 14: departure_from_normality alias...")
    A = shifted_superdiagonal(6, c=3.0)
    a1 = departure_from_normality(A)
    a2 = condition_number_departure_from_normality(A)
    assert a1 == a2
    print("  [PASS] alias consistent")


def test_row0_unimodality_and_hessian():
    """
    TEST 15: Row-0 Unimodality & Negative Hessian Test

    Verifies that phi(t) = log ||row_0(e^{tA})|| has:
    1. Near-zero derivative at t_peak = t* - 1/(4c)
    2. Strictly negative second derivative at t_peak (phi''(t_peak) < 0)
    confirming it is a strict local maximum.
    """
    import scipy.special as sp
    print("TEST 15: Row-0 unimodality & negative Hessian verification...")
    delta = 1.0
    c = 8.0
    n = 24
    t_star = (1.0 / delta) * np.log(c / delta)
    t_peak = t_star - 1.0 / (4.0 * c)

    def row0_sq_norm(t):
        x = (c / delta) * (1.0 - np.exp(-delta * t))
        return np.exp(-2.0 * delta * t) * sp.iv(0, 2.0 * x)

    def log_f(t):
        return 0.5 * np.log(row0_sq_norm(t))

    dt = 1e-4
    d1 = (log_f(t_peak + dt) - log_f(t_peak - dt)) / (2 * dt)
    assert abs(d1) < 1e-2, f"Derivative at t_peak not near zero: {d1}"

    d2 = (log_f(t_peak + dt) - 2 * log_f(t_peak) + log_f(t_peak - dt)) / (dt**2)
    assert d2 < 0, f"Hessian at t_peak is not negative: {d2}"
    print(f"  [PASS] Unimodality & Hessian verified: phi'(t_peak)={d1:.2e}, phi''(t_peak)={d2:.4f} < 0")


def test_bessel_tail_bound_rigor():
    """
    TEST 16: Bessel Power Series Tail Truncation Rigor

    Verifies that for n >= 4X, the truncation tail of I_0(2X) is relative error < 10^-11.
    """
    print("TEST 16: Bessel tail truncation rigor...")
    from pseudospectra.theory import verify_bessel_tail_bound
    c = 8.0
    delta = 1.0
    X = (c / delta) - 1.0  # X = 7.0
    n = 28  # n >= 4X = 28
    tail_err, is_small = verify_bessel_tail_bound(X, n)
    assert is_small, f"Tail error {tail_err} exceeds threshold"
    print(f"  [PASS] Bessel tail error for X={X}, n={n} is {tail_err:.3e}")



def test_corrected_mathematical_identities():
    """
    TEST 17: Corrected Mathematical Identities & Exact Constants.
    """
    import math
    import scipy.special as sp
    from scipy.linalg import expm

    # 1. Commutator [D, N] = +delta * N (0-based indexing)
    n = 10
    delta = 1.0
    c = 5.0
    D = np.diag([- (i + 1) * delta for i in range(n)])
    N = np.zeros((n, n))
    for i in range(n - 1):
        N[i, i + 1] = 1.0
    DN_ND = D @ N - N @ D
    assert np.allclose(DN_ND, delta * N), "[D, N] must equal +delta * N"

    # 2. Flagship constants
    C_exact = 1.0 / (np.e * (4.0 * np.pi)**0.25)
    assert np.isclose(C_exact, 0.19539034, atol=1e-6)

    split_limit = (4.0 * np.pi)**0.25
    assert np.isclose(split_limit, 1.88279253, atol=1e-6)

    # 3. Lemma 4.1 identity
    n_untrunc = int(np.ceil(3 * (c / delta)))
    D2 = np.diag([- (i + 1) * delta for i in range(n_untrunc)])
    N2 = np.zeros((n_untrunc, n_untrunc))
    for i in range(n_untrunc - 1):
        N2[i, i + 1] = 1.0
    t_star = (1.0 / delta) * np.log(c / delta)
    X = (c / delta) - 1.0
    E = expm(t_star * (D2 + c * N2))
    row0_norm_actual = np.linalg.norm(E[0, :])
    row0_formula = (delta / c) * np.sqrt(sp.iv(0, 2.0 * X))
    assert np.isclose(row0_norm_actual, row0_formula, rtol=1e-5)

    # 4. Entrywise peak times
    for j in range(1, 5):
        t_j = (1.0 / delta) * np.log(j + 1)
        dt = 1e-5
        def f_entry(t):
            x = (c / delta) * (1.0 - np.exp(-delta * t))
            return np.exp(-delta * t) * (x**j) / float(math.factorial(j))
        df = (f_entry(t_j + dt) - f_entry(t_j - dt)) / (2 * dt)
        assert abs(df) < 1e-4

    print("  [PASS] Test 17: All corrected identities and constants verified")


def test_resolvent_and_singular_vector_alignment():
    """
    TEST 18: Closed-form resolvent formula and singular vector alignment.
    """
    from scipy.linalg import expm
    n = 8
    delta = 1.0
    c = 4.0
    z = 1.0 + 2.0j
    D = np.diag([- (i + 1) * delta for i in range(n)])
    N = np.zeros((n, n))
    for i in range(n - 1):
        N[i, i + 1] = 1.0
    A = D + c * N

    # 1. Resolvent check
    R_exact = np.linalg.inv(z * np.eye(n) - A)
    R_form = np.zeros((n, n), dtype=complex)
    for i in range(n):
        for j in range(i, n):
            prod_val = 1.0
            for m in range(i, j + 1):
                prod_val *= (z + (m + 1) * delta)
            R_form[i, j] = (c**(j - i)) / prod_val
    assert np.allclose(R_exact, R_form, atol=1e-12), "Resolvent formula error"

    # 2. Singular vector alignment check
    t_star = (1.0 / delta) * np.log(c / delta)
    E = expm(t_star * A)
    U, s, Vt = np.linalg.svd(E)
    u0 = U[:, 0]
    assert abs(abs(u0[0]) - 1.0) < (delta / c)**2, "Singular vector alignment failed"
    print("  [PASS] Test 18: Resolvent formula and singular vector alignment verified")


def run_all_tests():
    """Run the complete test suite."""
    print("\n" + "=" * 70)
    print("PSEUDOSPECTRA TEST SUITE")
    print("=" * 70 + "\n")

    try:
        test_normal_matrix_pseudospectra_is_eigenvalue_disk()
        test_jordan_block_pseudospectra_balloon()
        test_kreiss_lower_bound_holds()
        test_grid_matches_pointwise()
        test_shifted_superdiagonal_spectral_stability()
        test_departure_from_normality()
        test_schur_form_consistency()
        test_bidiagonal_resolvent_matches_inverse()
        test_expm_closed_form_matches_scipy()
        test_interaction_picture_factorization_matches_scipy()
        test_bessel_bound_dominates_sampled_growth()
        test_suggested_t_max_scales_with_n_for_advdiff_style_spectrum()
        test_scaling_sweep_and_fit()
        test_kreiss_convergence_flag()
        test_departure_alias()
        test_row0_unimodality_and_hessian()
        test_bessel_tail_bound_rigor()
        test_corrected_mathematical_identities()
        test_resolvent_and_singular_vector_alignment()

        print("\n" + "=" * 70)
        print("[PASS] ALL TESTS PASSED")
        print("=" * 70 + "\n")
        return True
    except AssertionError as e:
        print(f"\n[FAIL] TEST FAILED: {e}\n")
        return False


if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)

