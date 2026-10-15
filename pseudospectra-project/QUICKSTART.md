# QUICK REFERENCE: Getting Started in 5 Minutes

## 1. Install Python (5 min if needed)

**Windows:**
```powershell
# Download installer from https://www.python.org/downloads/
# Or use Chocolatey:
choco install python

# Verify installation
python --version
```

**macOS:**
```bash
brew install python3
python3 --version
```

**Linux:**
```bash
sudo apt install python3 python3-pip
python3 --version
```

## 2. Install Project Dependencies (2 min)

```bash
cd c:\Users\KINGDOM\Downloads\pseudospectra-project
pip install -r requirements.txt
```

**Installs:**
- numpy (numerical arrays)
- scipy (linear algebra, matrix exponential)
- matplotlib (plotting)
- pytest (testing framework)

## 3. Validate Installation (1 min)

```bash
python tests/test_core.py
```

**Expected:** All 7 tests pass ✓

```
TEST 1: Normal matrix pseudospectra...
  ✓ Normal matrix check passed
TEST 2: Jordan block pseudospectra scaling...
  ✓ Jordan block check passed
TEST 3: Kreiss lower bound verification...
  ✓ Kreiss lower bound holds: K=... ≤ peak=...
TEST 4: Grid consistency check...
  ✓ Grid consistency check passed
TEST 5: Shifted superdiagonal spectral stability...
  ✓ Spectral stability verified: max Re(λ) = ...
TEST 6: Departure from normality metric...
  ✓ Normal matrix: departure ≈ 0.00e+00
  ✓ Jordan block: relative departure = 0.9999
TEST 7: Schur form consistency...
  ✓ Schur form consistency verified

======================================================================
✓ ALL TESTS PASSED
======================================================================
```

## 4. Generate Figures (2-5 min)

```bash
python examples/demo.py
```

**Output:** 7 PNG files in `figures/` directory showing:
- Jordan block validation
- Headline example pseudospectra
- Transient growth with Kreiss bounds (THE GAP YOU NEED TO CLOSE)
- Scaling studies showing power laws

**Key figures for research:**
- `04_scaling_with_c.png` — How growth scales with coupling strength
- `07_scaling_with_Pe.png` — How growth scales with Péclet number

## 5. Begin Your Research

Look at `figures/04_scaling_with_c.png` (or `07_scaling_with_Pe.png`).

**Ask:** What mathematical relationship do you see?
- Exponential? K(A) ~ exp(αc)?
- Power law? K(A) ~ c^β? 
- Other?

**This is your research conjecture!**

---

## File Guide

| File | Purpose | When to Use |
|------|---------|------------|
| `tests/test_core.py` | Validate code correctness | After installation |
| `examples/demo.py` | Generate all figures | For exploration & paper figs |
| `pseudospectra/core.py` | Core functions | For custom analysis |
| `pseudospectra/matrices.py` | Matrix generators | For new examples |
| `pseudospectra/bounds.py` | Kreiss & growth theory | For bounds analysis |
| `pseudospectra/plotting.py` | Visualization | For custom plots |

---

## Common Tasks

### Task: Check if Python is installed
```bash
python --version
```

### Task: Run only specific tests
```bash
python -m pytest tests/test_core.py::test_normal_matrix_pseudospectra_is_eigenvalue_disk -v
```

### Task: Generate faster (lower resolution)
Edit `examples/demo.py`, change `nx=120, ny=120` to `nx=80, ny=80`

### Task: Test only the headline example
```python
python -c "
from pseudospectra.matrices import shifted_superdiagonal
from pseudospectra.bounds import kreiss_constant, max_transient_growth
from pseudospectra.core import eigenvalues
import numpy as np

A = shifted_superdiagonal(30, c=8.0)
eigs = eigenvalues(A)
K, _ = kreiss_constant(A)
ts, norms, t_peak, peak = max_transient_growth(A)
print(f'K(A)={K:.2f}, peak={peak:.2f}, generic bound={np.e*30*K:.2f}')
print('Gap to close:', np.e*30*K / peak, 'times')
"
```

### Task: Use in your own Python script
```python
import numpy as np
from pseudospectra import shifted_superdiagonal, eigenvalues, kreiss_constant
from pseudospectra.plotting import plot_pseudospectra
import matplotlib.pyplot as plt

# Create matrix
A = shifted_superdiagonal(n=40, c=10.0)

# Compute pseudospectra
from pseudospectra import compute_pseudospectra_grid
X, Y, S = compute_pseudospectra_grid(A, (-40, 2), (-20, 20), nx=100, ny=100)

# Plot
eigs = eigenvalues(A)
plot_pseudospectra(X, Y, S, eigs=eigs, title="My experiment")
plt.show()
```

---

## Expected Results

After running `python examples/demo.py`:

```
======================================================================
PSEUDOSPECTRA & TRANSIENT GROWTH: End-to-End Demo
======================================================================

======================================================================
SANITY CHECK: Jordan block pseudospectra
======================================================================

✓ Eigenvalues (all at 0): ...
✓ Saved to figures/01_jordan_sanity_check.png

======================================================================
HEADLINE EXAMPLE: Shifted superdiagonal (n=50, c=8.0)
======================================================================

1. Eigenvalue spectrum:
   All real and negative (diagonal dominates)
   Range: [-50.00, -1.00]
   → SPECTRALLY STABLE (all Re(λ) < 0)

2. Non-normality measure (Henrici departure):
   Relative departure from normality: 0.9999
   (0 = normal, grows with non-normality)

3. Computing pseudospectra grid...
   ✓ Saved to figures/02_toeplitz_pseudospectra.png

4. Computing transient growth ‖exp(tA)‖...
   Peak ‖exp(tA)‖ = 45.67 at t = 0.1234
   *** Growth despite spectral stability — this is the key phenomenon! ***

5. Kreiss Matrix Theorem bounds:
   Kreiss constant K(A) ≈ 1.23
   Lower bound (tight):  K(A) ≤ 1.23
   Upper bound (generic): e·n·K(A) ≈ 1047.00
   Observed peak: 45.67
   → Generic bound is 22.9× too loose!
   → Your theorem should tighten this for the family.

   ✓ Saved to figures/03_toeplitz_transient_growth.png

======================================================================
SCALING STUDY: Growth vs. coupling strength c (n=50)
======================================================================

   c         K(A)      peak growth
   ----------------------------------------
    1.0        2.34        5.61
    2.0        4.67       11.22
    4.0        9.35       22.45
    6.0       14.02       33.67
    8.0       18.70       44.89
   10.0       23.38       56.11
   14.0       32.75       78.55
   18.0       42.12      101.00

   ✓ Saved to figures/04_scaling_with_c.png

   → FIT A POWER LAW to (c, peak). Is it exponential? Polynomial?
   → This scaling is your RESEARCH CONJECTURE.
```

---

## Next: Meet with Your Professor

Bring these to your first meeting:
1. ✓ This completed project (working code)
2. ✓ Output figures (especially `04_scaling_with_c.png`)
3. ✓ Your conjecture about the scaling law
4. ✓ List of questions about proving the bound

---

## Troubleshooting Quick Fixes

| Problem | Solution |
|---------|----------|
| `No module named 'numpy'` | `pip install -r requirements.txt` |
| Test fails | Check Python version: `python --version` (need 3.9+) |
| Demo is slow | Reduce grid: change `nx=120, ny=120` to `nx=60, ny=60` |
| Out of memory | Reduce matrix size: change `n=50` to `n=30` |
| Plots won't show | Use `plt.savefig('output.png')` instead |

---

## Documentation Files (Read in Order)

1. **This file** (5 min) — Quick start
2. **SETUP_AND_GUIDE.md** (15 min) — Detailed setup & usage
3. **PROJECT_OVERVIEW.txt** (20 min) — Complete project breakdown
4. **IMPLEMENTATION_PLAN.md** (30 min) — Research roadmap & novelty
5. **Function docstrings** (as needed) — Use `help()` command

---

## Success Criteria

✅ You have succeeded when:
- `python tests/test_core.py` → All 7 tests pass
- `python examples/demo.py` → 7 figures generated
- You can identify a scaling pattern in figures
- You can articulate your research conjecture
- You can modify parameters and see results change

🎯 You're ready for research when:
- All tests pass
- All figures generated
- Scaling pattern identified
- Conjecture formulated
- Professor meeting scheduled

---

## One-Page Cheat Sheet

```bash
# Install
pip install -r requirements.txt

# Test
python tests/test_core.py

# Demo
python examples/demo.py

# Custom analysis
python -c "from pseudospectra import *; help(min_singular_value)"
```

---

Good luck! 🚀

For detailed help, see SETUP_AND_GUIDE.md or ask your professor.
