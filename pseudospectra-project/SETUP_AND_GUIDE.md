# Setup and Quick Start Guide

## Installation

### 1. Install Python 3.9+ (if not already installed)

Download from [python.org](https://www.python.org/downloads/) or use a package manager:

**Windows:**
```bash
# Option A: Download from python.org
# https://www.python.org/downloads/

# Option B: Using Chocolatey
choco install python

# Option C: Using Anaconda/Miniconda
# https://docs.conda.io/projects/conda/en/latest/user-guide/install/
```

**macOS:**
```bash
# Using Homebrew
brew install python3

# Or download from python.org
```

**Linux (Ubuntu/Debian):**
```bash
sudo apt update
sudo apt install python3 python3-pip
```

### 2. Create a Virtual Environment (Recommended)

```bash
cd pseudospectra-project

# Windows
python -m venv venv
venv\Scripts\activate

# macOS/Linux
python3 -m venv venv
source venv/bin/activate
```

### 3. Install Dependencies

```bash
pip install --upgrade pip
pip install -r requirements.txt
```

**Dependencies:**
- `numpy>=1.24` - Numerical arrays and linear algebra
- `scipy>=1.10` - SVD, eigenvalue solvers, matrix exponential
- `matplotlib>=3.7` - Plotting
- `pytest>=7.0` - Testing framework

## Quick Start

### Run Tests (Validation)

Verify that the implementation is correct against known analytical results:

```bash
# Run all tests
python -m pytest tests/test_core.py -v

# Or run the test file directly
python tests/test_core.py
```

**Expected output:**
```
======================================================================
PSEUDOSPECTRA TEST SUITE
======================================================================

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
  ✓ Normal matrix: departure ≈ ...
  ✓ Jordan block: relative departure = ...
TEST 7: Schur form consistency...
  ✓ Schur form consistency verified

======================================================================
✓ ALL TESTS PASSED
======================================================================
```

### Run the Demo (Generate Figures)

Generate publication-quality figures demonstrating pseudospectra and transient growth:

```bash
python examples/demo.py
```

**Output:** Seven PNG files in `figures/`:
1. `01_jordan_sanity_check.png` - Validates pseudospectral scaling
2. `02_toeplitz_pseudospectra.png` - Headline example pseudospectra
3. `03_toeplitz_transient_growth.png` - Transient growth curve with Kreiss bounds
4. `04_scaling_with_c.png` - Growth scaling with coupling strength (research lead)
5. `05_advdiff_pseudospectra.png` - Advection-diffusion pseudospectra
6. `06_advdiff_transient_growth.png` - Advection-diffusion growth curve
7. `07_scaling_with_Pe.png` - Growth scaling with Péclet number

**Expected runtime:** ~2-5 minutes depending on system

## Project Structure

```
pseudospectra-project/
├── pseudospectra/              # Main package
│   ├── __init__.py
│   ├── core.py                 # Pseudospectra computation (min_singular_value, grids)
│   ├── matrices.py             # Matrix family generators
│   ├── bounds.py               # Kreiss constant & transient growth theory
│   └── plotting.py             # Visualization helpers
├── examples/
│   └── demo.py                 # End-to-end reproducible study
├── tests/
│   └── test_core.py            # Correctness tests vs. analytical results
├── docs/
│   ├── IMPLEMENTATION_PLAN.md   # Detailed research roadmap
│   └── literature_notes.md      # Reference materials
├── figures/                    # Generated plots (created by demo.py)
├── requirements.txt            # Python dependencies
├── README.md                   # Project overview
└── SETUP_AND_GUIDE.md         # This file
```

## Using the Library

### Basic Usage

```python
import numpy as np
from pseudospectra import (
    min_singular_value,
    compute_pseudospectra_grid,
    eigenvalues,
    shifted_superdiagonal,
    kreiss_constant,
    max_transient_growth,
    plot_pseudospectra,
    plot_transient_growth
)
import matplotlib.pyplot as plt

# 1. Create a matrix
n = 30
A = shifted_superdiagonal(n, c=8.0)

# 2. Compute pseudospectra
X, Y, SIGMA_MIN = compute_pseudospectra_grid(
    A, 
    x_range=(-30, 2), 
    y_range=(-15, 15),
    nx=120, ny=120
)

# 3. Overlay eigenvalues
eigs = eigenvalues(A)

# 4. Plot
plot_pseudospectra(X, Y, SIGMA_MIN, eigs=eigs)
plt.show()

# 5. Compute transient growth
ts, norms, t_peak, peak = max_transient_growth(A, t_max=10, n_t=200)

# 6. Compute Kreiss bounds
K, z_opt = kreiss_constant(A, x_max=5, n_real=60, n_imag=120)

# 7. Plot growth
plot_transient_growth(ts, norms, K=K, upper_bound=np.e * n * K)
plt.show()
```

### Understanding the Output

**Pseudospectra plot:**
- Eigenvalues appear as red × markers
- Contour levels show σ_min = ε (in log scale)
- For non-normal matrices, contours bulge far from eigenvalues
- The larger the bulge, the more "hidden" the instability

**Transient growth plot:**
- Blue curve: ‖exp(tA)‖₂ vs t
- Green dashed line: Kreiss constant K(A) (lower bound)
- Red dashed line: Generic upper bound e·n·K(A)
- Gap between green and red = looseness to improve

## Research Workflow

### Stage 1: Validation
```bash
python tests/test_core.py
```
Ensure all 7 tests pass. This verifies:
- ✓ Normal matrix check (most important)
- ✓ Jordan block pseudospectral scaling
- ✓ Kreiss lower bound holds
- ✓ Grid consistency
- ✓ Spectral stability for headline family
- ✓ Departure from normality metric
- ✓ Schur form equivalence

### Stage 2: Exploration
```bash
python examples/demo.py
```
Generate figures and identify scaling patterns (see `figures/04_scaling_with_c.png`).

### Stage 3: Conjecture
Fit the scaling data to a model:
- Is it exponential? K(A) ~ exp(αc)?
- Polynomial? K(A) ~ c^β?
- Logarithmic? K(A) ~ log(c)?

### Stage 4: Proof (with your professor)
Use the resolved formula to construct a rigorous upper bound.

## Troubleshooting

### "ModuleNotFoundError: No module named 'numpy'"
```bash
# Ensure virtual environment is activated
pip install -r requirements.txt
```

### "Plots not showing"
If using a Jupyter notebook or remote terminal:
```python
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
# ... rest of your code
```

### Slow demo execution
The computation time scales as O(n³ × grid_size²):
- n=30, grid=120×120: ~10 seconds
- n=50, grid=120×120: ~40 seconds
- n=70, grid=120×120: ~120 seconds

Reduce `nx` or `ny` in `examples/demo.py` for faster iteration:
```python
X, Y, S = compute_pseudospectra_grid(A, xr, yr, nx=80, ny=80)  # Faster
```

### Memory issues
For large matrices (n>100), the grid computation can use significant memory.
Options:
1. Reduce grid resolution (nx, ny)
2. Reduce matrix size (n)
3. Use `compute_pseudospectra_grid_fast()` (Schur form preprocessing)

## Next Steps

1. **Verify installation:** Run `python tests/test_core.py`
2. **Generate figures:** Run `python examples/demo.py`
3. **Read the literature:** Start with Trefethen & Embree, Chapters 1-4
4. **Explore:** Modify parameters in `examples/demo.py` and observe patterns
5. **Form a conjecture:** What scaling law does your chosen family follow?
6. **Collaborate with your professor:** Develop the mathematical proof

## Documentation

- **`docs/IMPLEMENTATION_PLAN.md`** - Detailed research plan and novelty checklist
- **`docs/literature_notes.md`** - Key references and background
- **`docs/results_log.md`** - Log of experimental results (to be filled as you research)
- **Function docstrings** - Type hints and detailed documentation in `.py` files

## Support

For questions or issues:
1. Check function docstrings: `help(pseudospectra.core.min_singular_value)`
2. Review test file: `tests/test_core.py` shows example usage
3. Consult Trefethen & Embree (Chapter 4 for grid method, Chapter 35 for Kreiss)
4. Ask your professor!

Good luck with your research! 🚀
