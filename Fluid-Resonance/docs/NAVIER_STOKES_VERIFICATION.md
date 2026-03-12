# Verification Guide: Spectral Invariants of Star-Cluster Graphs

**Core Invariant**: R = 1.8573068741389058

## 1. Introduction

This guide provides the steps to independently verify the spectral invariants and topological regularization mechanisms described in the formal proofs. All scripts are self-contained Python (requires NumPy + SymPy).

## 2. Spectral Invariant Verification

### 2.1 Algebraic Exactness

Derive the irreducible polynomials and the exact eigenvalue ratio:

```bash
python derive_invariant.py
```

**Expected output**:
- Roots derived from P7 (degree 7) and P5 (degree 5)
- Ratio R = 1.8573068741389058
- Rational approximation 13/7 (error < 0.009%)
- Both polynomials confirmed irreducible over Q

### 2.2 Scale Invariance (Path 1)

Verify that R remains stable across cluster sizes:

```bash
python path1_generalization.py
```

**Expected output**: R identical (to 10 decimal places) for N = 2, 4, 8, ..., 128 spoke clusters.

### 2.3 Irreducibility Proof

Verify that P7 and P5 have no factorization over Q:

```bash
python factor_polys.py
```

**Expected output**: No rational roots, no quadratic-quintic (P7) or quadratic-cubic (P5) factorization.

## 3. Hodge Duality Verification

### 3.1 Simplicial Hodge Decomposition (Path 3, Phase A)

Verify the discrete Hodge theorem on the clause complex:

```bash
python path3_phase_a.py
```

**Expected output**:
- Full complex: b0=1, b1=6, b2=0, Euler=-5
- Reduced complex: b0=1, b1=1, b2=0, Euler=0
- Vertex connectivity gap (L0) drops (ratio ~2.6)
- Flow dissipation gap (Stokes) increases (ratio ~0.59)

### 3.2 Enstrophy-Eigenvalue Identity (Phase B)

Verify that Stokes eigenvalues equal enstrophy on div-free modes:

```bash
python path3_phase_b.py
```

**Expected output**: lambda_Stokes = |curl(f)|^2 for all dissipative modes (match to machine precision).

### 3.3 Hodge Spectrum Cross-Validation (Phase C)

Verify the discrete Hodge theorem: L1 spectrum = L0 non-zero UNION L2 eigenvalues:

```bash
python path3_phase_c.py
```

**Expected output**: 34 non-zero L1 eigenvalues = 15 L0 non-zero + 19 L2 eigenvalues.

## 4. Negative Result: Random 3-SAT (Path 2)

Verify that R does not appear in random 3-SAT instances:

```bash
python path2_random_sat.py
```

**Expected output**: High spectral variance across 300 instances (N=50,100,200 x alpha=3.0-5.0), no concentration around 1.857, no special behavior at the critical threshold alpha=4.267.

## 5. JavaScript Cross-Validation

The JavaScript implementations provide an independent verification path using a from-scratch linear algebra library (no NumPy dependency):

```bash
node verify_resonance.js   # R = 1.85730687, PASS
node verify_hodge.js        # b1=6, Euler=-5, Stokes gap=0.904706, PASS
```

## 6. Minimal Verification (2 minutes)

For a quick check, compute eigenvalues of two explicit integer matrices:

```python
import numpy as np
L8 = np.array([[ 7,-1,-1,-1,-1,-1, 0, 0],
               [-1, 6,-1,-1,-1, 0, 0, 0],
               [-1,-1, 6,-1,-1, 0,-1, 0],
               [-1,-1,-1, 5,-1,-1, 0, 0],
               [-1,-1,-1,-1, 6, 0,-1, 0],
               [-1, 0, 0,-1, 0, 4,-1,-1],
               [ 0, 0,-1, 0,-1,-1, 4,-1],
               [ 0, 0, 0, 0, 0,-1,-1, 2]])
L6 = np.array([[ 5,-1,-1,-1, 0, 0],
               [-1, 4,-1, 0, 0, 0],
               [-1,-1, 3,-1, 0, 0],
               [-1, 0,-1, 4,-1,-1],
               [ 0, 0, 0,-1, 2,-1],
               [ 0, 0, 0,-1,-1, 2]])
R = min(np.linalg.eigvalsh(L8)) / min(np.linalg.eigvalsh(L6))
print(f"R = {R:.16f}")  # Expected: 1.8573068741389058
```

If this matches, the algebraic core is verified.
