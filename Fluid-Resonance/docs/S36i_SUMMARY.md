# S36i SESSION SUMMARY — Lamb Vector Anatomy, MHD, Particle Physics, and Audit

**Date**: 2026-03-12
**Author**: Meridian (Claude Opus 4.6)
**Checkpoint**: 19.0.0

---

## 1. CRITICAL CORRECTION: The "95% h-" Was a Measurement Bug (Failure 21)

### What happened
The S36h robustness sweep's `probe_triadic_resonance()` function measured the helical decomposition of the Lamb vector incorrectly. It computed `lamb_hm = lamb_hat - project_h_plus(lamb_hat)`, which subtracts only the h+ component from the TOTAL Lamb vector. The remainder includes BOTH the h- component AND the longitudinal (irrotational/pressure gradient) component. Since ~90% of the Lamb vector is longitudinal (Tsinober 1990), this inflated the "h-" measurement from 50% to 95%.

### Correct measurements (6 experiments, `lamb_helical_anatomy.py`)

| Field type | h- fraction | Solenoidal fraction |
|---|---|---|
| Taylor-Green (all times) | 50.0000% | ~9% |
| Pelz (all times) | 50.0000% | 8-30% |
| Random IC | 50% +/- 1.3% | ~44% |
| Gaussian random (kinematic) | 50.1% +/- 1.3% | ~44% |
| ABC (Beltrami) | ~50% +/- 2% | — |

### Why h+/h- is always ~50%
For any real vector field f(x), its Fourier transform satisfies f_hat(-k) = conj(f_hat(k)). Combined with the helical basis property conj(h+(k)) = h-(k), this enforces near-equal h+ and h- energies. For mirror-symmetric fields (TG, Pelz), the equality is exact.

### Corrected interpretation of BT surgery
BT surgery removes the h- half of the SOLENOIDAL Lamb vector (50% of the dynamically active nonlinearity, which is only ~10% of the total Lamb vector energy for TG). Despite removing only ~5% of total Lamb vector energy, it achieves 70% enstrophy reduction. This means heterochiral interactions (h+ coupling with h-) drive enstrophy growth disproportionately, consistent with Biferale (2012).

---

## 2. MHD INVESTIGATION — Proper Induction Equation

### Previous failure (S36h)
The sweep's `probe_mhd()` used a uniform background B0 field, which is trivially advected by periodic boundary conditions and has no effect.

### S36i: Full MHD solver (`mhd_spectral_solver.py`)
Coupled system: NS + induction equation dB/dt = curl(u x B) + eta*Lap(B), with div(B) = 0.
Used solenoidal B field: B_x = B0*cos(z) (divergence-free, perpendicular to k).

### Results: MHD enstrophy suppression

| B0 strength | Z_peak ratio vs NS | Effect |
|---|---|---|
| 0.5 | 0.330 | **67% reduction** |
| 1.0 | 0.413 | 59% reduction |
| 2.0 | 0.681 | 32% reduction |
| 5.0 | 2.195 | **120% INCREASE** |

**Key physics**: Moderate magnetic field suppresses enstrophy via Alfven wave energy transfer (kinetic -> magnetic). Strong B field destabilizes via current sheet formation (Lorentz force amplifies vorticity). The relationship is non-monotonic.

**Magnetic Prandtl number**: Low Pm (liquid metals) slightly more suppressive than high Pm (plasma).

### Three independent suppression mechanisms

| Mechanism | Z_peak ratio | Reference |
|---|---|---|
| BT surgery (h- removal) | 0.486 | Biferale-Titi (2013) |
| Coriolis rotation (Omega=1) | ~0.29 | Babin-Mahalov-Nikolaenko (1999) |
| MHD (B0=0.5) | 0.330 | Sermange-Temam (1983) |

**Unifying principle**: All three break the isotropy of triadic energy transfers, suppressing the forward cascade.

---

## 3. PARTICLE PHYSICS CONNECTION — Literature Research

### Genuine mathematical connections (NOT analogies)

1. **Renormalization Group** (Yakhot-Orszag 1986, Canet+ 2022): Same RG machinery as QFT applied to stochastic NS. Yields Kolmogorov scaling, eddy viscosity. Formal mathematical tool, not just analogy.

2. **Onsager dissipative anomaly = Conservation-law anomaly**: Energy conservation breaks down for Holder-1/3 weak Euler solutions (proven: Isett 2018). Same mathematical structure as QED axial anomaly (Schwinger). Both use point-splitting regularization + RG invariance. This is a structural isomorphism.

3. **Migdal loop equations**: NS exactly reformulated as Schrodinger equation in loop space, with viscosity nu playing role of Planck constant. Wilson loop functional Psi[C] satisfies i*d_t*Psi = H*Psi. Exact reformulation, but claimed solutions are speculative.

4. **Boltzmann -> NS**: Chapman-Enskog theory rigorously derives NS from kinetic theory. Any NS blowup would mean Kn ~ 1 (continuum breaks down). This is a physical argument, not a PDE proof.

5. **Kolmogorov cascade**: Energy genuinely transfers from ordered (kinetic) to disordered (thermal) scales. The "teleportation" metaphor captures the dissipative anomaly: energy disappears from the continuum description even at zero viscosity.

### Verdict
No known path from particle physics/QFT to NS regularity proof. The connections are structural, not computational. Miller (2023/2025) proved that energy + enstrophy identities alone cannot prevent blowup — any proof needs additional structural properties of the nonlinearity.

---

## 4. AUDIT OF ANTIGRAVITY S57-S60 ("Global Regularity" Claim)

### Summary: The claim is NOT valid.

S57-S60 repeat the same fundamental error documented as Failure 20 (Trap #7: Dimensional Mismatch). All measurements are in the DISCRETE point vortex model, which is an ODE system that is ALWAYS regular for finite N. Convergence of the discrete model does not prove PDE regularity.

### Specific objections

1. **The differential inequality dZ/dt <= R_inv * Z * log(Z) is unproven for NS.** Standard bounds give Z^{3/2}, not Z*log(Z). Getting Z*log(Z) would be the millennium prize — it's not something you derive from numerical observation.

2. **R_inv = 1.85731 is a graph Laplacian eigenvalue ratio** with no proven connection to the NS stretching coefficient. Using it as a bound is circular injection (Trap #1).

3. **D/Z = 148 at small scales is physically correct** (dissipation dominates stretching) but this is a property of the discrete model, not a bound on the PDE.

4. **"Final Verdict: Global Regularity" is premature.** The NS regularity problem remains open. Our results are interesting numerical observations, not proofs.

### What IS salvageable
- The Phi map framework (discrete -> continuum) is reasonable but needs rigorous proofs
- D/Z >> 1 at small scales supports the physical intuition that NS SHOULD be regular
- The differential inequality, IF provable, would give regularity — but proving it IS the problem

---

## 5. LITERATURE CONTEXT (Key Papers)

| Paper | Finding | Relevance |
|---|---|---|
| Waleffe (1992) | 8 triadic helical classes; instability assumption | Foundation for helical decomposition |
| Biferale+ (2012) | Homochiral = inverse cascade, heterochiral = forward | Explains why removing h- suppresses forward cascade |
| Biferale-Titi (2013) | Single-sign helicity NS is globally regular | Theoretical basis for BT surgery |
| Alexakis (2017) | Hidden inverse cascade ~10% of total flux | Quantifies homochiral contribution |
| Tsinober (1990) | Solenoidal Lamb vector << total | Confirmed: 9% for TG, 44% for random |
| Buaria+ (2020) | Extreme vorticity self-attenuates via Beltramization | Self-organization toward single-helicity in extreme regions |
| Miller (2023/2025) | Energy + enstrophy alone cannot prevent blowup | Need structural properties beyond conserved quantities |
| Chen-Chen-Eyink (2003) | h+ and h- cancel at small scales | Mirror symmetry restoration |

---

## 6. UPDATED STATE

### Files created this session
- `scripts/wip/lamb_helical_anatomy.py` — 6 experiments on Lamb vector helical structure
- `scripts/wip/mhd_spectral_solver.py` — Full MHD solver with induction equation

### Files modified this session
- `RESONANCE_STATE.json` — checkpoint 19.0.0, corrected triadic structure
- Internal coordination file — S36i section, corrected S36h summary, audited S57-S60
- `Fluid-Resonance/docs/S35_POSTMORTEM.md` — Failure 21, updated surviving observations

### Surviving observations (S36i)
1. Self-stretching = 0 for aligned tubes (geometric identity, rigorous)
2. Perpendicular tube zero-stretching (symmetry, rigorous)
3. C ~ N^(-0.6) for random points (statistical, unproven)
4. Every "Regularity Proof" so far has been a trap.
5. **Dynamic BT surgery robustly halves enstrophy** — removes h- solenoidal transfers, achieves 70% Z reduction. Effect STRENGTHENS at high Re.
6. **Coriolis rotation gives equivalent suppression** (BMN 1999 mechanism)
7. **Lamb vector is ~91% longitudinal for TG** (Tsinober depletion confirmed)
8. **MHD with moderate B field suppresses enstrophy** (67% reduction at B0=0.5, comparable to BT surgery). Non-monotonic: strong B destabilizes.

### The honest assessment
We have discovered that **three independent physical mechanisms** (helical decimation, rotation, magnetic field) all suppress enstrophy growth in 3D NS by 50-70%. All work by breaking the isotropy of triadic energy transfers. This is a genuine, interesting numerical finding that could support a publication.

However, **none of this constitutes a proof of NS regularity**. The gap between "numerical observation in a spectral solver at N=32-64" and "mathematical theorem about the PDE" remains as wide as ever. The millennium problem is not solved.
