# Current State — Fluid-Resonance Research (2026-03-13)

**Read this FIRST. Then CLAUDE_BRIDGE.md (active section only). Archive is in CLAUDE_BRIDGE_ARCHIVE.md.**

## Team
- **Architect** (Brendan) — human lead
- **Meridian 1** (Claude, original) — spectral analysis, triadic bounds, literature
- **Meridian 2** (Claude) — proof attempts, data audits, Antigravity oversight
- **Wanderer** (Claude) — independent verification, evolution tests, geometric analysis
- **Antigravity** (Gemini) — structural work, paper editing. **All numerical claims must be independently verified.**

## What's PROVED
- R = 1.8573068741389058 (ratio of λ_min of two integer Laplacians, 8×8 / 6×6)
- R < 2 for pure K_n cores, n≥5, w≥4 (algebraic: -16 < 16)
- P₇ and P₅ irreducible over Q
- L₁(K_n) = nI (Hodge 1-Laplacian on complete graph)
- Valve duality: connectivity weakens while dissipation strengthens
- Betti number: b₁ drops 6→1, Euler characteristic -5→0
- BT surgery: 51% enstrophy reduction, strengthens with Re

## Current Frontier: Leray Suppression Factor
**Best lead.** Cross-helical Lamb vector is 60-92% gradient (killed by Leray projection).
- **EXACT FORMULA (S93-W)**: α₊₋(θ, ρ) = 1 - (1+ρ)²(1+cosθ) / [(1+ρ²+2ρcosθ)(3-cosθ)]
- **Isotropic average (ρ=1)**: ⟨α₊₋⟩ = **1 - ln(2) = 0.30685...** (EXACT, proven analytically) — this is the **MINIMUM** over ρ
- **69.3% of cross-helical Lamb is gradient** (killed by Leray). Same-helical: only 13%.
- **CORRECTED (S94)**: ⟨α⟩(ρ) is minimized at ρ=1, increases to 2(1-ln2)≈0.614 as ρ→0 or ∞. Scale separation WEAKENS suppression.
- Meridian 2 (S93): Energy-weighted α stays BELOW geometric average (max 0.388). Cascade avoids dangerous triads.
- Meridian 1: Combined tightness 0.059–0.195 (gap = 5–17x, down from 600x).
- **Lamb-weighted avg (S94-W)**: α_E = **1/4** for ρ=1 incoherent (EXACT, stronger than 1-ln2). |P_sol(h⁺×h⁻)|² = sin²θ/4.
- **Remaining gap**: Phase coherence — per-triad bound ≤1/2 but multi-triad interference at same k3 can push α_E to ~0.9 at late times.
- **Key papers**: **Buaria, Lawson & Wilczek 2024** (Sci. Adv. eado1969) — anti-twist/self-regularization. **Buaria & Pumir 2020** (Nat. Comm.) — self-attenuation. *(S94-W misattributed to "Xiong & Yang" — corrected S94.)*

## Active Threads
1. **Enstrophy decomposition (S94)** — cross-helical = 65-81% of production, alpha reduces constant not exponent, C_eff ~ 10^{-3} << C_max = 200
2. **Leray suppression universality** — prove α < C for all solenoidal fields (Wanderer's spectral question)
3. **Buaria anti-twist integration** — Buaria, Lawson & Wilczek 2024 read (Meridian 2, S94). No formal bounds. Their Eq. 5 (stretching depends on twist ω_θ) is physical-space counterpart of our spectral α formula.
4. **Migdal loop-space** — his 2024 "No Explosion Theorem" + our Beltrami = ker(N_C) result

## Task Assignments (updated S94)
- **Wanderer**: ~~Spectral proof~~ DONE (S94) — per-triad bound proven, phase coherence gap identified
- **Meridian 2**: Read **Buaria, Lawson & Wilczek 2024** (not "Xiong & Yang") — DONE. No formal bounds. Phenomenological anti-twist via conditional averaging. Connection to our spectral α formula: amplification self-terminates because stretching creates alignment (small θ), where α → 0.
- **All**: Phase coherence control is the remaining gap — equivalent to regularity problem

## KILLED (do not re-explore)
- ❌ SAT connection (Path 2 negative — R doesn't predict satisfiability)
- ❌ s ~ r^alpha power law (Wanderer S92-W Task 2 — no correlation)
- ❌ Q2 local self-regulation (D > S reverses at extreme vorticity)
- ❌ Q3 cross-helicity boundary layer (depends on Q2, which failed)
- ❌ Q4 dissipative anomaly (not testable at our Re)
- ❌ Lyapunov functional F = solenoidal Lamb fraction (not monotone, S91)
- ❌ Star polygon = star cluster assumption
- ❌ Conjecture 9.1 for anchored clusters (counterexample: K4+8a, R=2.014)
- ❌ Triadic shell Bessel bound (worse than global C-S, Meridian 1)
- ❌ Antigravity's "Verified Shield" / Checkpoint 21.0.1 (overclaim, struck)
- ❌ Approach A r→s chain (helicity doesn't control solenoidal fraction)
- ❌ Claim 3 strong form for chiral flows (co-helical P_sol != 0, 31-53%)

## Last Handoff (Meridian 1, 2026-03-13)
- **Did**: Triadic C-S tightening script + literature survey (6 papers)
- **Result**: Shell bound failed; Leray suppression (5-17x gap) is the real mechanism
- **Next**: Integrate Buaria 2024 anti-twist with Wanderer's geometric α calculation
- **Blocked on**: Proving α stays bounded for all fields (the hard part)

## Key Files
- `scripts/wip/enstrophy_decomposition.py` — **Enstrophy production by helical sector** (Meridian 2, S94)
- `scripts/wip/spectral_weighting_proof.py` — **Lamb-weighted bound 1/4 + phase coherence analysis** (Wanderer S94)
- `scripts/wip/leray_analytical_formula.py` — **EXACT formula + verification + plots** (Wanderer S93-W)
- `scripts/wip/leray_suppression_geometry.py` — geometric α per-triad calculation (Wanderer S93)
- `scripts/wip/energy_weighted_leray.py` — energy-weighted α during NS evolution (Meridian 2)
- `scripts/wip/approach_a_triadic_bound.py` — triadic analysis (Meridian 1)
- `scripts/wip/approach_a_helicity_bound.py` — helicity sector analysis (Meridian 2)
- `scripts/wip/shared_algebraic_structure.py` — SpectralNS base solver
- `docs/LERAY_SUPPRESSION_NOTE.md` — **S88-S94 research note** (full write-up, 12 sections, summary table §12)
- `FORMAL_PROOFS.md` — proved theorems + gap analysis (4 RED gaps remain)
- `SPECTRAL_INVARIANT_RESULTS.md` — numerical verification archive
