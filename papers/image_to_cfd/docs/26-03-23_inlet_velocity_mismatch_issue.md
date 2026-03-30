# Inlet Velocity Mismatch Issue

> Date: 2026-03-23
> Severity: Critical — invalidates velocity magnitude comparison in 85/97 cases
> Status: Identified, fix in progress

## Problem

Phase 2 reference cases (0.18m mesh) and predicted cases (0.18m mesh) converge at **different solver presets**, each with a different inlet velocity:

| Preset | Inlet velocity | Ref cases | Pred cases |
|--------|---------------|-----------|------------|
| robust | 0.05 m/s | 10 | 90 |
| ultra_robust | 0.005 m/s | 75 | 0 |
| laminar_fallback | 0.01 m/s | 15 | 6 |
| conservative | 0.02 m/s | 0 | 1 |

The 0.18m mesh is finer than the 0.35m mesh used originally, causing more turbulence instability in reference cases (which have exact geometry) → cascading to ultra_robust preset with 10× lower inlet velocity.

Predicted cases have approximate geometry that happens to be more mesh-friendly → converge at robust preset.

## Impact

- 85/97 cases have mismatched inlet velocities
- Velocity magnitude similarity ≈ 0.023 (effectively zero)
- Overall CFD agreement score is suppressed by ~0.144 points

## Fix

For each predicted case, re-run the reference case with the **same preset** (same inlet velocity). This ensures fair pair-wise comparison.

- 18 reference geometries need re-running with `robust` preset (some also need `laminar_fallback` or `conservative`)
- Predicted cases remain unchanged
- CFD metrics re-computed on matched pairs only
