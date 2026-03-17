# Stabilization policy

## Framework assumption

Solver instability is a normal operating condition in this project.
The pipeline is designed to detect, classify, and recover from geometry/mesh/solver fragility.

## Current practical policy

The current implementation lives mainly in:
- `scripts/run_indoor_stabilized.py`

## Mesh-risk signals

Pay attention to:
- `defaultFaces`
- max non-orthogonality
- severely non-orthogonal face count
- aspect ratio
- skewness

These are combined into a simple risk score and then mapped to a risk level.

## Preset ordering

### Low-risk mesh
Try:
- baseline
- conservative
- robust
- ultra_robust
- laminar_fallback

### Medium-risk mesh
Try:
- conservative
- robust
- ultra_robust
- laminar_fallback
- baseline

### High-risk mesh
Try:
- robust
- ultra_robust
- laminar_fallback
- conservative
- baseline

## Repair policy

The default run flow attempts:
1. original scene
2. repaired scene

This is important because some scene abstractions are geometrically valid but still not mesh-friendly.

## Timeout policy

Each `simpleFoam` attempt should use an explicit timeout.
A timeout is not automatically the same thing as a physics failure; it may simply mean the case needs a longer solve window.

## Success criteria

A run is treated as successful when:
- solver log reaches end condition / target time
- return code is zero
- failure classifier does not flag a real failure condition

## What to record

Every summary should preserve:
- scene variant used
- mesh size
- preset name
- return code
- mesh risk level
- mesh risk reasons
- checkMesh metrics
- archived solver log path
