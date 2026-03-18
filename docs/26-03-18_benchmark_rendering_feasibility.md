# Benchmark Rendering Feasibility Update

> Date: 2026-03-18

## Summary

The frozen 12-scene benchmark bundle now exports the full 5-view input set required by the paper plan:

- `perspective`
- `birdseye`
- `floorplan`
- `wireframe`
- `section`

This closes the previous gap between the benchmark renderer implementation (4 views) and the paper evaluation matrix (5 views, including a section view).

## What changed

1. `scripts/render_benchmark_views.py`
   - added `section` as a first-class render mode
   - added automatic section-plane inference
   - for opposite-wall opening pairs, the section is aligned to the inferred flow corridor when possible
   - per-case render manifests now record the chosen section axis and plane coordinate

2. `scripts/run_benchmark_reference_batch.py`
   - default benchmark render-view list expanded to include `section`

3. Dataset artifacts refreshed
   - `benchmark/renderings/renderings_manifest.json`
   - `benchmark/renderings/<case>/manifest.json`
   - `benchmark/manifests/reference_batch_summary.json`

## Rendering contract

Canonical path pattern:

```text
benchmark/renderings/<case>/<view>/<case>_<view>.png
```

For section views:

```text
benchmark/renderings/<case>/section/<case>_section.png
```

## Current status

- Frozen subset size: 12 scenes
- 5-view render export: 12/12 scenes
- Reference CFD success remains: 12/12 scenes

## Later frozen-20 verification update

After the frozen bundle was expanded to 20 scenes, a dedicated integrity pass confirmed:

- 5-view render bundles: **20/20**
- Reference CFD links/artifacts: **20/20**
- Evaluation scaffold coverage: **100/100 tasks**
- Hard integrity issues: **0**
- Section-axis coverage across the 20-scene render bundle: **11 x-aligned / 9 y-aligned**

This is a useful feasibility signal for the paper: the 3D-to-2D rendering contract is no longer just producing files, it is now stable enough to support reproducible dataset packaging and downstream evaluation scaffolding.

## Notes

- The section renderer is intentionally geometry-only, like the other benchmark input views.
- It is a schematic vertical section rather than a photorealistic cutaway.
- This is appropriate for controlled benchmark input generation, but wild-image demonstrations remain a separate track.

## Suggested next step

Use the now-complete 5-view frozen-12 bundle to scaffold `benchmark/evaluations/` for image-conditioned scene generation and reference-vs-predicted comparison.
