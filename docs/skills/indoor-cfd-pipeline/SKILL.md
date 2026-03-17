---
name: indoor-cfd-pipeline
description: Orchestrate this project-local indoor CFD framework that converts a photo, rendering, text scene description, scene JSON, or partially prepared case into a simulation-ready indoor abstraction, Gmsh mesh, OpenFOAM run, stabilization retry loop, and paper-facing figures. Use when working inside this OpenFOAM project on image-to-scene abstraction, Gemini-backed scene generation, scene repair, meshing, solver stabilization, reruns, diagnostics, or result visualization/export.
---

# indoor-cfd-pipeline

Treat this as the **project-local orchestration skill** for `/home/yhjoo/projects/OpenFOAM`.

This skill exists so another agent can understand the framework quickly from the codebase:
- what the pipeline does
- which script is the canonical entrypoint
- how to route inputs to the correct stage
- which artifacts must exist after each stage
- how failures are handled and summarized

## Canonical entrypoint

For end-to-end runs, prefer:
- `scripts/run_indoor_stabilized.py`

This is the current **one-command pipeline entrypoint** for:
1. scene generation from text or image using Gemini
2. scene repair fallback
3. Gmsh geometry/mesh generation
4. OpenFOAM case creation
5. mesh-risk-aware solver preset selection
6. stabilization / retry loop
7. VTK export and comparison figure generation
8. 3D comparison figure and per-panel 3D figure generation

## Quick routing rule

Decide the starting stage before acting:
- **image / photo / rendering** → start at Gemini scene generation
- **text-only scenario** → start at Gemini scene generation
- **scene JSON** → skip Gemini and start at validation / repair / meshing
- **mesh (`.geo`, `.geo_unrolled`, `.msh`)** → skip scene generation and continue from mesh / case stage
- **existing OpenFOAM case** → inspect logs, mesh quality, and rerun stabilization only
- **completed run** → export or refine visualization / paper-facing outputs

## Core code map

Read these references for fast orientation:
- `references/pipeline-map.md`
- `references/stage-routing.md`
- `references/scene-abstraction-rules.md`
- `references/stabilization-policy.md`
- `references/output-contract.md`
- `references/paper-figure-standards.md`

Read these project docs when solver behavior or retry policy matters:
- `/home/yhjoo/projects/OpenFOAM/docs/26-03-14_indoor_cfd_stabilization_strategy.md`
- `/home/yhjoo/projects/OpenFOAM/docs/26-03-15_openfoam_solver_stabilization_compass.md`

## Operating principles

- Treat photo-derived geometry as **simulation-ready abstraction**, not literal reconstruction.
- Prefer the one-command pipeline unless the user explicitly wants a partial-stage run.
- Reuse existing artifacts instead of regenerating expensive upstream stages when possible.
- Treat divergence and mesh issues as normal failure modes, not exceptional surprises.
- Preserve provenance: requested model, actual model used, fallback usage, input image paths, mesh metrics, successful preset, and output artifact paths.
- If the workflow changes, update the references so future agents can recover context from the repository alone.

## Expected outputs from a successful end-to-end run

At minimum, the run should leave behind:
- scene JSON
- repaired scene JSON if used
- `.geo_unrolled`
- `.msh`
- OpenFOAM case directory
- solver log(s)
- stabilization summary JSON
- comparison figure (`comparison_1x2.png`)
- VTK export for post-processing

## Scope boundary

This is not a generic OpenFOAM skill.
It is the operating guide for **this repository’s indoor image/text → CFD → figure framework**.
he operating guide for **this repository’s indoor image/text → CFD → figure framework**.
