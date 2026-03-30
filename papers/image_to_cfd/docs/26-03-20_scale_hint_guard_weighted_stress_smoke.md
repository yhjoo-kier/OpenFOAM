# Guard-weighted wording smoke (stress lane)

> Date: 2026-03-20  
> Topic: SurfClaw / topic 2242  
> Goal: keep `layout-protected span-only` scale-hint pivot, but test a **view-guarded wording refinement** that weights the guard differently by view / scene density instead of globally hardening every case.

## 1. What changed

New setting added:

- `scale_hinted_longest_span_guard_weighted_v1`

Design intent:

- keep longest-span-only anchor
- keep layout / opening-wall / flow-path priority
- **do not** strengthen guard uniformly across all views
- emphasize:
  - `perspective`: suppress hidden depth / recessed-volume over-completion, especially for dense composite scenes
  - `section`: suppress unseen height / off-cut completion / opening-wall relocation
  - `dense composite`: prioritize opening topology + connected-room outline over obstacle detail
  - `layout-clear views` (`floorplan`, `birdseye`, `wireframe`): apply scale hint lightly, avoid over-guarding already legible layouts

## 2. Smoke subset

Stress lane fixed to 6 tasks:

1. `bench_a1_01 / wireframe`
2. `bench_a1_01 / section`
3. `bench_a3_03 / perspective`
4. `bench_a4_03 / floorplan`
5. `bench_a4_03 / perspective`
6. `bench_a4_03 / section`

All 6/6 completed successfully with CLI Gemini backend and CFD metrics.

## 3. Mean comparison over the 6-task stress lane

| Setting | Mean structural | Mean CFD | Mean rel Lx err | Mean rel Ly err | Mean rel Lz err | Opening-wall match | Room-kind match |
|---|---:|---:|---:|---:|---:|---:|---:|
| baseline (no scale hint) | 0.8542 | 0.5479 | 0.1800 | 0.1296 | 0.0819 | 6/6 | 6/6 |
| layout-protected | 0.7917 | 0.5085 | 0.0003 | 0.1791 | 0.1189 | 5/6 | 6/6 |
| view-guarded | 0.8319 | 0.4782 | 0.0003 | 0.2036 | 0.0946 | 5/6 | 6/6 |
| **guard-weighted** | **0.8542** | **0.4865** | **0.0003** | 0.1814 | **0.0734** | 5/6 | 5/6 |

Quick read:

- `guard-weighted` keeps the strong `Lx` scale gain.
- It recovers **mean structural score** back to baseline level.
- It improves **mean Lz error** vs both baseline and prior scale-hinted redesigns.
- But **mean CFD is still below baseline** (`0.4865 < 0.5479`).
- One residual regression remains in room-kind fidelity (`5/6`), so this is **not** full-rerun ready.

## 4. Per-task comparison highlights

### 4.1 `bench_a4_03 / perspective` — main win

| Setting | Structural | CFD | rel Ly err | rel Lz err | Opening-wall | Room-kind |
|---|---:|---:|---:|---:|---:|---:|
| baseline | 0.7083 | 0.5287 | 0.1139 | 0.0332 | ✓ | ✓ |
| layout-protected | 0.9167 | 0.4541 | 0.4462 | 0.1943 | ✓ | ✓ |
| view-guarded | 0.9167 | 0.5626 | 0.3354 | 0.0332 | ✓ | ✓ |
| **guard-weighted** | **0.9167** | 0.4535 | **0.1790** | **0.0134** | ✓ | ✓ |

Interpretation:

- The new perspective wording clearly reduces the hidden-depth / extra-back-span pathology.
- Bounding-box fidelity is much better than `layout-protected` and also tighter than `view-guarded`.
- However CFD **did not inherit** the `view-guarded` gain here; it dropped back near `layout-protected`.
- So perspective geometry improved, but flow-field agreement is still not stable enough.

### 4.2 `bench_a4_03 / section` — partial fix, not enough

| Setting | Structural | CFD | rel Lz err | Opening-wall | Room-kind |
|---|---:|---:|---:|---:|---:|
| baseline | 0.9167 | 0.5014 | 0.0332 | ✓ | ✓ |
| layout-protected | 0.8333 | 0.3801 | 0.0332 | ✓ | ✓ |
| view-guarded | **0.9500** | **0.4942** | 0.0332 | ✓ | ✓ |
| **guard-weighted** | 0.8334 | 0.4707 | **0.0010** | ✓ | **✗** |

Interpretation:

- The guard-weighted section wording suppresses unseen-height over-correction (`Lz` almost exact).
- But it still collapses composite → rectangular room-kind on this stress case.
- Opening walls remain correct, yet topology fidelity is still unstable.
- `view-guarded` remains the stronger section candidate on this case.

### 4.3 `bench_a1_01 / section` — residual opening-wall failure persists

| Setting | Structural | CFD | Opening-wall |
|---|---:|---:|---:|
| baseline | 1.0000 | 0.4427 | ✓ |
| layout-protected | 0.5000 | 0.3395 | ✗ |
| view-guarded | 0.5000 | 0.3745 | ✗ |
| **guard-weighted** | 0.6250 | 0.3956 | ✗ |

Interpretation:

- Section-specific wording helped somewhat on structure and CFD versus prior scale-hinted variants.
- But the simple section case still fails to restore baseline opening-wall fidelity.
- This means the section prompt is still over-correcting geometry around the cut/opening interpretation.

### 4.4 `bench_a4_03 / floorplan` — over-guard side-effect on layout-clear view

| Setting | Structural | CFD | rel Ly err | rel Lz err |
|---|---:|---:|---:|---:|
| baseline | 0.7500 | 0.5835 | 0.1077 | 0.1943 |
| layout-protected | 1.0000 | **0.6024** | 0.1139 | 0.0332 |
| view-guarded | 1.0000 | 0.4559 | **0.0031** | 0.0976 |
| **guard-weighted** | 1.0000 | 0.4321 | 0.1963 | 0.0332 |

Interpretation:

- The intended “apply lightly on layout-clear views” rule did **not** yet translate into stable behavior.
- Floorplan CFD regressed further, despite structural score staying perfect.
- This suggests the generic top-level wording is still globally biasing generation in a way that can harm already-clear plan views.

## 5. Verdict

`guard-weighted` is **useful as localization evidence**, but **not** a new frozen-20 main candidate yet.

### What it successfully demonstrated

- strong `Lx` scale correction remains reproducible
- perspective hidden-depth suppression can improve room-bbox fidelity without breaking opening-wall identity
- softer per-view weighting can recover some structural damage introduced by earlier globally guarded wording

### Why it is not promoted

- mean CFD on the stress lane remains below baseline
- `bench_a4_03/section` still loses composite room-kind
- `bench_a1_01/section` still misses opening-wall fidelity
- `bench_a4_03/floorplan` regresses despite being a layout-clear view

## 6. Next redesign direction

Do **not** full-rerun frozen-20 yet.

Recommended next move:

1. keep `layout-protected span-only` as the main parent design
2. split guard logic more sharply into:
   - `section`: explicitly prefer **opening-wall identity over section regularization**
   - `perspective`: preserve current hidden-depth suppression
   - `floorplan/birdseye/wireframe`: strip almost all extra guard text, keep only minimal topology-first sentence
3. run another **micro-smoke** first on:
   - `bench_a1_01/section`
   - `bench_a4_03/floorplan`
   - `bench_a4_03/section`
4. only reconsider frozen-20 rerun if these residual regressions clear without reintroducing the perspective failure

## 7. Files touched

- `scripts/run_benchmark_evaluation_task.py`
- `scripts/scaffold_benchmark_evaluations.py`
- `benchmark/evaluations_scale_hinted_guard_weighted/...`
