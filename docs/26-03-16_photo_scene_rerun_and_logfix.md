# 26-03-16 Photo scene rerun and log classification fix

## Summary
- Revisited the `photo_sample_43546_run_v2` result interpretation.
- Fixed a false-positive failure classification bug in `scripts/run_indoor_stabilized.py`.
- Added support for using an existing scene JSON directly via `--scenario <path>.json` without Gemini regeneration.
- Re-ran the full pipeline as `photo_sample_43546_run_v2_rerun` and confirmed successful solve.

## Root cause analysis
### 1. False SIGFPE detection
The solver log parser treated the standard OpenFOAM banner
`trapFpe: Floating point exception trapping enabled (FOAM_SIGFPE).`
as an actual floating point exception.

This caused successful runs to be labeled as `floating_point_exception`.

### 2. Scenario re-generation bug
`run_indoor_stabilized.py --scenario ...` was effectively interpreted as a natural-language scenario prompt only.
When a JSON path was passed, the script regenerated a new scene with Gemini instead of using the supplied JSON.
That broke reproducibility and caused reruns to diverge from the validated `photo_sample_43546_scene_v2.json` case.

## Code changes
### `scripts/run_indoor_stabilized.py`
- Refined failure classification to ignore the normal `trapFpe` banner.
- Prioritized successful completion detection (`End` / terminal target time with zero return code).
- Increased `simpleFoam` docker timeout from `180` to `300` seconds.
- If `--scenario` points to an existing `.json` file, the script now copies and uses it directly instead of regenerating via Gemini.

## Reevaluation of previous logs
Rechecked stored logs for `photo_sample_43546_run_v2`:
- `conservative`: completed
- `robust`: completed
- `laminar_fallback`: completed
- `ultra_robust`: actual floating point exception / blow-up

This means the previous interpretation of the photo-derived case as "mesh succeeds but all solver presets fail" was incorrect.
The case is in fact solvable under multiple presets.

## Full rerun result
### Run name
- `photo_sample_43546_run_v2_rerun`

### Outcome
- Scene JSON used directly from the validated photo-derived scene
- Meshing: success
- `gmshToFoam`: success
- `checkMesh`: success
- `simpleFoam`: success
- `foamToVTK -latestTime`: success
- `visualize_indoor_case.py`: success

### Successful configuration
- preset: `robust`
- mode: `RAS`
- mesh size: `0.35`
- scene variant: `original`

### Mesh quality notes
The run succeeded, but mesh risk remains high:
- `defaultFaces = 7862`
- `maxNonOrtho = 85.88`
- `severeNonOrthoFaces = 28`
- `maxAspectRatio = 30.61`
- `maxSkewness = 0.880`

Interpretation:
- Good enough for pipeline viability demonstration
- Not yet ideal for a robust/clean production-grade geometry-to-CFD workflow

## Output locations
### Local
- Result dir: `/home/yhjoo/projects/OpenFOAM/results/photo_sample_43546_run_v2_rerun`
- Summary: `/home/yhjoo/projects/OpenFOAM/results/photo_sample_43546_run_v2_rerun/stabilization_summary.json`
- Solver log: `/home/yhjoo/projects/OpenFOAM/results/photo_sample_43546_run_v2_rerun/log.simpleFoam_original_mesh0p35_robust.txt`
- Preview image: `/home/yhjoo/projects/OpenFOAM/results/photo_sample_43546_run_v2_rerun/comparison_1x2.png`

### Google Drive copy
- `/mnt/c/Users/User/GoogleDrive/photo_sample_43546_run_v2_rerun/`

## Recommended next steps
1. Reduce `defaultFaces` significantly during gmsh/OpenFOAM transfer.
2. Lower non-orthogonality and skewness for more stable and interpretable runs.
3. Add a stricter post-run success/quality report separating:
   - solver success
   - mesh-risk warning
   - field-instability hints
4. Preserve photo-derived scene JSON paths explicitly in future reruns to avoid accidental Gemini regeneration.
