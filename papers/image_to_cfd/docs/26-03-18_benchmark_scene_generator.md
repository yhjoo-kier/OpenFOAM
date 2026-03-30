# Benchmark Scene Generator Progress

> Date: 2026-03-18

## Summary

Implemented an initial rule-based benchmark scene generator for the Image-to-CFD paper workflow.

Created:
- `scripts/generate_benchmark_scenes.py`

Generated sample scenes:
- `benchmark/scenes/a1_01.json`
- `benchmark/scenes/a1_02.json`
- `benchmark/scenes/a2_01.json`
- `benchmark/scenes/a2_02.json`
- `benchmark/scenes/a3_01.json`
- `benchmark/scenes/a3_02.json`
- `benchmark/scenes/a4_01.json`
- `benchmark/scenes/a4_02.json`

## Supported benchmark matrix

- **A1**: rectangular room + simple obstacles (0-1)
- **A2**: rectangular room + complex obstacles (2-3)
- **A3**: composite/L-shaped room + simple obstacles (0-1)
- **A4**: composite/L-shaped room + complex obstacles (3-4)

## Current design choices

- Composite rooms are represented with `room.blocks` using **exactly 2 joined rectangular blocks**.
- Rectangular rooms use legacy `room.size` for backward compatibility.
- Openings are generated as one inlet and one outlet on distinct outer walls.
- Obstacles are placed conservatively with wall clearance and pairwise spacing constraints.
- All generated scenes are passed through `validate_indoor_scene.py` before being written.

## Validation status

All 8 sample scenes passed validation.

Notes:
- Some A1/A3 examples intentionally have 0 obstacles because the benchmark spec allows 0-1 for the simple categories.
- The validator still emits prototype-oriented warnings when obstacle count is outside the earlier prompt expectation of 3-5. These warnings are not benchmark-invalidating.

## Example command

```bash
python3 scripts/generate_benchmark_scenes.py --count 3 --overwrite
```

## Recommended next steps

1. Add category manifest CSV/JSON summarizing scene metadata.
2. Add optional minimum-obstacle forcing for A1/A3 when needed for balanced subsets.
3. Run `scene_to_gmsh.py` smoke tests over the generated set in an environment with `gmsh` Python module available.
4. Add rendering automation for the 5 planned view types.
5. Freeze a first benchmark subset (e.g. 3 scenes × 4 categories = 12 scenes).
