# Stage routing

## Decide the resume point from the user input

### A. Image / photo / rendering input
Use the **one-command default path** first.

Recommended entrypoint:
- `scripts/run_indoor_stabilized.py --backend api --image <path> --scenario <text> --name <run_name>`

Meaning:
1. Gemini creates a simulation-ready scene abstraction
2. scene is repaired if needed
3. Gmsh mesh is built
4. OpenFOAM case is created
5. solver stabilization is applied
6. figures and summary are exported

### B. Text-only scenario input
Route exactly like image input, but without `--image`.

### C. Existing scene JSON
Skip Gemini generation.
Use:
- `scripts/run_indoor_stabilized.py --scenario <path/to/json> --name <run_name>`

### D. Existing mesh (`.geo`, `.geo_unrolled`, `.msh`)
Do not force scene regeneration.
Continue from meshing/case generation or a lower-level script path.
Use project scripts directly when a mesh already exists.

### E. Existing OpenFOAM case
Inspect:
- `constant/polyMesh/boundary`
- `log.simpleFoam`
- `checkMesh` outputs
- prior summary JSON if present

Then decide whether to:
- rerun with a different stabilization preset
- regenerate upstream scene/mesh
- export figures only

### F. Completed run / figure request
Skip upstream regeneration unless the figure must be rebuilt deterministically.
Use existing summary, VTK, and case outputs first.

## Resume logic

Before restarting a full pipeline, check whether these already exist:
- scene JSON
- repaired scene JSON
- `.msh`
- case directory
- previous summary JSON

Prefer the earliest missing stage, not a full rerun by default.

## Failure-aware routing

Classify failure before choosing the next action:
- **scene/schema problem** → fix prompt / JSON / validation constraints
- **repair problem** → inspect repair rules and geometry assumptions
- **meshing problem** → geometry overlap, boundary issues, or invalid boolean result
- **mesh-import problem** → `gmshToFoam`, patch mapping, `defaultFaces`
- **solver-stability problem** → use mesh-risk-aware preset / retry ladder
- **visualization problem** → inspect VTK output and plotting script assumptions
