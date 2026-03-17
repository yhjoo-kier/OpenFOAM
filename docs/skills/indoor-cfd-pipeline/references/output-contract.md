# Output contract

## Minimum outputs for a successful end-to-end run

### Required intermediate artifacts
- scene JSON
- repaired scene JSON if used
- `.geo_unrolled`
- `.msh`
- OpenFOAM case directory
- mesh/import logs
- solver log(s)
- stabilization summary JSON

### Required final artifacts
- `comparison_1x2.png`
- VTK export directory
- concise machine-readable summary
- artifact paths for all deliverables

## Generation provenance

If Gemini is used, keep:
- requested backend
- requested model
- actual model used
- fallback usage flag
- fallback notice if relevant
- generation history / retry history
- input image paths
- original scenario text

## Solver provenance

Keep:
- mesh size
- selected preset
- scene variant used
- mesh risk level and reasons
- key `checkMesh` metrics
- archived solver log path

## Failure outputs

Even on failure, leave a useful summary:
- failure stage
- failure family
- retry attempts already performed
- log path(s)
- recommended next action if obvious

## Reproducibility rule

An agent reading only the summary JSON should be able to infer:
- what was requested
- what actually ran
- which artifacts were generated
- why fallback or repair happened
- where to continue from next time
