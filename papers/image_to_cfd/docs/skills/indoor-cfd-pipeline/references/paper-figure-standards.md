# Paper figure standards

## Purpose

The framework should not stop at solver completion.
It should produce outputs that are close to paper-ready and easy to refine.

## Minimum visualization output

For a successful run, generate at least:
- a 1x2 comparison figure showing scene/layout and flow result
- a 3D comparison figure showing geometry and flow field together

Current default artifacts:
- `results/<run_name>/comparison_1x2.png`
- `results/<run_name>/indoor_pipeline_3d_comparison.png`
- `results/<run_name>/indoor_pipeline_3d_comparison.pdf`
- `results/<run_name>/panel_geometry_3d.png`
- `results/<run_name>/panel_flow_3d.png`

## Preferred paper-facing outputs

When the request is publication-oriented, prefer producing:
- `PNG` for quick inspection and sharing
- `PDF` when vector/export quality is needed
- concise caption draft text
- explicit provenance linking figure → run → summary → source artifacts

## Consistency rules

Maintain consistent:
- font family and sizing
- colormap usage
- title / axis / legend formatting
- line widths and annotation density

## Reporting rule

A figure export should say:
- which run it came from
- whether original or repaired scene was used
- which solver preset succeeded
- where the figure file is located

## Future extension

A later paper-mode can standardize:
- comparison figure layout
- 3D render export
- PDF generation
- caption template generation
- manuscript-ready asset bundle creation
