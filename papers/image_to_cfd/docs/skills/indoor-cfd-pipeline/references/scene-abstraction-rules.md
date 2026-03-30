# Scene abstraction rules

## Purpose

Convert visual or textual indoor descriptions into a **solver-friendly abstraction** that preserves the important flow structure without attempting literal reconstruction.

## Primary rule

Do **not** treat a photo or rendering as a blueprint.
Treat it as evidence for a simplified CFD-ready layout.

## Abstraction priorities

Preserve these first:
- room envelope proportions
- major obstacle blocks
- likely airflow corridors / blockages
- plausible inlet / outlet placement
- gross obstruction density and arrangement

Simplify these aggressively:
- decorative details
- small objects
- fine geometry edges
- curved details that do not materially change the intended flow pattern

## Representation rule

Default scene representation should stay compatible with `indoor_cfd_scene_v1`:
- one rectangular room
- 3 to 5 box obstacles
- exactly one inlet and one outlet
- positive dimensions
- non-overlapping obstacles
- wall-aligned openings that stay within bounds

## Ambiguity rule

If image cues are ambiguous:
- choose the simpler layout
- prefer separation between obstacles
- avoid narrow passages unless clearly intended
- avoid opening placement that is likely to create impossible or unstable geometry

## Solver-awareness rule

Prefer scene layouts that are easier to mesh and solve:
- avoid obstacle overlap or near-touching faces
- preserve wall margin
- preserve obstacle clearance
- avoid openings that clip corners or intersect obstacle envelopes
- avoid over-detailed obstacle decomposition

## Prompting rule for Gemini

When using Gemini:
- explicitly ask for simulation-ready abstraction
- explicitly reject literal photoreal reconstruction
- request concise object IDs and schema compliance
- if needed, pair the visual cue with a short textual interpretation of the room intent
