# Appendix: VLM Prompt Template for Image-to-CFD Geometry Extraction

This document records the full prompt template sent to Google Gemini (as the VLM backbone) in the image-to-CFD pipeline. The prompt is defined in `scripts/generate_indoor_scene_with_gemini.py` as `PROMPT_TEMPLATE` and is populated at runtime by the `build_prompt()` function, which substitutes a scenario description and optionally appends an absolute-scale hint.

---

## Prompt Template

```text
You are generating a structured indoor CFD scene specification.
Return ONLY valid JSON. Do not include markdown, comments, or explanation.

Interpret any provided image(s), photos, or renderings as references for a **simulation-ready abstraction**,
not as a requirement for literal photoreal reconstruction.
Preserve the important flow-relevant layout qualitatively while simplifying geometry into solver-friendly boxes/openings.

Your first task is to decide the room topology from the image:
1. Is the visible room well represented by one rectangle?
2. Or does the image show a recessed alcove, interior corner, bent/L-shaped plan, or a clearly joined second room block?

If the image clearly shows a recessed or bent plan, DO NOT collapse it into a single rectangle.
In that case, use `room.blocks` with exactly 2 joined rectangular blocks.
Only use a single rectangular `room.size` when the visible layout is genuinely well-approximated by one rectangle.

Use exactly this schema:
{
  "schema_version": "indoor_cfd_scene_v1",
  "units": "m",
  "coordinate_system": {
    "origin": "room_min_corner",
    "x": "west_to_east",
    "y": "south_to_north",
    "z": "floor_to_ceiling"
  },
  "room": {
    EITHER "size": {
      "Lx": number,
      "Ly": number,
      "Lz": number
    }
    OR "blocks": [
      {
        "origin": {"x": number, "y": number, "z": number},
        "size": {"dx": number, "dy": number, "dz": number}
      }
    ]
  },
  "obstacles": [
    {
      "id": string,
      "type": "box",
      "min": {"x": number, "y": number, "z": number},
      "size": {"dx": number, "dy": number, "dz": number}
    }
  ],
  "openings": [
    {
      "id": string,
      "type": "inlet" | "outlet",
      "wall": "west" | "east" | "south" | "north" | "floor" | "ceiling",
      "center": {"u": number, "v": number},
      "size": {"du": number, "dv": number}
    }
  ],
  "meta": {
    "scene_type": string,
    "description": string
  }
}

Room rules:
- If the space is well-approximated by a single rectangular room, use `room.size`.
- If the layout clearly has an L-shaped or bent plan, use `room.blocks` with exactly 2 joined rectangular blocks.
- Strong cue examples for `room.blocks`: visible interior corner, recessed side alcove, one leg extending beyond another, or a clear non-rectangular floor perimeter.
- Do not use more than 2 room blocks.
- Do not create T-shaped, cross-shaped, or fragmented rooms.
- If the evidence is ambiguous, prefer a single rectangular room.
- For composite rooms, the 2 blocks must overlap or share a connected face segment so the room is one connected fluid domain.

Opening rules:
- Create exactly 1 inlet and 1 outlet.
- Each opening must lie fully on an exposed OUTER wall.
- For wall=`west` or `east`, `center.u` must be within the room extent along y and `center.v` within the room extent along z; do not place `center.u = 0` unless the room y-extent is actually centered there.
- For wall=`south` or `north`, `center.u` must be within the room extent along x and `center.v` within the room extent along z.
- Keep openings comfortably inside wall bounds; avoid edge-touching placements.
- Use moderate opening sizes (roughly du,dv around 0.3 to 0.8 m unless the image strongly suggests otherwise).

Obstacle rules:
- Create 3 to 5 non-overlapping box obstacles only when needed for the visible flow-relevant geometry.
- If the image is mostly empty architectural space, use the minimum stable obstacle set and do not hallucinate furniture.
- If obstacle detail is uncertain, preserve room topology and opening-wall identity first, then represent only the largest flow-relevant obstacles.
- In dense scenes, prefer fewer larger solver-friendly obstacles over many speculative small boxes.
- In dense scenes, keep clear separation between large obstacles; do not let boxes overlap or interpenetrate.
- All obstacles must lie fully inside the room; for composite rooms, each obstacle must lie inside one of the room blocks.
- Obstacles must not overlap each other.

General constraints:
- Openings must lie fully on an exposed outer wall and remain within room bounds.
- Use units of meters.
- All lengths must be positive.
- Use concise IDs like obs_001, inlet_001, outlet_001.
- Favor solver-friendly abstraction over visual detail.
- Favor the simplest stable geometry that preserves the main flow-relevant layout.
- If image cues are ambiguous, choose a conservative, simulation-stable layout.

Scenario:
{scenario}
```

## Scale-Hint Addendum

When a `--scale-hint` argument is provided, the `build_prompt()` function appends the following text to the scenario block before substitution:

```text
Absolute scale hint:
- <user-provided scale hint text>
- Treat this as a soft metric anchor for the dominant horizontal room span, not as an exact full bounding box.
- Preserve image-supported topology, opening-wall identity, and visible flow path before adjusting uncertain depth, height, or obstacle detail.
- If the image leaves some geometry ambiguous, keep those unsupported parts conservative instead of stretching them to satisfy the hint.
```

---

## Summary Statistics

| Metric | Value |
|---|---|
| Word count (base template, excluding scenario placeholder) | ~659 |
| Output schema | `indoor_cfd_scene_v1` (JSON) |

## Key Structural Elements Extracted by the Prompt

The prompt instructs the VLM to extract the following from input images:

1. **Room topology** -- single rectangular room (`room.size` with Lx/Ly/Lz) or composite L-shaped plan (`room.blocks` with up to 2 joined rectangular blocks, each with origin and size).
2. **Obstacles** -- 3 to 5 axis-aligned box obstacles representing flow-relevant furniture or equipment, each specified by min-corner position and size (dx/dy/dz).
3. **Openings** -- exactly 1 inlet and 1 outlet, each placed on a named outer wall with 2D center coordinates (u/v in the wall-local frame) and size (du/dv).
4. **Coordinate system** -- fixed convention: origin at room min corner, x = west-to-east, y = south-to-north, z = floor-to-ceiling, units in meters.
5. **Metadata** -- scene type label and free-text description.

The prompt enforces solver-stability constraints throughout: no overlapping obstacles, openings inside wall bounds, connected fluid domain for composite rooms, and a preference for conservative/simple geometry when image evidence is ambiguous.
