# Figure 9 split-asset rebuild path

> Date: 2026-03-21
> Scope: scale-calibrated wireframe branch only
> Status: tooling implemented, external QC not rerun yet

## Locked representative scope

- source manifest: `benchmark/manifests/evaluation_aggregate_summary_posthoc_scaled_longest_span.json`
- row 1: `bench_a3_01 / wireframe`
  - structural `0.750`
  - CFD `0.604`
  - obstacles `0 -> 3`
  - openings `N/W -> N/W`
- row 2: `bench_a3_03 / wireframe`
  - structural `0.875`
  - CFD `0.600`
  - obstacles `1 -> 3`
  - openings `E/S -> E/S`

Rationale is unchanged from the 2026-03-21 QC lock:
- both are scale-calibrated `opening_topology_preserved_with_obstacle_hallucination` cases,
- both preserve opening-wall topology,
- both keep the current wireframe-view rationale after birdseye was rejected on readability grounds.

## Concise diagnosis

- panel `(f)` is blocked by local crowding, not by a missing case or weak aggregate pattern,
- the old assembled figure was still trying to carry too much meaning inside the prediction panel or the adjacent cue lane,
- the GT cue was competing with the same geometry cluster it was supposed to clarify.

## Implemented rebuild path

Script:
- `scripts/paper_figures/make_figure9_obstacle_hallucination_split_assets.py`

Command:
```bash
python3 scripts/paper_figures/make_figure9_obstacle_hallucination_split_assets.py
```

Primary output directory:
- `results/paper_figures/figure9_split_assets/`

Implemented asset set:
- `figure9_row1_input`
- `figure9_row1_reference`
- `figure9_row1_prediction`
- `figure9_row1_evidence`
- `figure9_row2_input`
- `figure9_row2_reference`
- `figure9_row2_prediction`
- `figure9_row2_evidence`
- `figure9_row2_detail_pair`
- `figure9_legend`
- `figure9_geometry_grid`
- `figure9_split_asset_preview`
- `figure9_split_manifest.json`

## Message-allocation rule

- main geometry panels now show only the reference scene or the predicted scene plus openings,
- the main prediction panels intentionally drop the dashed GT overlay so `(f)` stops carrying both prediction evidence and GT-cue decoding at once,
- the crowded row-2 GT comparison is moved into `figure9_row2_detail_pair`, which uses the same crop for `pred.` and `GT`,
- row-level `limited CFD penalty` interpretation lives in standalone evidence assets, not in a dense right-side cue lane.

## Expected QC implications

Expected improvement:
- `(f)` should become materially less crowded because the GT overlay and source-box burden are removed from the main prediction panel,
- the evidence lane is no longer overloaded with both summary and blocker-resolution duties,
- the legend is simpler because the split path no longer needs a separate GT-overlay symbol in the main panel grammar.

Still required before any completion decision:
- local panel-level QC on `figure9_row2_prediction` and `figure9_row2_detail_pair`,
- local panel-level QC on `figure9_row1_evidence` and `figure9_row2_evidence`,
- local assembled QC on `figure9_split_asset_preview`,
- external subagent QC and Gemini QC on the exact same rebuilt preview revision.

Status consequence:
- Figure 9 remains **incomplete**,
- the current split-asset preview is a rebuild candidate, not a release figure,
- no Drive copy and no checklist `[x]` until one exact split-asset revision passes self + subagent + Gemini together.
