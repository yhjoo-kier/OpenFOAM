# Figure 6 QC Log — Aggregate performance across benchmark categories

- Date: 2026-03-19
- Figure: Figure 6
- Script: `/home/yhjoo/projects/OpenFOAM/scripts/paper_figures/make_figure6_category_aggregate.py`
- Outputs:
  - `/home/yhjoo/projects/OpenFOAM/results/paper_figures/figure6_category_aggregate_performance.pdf`
  - `/home/yhjoo/projects/OpenFOAM/results/paper_figures/figure6_category_aggregate_performance.png`
- Source artifacts:
  - `/home/yhjoo/projects/OpenFOAM/benchmark/manifests/evaluation_aggregate_summary.json`
  - `/home/yhjoo/projects/OpenFOAM/docs/26-03-19_cli_eval_aggregate_results.md`
  - `/home/yhjoo/projects/OpenFOAM/docs/26-03-19_cli_eval_case_comparison_100of100.md`

## Figure intent
Show category-level benchmark behavior for the frozen-20 image-conditioned evaluation, with the main category narrative that A1 acts as a positive control, A2 is dragged down by opening/blockage sensitivity, A3 shows a structure-vs-CFD decoupling regime, and A4 preserves relatively decent structure but weaker CFD fidelity.

## Source artifact / case lock
- Aggregate source fixed to the frozen-20 completed benchmark manifest (`100/100` tasks).
- Categories fixed to `A1`, `A2`, `A3`, `A4` in manuscript order.
- Metrics shown:
  - panel (a): mean structural score, mean CFD score
  - panel (b): room-kind match rate, opening-wall match rate

## Layout / typography decision
- Initial width candidate: `single-column`
- Final intended width: `double-column`
- Decision reason: single-column prototypes became too crowded once panel labels, values, and category semantics were preserved at readable journal scale.
- Panel layout: `1x2`
- Subfigure labels: yes — `(a)`, `(b)` inside figure
- Subcaption handling: no separate subcaption blocks inside image; main manuscript caption should explain panel roles
- Export size: `6.95 x 3.70 in`
- PNG dpi: `600`
- PDF backend: matplotlib PDF (vector)

## QC protocol completion
1. Source artifact/case fixed — PASS
2. Single/double-column decision + panel design fixed — PASS
3. Script written — PASS
4. PDF + PNG generated — PASS
5. Self visual inspection — PASS
6. Independent subagent visual QC — PASS
7. Gemini CLI visual QC — PASS

## Iteration notes
- First draft was attempted as single-column, but repeated QC flagged top crowding / annotation density.
- Figure was therefore promoted to double-column.
- Subsequent revisions removed nonessential explanatory text from the plot area, standardized bar-end annotation style, shortened category labels, and moved a compact legend to the bottom.

## Self visual QC
- Verdict: PASS
- Method: direct image inspection via vision tool on exported PNG
- Summary:
  - readable at double-column scale
  - no visible clipping
  - panel alignment clear
  - legend acceptable
  - main message visually obvious
- Minor polish applied after pass: left-panel numeric precision reduced from 3 decimals to 2 decimals.

## Subagent visual QC
- Verdict: PASS
- Session: `agent:main:subagent:f757141a-ffad-4f74-9204-3ceeddbc8084`
- Summary:
  - overall readability good at double-column scale
  - no clipping/cropping
  - category ordering and cross-panel comparison are clear
  - only minor recommendation was to consider slightly stronger grayscale separability for orange/purple in final manuscript assembly

## Gemini CLI visual QC
- Verdict: PASS
- Artifact: `/home/yhjoo/projects/OpenFOAM/.omx/artifacts/gemini-figure6-visual-qc-20260319-102116.md`
- Summary:
  - Gemini judged layout clean, unclipped, readable, and not overcrowded
  - Gemini explicitly confirmed the intended category-level headline:
    - A1 strongest overall / positive-control behavior
    - A2 clear opening-wall collapse
    - A3 decoupling behavior remains visible
    - A4 relatively decent structure but weakest CFD

## Final disposition
- Figure 6 is approved for current manuscript production use.
- Manuscript insertion recommendation: use the PDF version in LaTeX for the main paper, with this figure treated as a `double-column` result figure.
- If a later grayscale-print review is performed, only minor palette tuning may be considered; this is non-blocking.
