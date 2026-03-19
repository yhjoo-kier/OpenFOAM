# Gemini Figure QC Artifact

## Original user task
Figure 6 visual QC for aggregate performance across benchmark categories.

## Final prompt sent to Gemini CLI
```text
You are performing publication-readiness visual QC for a journal paper figure.
This figure uses TWO HORIZONTAL BAR CHART PANELS with shared category ordering.
Judge only MAJOR, publication-blocking issues. Ignore tiny polish issues.

Inspect the attached image and return PASS or FAIL first.
Check these items:
1) layout breakage
2) clipping or cutoff
3) unreadable text at normal paper scale
4) severe legend/label crowding or overlap
5) whether the core claim is visually obvious

Figure intent: aggregate performance across benchmark categories for the OpenFOAM Image-to-CFD frozen-20 benchmark.
Expected headline: A1 behaves as a positive control; A2 is opening/blockage-sensitive; A3 shows structure-vs-CFD decoupling; A4 keeps relatively decent structure but weaker CFD.

Keep the answer concise, bullet-based, and only fail if there is a clear blocking issue.

Reference image attachments:
@/home/yhjoo/projects/OpenFOAM/docs/figure_qc/_inputs/figure6_category_aggregate_performance.png
```

## Gemini output (raw)
```
**PASS**

*   **Layout:** Clean two-panel horizontal layout. Panels are well-aligned with a shared Y-axis (category labels).
*   **Clipping:** No clipping observed; all axis labels, titles, and data values are fully within the frame.
*   **Readability:** Text is crisp and sized appropriately for publication. Data labels next to bars are clear.
*   **Crowding:** Legend at the bottom is well-spaced. Data labels do not overlap grid lines or bars in a way that obscures information.
*   **Core Claim Verification:**
    *   **A1 (Control):** Shows the highest overall performance in both panels.
    *   **A2 (Opening-sensitive):** Shows a massive drop in "Opening-wall" match rate (48%) compared to others.
    *   **A3 (Decoupling):** Shows a relatively small gap between Structural/CFD scores compared to A4, but higher CFD than A2/A4.
    *   **A4 (Structure vs CFD):** Clearly shows "decent" structure (0.71) but the lowest CFD performance (0.43), supporting the "weaker CFD" claim.
```

## Concise summary
**PASS**

*   **Layout:** Clean two-panel horizontal layout. Panels are well-aligned with a shared Y-axis (category labels).
*   **Clipping:** No clipping observed; all axis labels, titles, and data values are fully within the frame.
*   **Readability:** Text is crisp and sized appropriately for publication. Data labels next to bars are clear.
*   **Crowding:** Legend at the bottom is well-spaced. Data labels do not overlap grid lines or bars in a way that obscures information.
*   **Core Claim Verification:**
    *   **A1 (Control):** Shows the highest overall performance in both panels.
    *   **A2 (Opening-sensitive):** Shows a massive drop in "Opening-wall" match rate (48%) compared to others.
    *   **A3 (Decoupling):** Shows a relatively small gap between Structural/CFD scores compared to A4, but higher CFD than A2/A4.
    *   **A4 (Structure vs CFD):** Clearly shows "decent" structure (0.71) but the lowest CFD performance (0.43), supporting the "weaker CFD" claim.

## Action items / next steps
- Use Gemini verdict as the third QC gate for Figure 6.
- If Gemini returns PASS, mark Figure 6 complete in the production plan.
- If Gemini returns FAIL, revise the script and rerun all QC steps.
