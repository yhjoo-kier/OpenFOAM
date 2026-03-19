# Original user task
Gemini CLI visual QC for OpenFOAM paper Figure 10.

# Final prompt sent to Gemini CLI

You are doing visual QC for a journal paper figure.

In the current working directory, inspect these two image files:
- ./figure10.png
- ./figure10_pdf_render.png

The first is the direct PNG export. The second is the PDF rendered back to PNG.

Evaluate:
1) layout breakage or panel misalignment
2) title/caption crowding
3) font readability at two-column journal scale
4) clipped annotations or callouts
5) label/legend crowding
6) PNG vs PDF-render mismatch
7) whether the main claim is clear from the figure alone: dense composite rooms can look structurally correct while CFD fidelity remains limited

Return exactly:
VERDICT: PASS or FAIL
- bullet findings
- if FAIL, exact fixes to make
Keep it concise.

# Gemini output (raw)

I will read both image files to compare the direct PNG export and the PDF-rendered version.
VERDICT: PASS

- **Layout & Alignment**: All four panels (a–d) are perfectly aligned in a 2x2 grid. The flow comparison sub-panels in (b) and (d) are symmetrically divided by a clear vertical separator.
- **Font & Readability**: Headers and main axis labels use a bold, high-contrast sans-serif font that will remain legible at two-column journal scale. The white-boxed annotations ("room F1", etc.) use a smaller but clear typeface.
- **PNG vs. PDF Mismatch**: No significant discrepancies found. Transparency in the 3D flow volumes and the dashed vector overlays are correctly preserved in the PDF render.
- **Main Claim Clarity**: The figure successfully demonstrates the claim. Panels (a) and (c) show high "Structural" scores (dark grey bars) contrasting with low "CFD" scores (teal bars), visually proving that structural correctness does not guarantee CFD fidelity.
- **Crowding & Clipping**: Annotations are well-contained within panel boundaries. The orange callout boxes are placed in "dead space" to avoid obscuring geometry, and white backgrounds on labels prevent interference with grid lines.

# Concise summary
TBD

# Action items / next steps
TBD
