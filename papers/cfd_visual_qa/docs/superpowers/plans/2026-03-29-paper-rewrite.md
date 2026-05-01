# CFD-VisQA Paper Full Rewrite Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Rewrite the CFD-VisQA paper from scratch using verified API-isolated data, with proper section depth for a full journal paper (~5,500 words body text).

**Architecture:** LaTeX rewrite of `papers/cfd_visual_qa/latex/main.tex` using elsarticle template. All numerical claims must be verified against the API-isolated JSON files in `benchmark/labels/vlm_eval_*_api_*`. Introduction and Related Work are placeholder-level (to be expanded later with PaperSearch). All other sections are full depth.

**Tech Stack:** LaTeX (elsarticle, natbib), pdflatex build chain, existing matplotlib/VTK figures

---

## Critical Data Verification

**MUST-USE data source:** `benchmark/labels/vlm_eval_{model}_api_t{1,2,3}_{setup,imageonly}.json`

**Verified numbers from API-isolated files (2026-03-29):**

| Evaluator | Condition | T1 | T2 | T3 | Mean | Mean FP | Mean FN |
|-----------|-----------|-----|-----|-----|------|---------|---------|
| Claude Opus 4.6 | Setup | 90.0 | 90.0 | 86.7 | 88.9% | 3.3 | 0.0 |
| Claude Opus 4.6 | Image-only | 33.3 | 33.3 | 33.3 | 33.3% | 0.0 | 20.0 |
| GPT-5.4 | Setup | 56.7 | 70.0 | 76.7 | **67.8%** | 0.7 | 9.0 |
| GPT-5.4 | Image-only | 43.3 | 56.7 | 46.7 | **48.9%** | 0.7 | 14.7 |
| Gemini 3.1 Pro | Setup | 56.7 | 53.3 | 63.3 | 57.8% | 5.3 | 7.3 |
| Gemini 3.1 Pro | Image-only | 63.3 | 56.7 | 40.0 | 53.3% | 4.3 | 9.7 |
| Expert | Setup | 73.8 | — | — | 73.8% | 2.0 | 19.0* |
| Expert | Image-only | 66.7 | — | — | 66.7% | 2.0 | 8.0 |

*Expert FN count needs verification from expert label files.

**CRITICAL CORRECTION:** The previous paper reported GPT-5.4 setup = 80.0% (all trials identical). The actual API-isolated data shows 56.7%/70.0%/76.7% → mean 67.8%. This changes the ranking and narrative.

**Corrected setup-conditioned ranking:**
1. Claude 88.9% (FN=0)
2. Expert 73.8%
3. GPT-5.4 67.8%
4. Gemini 57.8%

**Corrected image-only ranking:**
1. Expert 66.7%
2. Gemini 53.3%
3. GPT-5.4 48.9%
4. Claude 33.3%

**Rank reversal still holds** — Claude goes from #1 to #4, Expert from #2 to #1. But now GPT-5.4 no longer outperforms the expert in setup-conditioned evaluation. The expert is #2 in both conditions, making the expert's stability even more striking.

**Updated delta (pp):**
- Claude: -55.6 (88.9 → 33.3)
- GPT-5.4: -18.9 (67.8 → 48.9) — was -36.7
- Gemini: -4.4 (57.8 → 53.3)
- Expert: -7.1 (73.8 → 66.7)

---

## File Map

| File | Action | Responsibility |
|------|--------|---------------|
| `latex/main.tex` | **REWRITE** | Full paper LaTeX source |
| `latex/main.pdf` | Rebuild | pdflatex output |
| `results/fig_*.png` | Keep/Regenerate | Figures (some need number updates) |
| `scripts/gen_fig_*.py` | Modify | Update hardcoded numbers where needed |

Existing figures to keep as-is (no number changes):
- `fig_benchmark_overview.png` — benchmark structure diagram
- `fig_error_examples.png` — correct vs error flow fields
- `fig_eval_protocol.png` — evaluation protocol flowchart
- `fig_result_bc_swap_illusion.png` — BC swap qualitative example
- `fig_result_setup_vs_imageonly.png` — wrong-viscosity qualitative example
- `fig_result_three_evaluators.png` — mesh comparison example
- `fig_result_shared_fp.png` — shared false positive example

Figures that may need regeneration (contain hardcoded accuracy numbers):
- `fig_result_rank_reversal.png` — GPT-5.4 numbers changed
- `fig_result_main_accuracy.png` — GPT-5.4 numbers changed
- `fig_result_fp_fn_profile.png` — GPT-5.4 FP/FN changed
- `fig_result_trial_consistency.png` — GPT-5.4 trial values changed

---

### Task 1: Verify Expert Baseline Data

**Files:**
- Read: `benchmark/labels/expert_image_only.json`
- Read: `benchmark/labels/vlm_eval_*_api_*.json` (all 18 trial files)

- [ ] **Step 1: Extract and verify all numbers**

```bash
cd papers/cfd_visual_qa
python3 << 'PYEOF'
import json, os

print("=== API-Isolated VLM Results ===")
for model in ['claude', 'gemini', 'gpt5']:
    for cond in ['setup', 'imageonly']:
        for trial in [1, 2, 3]:
            path = f'benchmark/labels/vlm_eval_{model}_api_t{trial}_{cond}.json'
            d = json.load(open(path))
            print(f"{model} t{trial} {cond}: acc={d['accuracy']*100:.1f}% FP={d['fp']} FN={d['fn']} n={d['total']}")

print("\n=== Expert Image-Only ===")
d = json.load(open('benchmark/labels/expert_image_only.json'))
print(f"expert image_only: acc={d['accuracy']*100:.1f}% FP={d['fp']} FN={d['fn']} n={d['total']}")
PYEOF
```

Expected: All numbers match the table above. If any discrepancy, update the table before proceeding.

- [ ] **Step 2: Locate expert setup-conditioned data**

The expert setup-conditioned accuracy (73.8%) must come from a file. Search for it:
```bash
find benchmark/labels/ -name "expert*" -o -name "*expert*" | head -20
```

Verify the 73.8% figure and its FP/FN breakdown.

- [ ] **Step 3: Create a single canonical verified data file**

Write `benchmark/labels/paper_canonical_apiisolated.json` with all verified numbers in one place. This is the single source of truth for the paper.

---

### Task 2: Regenerate Affected Figures

**Files:**
- Modify: `scripts/gen_fig_result_rank_reversal.py`
- Modify: `scripts/gen_fig_result_main_accuracy.py`
- Modify: `scripts/gen_fig_result_fp_fn_profile.py`
- Modify: `scripts/gen_fig_result_trial_consistency.py`
- Output: `results/fig_result_*.png` and `results/fig_result_*.pdf`

- [ ] **Step 1: Read each gen_fig script to find hardcoded GPT-5.4 numbers**

```bash
grep -n "80.0\|GPT\|gpt" scripts/gen_fig_result_rank_reversal.py
grep -n "80.0\|GPT\|gpt" scripts/gen_fig_result_main_accuracy.py
grep -n "80.0\|GPT\|gpt" scripts/gen_fig_result_fp_fn_profile.py
grep -n "80.0\|GPT\|gpt" scripts/gen_fig_result_trial_consistency.py
```

- [ ] **Step 2: Update GPT-5.4 numbers in each script**

Replace old values with verified values:
- Setup accuracy: 56.7, 70.0, 76.7 (mean 67.8)
- Image-only accuracy: 43.3, 56.7, 46.7 (mean 48.9)
- FP/FN: per verified data

- [ ] **Step 3: Regenerate all four figures**

```bash
cd papers/cfd_visual_qa
source ../../.venv/bin/activate  # or appropriate venv
python scripts/gen_fig_result_rank_reversal.py
python scripts/gen_fig_result_main_accuracy.py
python scripts/gen_fig_result_fp_fn_profile.py
python scripts/gen_fig_result_trial_consistency.py
```

- [ ] **Step 4: Visual QC — verify figures show corrected GPT-5.4 numbers**

Read each output PNG and confirm the bars/markers reflect the new values.

---

### Task 3: Write LaTeX — Frontmatter + Abstract

**Files:**
- Modify: `latex/main.tex`

- [ ] **Step 1: Write preamble, author block, and abstract**

The abstract must reflect the corrected numbers. Key changes:
- GPT-5.4 setup: 80.0% → 67.8%
- Ranking: Claude > Expert > GPT-5.4 > Gemini (not Claude > GPT > Expert > Gemini)
- Expert is #2 in BOTH conditions — narrative emphasis shift
- GPT-5.4 drop: 36.7pp → 18.9pp

Abstract target: ~230 words. Core claims:
1. CFD-VisQA benchmark (60 cases, 258 images, 279 questions)
2. API-isolated protocol with contamination discovery
3. Rank reversal: Claude #1 with setup → #4 without
4. Expert is most robust evaluator across conditions
5. Setup-image cross-referencing is model-specific, not general VLM capability

---

### Task 4: Write LaTeX — Introduction (Placeholder ~400w)

**Files:**
- Modify: `latex/main.tex`

- [ ] **Step 1: Write Introduction with placeholder citations**

Structure:
- P1: CFD in automated pipelines, validation bottleneck (~100w)
- P2: VLM potential for visual validation, unknown capability (~100w)
- P3: Gap — no CFD-specific VLM benchmark exists (~80w)
- P4: This paper introduces CFD-VisQA + key findings summary (~80w)
- Contribution list (4 items, same as current but with corrected numbers)

Use `[CITE:keyword]` markers for citations to be filled later with PaperSearch.

---

### Task 5: Write LaTeX — Related Work (Placeholder ~300w)

**Files:**
- Modify: `latex/main.tex`

- [ ] **Step 1: Write Related Work skeleton**

Three subsections with placeholder references:
- 2.1 VLM Benchmarks for Scientific Visuals (~120w)
- 2.2 AI for CFD Result Assessment (~120w)
- 2.3 Research Gap (~60w)

Mark all citations as `[CITE:keyword]`. Keep existing 8 bibitems but note they need Zotero replacement.

---

### Task 6: Write LaTeX — Benchmark Design (~1,200w)

**Files:**
- Modify: `latex/main.tex`

- [ ] **Step 1: Write Section 3.1 Design Philosophy (~150w)**

Three principles: (1) spectrum of visual understanding, (2) correct + erroneous pairs, (3) full reproducibility. Include Fig reference to benchmark overview.

- [ ] **Step 2: Write Section 3.2 Flow Scenarios + Table (~200w)**

Keep existing Table 1 (scenarios). Add text describing:
- Reynolds/Rayleigh number ranges
- Solver selection rationale (simpleFoam vs buoyantBoussinesq)
- 2D steady-state scope and rationale

- [ ] **Step 3: Write Section 3.3 Error Taxonomy (~250w)**

Expand from current ~1 sentence per error type to ~40 words each:
- E1: Under-convergence — stopped at 30 iterations, residuals not converged
- E2: BC swap — inlet/outlet or hot/cold wall exchanged
- E3: Wrong viscosity — 10-100x perturbation, changes Re/Ra
- E4: Wrong turbulence model — laminar solver at turbulent Re
- E5: Coarse mesh — ~10x cell reduction, discretization artifacts
- E8: Gravity/direction reversal — buoyancy direction flipped

Explain the design rationale: errors span a spectrum from visually obvious (E1, E5) to visually subtle but setup-inconsistent (E2, E3).

- [ ] **Step 4: Write Section 3.4 Visualization Standardization (NEW, ~150w)**

Describe:
- Rendering pipeline (PyVista for VTK, matplotlib for 2D)
- Fixed colormap (jet/viridis), fixed resolution (1920x1080)
- Standardized axis labels, colorbar placement
- 6 visualization types: velocity contour, vector field, streamlines, pressure contour, temperature contour, composite overlay

- [ ] **Step 5: Write Section 3.5 Question Hierarchy (NEW, ~150w)**

Five difficulty levels:
- L1: Visual reading ("What is the maximum velocity?")
- L2: Pattern recognition ("Is there a recirculation zone?")
- L3: Physics judgment ("Is this flow physically plausible?")
- L4: Quantitative reasoning ("Is the boundary layer thickness consistent with Re?")
- L5: Cross-field reasoning ("Are velocity and temperature fields mutually consistent?")

Explain that the main evaluation in this paper uses L3 (binary plausibility assessment) as the primary task.

- [ ] **Step 6: Write Section 3.6 Expert Annotation Protocol (NEW, ~150w)**

Describe:
- Single domain expert (thermal-fluid CFD, 10+ years)
- Blind evaluation via Google Drive DOCX
- Setup-conditioned: 80 items evaluated
- Image-only: 30 items (subset), completed separately
- Binary verdict (OK/Anomaly) + free-text reasoning
- Acknowledge single-expert limitation

- [ ] **Step 7: Add Fig and Table references**

Ensure all figure/table refs are correct. Benchmark overview figure + scenarios table + error examples figure.

---

### Task 7: Write LaTeX — Evaluation Protocol as Independent Section (~600w)

**Files:**
- Modify: `latex/main.tex`

- [ ] **Step 1: Write Section 4.1 API-Isolated Design (~200w)**

Describe:
- Base64 image encoding, direct API calls
- Zero filesystem access — no project context, no shared state
- Each item is independent API call
- Models: Claude Opus 4.6 (Anthropic API), GPT-5.4 (OpenAI API), Gemini 3.1 Pro (Google GenAI API)
- Default temperature settings

- [ ] **Step 2: Write Section 4.2 Two Evaluation Conditions (~150w)**

Setup-conditioned:
- Image + problem setup text (geometry, Re/Ra, BCs, solver, expected regime)
- Tests: can VLM cross-reference stated parameters against visual features?

Image-only:
- Image alone, no textual context
- Tests: can VLM assess plausibility from visual features alone?

Explain this is a deliberate ablation, not two separate benchmarks.

- [ ] **Step 3: Write Section 4.3 Multi-Trial Protocol (~100w)**

- 3 independent trials per model per condition
- Default API temperature (non-zero for most models)
- Report per-trial accuracy, mean, and variance
- 30 items per evaluator per condition (common subset)

- [ ] **Step 4: Write Section 4.4 Contamination Discovery (~150w)**

- Subagent evaluation: Claude 99.6% (with filesystem access)
- API-isolated evaluation: Claude 88.9% (zero filesystem access)
- 10.7pp inflation from filesystem contamination
- This discovery motivated the API-isolated protocol
- Establishes API isolation as methodological requirement

Include eval protocol figure reference.

---

### Task 8: Write LaTeX — Results (~800w)

**Files:**
- Modify: `latex/main.tex`

All numbers MUST come from the verified data in Task 1.

- [ ] **Step 1: Write Section 5.1 Overall Performance (~200w)**

Table 2 (corrected): setup-conditioned accuracy across 3 trials.

| Evaluator | T1 | T2 | T3 | Mean | FP | FN |
|-----------|-----|-----|-----|------|-----|-----|
| Claude | 90.0 | 90.0 | 86.7 | 88.9% | 3.3 | 0.0 |
| Expert | 73.8 | — | — | 73.8% | 2 | 19* |
| GPT-5.4 | 56.7 | 70.0 | 76.7 | 67.8% | 0.7 | 9.0 |
| Gemini | 56.7 | 53.3 | 63.3 | 57.8% | 5.3 | 7.3 |

Key text points:
- Claude #1: FN=0 across all trials (catches every error)
- Expert #2: low FP but high FN (misses subtle errors)
- GPT-5.4 #3: large trial-to-trial variance (56.7-76.7%)
- Gemini #4: moderate accuracy, highest FP rate

- [ ] **Step 2: Write Section 5.2 Rank Reversal (~200w)**

Table 3 (corrected): setup vs image-only.

| Evaluator | Setup | Image-only | Delta |
|-----------|-------|------------|-------|
| Claude | 88.9% | 33.3% | -55.6 |
| GPT-5.4 | 67.8% | 48.9% | -18.9 |
| Gemini | 57.8% | 53.3% | -4.4 |
| Expert | 73.8% | 66.7% | -7.1 |

Narrative: Complete rank reversal between conditions. Claude goes from best to worst. Expert is #2 in setup but #1 in image-only — the most robust evaluator. Setup dependency magnitude is diagnostic of assessment strategy.

Reference rank reversal figure.

- [ ] **Step 3: Write Section 5.3 Per-Error-Type Analysis (~150w)**

Reference FP/FN profile figure. Key findings:
- Claude: 100% recall on E3, E4, E5; lower on E2 (BC swap), E8 (gravity)
- Expert: blind spots on E2 (20%), E3 (29%) — subtle errors
- BC swap illusion: mirror image is physically self-consistent

Reference BC swap figure.

- [ ] **Step 4: Write Section 5.4 Trial Consistency (~100w)**

Reference trial consistency figure. Key findings:
- Claude: most consistent (86.7-90% setup, identical 33.3% image-only)
- GPT-5.4: high variance (56.7-76.7% setup) — notable change from old narrative
- Gemini: highest variance both conditions

- [ ] **Step 5: Write Section 5.5 Qualitative Case Studies (~150w)**

Reference setup-vs-imageonly figure, three-evaluators figure, shared-FP figure.
Brief narrative for each: wrong-viscosity mechanism, coarse-mesh detection, complex-flow false positive.

---

### Task 9: Write LaTeX — Discussion (~900w)

**Files:**
- Modify: `latex/main.tex`

- [ ] **Step 1: Write Section 6.1 Setup Text as Critical Variable (~300w)**

Core argument:
- Rank reversal demonstrates VLMs don't "understand" CFD physics
- Claude performs setup-image cross-referencing (explicit parameter checking)
- Expert uses gestalt assessment (internalized physical intuition)
- Expert is robust to setup removal because strategy doesn't depend on text
- Expert vulnerable to subtle errors (BC swap, wrong viscosity) that are qualitatively correct
- Updated: GPT-5.4 shows moderate cross-referencing but with high variance

- [ ] **Step 2: Write Section 6.2 Model-Specific Capabilities (~200w)**

Three profiles:
- Claude: strongest cross-referencer, FN=0 with setup, worst without (33.3%)
- GPT-5.4: moderate capability, high trial variance suggests less deterministic reasoning
- Gemini: weakest setup dependency (-4.4pp), suggests limited cross-referencing; image-only performance second only to expert

- [ ] **Step 3: Write Section 6.3 Implications for CFD Automation (~150w)**

Four practical points:
1. Setup text is mandatory for VLM-based validation
2. Model selection matters: only Claude achieves FN=0
3. False positives are the trade-off: flag for expert review
4. API isolation is a protocol requirement for benchmarking

- [ ] **Step 4: Write Section 6.4 Limitations (~250w)**

Expanded from current 1 paragraph to thorough treatment:
- Sample size: 30 items for cross-model comparison
- Single expert baseline — multi-expert needed
- 2D steady-state RANS only — 3D, transient, multi-physics not tested
- Evaluation asymmetry: expert has 80+30 items, VLMs have 30+30
- Default temperature settings — deterministic decoding not tested
- No open-source VLMs tested (LLaVA, Qwen2.5-VL, InternVL)
- Error taxonomy is controlled/templated — real-world errors may differ
- Setup-conditioned task rewards text-image consistency checking more than deep physical understanding (acknowledge this explicitly)

---

### Task 10: Write LaTeX — Conclusion + Back Matter (~300w)

**Files:**
- Modify: `latex/main.tex`

- [ ] **Step 1: Write Conclusion (~200w)**

Three paragraphs:
- P1: What we did (benchmark + protocol + evaluation)
- P2: What we found (rank reversal, corrected numbers)
- P3: Practical implication (Claude + setup text as first-pass validator, expert review for FP)

- [ ] **Step 2: Write Data Availability + Declarations**

Keep existing text, verify GitHub URL.

- [ ] **Step 3: Write Bibliography (placeholder bibitems)**

Keep existing 8 bibitems. Add `[CITE:keyword]` markers in text for future PaperSearch expansion.

---

### Task 11: Build and Verify PDF

**Files:**
- Build: `latex/main.tex` → `latex/main.pdf`

- [ ] **Step 1: Run pdflatex twice**

```bash
cd papers/cfd_visual_qa/latex
pdflatex -interaction=nonstopmode main.tex
pdflatex -interaction=nonstopmode main.tex
```

Expected: No errors, PDF generated.

- [ ] **Step 2: Check page count**

Target: 16-20 pages (elsarticle preprint 12pt with 11 figures).

- [ ] **Step 3: Verify all figures render**

Read the PDF and confirm all 11 figures and 3 tables are present and correctly placed.

- [ ] **Step 4: Verify no stale numbers remain**

Search for old incorrect values in the tex source:
```bash
grep -n "80\.0" latex/main.tex    # old GPT-5.4 number
grep -n "36\.7" latex/main.tex    # old GPT-5.4 delta
grep -n "99\.6" latex/main.tex    # should only appear in contamination section
```

- [ ] **Step 5: Send PDF to user via Telegram**

---

### Task 12: Word Count Verification

- [ ] **Step 1: Count words per section**

Run the word-counting script from earlier against the new main.tex. Verify:
- Total body text: ~5,000-5,500 words
- Benchmark Design: ~1,200 words
- Evaluation Protocol: ~600 words
- Results: ~800 words
- Discussion: ~900 words

If any section is significantly under target, flag for expansion.

---

## Execution Notes

- **Tasks 1-2** must complete before Task 3 (verified data + updated figures needed for writing)
- **Tasks 3-10** are sequential (writing order matters for narrative flow)
- **Tasks 11-12** are post-writing verification
- All LaTeX writing goes into a single `main.tex` file (overwrite existing)
- Introduction and Related Work are placeholder-level — will be expanded later with PaperSearch in a research-library-accessible environment
