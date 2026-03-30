# CFD Visual QA — Paper Outline

> Date: 2026-03-26
> Target: ~6,000-8,000 words
> Target venue: TBD (Computers & Fluids, Building and Environment, or NeurIPS Datasets)

---

## Title (Working)

**"Can Vision-Language Models Evaluate CFD? A Benchmark for Physical Plausibility Assessment of Flow Field Visualizations"**

Alternative:
- "CFD-VisQA: Benchmarking VLM Understanding of Computational Fluid Dynamics Visualizations"
- "Beyond Charts: Evaluating Vision-Language Models on Scientific Flow Field Imagery"

---

## Abstract (~250 words)

- Motivation: CFD automation pipelines need automated result validation
- Gap: No benchmark exists for VLM evaluation on CFD flow field images
- Contribution: CFD-VisQA benchmark (10 scenarios, 60 cases, 258 images, 279 QA)
- Method: Correct + intentionally erroneous flow fields, multi-level questions, blind evaluation
- Results: Claude Opus 100%, Gemini 80% (valid), Expert 87% on 30-item pilot
- Key finding: VLMs detect subtle errors experts miss (setup-conditioned cross-referencing)
- Implication: VLMs can serve as automated CFD validators, but with different blind spots

---

## 1. Introduction (~800 words)

### 1.1 Motivation
- CFD automation growing (image-to-CFD pipelines, AI-assisted simulation)
- Result validation remains manual bottleneck — expert visually inspects flow fields
- VLMs show promise in scientific image understanding but untested on CFD

### 1.2 Problem Statement
- Can current frontier VLMs assess physical plausibility of CFD flow field visualizations?
- What types of errors can they detect vs miss?
- How do they compare to domain experts?

### 1.3 Contributions
1. **CFD-VisQA benchmark**: First dedicated dataset for VLM evaluation on CFD visualizations
   - 10 thermal-fluid scenarios × 6 cases (2 correct + 4 error) = 60 CFD cases
   - 258 images across 5 visualization types
   - 279 QA questions at 5 difficulty levels
   - Expert annotations (35+ labeled items)
2. **Systematic error taxonomy**: 8 error types with 3 severity levels
3. **Blind evaluation protocol**: Setup-conditioned vs image-only tasks
4. **Comparative analysis**: Claude Opus 4.6, Gemini 3.1, domain expert
5. **Key findings**: Contamination paradox, differential blind spots, setup-conditioned advantage

---

## 2. Related Work (~600 words)

### 2.1 VLM Benchmarks for Scientific Visuals
- MMMU, SciFIBench, PhysBench — no CFD coverage
- ClimateIQA — closest analog (meteorological heatmaps)

### 2.2 AI for CFD
- CFDLLMBench — text/code only, no visual
- ML4CFD — surrogate models, no VLM
- MDPI VLM-CFD (2026) — only prior work, small scale

### 2.3 Flow Visualization Understanding
- Classical CV: VortexNet, Shock-Net — no language interface
- Kashefi (2024) — generative AI fails on fluid mechanics imagery

### 2.4 Gap
- No benchmark for VLM evaluation on CFD output images with physics-grounded QA

---

## 3. Benchmark Design (~1,500 words)

### 3.1 Design Philosophy
- 3-axis: scenario × visualization × question level
- Correct + error pairs for each scenario
- Reproducible (OpenFOAM, scripts public)

### 3.2 Flow Scenarios (Table)
| ID | Scenario | Physics | Solver |
|----|----------|---------|--------|
| S1 | Heated vertical plate | Natural convection | buoyantBoussinesqSimpleFoam |
| S2 | Channel flow (laminar) | Forced convection | simpleFoam |
| S3 | Channel flow (turbulent) | Forced convection | simpleFoam + k-ω SST |
| S4 | Backward-facing step | Separation/reattachment | simpleFoam |
| S5 | Lid-driven cavity | Recirculation | simpleFoam |
| S6 | Differentially heated cavity | Natural convection | buoyantBoussinesqSimpleFoam |
| S7 | Mixed convection channel | Mixed convection | buoyantBoussinesqSimpleFoam |
| S8 | Heated obstacle | Wake + thermal | buoyantBoussinesqSimpleFoam |
| S9 | Rayleigh-Bénard | Convection cells | buoyantBoussinesqSimpleFoam |
| S10 | Ventilated room | Indoor airflow | simpleFoam |

### 3.3 Error Types (Table)
| ID | Error Type | Generation Method | Severity |
|----|-----------|------------------|----------|
| E1 | Under-converged | Early stop (30 iters) | Moderate-Severe |
| E2 | BC swap | Inlet↔outlet or hot↔cold | Subtle-Severe |
| E3 | Wrong viscosity | 10-100× increase | Subtle |
| E4 | Wrong turbulence model | Laminar at turbulent Re | Subtle-Moderate |
| E5 | Coarse mesh | 10×reduced resolution | Subtle-Moderate |
| E6 | Artificial modification | Post-solve field edit | Subtle |
| E7 | Symmetry violation | (reserved) | — |
| E8 | Gravity/direction flip | Reverse buoyancy | Moderate-Severe |

### 3.4 Visualization Types
- V1: Velocity magnitude contour
- V2: Temperature contour
- V3: Pressure contour
- V4: Velocity vector field
- V5: Streamlines
- V6: Composite (contour + streamlines)

### 3.5 Question Levels
- L1: Visual reading (colorbar values)
- L2: Qualitative pattern (recirculation presence)
- L3: Physical reasoning (plausibility judgment) ← main benchmark task
- L4: Quantitative estimation (BL thickness)
- L5: Comparative judgment (vs reference)

### 3.6 Annotation Schema
- Multi-layer labels: binary (correct/error) + error type + severity + region + explanation
- Expert + VLM labels in same format for direct comparison

---

## 4. Evaluation Methodology (~800 words)

### 4.1 Task Definition
- **Setup-conditioned**: Problem setup text + image → plausibility judgment
- **Image-only**: Image alone → anomaly detection (harder variant)

### 4.2 Blind Protocol
- Anonymized case codes (CFD-XXXX)
- No ground truth in evaluation context
- Independent evaluation sessions (subagent or CLI)

### 4.3 Models Evaluated
- Claude Opus 4.6 (Anthropic)
- Gemini 3.1 (Google)
- Domain expert (thermal-fluid PhD)

### 4.4 Metrics
- Overall accuracy
- Per-error-type recall
- Per-scenario accuracy
- False positive / false negative rates
- No-response rate (API reliability)

---

## 5. Results (~1,200 words)

### 5.1 Overall Performance (Table — main result)

| Evaluator | Accuracy | FN Rate | FP Rate | No-Response |
|-----------|----------|---------|---------|-------------|
| Claude Opus 4.6 | 100% (30/30) | 0% | 0% | 0% |
| Gemini 3.1 | 80% (16/20 valid) | 20% | 0% | 33% |
| Expert | 87% (26/30) | 13% | 0% | 0% |

### 5.2 Per-Error-Type Analysis (Table)

| Error Type | Expert | Claude | Gemini |
|-----------|--------|--------|--------|
| E1 Under-converged | 100% | 100% | ~100% |
| E2 BC swap | 33% | 100% | ~50% |
| E3 Wrong viscosity | 50% | 100% | ~67% |
| E4 Wrong turb model | 50% | 100% | ~100% |
| E5 Coarse mesh | 67% | 100% | ~50% |
| E8 Gravity/dir flip | 100% | 100% | ~80% |

### 5.3 Contamination Paradox
- Claude with case names (contaminated): 87%
- Claude blind: 100%
- Hypothesis: Setup-text cross-referencing is more effective without name-based priors

### 5.4 Differential Blind Spots
- Expert misses: BC swap, wrong viscosity — gestalt judgment
- Gemini misses: coarse mesh, some gravity flips — image processing issues
- Claude: no misses at current sample size — systematic cross-referencing

### 5.5 Setup-Conditioned vs Image-Only
- (Planned: run image-only variant and compare)

---

## 6. Discussion (~1,000 words)

### 6.1 Why VLMs Outperform on Subtle Errors
- Systematic: reads colorbar values, compares to stated BCs
- Cross-referential: setup text provides checkable constraints
- No fatigue: consistent attention to detail across all items

### 6.2 Why Experts Outperform on Certain Tasks
- Holistic: "feels wrong" before articulating why
- Domain transfer: applies knowledge from similar problems
- (Not observed in current data — expert may outperform on larger/noisier datasets)

### 6.3 Limitations
- Sample size (30 items — "initial finding" level)
- Single expert
- Two VLMs (GPT-5.4 planned but CLI issue)
- 2D only (no 3D rendering evaluation)
- Steady-state only (no transient)

### 6.4 Implications for CFD Automation
- VLMs as automated validators in CFD pipelines
- Setup-conditioned prompting essential (not image-only)
- Ensemble approach: VLM + expert for coverage of different error types

---

## 7. Conclusion (~300 words)

- First benchmark for VLM evaluation on CFD flow field visualizations
- 10 scenarios, 60 cases, 258 images, 279 QA — publicly available
- VLMs (Claude Opus) achieve 100% on 30-item blind evaluation
- Key insight: setup-conditioned cross-referencing detects subtle errors experts miss
- Future: scale to 3D, transient, multi-physics; fine-tuned domain VLMs

---

## Figures (Planned)

1. **Benchmark overview**: 3-axis diagram (scenario × viz × question level)
2. **Error type examples**: 2×4 grid (correct vs error for 4 error types)
3. **Evaluation protocol**: Blind vs contaminated workflow diagram
4. **Main results**: Bar chart — accuracy by evaluator
5. **Per-error heatmap**: Error type × evaluator detection rate
6. **Example VLM reasoning**: Side-by-side correct vs error with VLM explanations
7. **Contamination paradox**: Bar chart showing blind > contaminated

## Tables

1. Scenario summary (ID, name, physics, solver, cases)
2. Error type taxonomy (ID, generation method, severity)
3. Main results (evaluator × metrics)
4. Per-error recall (error type × evaluator)
5. Dataset statistics (total cases, images, QA by scenario)
