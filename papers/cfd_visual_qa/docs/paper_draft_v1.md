# Can Vision-Language Models Evaluate CFD? A Benchmark for Physical Plausibility Assessment of Flow Field Visualizations

## Abstract

Validating computational fluid dynamics (CFD) simulation results remains a manual bottleneck in automated pipelines, requiring expert visual inspection of flow field images. We introduce CFD-VisQA, the first benchmark for evaluating vision-language model (VLM) capabilities in assessing physical plausibility of CFD visualizations. The benchmark comprises 10 canonical thermal-fluid scenarios with 60 OpenFOAM cases, 258 flow field images, and 279 physics-grounded questions at 5 difficulty levels. Each scenario includes correct solutions and systematically generated errors spanning 6 types at varying severity levels. Using a blind, setup-conditioned protocol with anonymized case identifiers, we evaluate Claude Opus 4.6 (75 items, 3 independent trials), Gemini 3.1 (30 items, 3 trials), and a domain expert (80 items). Claude Opus achieves 99.6% accuracy across 225 judgments (majority vote: 100%), compared to 73.8% for the domain expert and 87.5% mean accuracy for Gemini 3.1. The expert primarily misses "subtle" errors that produce plausible-looking flows (boundary condition swaps: 20% recall, wrong viscosity: 29% recall), while Claude detects these through systematic cross-referencing of setup text against visual features. We document a contamination paradox where providing case-name metadata decreases VLM accuracy (87% vs. 100% blind), and show that per-item consistency across 3 trials reaches 98.7% for Claude. These findings suggest that frontier VLMs show promise as automated CFD validators with error detection patterns complementary to human experts, though evaluation scope and model coverage require further expansion. The benchmark dataset and evaluation pipeline are publicly available.

---

# 1. Introduction

Computational fluid dynamics (CFD) is increasingly integrated into automated design and analysis pipelines, from building ventilation assessment to industrial heat exchanger optimization. Recent advances in vision-language models (VLMs) have enabled new paradigms such as image-to-CFD frameworks, where a VLM extracts geometric information from architectural images to automatically generate CFD simulation domains [Joo 2026]. As these automated pipelines proliferate, a critical bottleneck emerges: *who validates the results?*

In traditional CFD workflows, an experienced engineer inspects velocity contours, temperature fields, and streamline plots to confirm that the solution is physically reasonable—checking for proper boundary layer development, expected recirculation patterns, and absence of numerical artifacts. This visual inspection requires domain expertise accumulated over years and cannot scale with the throughput of automated pipelines. If a VLM could perform this validation step—assessing whether a CFD flow field visualization is physically plausible—the automation loop would close entirely.

Yet it remains unknown whether current VLMs possess the physical reasoning capability to evaluate CFD output images. Existing VLM benchmarks test scientific knowledge through multiple-choice questions on textbook diagrams [Yue 2024 MMMU], chart comprehension [Pandey 2025], or everyday physics understanding [Chow 2025 PhysBench], but none evaluate the ability to assess *simulation output images* against known physical laws. The closest work evaluates VLMs on meteorological heatmaps [Chen 2024 ClimateIQA] or tests generative AI on fluid mechanics prompts [Kashefi 2024], but neither addresses the specific task of CFD result validation.

This paper introduces **CFD-VisQA**, the first benchmark designed to evaluate VLM capabilities in assessing physical plausibility of CFD flow field visualizations. The benchmark comprises 10 canonical thermal-fluid scenarios (natural convection, forced convection, mixed convection), each with correct solutions and systematically generated erroneous solutions spanning 6 error types. We evaluate two frontier VLMs (Claude Opus 4.6 and Gemini 3.1) against a domain expert using a blind, multi-trial protocol on 75–80 evaluation items, revealing distinct error detection patterns across evaluator types.

Our contributions are:

1. **CFD-VisQA benchmark dataset**: 60 OpenFOAM CFD cases across 10 thermal-fluid scenarios, producing 258 flow field images with 279 physics-grounded questions and expert annotations. The dataset includes both correct solutions and 6 types of intentionally generated errors at varying severity levels.

2. **Setup-conditioned evaluation protocol**: We demonstrate that providing problem setup text alongside the image (boundary conditions, Reynolds number, solver type) dramatically improves VLM performance compared to image-only evaluation, and introduce a blind evaluation protocol that prevents bias from case-name leakage.

3. **Comparative analysis revealing differential blind spots**: Claude Opus achieves 98.7% accuracy on 75 items (99.6% across 3 trials, majority vote 100%), the domain expert 73.8% on 80 items, and Gemini 3.1 87.5% mean on 30 items. Each evaluator misses different error types, suggesting complementary strengths.

4. **Multi-trial consistency and contamination paradox**: Three independent blind trials show 98.7% per-item consistency for Claude. Providing case-name metadata decreases VLM accuracy (87% contaminated vs. 100% blind), demonstrating that evaluation protocol design materially affects results.

# 2. Related Work

## 2.1 VLM Benchmarks for Scientific Visuals

Several benchmarks evaluate VLM performance on scientific and technical imagery, though none address CFD specifically. MMMU [Yue 2024] tests college-level reasoning across 30 disciplines with 11,500 questions, including engineering diagrams, but contains no flow field images. SciFIBench [Roberts 2024] evaluates 28 large multimodal models on 2,000 scientific figure interpretation tasks from arXiv papers, finding that VLM capabilities in scientific figure understanding are "not well characterised." PhysBench [Chow 2025] tests 75 VLMs on physical world understanding across 10,002 entries, finding that the best model (GPT-4o) achieves only 49.5% vs. 95.9% human accuracy; however, its "fluid dynamics" category covers everyday scenarios (water pouring) rather than simulation outputs.

Most relevant to our work is ClimateIQA [Chen 2024], which constructs 26,280 meteorological heatmap images with 762,120 instruction samples for VLM evaluation. Like CFD contour plots, meteorological heatmaps require colormap interpretation and spatial pattern recognition. ClimateIQA finds that GPT-4o and Qwen-VL "struggle with precise color identification and spatial localization"—a finding directly applicable to CFD visualization understanding. We extend this line of work from meteorological to engineering simulation imagery, with the additional challenge of physics-grounded plausibility assessment.

Pandey and Ottley [2025] benchmark VLMs on standardized visualization literacy tests, finding accuracy drops from 76–96% on simple line charts to 18–61% on multi-encoding visualizations and 25–30% on anomaly detection. CFD colormapped fields are inherently multi-encoding visualizations, suggesting VLMs will struggle with quantitative reading of flow field images.

## 2.2 AI for CFD Result Assessment

The intersection of AI and CFD has been explored primarily through surrogate modeling and simulation automation, not result validation. CFDLLMBench [Somasekharan 2025] evaluates LLMs on CFD knowledge (92% on conceptual MCQs) and OpenFOAM case generation (25–34% success), but includes no visual component. A fine-tuned LLM for OpenFOAM automation [arXiv 2504.09602] achieves 88.7% solution accuracy on case setup but does not evaluate simulation results.

The only published work directly evaluating VLMs on CFD visualizations is a 2026 study on mixed-reality indoor CFD images [MDPI Technologies 2026], which fine-tunes Qwen2.5-VL from <30% to >60% accuracy on domain-specific questions. However, this work uses a small dataset, a single VLM architecture, and focuses on mixed-reality overlays rather than standard CFD post-processing output.

Kashefi [2024] demonstrates that generative AI systems (Midjourney, DALL-E, Gemini) produce physically incorrect fluid mechanics imagery when prompted with classical phenomena (von Kármán vortex streets, Kelvin-Helmholtz instabilities), attributing this to limited training data from copyright-protected journal figures. This negative result motivates the need for a public CFD visualization benchmark.

## 2.3 Flow Visualization with Deep Learning

Classical deep learning approaches to flow visualization focus on specific feature detection tasks: VortexNet and Vortex-U-Net for vortex identification [Liu 2022 survey], Shock-Net for shock wave detection, and CNNs for Reynolds number classification from flow images [arXiv 2305.11884]. These methods produce categorical outputs (vortex/no-vortex, Re class) without natural language understanding or explanatory capability.

Physics-informed computer vision [Banerjee 2024] reviews 250+ papers on incorporating physical laws into vision tasks, covering turbulent flow super-resolution and PINN-based reconstruction. However, this literature predates the VLM era and does not address language-based understanding or natural language queries about flow field images.

## 2.4 Research Gap

No existing benchmark evaluates VLM capability to assess physical plausibility of CFD simulation output images. The gap sits at the intersection of scientific VLM benchmarks (which lack CFD imagery), CFD AI tools (which lack visual understanding), and flow visualization ML (which lacks language interfaces). CFD-VisQA fills this gap with a dedicated benchmark combining physics-grounded questions, systematically generated errors, and a blind evaluation protocol enabling fair comparison between VLMs and domain experts.



# 3. Benchmark Design

## 3.1 Design Philosophy

The CFD-VisQA benchmark is built on three principles. First, it must measure a *spectrum* of visual understanding capabilities, from basic colorbar reading to physics-grounded plausibility assessment. Second, it requires both correct and intentionally erroneous flow fields—evaluating only correct cases cannot reveal whether a model detects anomalies or merely describes what it sees. Third, all data must be reproducible: every flow field is generated from OpenFOAM simulations with published scripts, enabling independent verification and extension.

The benchmark is organized along three orthogonal axes: flow scenarios (10 canonical thermal-fluid problems), visualization types (6 rendering modes), and question difficulty levels (5 tiers from visual reading to quantitative estimation). This factored design enables fine-grained analysis of where VLM capabilities break down.

## 3.2 Flow Scenarios

We select ten canonical thermal-fluid scenarios spanning natural convection, forced convection, and mixed convection regimes [Table:scenarios]. Each scenario represents a well-characterized problem with known analytical or experimental reference solutions, ensuring unambiguous ground truth for the correct cases.

The scenarios range from fundamental flows (laminar channel, lid-driven cavity) to applied configurations (ventilated room, heated obstacle). Reynolds numbers span $10^2$–$10^4$ and Rayleigh numbers span $10^3$–$3 \times 10^7$, covering laminar and turbulent regimes. All simulations are two-dimensional, using OpenFOAM's `simpleFoam` for isothermal flows and `buoyantBoussinesqSimpleFoam` for buoyancy-driven flows with the Boussinesq approximation.

For each scenario, we generate two correct cases (baseline and a variant, e.g., different Re or finer mesh) and four error cases, yielding 60 CFD simulations in total.

## 3.3 Error Taxonomy

A key contribution of this benchmark is the systematic generation of physically incorrect flow fields. We define six error types [Table:errors], each produced by a controlled perturbation to the correct simulation setup:

**E1 — Under-convergence.** The solver is stopped after only 30 iterations (vs. thousands for convergence). This produces partially developed flow patterns with residual artifacts, representing a common failure mode in automated CFD pipelines.

**E2 — Boundary condition swap.** Inlet and outlet (or hot and cold wall) boundary conditions are exchanged. In some scenarios (e.g., backward-facing step), this produces a subtly different but plausible-looking flow; in others (e.g., lid-driven cavity with lid on bottom), the error is visually obvious.

**E3 — Wrong viscosity.** Kinematic viscosity is increased by 10–100×, effectively reducing the Reynolds or Rayleigh number well below the intended value. The resulting flow is physically self-consistent but incorrect for the stated conditions—a "subtle" error detectable only by quantitative comparison with the problem setup.

**E4 — Wrong turbulence model.** A laminar solver is applied to a turbulent-regime flow (Re = 10,000). The result lacks turbulent mixing, producing unrealistic velocity profiles, but may appear superficially smooth.

**E5 — Coarse mesh.** Grid resolution is reduced by approximately 10× in each direction, eliminating boundary layer resolution and fine-scale features. Detectability depends on whether the flow has small-scale structures (e.g., corner vortices in lid-driven cavity, convection cells in Rayleigh-Bénard).

**E8 — Gravity/direction reversal.** The gravity vector is flipped (natural convection) or the lid direction is reversed (cavity flow). This inverts the buoyancy-driven circulation, producing a flow that violates physical expectations for the stated boundary conditions.

Each error type is assigned a severity level based on pilot expert labeling: *severe* (detectable at a glance), *moderate* (requires careful inspection), or *subtle* (requires quantitative comparison with the problem setup). Severity varies by scenario—e.g., E2 (BC swap) is severe in lid-driven cavity but subtle in differentially heated cavity.

## 3.4 Visualization Pipeline

All images are generated by a standardized rendering pipeline to eliminate visualization artifacts as a confounding variable. Flow fields are sampled on a uniform 2D grid using PyVista's point interpolation, then rendered with Matplotlib at 300 DPI resolution.

Six visualization types are produced:
- **V1** — Velocity magnitude contour (colormap: viridis)
- **V2** — Temperature contour (colormap: coolwarm) [thermal scenarios only]
- **V3** — Pressure contour (colormap: RdBu_r)
- **V4** — Velocity vector field (quiver overlay)
- **V5** — Streamlines (colored by speed)
- **V6** — Composite: velocity contour + white streamlines

All images include labeled axes (SI units), colorbar with physical quantity and units, and domain outline. Font, colormap, and layout are fixed across all cases to prevent evaluators from using rendering style as a cue.

## 3.5 Question Design

Questions are organized in five difficulty levels:

- **L1 (Visual Reading):** Direct extraction from the image, e.g., "What is the approximate maximum velocity?"
- **L2 (Qualitative Pattern):** Pattern identification, e.g., "Is a recirculation zone present?"
- **L3 (Physical Reasoning):** Plausibility assessment, e.g., "Is this flow field physically plausible for the stated conditions?" — This is the **primary benchmark task**.
- **L4 (Quantitative Estimation):** Numerical estimation, e.g., "Estimate the boundary layer thickness."
- **L5 (Comparative Judgment):** Cross-image comparison, e.g., "How does this differ from the reference case?"

The benchmark focuses on L3 questions for the main evaluation, as these most directly test whether VLMs can serve as automated CFD validators. L1–L2 questions establish baseline visual capability, while L4–L5 questions probe deeper reasoning.

# 4. Evaluation Methodology

## 4.1 Task Definition

We define two evaluation tasks with differing information levels:

**Setup-Conditioned Visual QA (primary).** The evaluator receives both a problem setup description and the flow field image. The setup specifies geometry, boundary conditions (including heated wall positions, inlet/outlet locations, lid direction), Reynolds/Rayleigh number, and solver type. The evaluator judges whether the visualized flow field is physically plausible given these conditions. This reflects the practical use case: a CFD pipeline generates both the setup and the result, and the validator has access to both.

**Image-Only Visual QA (secondary).** The evaluator receives only the image, without problem setup information, and judges whether any non-physical patterns are present. This is a harder task requiring the evaluator to infer the problem type from visual cues alone.

## 4.2 Blind Evaluation Protocol

To prevent evaluation bias, we implement a strict blind protocol:

1. **Anonymized identifiers.** Case names (which may encode error types, e.g., "E1_underconverged") are replaced with opaque blind codes (e.g., "CFD-3E91").
2. **No ground truth exposure.** Evaluators have no access to correct/error labels during evaluation.
3. **Independent context.** Each VLM evaluation runs in a separate session (subagent or CLI invocation) with no memory of previous items or results.
4. **Standardized prompts.** All evaluators receive identical problem setup text and the same question format.

We discovered the importance of this protocol empirically: an initial non-blind evaluation of Claude Opus (with case names visible in file paths and ground truth loaded in context) achieved 87% accuracy, while the subsequent blind evaluation achieved 100%—a phenomenon we term the *contamination paradox* (Section 5.3).

## 4.3 Models and Expert

We evaluate three frontier VLMs and one domain expert:

- **Claude Opus 4.6** (Anthropic): Evaluated via independent subagent invocations with blind codes. Each item evaluated in a fresh context with no prior results.
- **Gemini 3.1** (Google): Evaluated via Gemini CLI in headless mode, with image paths provided in the prompt text. The model reads images via its file access capability.
- **Domain expert**: A thermal-fluid dynamics researcher with expertise in CFD simulation and post-processing. Labeling conducted via structured DOCX forms with embedded images and problem setup descriptions.

## 4.4 Metrics

We report the following metrics for each evaluator:

- **Overall accuracy**: Fraction of items correctly classified as plausible (correct cases → OK) or implausible (error cases → Anomaly).
- **Per-error-type recall**: Detection rate for each error type (E1–E8), measuring which error categories each evaluator can and cannot detect.
- **False negative rate**: Fraction of error cases judged as plausible (missed errors—the critical failure mode for a CFD validator).
- **False positive rate**: Fraction of correct cases judged as implausible.
- **No-response rate**: Fraction of items where the evaluator failed to produce a valid response (relevant for API-based VLM evaluation).



# 5. Results

## 5.1 Overall Performance

[Table:main-results] summarizes the overall accuracy of each evaluator. Claude Opus 4.6 was evaluated on 75 blind items, Gemini 3.1 on 29 valid items (from 30 attempted, 1 no-response), and the domain expert on 80 items. All items span the full 10 scenarios, 6 error types, and multiple visualization modes.

Claude Opus 4.6 achieves 98.7% accuracy (74/75), with a single false positive on one correct case (S8 heated obstacle, whose complex multi-vortex wake was judged as non-physical). The domain expert achieves 73.8% accuracy (59/80), with 19 false negatives and 2 false positives. Gemini 3.1 achieves 86.2% accuracy on valid responses (25/29), with 2 false positives and 2 false negatives.

The critical failure mode for a CFD validator is the false negative: accepting an erroneous flow field as correct. On this metric, Claude achieves a 0% false negative rate (all errors detected), Gemini 6.9% (2/29), and the expert 23.8% (19/80). Claude's single error is a false positive—flagging a physically correct but visually complex flow as suspicious—which is the *safer* failure mode for a validator.

Notably, all three evaluators produced false positives on the same category of case: complex wake/convection patterns with multiple vortices that appear "excessive" despite being physically correct. This suggests an inherent difficulty boundary in the benchmark where visual complexity triggers false alarms regardless of evaluator type.

## 5.2 Per-Error-Type Analysis

[Table:per-error-recall] reveals that error detection capability is strongly dependent on error type, and different evaluators have different blind spots.

**Consistently detectable (recall ≥ 80% across all evaluators):**
- *E1 Under-convergence*: Expert 73%, Claude 100%, Gemini ~100%. Produces visually obvious artifacts (non-physical velocity peaks, incomplete flow development), though detectability decreases in complex geometries (e.g., expert missed S8 under-convergence).
- *E8 Gravity/direction reversal*: Expert 100%, Claude 100%, Gemini ~80%. Inverts expected flow patterns. Detectable through physics reasoning about buoyancy direction.
- *E6 Artificial modification*: Expert 100%, Claude 100%. Suppressed recirculation behind backward-facing step correctly identified.

**Expert-challenging, VLM-detectable:**
- *E2 Boundary condition swap*: Expert 20%, Claude 100%, Gemini ~50%. The most consistent expert blind spot. BC swaps that produce mirror images (hot↔cold wall) appear plausible on visual inspection. Claude detects all BC swaps by cross-referencing stated wall temperatures against colorbar values.
- *E3 Wrong viscosity*: Expert 29%, Claude 100%, Gemini ~67%. Increased viscosity produces qualitatively correct but quantitatively wrong flows. Claude detects via development length and thermal penetration inconsistencies with stated Re/Ra.

**Context-dependent:**
- *E5 Coarse mesh*: Expert 44%, Claude 100%, Gemini ~50%. Detection depends on whether fine-scale features (corner vortices, convection cells) are expected. In scenarios without such features, coarse mesh may appear acceptable.
- *E4 Wrong turbulence model*: Expert 50%, Claude 100%, Gemini ~100%. Laminar solver at turbulent Re produces smooth flows that experts may judge as "reasonable."
- *E3 Wrong viscosity*: The expert misses this in laminar channel flow (S2) where the flow looks qualitatively correct at a different Re, but detects it in natural convection (S1) where the boundary layer is completely suppressed. Claude detects both by noting quantitative inconsistencies (development length too short for stated Re, thermal penetration too deep).

**Context-dependent detectability:**
- *E5 Coarse mesh*: Detection depends on whether the scenario has fine-scale structures that become visibly degraded. In lid-driven cavity (S5), the expert misses it ("looks pretty"), but in heated plate (S1) and Rayleigh-Bénard (S9), coarse mesh artifacts are more apparent. Claude detects all except the S5 case in the initial 15-item evaluation but detects it in the expanded 30-item evaluation.
- *E4 Wrong turbulence model*: Using a laminar solver at Re=10,000 produces smooth, seemingly plausible flow that the expert judges as "reasonable-looking." Claude detects it by reasoning that a flat turbulent profile from a laminar solver is non-physical at this Re.

## 5.3 Multi-Trial Consistency

To assess the reliability of VLM judgments, we conducted three independent blind evaluations of Claude Opus on all 75 items. Each trial used fresh blind codes, independent evaluation sessions, and no access to previous trial results.

Claude achieved 98.7% (74/75) on Trial 1, 100% (75/75) on Trial 2, and 100% (75/75) on Trial 3. Across 225 total judgments, 224 were correct (99.6%). The single error — a false positive on S8's complex multi-vortex wake — appeared only in Trial 1 and was not reproduced in Trials 2 or 3, confirming it as a stochastic single-trial event (probability ≈ 0.4%).

Using majority vote (≥2/3 correct), Claude achieves 100% accuracy on all 75 items. Per-item consistency shows 74/75 items (98.7%) with perfect 3/3 agreement, and 1 item with 2/3 agreement. No false negatives appeared in any trial, indicating that Claude's error detection capability is highly reliable.

## 5.4 The Contamination Paradox

An unexpected finding emerged from our protocol development. An initial, non-blind evaluation of Claude Opus—where case names (e.g., "S5_E1_underconverged") were visible in file paths and ground truth labels were loaded in the evaluation context—achieved 87% accuracy with 3 false negatives (S5_E5, S6_E2, S2_E3). The subsequent blind evaluation, with anonymized codes and no ground truth, achieved 100%.

We hypothesize two mechanisms: (1) case names containing error type labels ("underconverged," "bc_swap") may have triggered *anchoring bias*, causing the model to match its response to the expected label rather than performing independent visual analysis; (2) without name-based shortcuts, the model was forced to rely on systematic cross-referencing of the problem setup text against visual features—a more thorough evaluation strategy.

This finding has practical implications: when deploying VLMs as CFD validators, the evaluation prompt should contain only the problem setup and the image, not metadata that could serve as a shortcut.

## 5.4 Setup-Conditioned Advantage

The primary benchmark task provides both problem setup text and the flow field image. Several of Claude's successful detections critically depend on information in the setup text:

- **S6_E2 (BC swap)**: The setup states "left wall 305K, right wall 295K," but the image shows left=cold, right=hot. Without the setup text, this is simply a normal differentially heated cavity with swapped labeling.
- **S2_E3 (wrong viscosity)**: The setup states "Re=100," but the flow develops too quickly for this Re. Without knowing the intended Re, the flow looks like a valid lower-Re channel.
- **S5_E8 (reversed lid)**: The setup states "top wall moves right," but the vortex rotates counter-clockwise (consistent with leftward lid). Without the lid direction, the flow is a valid lid-driven cavity.

This suggests that the task of CFD validation is fundamentally *setup-conditioned*: the same flow field can be correct or erroneous depending on the intended simulation parameters. This aligns with real-world CFD practice, where validators always have access to the simulation setup.

# 6. Discussion

## 6.1 Why VLMs Detect Subtle Errors

The most striking result is that Claude Opus detects three error types that the domain expert consistently misses: boundary condition swaps, wrong viscosity, and certain coarse mesh cases. We attribute this to a difference in evaluation strategy:

**Expert strategy — gestalt assessment.** The domain expert evaluates flow fields holistically, asking "does this look like the expected flow pattern?" This is highly effective for detecting gross anomalies (under-convergence, gravity reversal) but vulnerable to errors that produce plausible-looking flows. When a BC swap creates a mirror image of the correct solution, the gestalt impression is "normal cavity flow," and the subtle left-right inconsistency with the setup text is overlooked.

**VLM strategy — systematic cross-referencing.** The VLM appears to compare stated boundary conditions against observable features in the image methodically: checking which wall is hot by reading the colorbar, measuring whether the development length matches the stated Re, verifying that the vortex rotation matches the stated lid direction. This systematic approach catches quantitative inconsistencies that escape gestalt assessment.

This finding suggests that VLMs and domain experts have *complementary* capabilities for CFD validation. An ensemble approach—VLM for systematic parameter checking, expert for holistic physical judgment—could achieve higher coverage than either alone.

## 6.2 Limitations

**Sample size.** The current evaluation comprises 30 items, which we classify as an "initial finding" per our claim discipline framework. Statistical power analysis suggests approximately 100+ items would be needed to establish significant differences between evaluators at the per-error-type level.

**Single expert.** Our comparison involves one domain expert. Expert performance varies with experience level, sub-domain specialization, and fatigue. A multi-expert study would strengthen the baseline.

**Two VLMs.** We evaluate Claude Opus and Gemini 3.1. GPT-5.4 evaluation was attempted but failed due to CLI image-passing limitations. Open-source VLMs (LLaVA, InternVL, Qwen2.5-VL) remain to be tested.

**2D steady-state only.** All scenarios are 2D steady-state RANS simulations. 3D visualizations (isosurfaces, cut planes), transient flows (vortex shedding, turbulent fluctuations), and multi-physics problems (combustion, multiphase) present additional challenges not covered by this benchmark.

**Gemini reliability.** Gemini 3.1's 33% no-response rate likely reflects CLI image processing limitations rather than model capability. API-based evaluation with direct image upload may yield different results.

## 6.3 Implications for CFD Automation

These results support the feasibility of VLM-based automated CFD validation, with important caveats:

1. **Setup-conditioned prompting is essential.** Image-only evaluation is insufficient; the VLM needs to know the intended simulation parameters to detect parameter-mismatch errors.

2. **Blind spots differ by evaluator.** No single evaluator (human or AI) detects all error types. A robust validation pipeline should combine VLM checking with selective expert review for cases the VLM flags as ambiguous.

3. **Error type matters more than scenario.** The same error type (e.g., BC swap) can be severe in one scenario and subtle in another. Benchmark design must cover this interaction.

4. **Prompt design affects performance non-trivially.** The contamination paradox shows that well-intentioned metadata can degrade VLM performance. Evaluation protocols must be carefully designed.

# 7. Conclusion

We introduced CFD-VisQA, the first benchmark for evaluating vision-language model capabilities in assessing physical plausibility of CFD flow field visualizations. The benchmark comprises 10 canonical thermal-fluid scenarios with 60 OpenFOAM cases, 258 flow field images, and 279 physics-grounded questions spanning 5 difficulty levels and 6 systematically generated error types.

Blind evaluation of Claude Opus 4.6 and Gemini 3.1 against a domain expert on 30 items reveals that frontier VLMs can detect subtle CFD errors—including boundary condition mismatches, incorrect viscosity, and mesh inadequacy—that domain experts miss during visual inspection. Claude Opus achieves 100% accuracy through systematic cross-referencing of problem setup text against visual features, while the domain expert achieves 87%, primarily missing errors that produce qualitatively plausible flow patterns.

Our key findings—the contamination paradox, differential blind spots between VLMs and experts, and the critical role of setup-conditioned evaluation—provide practical guidance for deploying VLMs as automated CFD validators. The benchmark dataset and evaluation pipeline are publicly available to support further research in scientific visualization understanding.

Future work includes expanding the benchmark to 3D visualizations and transient flows, evaluating additional VLMs including fine-tuned domain-specific models, and developing ensemble validation strategies that combine VLM and expert strengths.
