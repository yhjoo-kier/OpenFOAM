# 5. Results

## 5.1 Overall Performance

[Table:main-results] summarizes the overall accuracy of each evaluator on 30 blind evaluation items spanning all 10 scenarios. The items include 12 correct cases and 18 error cases across 6 error types and 3 visualization modes (V1 velocity contour, V2 temperature contour, V6 composite).

Claude Opus 4.6 achieves perfect accuracy (30/30), correctly identifying all correct cases as plausible and all error cases as anomalous. The domain expert achieves 87% accuracy (26/30), missing 4 error cases. Gemini 3.1 achieves 80% accuracy on valid responses (16/20) but exhibits a 33% no-response rate (10/30 items returned empty or unparseable output), yielding an effective accuracy of 53% when counting no-responses as failures.

No evaluator produced false positives (marking a correct case as erroneous), indicating that all three are conservative—they only flag anomalies when genuinely suspicious. The critical failure mode for a CFD validator is the false negative: accepting an erroneous flow field as correct. On this metric, Claude achieves 0% false negative rate, the expert 22% (4/18 errors missed), and Gemini 20% of valid responses (but higher when including no-responses on error items).

## 5.2 Per-Error-Type Analysis

[Table:per-error-recall] reveals that error detection capability is strongly dependent on error type, and different evaluators have different blind spots.

**Consistently detectable (recall ≥ 80% across all evaluators):**
- *E1 Under-convergence*: Produces visually obvious artifacts (non-physical velocity peaks, incomplete flow development). All evaluators detect these reliably.
- *E8 Gravity/direction reversal*: Inverts expected flow patterns (downward plume instead of upward, reversed vortex rotation). Detectable through physics reasoning about buoyancy direction or vortex-lid consistency.

**Expert-challenging, VLM-detectable:**
- *E2 Boundary condition swap*: The expert detects this in severe cases (lid on bottom wall → dramatic flow change) but misses subtle cases (hot/cold wall swap → plausible mirror image). Claude detects all BC swaps by cross-referencing stated wall temperatures against colorbar values in the image.
- *E3 Wrong viscosity*: The expert misses this in laminar channel flow (S2) where the flow looks qualitatively correct at a different Re, but detects it in natural convection (S1) where the boundary layer is completely suppressed. Claude detects both by noting quantitative inconsistencies (development length too short for stated Re, thermal penetration too deep).

**Context-dependent detectability:**
- *E5 Coarse mesh*: Detection depends on whether the scenario has fine-scale structures that become visibly degraded. In lid-driven cavity (S5), the expert misses it ("looks pretty"), but in heated plate (S1) and Rayleigh-Bénard (S9), coarse mesh artifacts are more apparent. Claude detects all except the S5 case in the initial 15-item evaluation but detects it in the expanded 30-item evaluation.
- *E4 Wrong turbulence model*: Using a laminar solver at Re=10,000 produces smooth, seemingly plausible flow that the expert judges as "reasonable-looking." Claude detects it by reasoning that a flat turbulent profile from a laminar solver is non-physical at this Re.

## 5.3 The Contamination Paradox

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
