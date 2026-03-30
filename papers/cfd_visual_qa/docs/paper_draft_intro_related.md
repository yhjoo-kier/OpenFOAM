# 1. Introduction

Computational fluid dynamics (CFD) is increasingly integrated into automated design and analysis pipelines, from building ventilation assessment to industrial heat exchanger optimization. Recent advances in vision-language models (VLMs) have enabled new paradigms such as image-to-CFD frameworks, where a VLM extracts geometric information from architectural images to automatically generate CFD simulation domains [Joo 2026]. As these automated pipelines proliferate, a critical bottleneck emerges: *who validates the results?*

In traditional CFD workflows, an experienced engineer inspects velocity contours, temperature fields, and streamline plots to confirm that the solution is physically reasonable—checking for proper boundary layer development, expected recirculation patterns, and absence of numerical artifacts. This visual inspection requires domain expertise accumulated over years and cannot scale with the throughput of automated pipelines. If a VLM could perform this validation step—assessing whether a CFD flow field visualization is physically plausible—the automation loop would close entirely.

Yet it remains unknown whether current VLMs possess the physical reasoning capability to evaluate CFD output images. Existing VLM benchmarks test scientific knowledge through multiple-choice questions on textbook diagrams [Yue 2024 MMMU], chart comprehension [Pandey 2025], or everyday physics understanding [Chow 2025 PhysBench], but none evaluate the ability to assess *simulation output images* against known physical laws. The closest work evaluates VLMs on meteorological heatmaps [Chen 2024 ClimateIQA] or tests generative AI on fluid mechanics prompts [Kashefi 2024], but neither addresses the specific task of CFD result validation.

This paper introduces **CFD-VisQA**, the first benchmark designed to evaluate VLM capabilities in assessing physical plausibility of CFD flow field visualizations. The benchmark comprises 10 canonical thermal-fluid scenarios (natural convection, forced convection, mixed convection), each with correct solutions and systematically generated erroneous solutions spanning 6 error types. We evaluate two frontier VLMs—Claude Opus 4.6 and Gemini 3.1—against a domain expert on 30 blind evaluation items, revealing surprising performance patterns.

Our contributions are:

1. **CFD-VisQA benchmark dataset**: 60 OpenFOAM CFD cases across 10 thermal-fluid scenarios, producing 258 flow field images with 279 physics-grounded questions and expert annotations. The dataset includes both correct solutions and 6 types of intentionally generated errors at varying severity levels.

2. **Setup-conditioned evaluation protocol**: We demonstrate that providing problem setup text alongside the image (boundary conditions, Reynolds number, solver type) dramatically improves VLM performance compared to image-only evaluation, and introduce a blind evaluation protocol that prevents bias from case-name leakage.

3. **Comparative analysis revealing differential blind spots**: On 30 evaluation items, Claude Opus achieves 100% accuracy while the domain expert achieves 87%, with each missing different error types—the expert misses "subtle" errors (boundary condition swaps, wrong viscosity) while Gemini shows higher no-response rates on certain image types. This suggests VLMs and experts have complementary strengths for CFD validation.

4. **Contamination paradox**: We document that providing case-name metadata to VLMs *decreases* performance (87% contaminated vs. 100% blind), hypothesizing that name-based priors interfere with systematic visual-textual cross-referencing.

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
