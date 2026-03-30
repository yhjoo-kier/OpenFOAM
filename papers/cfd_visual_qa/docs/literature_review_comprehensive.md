# Comprehensive Literature Review: CFD-VisQA

> **Purpose**: Foundation document for Introduction and Related Work sections
> **Paper**: CFD-VisQA: Benchmarking VLM Understanding of CFD Flow Field Visualizations
> **Date**: 2026-03-30
> **Citation format**: [Author Year] (to be converted to BibTeX citekeys later)

---

## 1. VLM Benchmarks for Scientific and Engineering Imagery

The rapid proliferation of vision-language models (VLMs) has generated a parallel effort to benchmark their capabilities on domain-specific tasks. While general-purpose benchmarks such as VQAv2 and GQA established foundational evaluation paradigms, the scientific and engineering communities have increasingly recognized that general visual question-answering performance does not predict domain-specific competence. This section surveys benchmarks that evaluate VLMs on scientific, technical, and engineering imagery, identifying a persistent gap in computational fluid dynamics (CFD) visualization assessment.

**Broad scientific benchmarks.** MMMU [Yue 2024] represents the most comprehensive multidisciplinary effort, assembling 11,500 college-level questions across 30 disciplines and 183 subfields, drawn from textbooks and lecture materials. With over 403 citations since its CVPR 2024 publication, MMMU has become a de facto standard for measuring VLM reasoning on expert-level content. Critically, however, engineering and science subjects yield the lowest model accuracy, and the benchmark contains no CFD flow field imagery—its engineering questions feature circuit diagrams, structural mechanics sketches, and thermodynamic cycles rather than simulation output. MathVista [Lu 2024], presented at ICLR 2024 with 184 citations, extends visual reasoning evaluation to mathematical and scientific plots, testing competencies in chart deduction, geometry problem-solving, and function understanding. While it includes some scientific visualizations, the focus remains on data plots and geometric diagrams rather than physics simulation imagery.

MMBench [Liu 2024] introduces CircularEval, a systematic protocol that rotates answer choices to reduce positional bias in VLM evaluation. Published at ECCV 2024 with 182 citations, its methodological contribution to robust evaluation design is relevant to CFD-VisQA's blind protocol development, but its content spans generic visual understanding tasks without scientific specialization.

**Scientific figure understanding.** A more targeted line of work evaluates VLMs specifically on scientific figures extracted from research publications. SciFIBench [Roberts 2024] constructs 2,000 questions on eight categories of scientific figures from arXiv papers, evaluating 28 large multimodal models. The authors conclude that VLM capabilities in scientific figure understanding are "poorly characterized," a finding that motivated subsequent domain-specific benchmarks. Multimodal ArXiv and ArXivCap [Li 2024], presented at ACL 2024, compile 6.4 million scientific figures with captions, providing a large-scale pretraining resource. While fluid mechanics figures appear incidentally, no CFD-specific curation or labeling exists. SPIQA [Pramanick 2024] focuses on question-answering about scientific paper figures from NeurIPS 2024, and CharXiv [Wang 2024] specifically targets chart understanding from arXiv papers. Both advance the methodological frontier of scientific figure QA but do not extend to simulation output imagery.

**Physics and physical reasoning.** PhysBench [Chow 2025], an ICLR 2025 oral paper, represents the most ambitious effort to evaluate VLM physical reasoning, testing 75 VLMs on 10,002 entries across four physical domains. The headline finding—GPT-4o achieves only 49.5% versus 95.9% human accuracy—underscores a substantial gap between VLM and human physical understanding. However, PhysBench's "fluid dynamics" category comprises everyday scenarios (water pouring, oil spreading) captured from natural images and videos, not CFD simulation output. The physical reasoning required is intuitive and qualitative (will the water spill?) rather than the quantitative cross-referencing of simulation parameters against flow field features that CFD validation demands.

**Visualization literacy.** Pandey and Ottley [2025], presented at EuroVis, benchmark VLMs on standardized visualization literacy assessments (VLAT and CALVI). Their finding that accuracy drops from 76--96% on simple line charts to 18--61% on multi-encoding visualizations and 25--30% on anomaly detection is directly pertinent to CFD visualization understanding. CFD colormapped fields are inherently multi-encoding: a single contour plot simultaneously encodes scalar magnitude through color, spatial structure through contour geometry, and boundary conditions through edge behavior. The measured degradation on multi-encoding tasks predicts that VLMs will struggle with quantitative interpretation of CFD imagery.

**Domain-specific scientific benchmarks.** Several recent benchmarks target specific scientific domains, providing models for CFD-VisQA's design. ClimateIQA [Chen 2024] is the most analogous, constructing 26,280 meteorological heatmap images with 762,120 instruction samples. Like CFD contour plots, meteorological heatmaps require colormap interpretation and spatial pattern recognition. ClimateIQA finds that GPT-4o and Qwen-VL "struggle with precise color identification and spatial localization"—a finding directly transferable to CFD visualization. The benchmark further develops Climate-Zoo, a fine-tuned VLM, and a SPOT algorithm for spatial reasoning, demonstrating that domain-specific training substantially improves performance on colormapped imagery. MicroVQA [Burgess 2025] targets microscopy visual question-answering at CVPR 2025, reporting peak accuracy of only 53%, reinforcing the pattern that domain-specific scientific imagery remains challenging for current VLMs. MaCBench [Alampara 2025], published in Nature Computational Science, benchmarks VLMs on chemistry and materials science imagery, extending the multi-domain evaluation landscape.

**Engineering-specific benchmarks.** DesignQA [Doris 2025], published in ASME Journal of Computing and Information Science in Engineering, evaluates VLMs on engineering documentation including CAD images and technical drawings. Picard et al. [2025], in an AI Review survey titled "From concept to manufacturing," systematically assess VLM performance across the engineering design pipeline. Both demonstrate that engineering imagery poses distinct challenges—tolerance interpretation, assembly reasoning, material identification—that differ from but complement the physical plausibility assessment required for CFD.

**Gap.** Despite this substantial body of work, no existing benchmark evaluates VLM capabilities on CFD simulation output images. The benchmarks above collectively cover textbook diagrams, data charts, scientific paper figures, natural physics scenes, meteorological heatmaps, microscopy, chemistry, and engineering drawings—but not the velocity contours, temperature fields, streamline plots, and pressure distributions that constitute the daily visual language of CFD practice. CFD-VisQA fills this gap with the first dedicated benchmark for VLM evaluation on physics simulation imagery, introducing the novel task of plausibility judgment (detecting errors in simulations) rather than comprehension of correct visualizations.

---

## 2. AI for CFD: Automation, Surrogate Modeling, and Validation

The application of artificial intelligence to computational fluid dynamics has accelerated dramatically since 2023, spanning simulation setup automation, surrogate modeling, and—most recently—natural language interfaces. This section traces the evolution from text-only LLM tools to emerging multimodal approaches, identifying the persistent absence of visual result validation.

**LLM-based CFD knowledge and automation.** CFDLLMBench [Somasekharan 2025] provides the most systematic evaluation of LLM capabilities in the CFD domain. Across a suite of conceptual multiple-choice questions, LLMs achieve 92% accuracy on CFD knowledge, demonstrating strong grasp of fluid mechanics fundamentals. However, when tasked with generating OpenFOAM case files from natural language descriptions, success rates drop to 25--34%, revealing a gap between declarative knowledge and procedural competence. Critically, CFDLLMBench is entirely text-based—no images are presented or evaluated. OpenFOAMGPT [Pandey 2025], published in Physics of Fluids, implements a retrieval-augmented generation (RAG) agent specifically for OpenFOAM simulation setup, enabling natural language queries about solver selection, boundary condition configuration, and mesh parameters. Dong et al. [2025] fine-tune Qwen2.5 on a natural-language-to-OpenFOAM (NL2FOAM) dataset, achieving 88.7% accuracy on case generation tasks, as reported in Theoretical and Applied Mechanics Letters (TAML). Wang et al. [2025], also in TAML, provide a broader systematic evaluation of LLM capabilities on CFD problems spanning formulation, discretization, and solver configuration.

**LLM-driven CFD innovation.** Beyond routine automation, researchers have explored LLMs as creative tools for CFD methodology development. AutoTurb [Zhang 2025], published in Physics of Fluids, combines LLM-generated candidate expressions with evolutionary search to discover novel turbulence model closures, demonstrating that LLMs can contribute to CFD methodology rather than merely automating existing workflows. Wang and Cheng [2025], in Engineering Applications of Artificial Intelligence (EAAI), investigate whether foundation LLMs (specifically LLaMA-3) can directly predict fluid dynamics quantities from textual problem descriptions, testing the boundary between language modeling and physics simulation. Alexiadis and Ghiassi [2024], in Results in Engineering, propose a "Text to Tech" pipeline connecting LLM-generated specifications to CFD mesh generation and solver configuration.

**Multimodal and vision-based CFD applications.** The extension from text-only to multimodal AI for CFD remains nascent. The most directly relevant work is a 2026 study in MDPI Technologies that fine-tunes Qwen2.5-VL on indoor CFD mixed-reality (MR) images, improving accuracy from below 30% to above 60% on domain-specific questions. While this demonstrates that domain-specific fine-tuning can substantially improve VLM performance on CFD imagery, the study is limited to a single VLM architecture, a small dataset, mixed-reality overlays (rather than standard post-processing output), and comprehension questions (rather than plausibility judgment). DrivaerNet++ [Elrefaie 2025], presented at NeurIPS 2025, contributes a large-scale multimodal automotive aerodynamics dataset with CFD simulation results, establishing infrastructure for training and evaluating models on engineering simulation data, though without VLM evaluation.

**Domain-specific CFD applications with AI.** Several application-focused studies demonstrate the growing integration of AI with CFD in specific domains. Dai et al. [2025], in Building and Environment, explore GPT-4o for indoor HVAC system design, testing whether VLMs can reason about thermal comfort and airflow patterns in architectural contexts. Alanis Ruiz et al. [2025], also in Building and Environment, employ conditional deep convolutional generative adversarial networks (cDCGANs) for indoor pollutant dispersion CFD, using generative models as surrogate predictors. Rianto et al. [2025], in Fire Safety Journal, apply generative AI to fire dynamics CFD simulation, and Calix et al. [2024] construct a CFD furnace image dataset specifically designed for GAN training. Calzolari and Liu [2021], in a widely cited Building and Environment review, survey the use of deep learning to replace or improve CFD simulations for building applications, establishing the conceptual framework within which VLM-based validation would operate.

**Gap.** The AI-for-CFD literature has addressed simulation setup, surrogate prediction, and domain-specific applications, but has not addressed visual result validation. LLM tools can generate OpenFOAM cases but cannot inspect whether the resulting flow fields are physically reasonable. The emerging multimodal work (MDPI 2026, DrivaerNet++) establishes that VLMs can process CFD imagery, but no study has systematically evaluated whether VLMs can detect errors in simulation outputs—the specific capability required to close automated CFD pipelines. CFD-VisQA addresses this gap by defining plausibility assessment as a distinct evaluation task and providing the first benchmark for measuring it.

---

## 3. VLM Physical Reasoning

A growing body of work investigates whether VLMs possess genuine physical understanding or merely pattern-match visual features. This question is central to CFD-VisQA: assessing flow field plausibility requires reasoning about conservation laws, boundary layer physics, and buoyancy-driven phenomena—not just recognizing visual patterns.

**Intuitive physics in VLMs.** Schulze Buschoff et al. [2025], in Nature Machine Intelligence, conduct a systematic comparison of VLM and human performance on intuitive physics tasks, finding that while VLMs can solve some physical reasoning problems, their error patterns diverge substantially from human cognition. This divergence is directly relevant to CFD-VisQA's finding that VLMs and human experts exhibit complementary blind spots—with Claude detecting subtle parameter-image inconsistencies that experts miss, and experts excelling at holistic physical judgment in image-only conditions.

**Structured physical reasoning benchmarks.** GRASP [Jassim 2024], presented at IJCAI, constructs a Unity simulation-based benchmark for physical reasoning, finding VLMs perform at or below 50%—near chance level—on tasks requiring understanding of gravity, collision, and object permanence. DynSuperCLEVR [Wang 2025], at ICLR 2025, extends to 4D physics dynamics VQA, testing whether models can reason about temporal evolution of physical systems. Ghaffari and Krishnaswamy [2024], at AAAI Spring Symposium, systematically explore failure cases in multimodal reasoning about physical dynamics, cataloging the specific failure modes that VLMs exhibit. Dewantoro et al. [2025], at IEEE Conference on Games, investigate whether multimodal LLMs can reason about structural stability—a question analogous to CFD-VisQA's plausibility assessment but in the solid mechanics domain.

**Applied physical reasoning.** Gao et al. [2024], at ICRA, develop physically grounded VLMs for robotic manipulation, demonstrating that physical reasoning can be elicited through appropriate prompting and fine-tuning for specific application domains. Duan et al. [2025], at IEEE VR Workshop, evaluate VLMs for augmented reality scene evaluation, testing spatial and physical reasoning in mixed-reality contexts. Lai et al. [2024], at IROS, specifically address VLM capabilities for liquid perception—the closest existing work to fluid dynamics visual understanding, though focused on everyday liquids in containers rather than CFD simulation output.

**Implications for CFD-VisQA.** The physical reasoning literature establishes two key findings. First, current VLMs possess limited intuitive physics capabilities, performing near chance on many structured reasoning tasks. Second, domain-specific prompting and setup information can substantially improve physical reasoning performance. CFD-VisQA's setup-conditioned protocol—providing boundary conditions, Reynolds number, and solver type alongside the image—aligns with this finding, and our results confirm it: in setup-conditioned evaluation, Claude Opus achieves 88.9% accuracy, but in image-only evaluation (requiring purely visual physical reasoning), accuracy drops to 33.3%. This rank reversal—where the domain expert leads at 66.7% in image-only conditions—suggests that current VLMs rely heavily on text-image cross-referencing rather than genuine visual physical understanding.

---

## 4. Flow Visualization and Deep Learning (Pre-VLM Era)

Before the advent of VLMs, the intersection of deep learning and flow visualization focused on task-specific feature detection rather than language-grounded understanding. This body of work establishes what was possible with classical computer vision approaches and highlights the qualitative leap that VLMs represent.

Liu et al. [2022], in a comprehensive survey in Advances in Aerodynamics, review deep learning methods for flow visualization analysis, covering vortex identification (VortexNet, Vortex-U-Net), shock wave detection (Shock-Net), turbulence structure classification, and flow regime identification from simulation images. These methods achieve high accuracy on their specific tasks but produce categorical outputs without natural language explanations—a VortexNet can detect a vortex but cannot explain why its position is physically inconsistent with the stated boundary conditions.

Kashefi [2024] provides a particularly relevant negative result, demonstrating that generative AI systems (Midjourney, DALL-E 3, Gemini) produce physically incorrect fluid mechanics imagery when prompted with classical phenomena such as von Karman vortex streets, Kelvin-Helmholtz instabilities, and Taylor-Couette flow. The generated images exhibit plausible visual aesthetics but violate fundamental physical constraints—symmetric vortex streets, wrong instability wavelengths, non-physical streamline topologies. Kashefi attributes this to the scarcity of fluid mechanics imagery in training data, much of which is locked behind journal copyright. This finding motivates the need for a public benchmark of authentic CFD visualizations.

Banerjee et al. [2024], in ACM Computing Surveys, provide an exhaustive review of over 250 papers on physics-informed computer vision, covering turbulent flow super-resolution, physics-informed neural network (PINN) reconstruction of flow fields from sparse measurements, and data-driven turbulence closure models. While this literature demonstrates that physical constraints can improve computer vision for fluid mechanics, it predates the VLM era and does not address natural language understanding of or queries about flow field images.

**Gap.** Pre-VLM approaches to flow visualization understanding are task-specific (vortex detection, regime classification), lack natural language interfaces, and cannot perform the open-ended plausibility assessment that CFD validation requires. The transition from categorical classification to language-grounded judgment represents a qualitative capability gap that VLMs could potentially fill—but that has not been evaluated until CFD-VisQA.

---

## 5. Direct Competitor and Concurrent Work

The most directly comparable work to CFD-VisQA is Ezemba et al. [2025], "Simulation vs. Hallucination: Evaluating Vision Language Models on Physical Simulations," presented at DESTION 2025. This study evaluates seven VLMs on 12 Ansys simulation images (6 structural mechanics, 6 fluid dynamics) using 40 true/false comprehension questions. The headline result—VLMs achieve approximately 69.2% accuracy on text-based questions but only approximately 53% on image-based questions, near random chance for binary classification—resonates with CFD-VisQA's finding that setup text is critical for VLM performance on simulation imagery.

However, critical methodological differences distinguish the two works.

**Task definition.** Ezemba et al. test *comprehension* of correct simulations: "What is the maximum temperature in this image?" or "Is the flow laminar or turbulent?" All 12 simulations are physically correct; the evaluation measures whether VLMs can read and understand simulation output. CFD-VisQA tests *plausibility judgment*: given a flow field that may contain intentional errors, can the evaluator determine whether the result is physically reasonable? This requires not just reading the visualization but reasoning about whether the observed patterns are consistent with the stated physics—a fundamentally different and more demanding cognitive task.

**Error generation.** Because Ezemba et al. use only correct simulations, their benchmark cannot distinguish between a VLM that genuinely understands the physics and one that defaults to "looks correct" for all inputs. CFD-VisQA's systematic error generation across 6 types and 3 severity levels enables fine-grained analysis of which physical inconsistencies VLMs can and cannot detect.

**Scale and reproducibility.** Ezemba et al. use 12 Ansys simulations with 40 questions, while CFD-VisQA provides 60 OpenFOAM cases with 258 images and 279 questions. The use of commercial Ansys software limits reproducibility, whereas CFD-VisQA's OpenFOAM-based pipeline with published scripts enables independent verification and extension.

**Expert comparison.** CFD-VisQA includes a domain expert as a calibration evaluator, enabling direct comparison of VLM and human performance. Ezemba et al. do not include expert evaluation, making it impossible to contextualize VLM performance against human capability.

**Evaluation protocol.** CFD-VisQA employs an API-isolated blind protocol with anonymized case identifiers and three independent trials per model, preventing contamination from case metadata and enabling consistency measurement. Ezemba et al. do not describe isolation or blinding procedures.

**Model selection.** Ezemba et al. evaluate small-scale open-source VLMs (LLaVA, CogVLM, InternVL), while CFD-VisQA evaluates frontier models (Claude Opus 4.6, GPT-5.4, Gemini) alongside the domain expert. The choice of frontier models yields substantially different performance profiles; CFD-VisQA's setup-conditioned accuracies (Claude 88.9%, GPT-5.4 67.8%, Gemini 57.8%) far exceed the approximately 53% image-based accuracy reported by Ezemba et al.

**Complementary finding.** Despite these differences, the Ezemba et al. result that image-based VLM accuracy approaches random chance on simulation imagery actually *strengthens* CFD-VisQA's motivation. It independently confirms that simulation output images pose a distinct challenge for VLMs and that the gap between text-based and image-based performance is a robust phenomenon across different model families, simulation software, and task definitions. CFD-VisQA advances beyond this observation by characterizing *which* errors are detectable, under what conditions, and how VLM detection patterns compare to human expert judgment.

---

## 6. Research Gap and Positioning of CFD-VisQA

The literature surveyed above converges on a well-defined research gap at the intersection of four fields.

**From scientific VLM benchmarks:** Extensive evaluation infrastructure exists for VLMs on scientific imagery—textbook diagrams, data charts, paper figures, meteorological heatmaps, microscopy, chemistry—but CFD simulation output imagery is entirely absent. The closest analog, ClimateIQA's meteorological heatmaps, shares the colormap interpretation challenge but not the physics-grounded plausibility judgment that CFD validation requires.

**From AI for CFD:** LLMs have been evaluated on CFD knowledge (92% MCQ accuracy), case generation (25--34% success), and even turbulence model discovery, but all evaluations are text-based or code-based. The visual modality—inspecting whether the resulting flow field is physically correct—remains unevaluated. The single VLM-on-CFD study (MDPI 2026) addresses comprehension rather than error detection and is limited to a single model architecture on mixed-reality imagery.

**From VLM physical reasoning:** VLMs perform near chance on structured physical reasoning benchmarks (GRASP: approximately 50%, PhysBench: 49.5%), and their physical understanding diverges qualitatively from human cognition. Domain-specific prompting improves performance, but no study has tested physical reasoning on engineering simulation output—a domain requiring quantitative reasoning about conservation laws, boundary conditions, and numerical method artifacts.

**From flow visualization ML:** Pre-VLM deep learning achieved task-specific successes (vortex detection, regime classification) but lacked natural language interfaces and open-ended reasoning capability. Generative AI fails to produce physically correct fluid mechanics imagery, and no public benchmark of authentic CFD visualizations exists for evaluation purposes.

**From concurrent work:** Ezemba et al. [2025] independently confirm that VLMs struggle with simulation imagery (approximately 53% image-based accuracy), but their comprehension-only task on correct simulations cannot reveal error detection capabilities, lacks expert comparison, and uses commercial software limiting reproducibility.

**CFD-VisQA fills this gap along six dimensions:**

1. **First CFD-specific VLM benchmark.** 10 canonical thermal-fluid scenarios, 60 OpenFOAM cases, 258 images, 279 physics-grounded questions—purpose-built for evaluating VLM understanding of CFD flow field visualizations.

2. **Error detection as evaluation task.** Six systematic error types (under-convergence, boundary condition swap, wrong viscosity, wrong turbulence model, coarse mesh, gravity flip) at varying severity levels enable fine-grained analysis of VLM error detection capabilities, going beyond comprehension to plausibility judgment.

3. **Setup-conditioned versus image-only evaluation.** The dual-protocol design reveals the critical role of simulation context: setup-conditioned evaluation yields Claude 88.9%, Expert 73.8%, GPT-5.4 67.8%, Gemini 57.8%, while image-only evaluation produces a rank reversal with Expert leading at 66.7% and Claude dropping to 33.3%. This finding quantifies the contribution of text-image cross-referencing versus purely visual physical reasoning.

4. **Expert calibration.** Including a domain expert as a baseline evaluator enables contextualizing VLM performance against human capability, revealing complementary blind spots: experts miss subtle parameter-image inconsistencies (boundary condition swaps: 20% recall, wrong viscosity: 29% recall) while excelling at holistic visual assessment.

5. **API-isolated blind protocol.** Anonymized case identifiers, three independent trials per model, and strict API isolation prevent contamination. The contamination paradox—where providing case metadata *decreases* accuracy—demonstrates that protocol design materially affects VLM evaluation outcomes, a methodological contribution applicable beyond CFD.

6. **Open and reproducible.** All CFD cases are generated with OpenFOAM using published scripts, all evaluation code is provided, and all results are independently reproducible—contrasting with commercial-software-based evaluations that cannot be verified or extended.

**Positioning statement.** CFD-VisQA is situated at the convergence of scientific VLM benchmarking and AI-assisted CFD validation. It contributes the first evidence that frontier VLMs can detect specific categories of CFD simulation errors through systematic cross-referencing of setup parameters against visual features—a capability not predicted by existing physical reasoning benchmarks—while simultaneously revealing that this capability collapses without textual context, exposing a fundamental limitation in current VLM visual physics understanding. The benchmark provides both a diagnostic tool for measuring progress toward visually grounded physical reasoning and a practical evaluation framework for deploying VLMs as automated CFD validators.

---

## References (to be converted to BibTeX citekeys)

### Category 1: VLM Benchmarks
- [Yue 2024] Yue et al., "MMMU: A Massive Multi-discipline Multimodal Understanding and Reasoning Benchmark," CVPR 2024.
- [Roberts 2024] Roberts et al., "SciFIBench: Benchmarking Large Multimodal Models for Scientific Figure Interpretation," NeurIPS 2024.
- [Chow 2025] Chow et al., "PhysBench: Benchmarking and Enhancing Vision-Language Models for Physical World Understanding," ICLR 2025.
- [Lu 2024] Lu et al., "MathVista: Evaluating Mathematical Reasoning of Foundation Models in Visual Contexts," ICLR 2024.
- [Liu 2024] Liu et al., "MMBench: Is Your Multi-modal Model an All-around Player?" ECCV 2024.
- [Chen 2024] Chen et al., "ClimateIQA: A Multimodal Benchmark for Climate-Related Visual Question Answering," 2024.
- [Doris 2025] Doris et al., "DesignQA: A Multimodal Benchmark for Evaluating Large Language Models' Understanding of Engineering Documentation," ASME JCISE 2025.
- [Picard 2025] Picard et al., "From concept to manufacturing: Evaluating vision-language models for engineering design," AI Review 2025.
- [Alampara 2025] Alampara et al., "MaCBench: Benchmarking Vision-Language Models on Chemistry and Materials Science," Nature Computational Science 2025.
- [Li 2024] Li et al., "Multimodal ArXiv: A Dataset for Improving Scientific Comprehension of Large Vision-Language Models," ACL 2024.
- [Pramanick 2024] Pramanick et al., "SPIQA: A Dataset for Multimodal Question Answering on Scientific Papers," NeurIPS 2024.
- [Burgess 2025] Burgess et al., "MicroVQA: A Microscopy Visual Question Answering Benchmark," CVPR 2025.
- [Wang 2024 CharXiv] Wang et al., "CharXiv: Charting Gaps in Realistic Chart Understanding in Multimodal LLMs," NeurIPS 2024.
- [Pandey 2025] Pandey and Ottley, "Visualization Literacy of Vision-Language Models," EuroVis 2025.

### Category 2: AI for CFD
- [Somasekharan 2025] Somasekharan et al., "CFDLLMBench: A Benchmark for Evaluating LLMs on CFD," 2025.
- [Pandey 2025 OpenFOAMGPT] Pandey et al., "OpenFOAMGPT: A RAG-based LLM Agent for OpenFOAM," Physics of Fluids 2025.
- [Dong 2025] Dong et al., "Fine-tuning LLM for CFD: NL2FOAM," Theoretical and Applied Mechanics Letters 2025.
- [Wang 2025 TAML] Wang et al., "LLM evaluations in CFD," Theoretical and Applied Mechanics Letters 2025.
- [Zhang 2025] Zhang et al., "AutoTurb: LLM + Evolutionary Search for Turbulence Models," Physics of Fluids 2025.
- [Wang 2025 EAAI] Wang and Cheng, "Can Foundation LLMs Predict Fluid Dynamics?" Engineering Applications of AI 2025.
- [Alexiadis 2024] Alexiadis and Ghiassi, "Text to Tech: LLM + CFD Pipeline," Results in Engineering 2024.
- [MDPI 2026] Anonymous, "Qwen2.5-VL Fine-tuned on Indoor CFD MR Images," MDPI Technologies 2026.
- [Calzolari 2021] Calzolari and Liu, "Deep learning to replace CFD for building applications," Building and Environment 2021.
- [Elrefaie 2025] Elrefaie et al., "DrivaerNet++: A Large-Scale Multimodal Car Dataset with CFD Simulations," NeurIPS 2025.
- [Dai 2025] Dai et al., "GPT-4o for indoor HVAC design," Building and Environment 2025.
- [Alanis Ruiz 2025] Alanis Ruiz et al., "cDCGAN for indoor pollutant dispersion CFD," Building and Environment 2025.
- [Rianto 2025] Rianto et al., "Generative AI for fire CFD simulation," Fire Safety Journal 2025.
- [Calix 2024] Calix et al., "CFD furnace image dataset for GANs," 2024.

### Category 3: VLM Physical Reasoning
- [Schulze Buschoff 2025] Schulze Buschoff et al., "Visual cognition in multimodal large language models," Nature Machine Intelligence 2025.
- [Gao 2024] Gao et al., "Physically Grounded Vision-Language Models for Robotic Manipulation," ICRA 2024.
- [Jassim 2024] Jassim et al., "GRASP: A Grid-Based Benchmark for Evaluating Commonsense Spatial Reasoning," IJCAI 2024.
- [Wang 2025 DynSuperCLEVR] Wang et al., "DynSuperCLEVR: Dynamic Scene Understanding via 4D Physics Simulation," ICLR 2025.
- [Ghaffari 2024] Ghaffari and Krishnaswamy, "Exploring Failure Cases in Multimodal Reasoning about Physical Dynamics," AAAI Spring Symposium 2024.
- [Dewantoro 2025] Dewantoro et al., "Can Multimodal LLMs Reason About Stability?" IEEE Conference on Games 2025.
- [Duan 2025] Duan et al., "VLM for AR Scene Evaluation," IEEE VR Workshop 2025.
- [Lai 2024] Lai et al., "VLM for Liquid Perception," IROS 2024.

### Category 4: Flow Visualization + Deep Learning
- [Liu 2022] Liu et al., "Deep learning for flow visualization," Advances in Aerodynamics 2022.
- [Kashefi 2024] Kashefi, "Generative AI and fluid mechanics imagery," 2024.
- [Banerjee 2024] Banerjee et al., "Physics-Informed Computer Vision: A Review and Perspectives," ACM Computing Surveys 2024.

### Category 5: Direct Competitor
- [Ezemba 2025] Ezemba et al., "Simulation vs. Hallucination: Evaluating Vision Language Models on Physical Simulations," DESTION 2025.
