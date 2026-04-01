# Title

Can Vision-Language Models Evaluate Computational Fluid Dynamics? A Benchmark for Physical Plausibility Assessment of Flow Field Visualizations

# Abstract

Validating computational fluid dynamics (CFD) simulation results remains a manual bottleneck in automated pipelines, requiring expert visual inspection of flow field images. Vision-language models (VLMs) could close this automation gap, yet their capability for physics-grounded assessment of CFD visualizations is unknown. We introduce CFD-VisQA, the first benchmark for this task, comprising 10 canonical thermal-fluid scenarios with 60 OpenFOAM cases, 258 flow field images, and 279 physics-grounded questions spanning 6 systematically generated error types. Using an API-isolated, blind evaluation protocol with base64-encoded images and three independent trials per model, we evaluate Claude Opus 4.6, GPT-5.4, Gemini 3.1 Pro, and a domain expert on 30 common items under two conditions: setup-conditioned (with problem description) and image-only. The results reveal a striking rank reversal: with setup text, Claude leads at 88.9% (zero false negatives across all trials), followed by the expert at 73.8%, GPT-5.4 at 67.8%, and Gemini at 57.8%; without setup text, the expert leads at 66.7% while Claude drops to 33.3%, the lowest among all evaluators. Claude's accuracy drop of 55.6 percentage points confirms that its strength derives from systematic cross-referencing of stated parameters against visual features, a capability that collapses without textual context. The expert's stability (-7.1 pp) reflects gestalt assessment grounded in physical intuition. We further demonstrate that subagent-based evaluation inflates VLM accuracy from 88.9% to 99.6% due to potential filesystem contamination, establishing API isolation as a methodological requirement for trustworthy VLM benchmarking. The benchmark dataset and evaluation pipeline are publicly available.

# Keywords

vision-language model; computational fluid dynamics; benchmark; flow visualization; physical plausibility; automated validation

# Target Journal

Engineering Applications of Artificial Intelligence (Elsevier)

# Author

Younghwan Joo, Energy Efficiency Research Division, Korea Institute of Energy Research (KIER); Energy Engineering, University of Science & Technology (UST)

# Corresponding Author

Younghwan Joo, yhjoo@kier.re.kr

# Data Availability

https://github.com/yhjoo-kier/OpenFOAM (directory: papers/cfd_visual_qa/)
