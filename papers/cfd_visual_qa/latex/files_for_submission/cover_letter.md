2026-04-01

Dear Editor:

I wish to submit an original article entitled "Can Vision-Language Models Evaluate CFD? A Benchmark for Physical Plausibility Assessment of Flow Field Visualizations" for consideration in Engineering Applications of Artificial Intelligence.

Automated CFD pipelines increasingly rely on AI for geometry extraction and simulation setup, yet no systematic method exists to validate the resulting flow fields without expert visual inspection. This manuscript introduces CFD-VisQA, the first benchmark for evaluating whether vision-language models can assess the physical plausibility of CFD flow field visualizations. To my knowledge, no prior work has evaluated VLM capability for CFD result validation with a controlled benchmark and an API-isolated evaluation protocol.

The study uses a 60-case benchmark (10 canonical thermal-fluid scenarios with 6 systematically generated error types), evaluating three frontier VLMs (Claude Opus 4.6, GPT-5.4, Gemini 3.1 Pro) and a domain expert under two conditions: setup-conditioned and image-only. The central finding is a rank reversal: with problem setup text, Claude leads at 88.9% accuracy with zero false negatives; without setup text, the expert leads at 66.7% while Claude drops to 33.3%. This reveals that VLM performance derives from setup-image cross-referencing rather than visual understanding of physics. We further demonstrate that subagent-based evaluation inflates accuracy from 88.9% to 99.6% due to filesystem contamination, establishing API isolation as a methodological requirement for VLM benchmarking.

I believe this work fits well within the scope of Engineering Applications of Artificial Intelligence. The study addresses a practical engineering problem --- automated validation in CFD pipelines --- using frontier AI models, and the benchmark provides a reusable resource for the research community. The findings have direct implications for deploying VLMs as automated validators in engineering simulation workflows.

This manuscript has not been published or presented elsewhere in part or in entirety and is not under consideration by another journal. I have no conflicts of interest to declare. I have read and understood the journal's policies and submit this manuscript in accordance with them.

Thank you for your consideration. I look forward to hearing from you.

Sincerely,

Younghwan Joo

Energy Efficiency Research Division

Korea Institute of Energy Research

152 Gajeong-ro, Yuseong-gu, Daejeon 34129, Republic of Korea

E-mail: yhjoo@kier.re.kr
