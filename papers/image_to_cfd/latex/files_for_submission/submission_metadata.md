# Title

Automated Indoor CFD Analysis from a Single 2D Image via Vision-Language Model Abstraction

# Abstract

CFD simulation of indoor environments requires three-dimensional geometric models that take hours to days of manual effort to construct. This paper presents a framework that converts a single 2D image of an indoor space into a steady-state CFD solution by using a vision-language model (VLM) as the geometric abstraction engine. Given an architectural drawing or rendered image, the VLM (Gemini 3.1 Pro) extracts a structured 3D scene description — room geometry, internal obstacles, and ventilation openings — which is then scale-calibrated using one reference dimension and automatically meshed and solved with OpenFOAM. A rule-based benchmark of 20 indoor geometries at four complexity levels and five input view types (100 cases total) provides independently computed reference CFD solutions with no VLM involvement in the ground truth. Of 100 cases, 97 converge; the framework achieves a mean structural score of 0.781 and a mean CFD agreement score of 0.477. Floor-plan inputs perform best (structural 0.884, CFD 0.572); section views perform worst (0.679, 0.396). Room topology is correctly identified in 95% of cases, and opening-wall assignment matches the reference in 70% overall (100% for floor plans). Three failure modes are identified: composite-room collapse under section views, obstacle hallucination with limited CFD penalty, and a structure-versus-fidelity gap in dense configurations. The pipeline, benchmark, and evaluation code are publicly available.

# Keywords

indoor CFD; vision-language model; geometric abstraction; automated simulation; benchmark
