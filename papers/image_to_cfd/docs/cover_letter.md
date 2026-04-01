2026-04-01

Dear Editor:

I wish to submit an original article entitled "Automated Indoor CFD Analysis from a Single 2D Image via Vision-Language Model Abstraction" for consideration in Engineering Applications of Artificial Intelligence.

CFD is a standard tool for indoor airflow diagnostics, yet it is rarely used in routine energy audits because the 3D geometry has to be built by hand from drawings or photos, which takes hours of CAD work. This manuscript shows that a vision-language model can extract a 3D scene description from a single image accurately enough to produce a physically meaningful CFD solution. To my knowledge, no prior work has tested VLM spatial perception against downstream CFD accuracy with a controlled benchmark.

The study uses a 100-case benchmark (20 geometries at four complexity levels, five view types each) with reference CFD solutions computed independently of the VLM. Geometric fidelity and flow-field accuracy are scored separately, which revealed that a geometrically accurate prediction does not guarantee a good CFD match --- an important finding for setting realistic expectations. Three systematic failure modes define where the approach works and where it breaks down. The solver is validated against the Nielsen et al. (1978) ventilated room experiment, and the full pipeline and benchmark are publicly available.

I believe this work fits well within the scope of Engineering Applications of Artificial Intelligence. The study demonstrates a practical engineering application of frontier AI models --- converting architectural images into CFD simulations --- with a systematic evaluation framework. The VLM-based geometric abstraction addresses a real bottleneck in building performance assessment, and the benchmark methodology is applicable to other engineering domains where AI-generated inputs feed into physics-based solvers.

This manuscript has not been published or presented elsewhere in part or in entirety and is not under consideration by another journal. I have no conflicts of interest to declare. I have read and understood the journal's policies and submit this manuscript in accordance with them.

Thank you for your consideration. I look forward to hearing from you.

Sincerely,

Younghwan Joo

Energy Efficiency Research Division

Korea Institute of Energy Research

152 Gajeong-ro, Yuseong-gu, Daejeon 34129, Republic of Korea

E-mail: yhjoo@kier.re.kr
