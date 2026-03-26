# Cover Letter — Building and Environment

Dear Editor,

I submit the manuscript "Automated Indoor CFD Analysis from a Single 2D Image via Vision-Language Model Abstraction" for consideration in Building and Environment.

CFD is a standard tool for indoor airflow diagnostics, yet it is rarely used in routine energy audits because the 3D geometry has to be built by hand from drawings or photos, which takes hours of CAD work. This manuscript shows that a vision-language model can extract a 3D scene description from a single image accurately enough to produce a physically meaningful CFD solution. To my knowledge, no prior work has tested VLM spatial perception against downstream CFD accuracy with a controlled benchmark.

The study uses a 100-case benchmark (20 geometries at four complexity levels, five view types each) with reference CFD solutions computed independently of the VLM. Geometric fidelity and flow-field accuracy are scored separately, which revealed that a geometrically accurate prediction does not guarantee a good CFD match — an important finding for setting realistic expectations. Three systematic failure modes define where the approach works and where it breaks down. The solver is validated against the Nielsen et al. (1978) ventilated room experiment, and the full pipeline and benchmark are publicly available.

I believe this work fits Building and Environment because it sits at the intersection of indoor CFD simulation and computational tools for the built environment. The approach could make airflow screening accessible to practitioners who currently lack the CAD skills or time to prepare geometry models by hand.

This manuscript has not been published previously and is not under consideration elsewhere. I declare no competing interests.

Sincerely,

Younghwan Joo
Energy Efficiency Research Division
Korea Institute of Energy Research
152 Gajeong-ro, Yuseong-gu, Daejeon 34129, Republic of Korea
Email: yhjoo@kier.re.kr
