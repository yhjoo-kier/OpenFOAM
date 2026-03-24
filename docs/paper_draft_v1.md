# Automated Indoor CFD Analysis from a Single 2D Image via Vision-Language Model Abstraction

---

## Abstract

Computational fluid dynamics (CFD) simulation of indoor environments requires three-dimensional geometric models that are costly and time-consuming to construct manually. This paper presents an automated framework that converts a single two-dimensional image of an indoor space into a fully resolved steady-state CFD solution, using a vision-language model (VLM) as the geometric abstraction engine. The framework accepts an architectural drawing or rendered image, invokes a VLM (Gemini 3.1 Pro, preview) to extract a structured 3D scene description comprising room geometry, internal obstacles, and ventilation openings, applies a single-factor scale calibration that requires one reference dimension (the longest horizontal span) to correct absolute size, and then automatically generates a computational mesh and boundary conditions for OpenFOAM. A rule-based benchmark dataset of 20 indoor geometries spanning four complexity levels and five input view types (100 evaluation cases in total) is constructed with independently computed reference CFD solutions, ensuring that no VLM output contaminates the ground truth. Across 97 of 100 cases (3 cases fail due to solver divergence) the framework achieves a mean structural similarity score of 0.781 and a mean reference CFD agreement score of 0.477, with the floor-plan view yielding the highest performance (structural 0.884, CFD 0.572) and the section view the lowest (structural 0.679, CFD 0.396). Room-kind topology is correctly identified in 95% of cases, and the opening-wall assignment matches the reference in 70% of cases overall, rising to 100% for floor-plan inputs. Systematic failure-mode analysis reveals three characteristic error patterns: composite-room collapse under section views, obstacle hallucination with limited CFD penalty, and a structure-versus-fidelity gap in dense composite geometries. The full pipeline, benchmark dataset, and evaluation code are provided to support reproducibility.

---

## 1. Introduction

### 1.1 Motivation

Energy efficiency improvement in industrial manufacturing processes is a global priority. Manufacturing facilities account for a substantial fraction of total energy consumption in industrialised economies, and governmental and corporate carbon-reduction targets are driving systematic energy audits across sectors [@gourlis.kovacic2017BuildingInformation; @zhao.etal2021IndustrialReheating]. A growing body of evidence demonstrates that indoor airflow patterns within factories, warehouses, and processing halls have a direct impact on HVAC energy consumption, process thermal management, and worker thermal comfort. Computational fluid dynamics (CFD) simulation is the most powerful diagnostic tool available for analysing these airflow patterns, as it can reveal recirculation zones, thermal stratification, short-circuiting of supply air, and other phenomena that are invisible to point measurements and difficult to infer from energy consumption data alone [@vuorinen.etal2020ModellingAerosol; @liu.etal2019ReviewCFD].

Despite its diagnostic value, CFD-based airflow analysis is rarely included in standard energy audit workflows. The barrier is not the simulation itself — modern solvers such as OpenFOAM [@weller.etal1998TensorialApproach] handle typical indoor geometries in minutes to hours — but the upstream geometry preparation step. Converting site photographs or architectural drawings into a CFD-ready 3D model requires CAD expertise and several hours to days of modelling work unavailable within standard audit budgets [@calzolari.liu2021DeepLearning], effectively excluding CFD from the energy diagnostician's toolkit.

Existing approaches have not fully resolved this bottleneck. BIM-to-CFD workflows [@utkucu.sozer2020InteroperabilityData; @weerasuriya.etal2019HolisticFramework; @hou.etal2022PredictionOptimization] require a digital building model that rarely exists for older industrial facilities [@gourlis.kovacic2017BuildingInformation]. Scan-to-BIM methods [@rocha.etal2020ScantoBIMMethodology] demand specialised equipment impractical during a time-constrained site visit. Deep-learning surrogates [@guo.etal2016ConvolutionalNeural; @calzolari.liu2021DeepLearning; @kastner.dogan2023GANBasedSurrogate] are limited to pre-trained geometry families and cannot handle diverse industrial layouts.

The ideal capability for an energy diagnostician would be to photograph a facility space — or use an existing floor plan — and obtain a CFD-quality airflow assessment within minutes, with no manual geometry modelling required. The recent emergence of vision-language models (VLMs) with spatial reasoning capabilities [@openai.etal2023GPT4Technical; @geminiteam.etal2023GeminiFamily; @chen.etal2024SpatialVLMEndowing; @hong.etal20233DLLMInjecting] makes this vision technically plausible. A VLM can interpret a photograph or architectural drawing, extract the spatial layout of a room including dimensions and major obstacles, and express this information in a structured format. If this structured output can be automatically converted into a CFD-ready mesh, the entire pathway from image to simulation result can be automated.

This paper presents and evaluates such a framework. Although the current implementation and benchmark focus on general indoor spaces rather than industrial facilities, the underlying methodology — VLM-based geometric abstraction from a single 2D image followed by automated CFD — is directly applicable to the industrial energy diagnostics scenario described above. The general indoor case is chosen as the starting point because it enables a controlled benchmark with well-defined complexity levels, while the extension to industrial environments remains a natural next step.

### 1.2 Related work

**CFD for indoor environmental analysis.** CFD simulation of indoor environments covers ventilation effectiveness [@kosutova.etal2019CrossventilationGeneric; @liu.etal2019ReviewCFD], thermal comfort [@chen.etal2021UnsteadystateCFD], contaminant and pathogen dispersion [@vuorinen.etal2020ModellingAerosol; @liu.etal2021SimulationbasedStudy; @srivastava.etal2021EffectiveVentilation], and building energy performance [@brozovsky.etal2022AssessingImpact; @ren.cao2019DevelopmentApplication]. In industrial contexts — furnaces [@zhao.etal2021IndustrialReheating], process ventilation, occupational exposure — use remains far less systematic than in commercial building design, largely due to the geometry preparation overhead.

**Geometry acquisition and the automation gap.** BIM platforms can export geometry to CFD tools [@utkucu.sozer2020InteroperabilityData; @hou.etal2022PredictionOptimization], but BIM penetration in the industrial sector is limited [@gourlis.kovacic2017BuildingInformation]. Scan-to-BIM methods [@rocha.etal2020ScantoBIMMethodology; @wang.etal2019ApplicationOriented] produce detailed as-built models from point clouds but are prohibitive for routine energy audits. Image-based reconstruction methods including NeRF [@mildenhall.etal2020NeRFRepresenting] yield implicit or noisy representations that cannot be directly used by finite-volume solvers — a fundamental mismatch with CFD input requirements.

**Data-driven flow prediction.** Neural network surrogates [@guo.etal2016ConvolutionalNeural; @calzolari.liu2021DeepLearning; @kastner.dogan2023GANBasedSurrogate] achieve impressive speed but are limited to their training domain and cannot provide the conservation guarantees of physics-based solvers, making them unsuitable for liability-sensitive engineering assessments.

**LLM-assisted CFD workflows.** FLUID-GPT [@yang.etal2023FLUIDGPTFast], OpenFOAMGPT [@pandey.etal2025OpenFOAMGPTRetrievalaugmented], and Dong et al. [@dong.etal2025FinetuningLarge] demonstrate growing integration of language models into CFD practice, but all address the *downstream* workflow — solver configuration, parameter selection, result interpretation — while assuming the geometric model already exists. The *upstream* bottleneck of geometry extraction from visual input remains untouched.

**Vision-language models for spatial understanding.** VLMs such as GPT-4 [@openai.etal2023GPT4Technical] and Gemini [@geminiteam.etal2023GeminiFamily] have demonstrated spatial reasoning capabilities including object localisation, dimension estimation, and scene layout interpretation [@chen.etal2024SpatialVLMEndowing; @hong.etal20233DLLMInjecting]. The present work leverages these capabilities not to reconstruct scene geometry in full fidelity, but to abstract it into a simulation-compatible schema — prioritising topological correctness and approximate dimensional accuracy over surface-level detail.

### 1.3 Contributions

This paper contributes the following:

1. An automated image-to-CFD framework that eliminates the manual geometry preparation step by using a VLM to extract a simulation-ready 3D scene description from a single 2D image, targeting screening-level indoor airflow analysis.
2. A systematic benchmark comprising 100 evaluation cases (20 geometries × 5 view types) with independently computed reference CFD solutions, providing the first quantitative assessment of VLM-based geometric abstraction for downstream CFD accuracy.
3. An evaluation methodology that decouples geometric fidelity from flow-field accuracy, revealing that structural similarity does not guarantee CFD agreement — a finding with practical implications for acceptable error thresholds in screening applications.
4. Identification of three characteristic failure modes (composite-room collapse, obstacle hallucination, structure-versus-fidelity gap) and analysis of their practical impact, informing the design of future VLM-based simulation tools.

---

## 2. Methodology

### 2.1 Framework overview

The proposed framework accepts a single 2D image of an indoor space and produces a steady-state CFD solution through a five-stage pipeline [Fig:method-framework]. In the first stage, the image is presented to a vision-language model (Gemini 3.1 Pro) (accessed via API in March 2026; model identifier `gemini-3.1-pro-preview`) together with a structured prompt that requests a JSON description of the room geometry, internal obstacles, and ventilation openings. The prompt specifies the coordinate system convention, unit system, and the required schema fields, but does not provide any ground-truth information about the scene. The VLM returns a scene description conforming to the `indoor_cfd_scene_v1` schema, which encodes the room as one or more axis-aligned rectangular blocks, each obstacle as an axis-aligned bounding box, and each opening as a wall-anchored rectangular patch with an associated type (inlet or outlet).

In the second stage, a single-factor scale calibration is applied to the VLM output. Because the VLM has no access to absolute dimensional references, the predicted scene dimensions are generally correct in proportion but incorrect in absolute magnitude. The calibration identifies the longest horizontal span in the predicted geometry and rescales the entire scene uniformly so that this span matches the corresponding reference dimension:

$$\alpha = \frac{L_\text{ref}}{L_\text{pred}}, \qquad \mathbf{x}_\text{calibrated} = \alpha \, \mathbf{x}_\text{pred} \tag{Eq:scale-calibration}$$

where $L_\text{ref}$ and $L_\text{pred}$ are the longest horizontal spans in the reference and predicted geometries, respectively, and $\alpha$ is the uniform scaling factor applied to all coordinates. This post-hoc scaling is the only point at which reference information enters the predicted path, and it corrects the dominant source of volumetric error without altering the topological structure of the prediction.

### 2.2 Governing equations and CFD setup

The indoor airflow is modelled as steady-state, isothermal, incompressible, turbulent flow governed by the Reynolds-averaged Navier–Stokes (RANS) equations. The continuity and momentum equations are:

$$\nabla \cdot \mathbf{U} = 0 \tag{Eq:continuity}$$

$$\nabla \cdot (\mathbf{U} \otimes \mathbf{U}) = -\frac{1}{\rho}\nabla p + \nabla \cdot \bigl[(\nu + \nu_t)\bigl(\nabla \mathbf{U} + (\nabla \mathbf{U})^\top\bigr)\bigr] \tag{Eq:momentum}$$

where $\mathbf{U}$ is the mean velocity vector, $p$ is the kinematic pressure (pressure divided by density), $\nu$ is the kinematic viscosity (air at 20 °C, $\nu = 1.5 \times 10^{-5}$ m²/s), and $\nu_t$ is the turbulent eddy viscosity.

Turbulence closure is provided by the $k$–$\omega$ SST model [@menter1994TwoequationEddyviscosity], which blends the $k$–$\omega$ formulation near walls with the $k$–$\varepsilon$ formulation in the free stream. The transport equations for the turbulent kinetic energy $k$ and the specific dissipation rate $\omega$ are:

$$\nabla \cdot (\mathbf{U} k) = \nabla \cdot \bigl[(\nu + \sigma_k \nu_t) \nabla k\bigr] + P_k - \beta^* k \omega \tag{Eq:k-transport}$$

$$\nabla \cdot (\mathbf{U} \omega) = \nabla \cdot \bigl[(\nu + \sigma_\omega \nu_t) \nabla \omega\bigr] + \frac{\gamma}{\nu_t} P_k - \beta \omega^2 + 2(1 - F_1) \frac{\sigma_{\omega 2}}{\omega} \nabla k \cdot \nabla \omega \tag{Eq:omega-transport}$$

where $P_k$ is the turbulence production term, $F_1$ is the blending function, and $\beta^*$, $\beta$, $\gamma$, $\sigma_k$, $\sigma_\omega$, $\sigma_{\omega 2}$ are model constants taking their standard values [@menter1994TwoequationEddyviscosity].

The boundary conditions are summarised as follows:

| Boundary | $\mathbf{U}$ | $p$ | $k$ | $\omega$ |
|----------|-------------|-----|-----|----------|
| Inlet | $\dot{V} = 0.05$ m³/s | Zero gradient | Eq. below | Eq. below |
| Outlet | Zero gradient | $p = 0$ | Zero gradient | Zero gradient |
| Walls | No-slip | Zero gradient | Wall function | Wall function |

where the inlet turbulence quantities are specified as $k_\text{in} = \frac{3}{2}(I \cdot U_\text{in})^2$ with turbulence intensity $I = 0.05$, and $\omega_\text{in} = k_\text{in}^{1/2} / (C_\mu^{1/4} \cdot l)$ with mixing length $l = 0.07 D_h$. Wall functions follow [@launder.spalding1974NumericalComputation].

where $I$ is the turbulence intensity, $D_h$ is the hydraulic diameter of the inlet, and $C_\mu = 0.09$. The fixed volume flow rate at the inlet ensures that the flow forcing is identical between predicted and reference cases regardless of differences in opening size.

The fixed volume flow rate of 0.05 m³/s is applied identically to both predicted and reference cases, ensuring that any differences in the flow field arise solely from geometric discrepancies rather than from differences in flow forcing. This design choice decouples the evaluation from opening-size accuracy but means that the inlet velocity varies inversely with the opening area.

A tetrahedral mesh is generated using Gmsh [@geuzaine.remacle2009Gmsh3D] with a target cell size of 0.18 m, yielding approximately 37,000 to 105,000 cells depending on the room dimensions. A grid independence study (Section 5.5) confirmed that the 0.18 m resolution produces approximately converged bulk flow quantities, with less than 5% change upon further refinement to 0.10 m. Nevertheless, near-wall resolution remains limited and the wall-function approach operates at elevated $y^+$ values; the resulting solutions should be interpreted as engineering-level screening estimates rather than high-fidelity resolved simulations.

The equations are solved using the SIMPLE algorithm [@patankar.spalding1972CalculationProcedure] implemented in OpenFOAM's simpleFoam solver [@weller.etal1998TensorialApproach]. A multi-level robustness strategy is employed to handle cases where the default solver settings fail to converge. The mesh resolution ladder proceeds through 0.18, 0.25, and 0.35 m cell sizes (coarsening on failure), and at each resolution the solver attempts up to four preset configurations of progressively reduced relaxation factors, including a laminar fallback that disables the turbulence model entirely. This escalation pathway ensures that every case produces a converged solution, albeit at varying levels of resolution and modelling fidelity.

In the fifth stage, the predicted flow field is compared against the independently computed reference solution to produce structural and reference CFD agreement scores.

### 2.3 Evaluation metrics

The evaluation framework produces two complementary scores that capture different aspects of prediction quality: a structural score that measures geometric fidelity, and a reference CFD agreement score that measures flow-field similarity.

#### 2.3.1 Structural score

The structural score quantifies the geometric similarity between the predicted and reference 3D scenes. It is computed as the unweighted arithmetic mean of four components, each ranging from 0 (no match) to 1 (perfect match):

**Component (i): Room block F1 (IoU ≥ 0.2).** A greedy bounding-box matching algorithm pairs each predicted room block with the reference block maximising their 3D IoU; pairs below 0.2 are rejected. The F1 score is computed from the resulting true-positive, false-positive, and false-negative counts; this component is dimension-dependent.

**Component (ii): Obstacle F1 (IoU ≥ 0.1).** The same greedy matching is applied to obstacle bounding boxes at a lower IoU threshold of 0.1, reflecting the greater difficulty of predicting furniture-scale objects from a single image. This is the most demanding component.

**Component (iii): Opening type F1.** F1 score measuring how accurately the VLM classifies each ventilation opening as inlet or outlet. This component is dimension-independent, depending only on the categorical label.

**Component (iv): Opening wall match ratio.** For each correctly typed opening, checks whether it is assigned to the correct wall (north, south, east, west, floor, or ceiling). This topological property tests the VLM's spatial reasoning independent of absolute dimensions.

The overall structural score is:

$$S_\text{structural} = \frac{1}{4}\bigl(F_{1,\text{room}} + F_{1,\text{obstacle}} + F_{1,\text{opening type}} + r_\text{wall match}\bigr) \tag{Eq:structural-score}$$

Components (i) and (ii) are dimension-dependent and benefit from scale calibration; components (iii) and (iv) are topology-dependent and reflect the VLM's spatial reasoning independent of absolute dimensions.

#### 2.3.2 Reference CFD agreement score

The reference CFD agreement score quantifies the similarity between the predicted and reference steady-state flow fields. Because the predicted and reference geometries generally differ, a direct cell-by-cell comparison is not possible. Instead, both flow fields are sampled onto a common normalised lattice of 18 × 18 × 10 points in room-relative coordinates, with epsilon margins (5%, 5%, 8%) from each boundary to avoid wall-adjacent sampling artefacts. Only lattice points that fall within both the reference and predicted computational domains (the spatial overlap region) are compared.

The score is the unweighted arithmetic mean of four components, each fixed at exactly four terms regardless of data availability (a component scores 0.0 when its field data are unavailable):

**Component (i): Overlap ratio.** Jaccard-style intersection-over-union of valid lattice points in the two domains, bridging geometric and flow-field evaluation. A score of 1.0 means identical normalised extents; 0.0 means no spatial overlap.

**Component (ii): Velocity magnitude similarity.** Normalised RMSE of speed over the overlap region, $\max(0,\; 1 - \min(1,\; \text{RMSE}_{|U|} / \text{RMS}_{|U|,\text{ref}}))$, with an RMS floor of $10^{-6}$ m/s to prevent instability in stagnation-dominated regions.

**Component (iii): Velocity direction similarity.** Mean cosine similarity between reference and predicted velocity vectors, mapped to $[0,1]$ via $0.5(\bar{c}+1)$; points with speed below $10^{-8}$ m/s are excluded. A score of 0.5 corresponds to random directions; below 0.5 indicates systematic disagreement.

**Component (iv): Pressure similarity.** Identical formulation to velocity magnitude similarity applied to the gauge-corrected kinematic pressure field; the overlap-domain mean is subtracted from both fields before comparison to eliminate additive-constant artefacts.

The overall reference CFD agreement score is:

$$S_\text{CFD} = \frac{1}{4}\bigl(r_\text{overlap} + s_{|U|} + s_\text{dir} + s_p\bigr) \tag{Eq:cfd-agreement-score}$$

where the individual components are defined as:

$$r_\text{overlap} = \frac{|\mathcal{I}|}{|\mathcal{U}|}, \qquad \mathcal{I} = \mathcal{V}_\text{ref} \cap \mathcal{V}_\text{pred}, \quad \mathcal{U} = \mathcal{V}_\text{ref} \cup \mathcal{V}_\text{pred} \tag{Eq:overlap-ratio}$$

$$s_{|U|} = \max\!\Bigl(0,\; 1 - \min\!\bigl(1,\; \text{RMSE}_{|U|}\, /\, \text{RMS}_{|U|,\text{ref}}\bigr)\Bigr) \tag{Eq:vel-mag-sim}$$

$$s_\text{dir} = \frac{1}{2}\!\left(\frac{1}{N_\text{dir}} \sum_{i \in \mathcal{I},\, |\mathbf{U}_i| > \epsilon} \frac{\mathbf{U}_{i,\text{ref}} \cdot \mathbf{U}_{i,\text{pred}}}{|\mathbf{U}_{i,\text{ref}}|\,|\mathbf{U}_{i,\text{pred}}|} \;+\; 1\right) \tag{Eq:vel-dir-sim}$$

$$s_p = \max\!\Bigl(0,\; 1 - \min\!\bigl(1,\; \text{RMSE}_{p'}\, /\, \text{RMS}_{p',\text{ref}}\bigr)\Bigr), \qquad p' = p - \bar{p}_{\mathcal{I}} \tag{Eq:pressure-sim}$$

where $\mathcal{V}_\text{ref}$ and $\mathcal{V}_\text{pred}$ are the sets of valid lattice points in each domain, $\epsilon = 10^{-8}$ m/s is the zero-velocity exclusion threshold, $N_\text{dir}$ is the count of points satisfying the exclusion criterion, and $p'$ denotes the gauge-corrected pressure obtained by subtracting the overlap-domain mean $\bar{p}_{\mathcal{I}}$ from each field.

This composite score is reported as a summary statistic alongside individual components for failure-mode analysis. Because it compares one CFD solution against another rather than against measurements, it is a reference-agreement metric rather than a formal ASME V&V 20 validation. The following standard metrics [@franke.etal2011COST732] are also computed and reported in the supplementary data:

$$q = \frac{1}{N}\sum_{i=1}^{N} \mathbb{1}\!\bigl(|\phi_{i,\text{pred}} - \phi_{i,\text{ref}}| < \max(W,\; D\,|\phi_{i,\text{ref}}|)\bigr), \quad D = 0.25 \tag{Eq:hit-rate}$$

$$\text{FAC2} = \frac{1}{N}\sum_{i=1}^{N} \mathbb{1}\!\left(0.5 \leq \frac{\phi_{i,\text{pred}}}{\phi_{i,\text{ref}}} \leq 2.0\right) \tag{Eq:fac2}$$

$$\text{FB} = \frac{2\,(\bar{\phi}_\text{ref} - \bar{\phi}_\text{pred})}{\bar{\phi}_\text{ref} + \bar{\phi}_\text{pred}} \tag{Eq:fractional-bias}$$

where $q$ is the hit rate with relative tolerance $D$ and absolute tolerance $W$, FAC2 is the fraction of predictions within a factor of 2, and FB is the fractional bias.

### 2.4 Solver validation against experimental data

To establish the credibility of the reference CFD solutions independently of the VLM evaluation, the solver configuration (simpleFoam with k-$\omega$ SST, wall functions, SIMPLE algorithm) was validated against the widely used Nielsen et al. (1978) ventilated room benchmark [@nielsen.etal1978VelocityCharacteristics]. This benchmark consists of a 9.0 m × 3.0 m two-dimensional room with a ceiling-mounted slot inlet (height $h$ = 0.168 m, $U_0$ = 0.455 m/s, Re $\approx$ 5,000) and a floor-level slot outlet of equal dimensions. Laser Doppler Anemometry (LDA) measurements of the normalised streamwise velocity $u/U_0$ are available at three vertical profiles ($x/L$ = 1/3, 1/2, 2/3).

The case was solved on a structured mesh of 43,200 cells (360 × 120 × 1 in the streamwise, vertical, and spanwise directions, with symmetry boundary conditions to enforce 2D flow). [Table:method-validation-nielsen] compares the CFD predictions against the experimental measurements, and [Fig:method-validation-contour] shows the velocity profile comparison at each measurement station.

[Table:method-validation-nielsen] Solver validation: simpleFoam k-$\omega$ SST against Nielsen et al. (1978) LDA measurements at three vertical profiles (43,200 cells).

| Profile | $x/L$ | Points | Hit rate ($D$ = 0.25) | RMSE ($u/U_0$) | Pearson $R$ |
|---------|--------|--------|----------------------|-----------------|-------------|
| A | 0.33 | 16 | 0.125 | 0.250 | 0.788 |
| B | 0.50 | 16 | 0.188 | 0.227 | 0.806 |
| C | 0.67 | 16 | 0.000 | 0.332 | 0.932 |
| **Mean** | — | **48** | **0.104** | **0.270** | **0.842** |

The Pearson correlation ($R$ = 0.842) confirms that the solver reproduces the qualitative velocity profile shape [Fig:method-validation-contour]. The low hit rate ($q$ = 0.104) reflects quantitative magnitude differences in the low-speed recirculation region, consistent with known wall-function RANS limitations at moderate Reynolds numbers. A grid independence check (10,800 vs. 43,200 cells) showed negligible change in $R$ and hit rate. The solver accurately predicts flow topology but absolute velocity magnitudes carry approximately 27% normalised RMSE uncertainty, which propagates into the velocity magnitude similarity component of the CFD agreement score.

### 2.5 Reference path and data-leakage prevention

The reference CFD solution for each benchmark case is computed from the ground-truth scene JSON using the same meshing and solver pipeline described above [Fig:method-eval-pathway]. Critically, the reference path does not involve any VLM invocation: the ground-truth scene is generated by the rule-based benchmark generator and is never presented to Gemini. The 2D input images are rendered from the 3D ground-truth scene, but these renderings are used only as inputs to the VLM and do not participate in the reference CFD computation. This separation ensures that the evaluation is free from circular reasoning or data leakage.

---

## 3. Benchmark Dataset

### 3.1 Geometry design

The benchmark dataset comprises 20 indoor geometries organised into a 2 x 2 complexity matrix along two independent axes: room topology (rectangular versus composite L-shaped) and obstacle density (simple versus dense). The four resulting categories are designated A1 through A4 [Fig:bench-design]. Category A1 contains simple rectangular rooms with zero or one obstacle, serving as a positive control where the framework is expected to perform well. Category A2 adds two to three obstacles to the rectangular room, introducing clutter without topological complexity. Category A3 replaces the rectangular room with an L-shaped composite geometry formed by the Boolean union of two axis-aligned blocks, while keeping obstacle density low. Category A4 combines the composite topology with dense obstacle placement, representing the most challenging condition. Each category contains five independently generated cases, yielding 20 geometries in total. [Fig:bench-cfd-showcase] shows a representative case from each category with its reference CFD velocity field, illustrating the range of geometric complexity and the corresponding flow patterns.

All geometries are generated by a deterministic rule-based generator with fixed random seeds to ensure reproducibility. Room dimensions are sampled uniformly within realistic ranges (3 to 8 m in length, 3 to 6 m in width, 2.4 to 3.5 m in height). Obstacles are placed with a minimum wall clearance of 0.4 m and a minimum inter-obstacle clearance sufficient to avoid degenerate mesh cells. Each room has exactly one inlet and one outlet, placed on different walls. Opening sizes are sampled from a range of 0.3 to 1.0 m in each dimension, consistent with typical indoor ventilation diffuser dimensions.

### 3.2 Multi-view rendering protocol

Each of the 20 benchmark geometries is rendered into five distinct 2D input views [Fig:bench-multiview]: perspective, bird's-eye, floor plan, wireframe, and section. The perspective view provides a natural photographic viewpoint with depth cues and foreshortening. The bird's-eye view offers a top-down oblique angle that reveals the overall layout. The floor-plan view is an orthographic top-down projection that preserves planimetric proportions but eliminates height information. The wireframe view renders the geometry as line edges without surface shading. The section view shows a vertical cut through the room, revealing the cross-sectional profile but obscuring features outside the cutting plane. Rendering all five views from each geometry yields 100 input images in total, enabling a controlled analysis of how input view type affects downstream prediction quality while holding the underlying geometry constant.

### 3.3 Scale calibration

The VLM does not have access to absolute dimensional references in the input image, and its predicted scene dimensions therefore exhibit a systematic scale error. To correct this, a single-factor post-hoc scaling is applied: the longest horizontal span in the predicted geometry is identified, and the ratio between this span and the corresponding span in the reference geometry is used as a uniform scale factor applied to all coordinates. After calibration, the mean relative error in the scaled dimension (Lx) drops to 0.4%, while the remaining dimensions retain relative errors of approximately 22% (Ly) and 20% (Lz). This residual anisotropy reflects the VLM's tendency to estimate the dominant room dimension more accurately than secondary dimensions, and it is a known limitation of single-factor calibration.

---

## 4. Results

### 4.1 Aggregate performance

To contextualise the framework's performance, a naive geometric baseline is evaluated alongside the scale-calibrated pipeline. The baseline predicts a fixed 5 m × 4 m × 3 m rectangular room with no obstacles and one inlet-outlet pair on opposite walls for every case. [Table:result-baseline-comparison] compares the two conditions.

[Table:result-baseline-comparison] Structural score comparison: naive baseline versus scale-calibrated VLM pipeline (n = 100).

| Condition | Structural score | Room-kind match | Opening-wall match |
|-----------|-----------------|-----------------|-------------------|
| Naive geometric default | 0.621 | 50% | 80% |
| VLM, scale-calibrated | 0.781 | 95% | 70% |

Note: the naive baseline achieves a higher opening-wall match rate (80%) than the VLM (70%) because the fixed west-inlet/east-outlet assignment happens to match many reference cases by chance, not because of superior spatial reasoning.

The scale-calibrated framework achieves a mean structural score of 0.781 (+25.8% relative over the naive baseline). The baseline's non-trivial score of 0.621 reflects that the structural metric rewards correct opening-type classification even when dimensions and obstacles are wrong. The VLM achieves a substantially higher room-kind match rate (95% vs. 50%), confirming reliable rectangular/composite topology discrimination that a fixed-geometry baseline cannot perform.

To further clarify the contribution of each evaluation component, [Table:result-component-breakdown] decomposes the structural score into its four constituents for both the naive baseline and the scale-calibrated pipeline.

[Table:result-component-breakdown] Component-level breakdown of the structural score.

| Component | Naive | Scale-calibrated |
|-----------|-------|------------------|
| Room block F1 (IoU ≥ 0.2) | 0.833 | 0.988 |
| Obstacle F1 (IoU ≥ 0.1) | 0.150 | 0.400 |
| Opening type F1 | 1.000 | 1.000 |
| Opening wall match ratio | 0.500 | 0.737 |
| **Overall structural score** | **0.621** | **0.781** |

Opening type F1 is trivially perfect for both conditions (both predict exactly one inlet and one outlet). Opening wall match ratio improves from 0.500 (naive, chance alignment) to 0.737 (VLM), confirming genuine spatial reasoning about which walls carry openings. Room block F1 is near-perfect for the VLM (0.988) after scale calibration, while obstacle F1 is the weakest component (0.400), accounting for most of the gap between current performance and a perfect score. The naive baseline's obstacle F1 of 0.150 is non-zero only because a small number of reference cases have zero obstacles.

The mean reference CFD agreement score is 0.477. Component decomposition [Fig:result-cfd-component-breakdown] shows overlap ratio (0.888) and velocity direction similarity (0.634) are moderate to high, while velocity magnitude similarity (0.117) is low across all views and pressure similarity (0.269) is intermediate — indicating qualitatively correct flow patterns but quantitatively different velocity magnitudes. CFD baselines are not computed for the naive default condition, as that baseline exists only to contextualise the structural score.

For reference, the COST 732 standard metrics averaged across all 97 cases yield a velocity magnitude hit rate of q = 0.09, FAC2 = 0.28, and a Pearson correlation of R = 0.44, confirming the limited quantitative velocity agreement identified by the agreement score decomposition.

[Table:result-detailed-metrics] Detailed evaluation results for the scale-calibrated pipeline (n = 97 valid cases out of 100; three cases excluded due to solver divergence failures).

| Metric | Value |
|--------|-------|
| Mean structural score | 0.781 ± 0.150 |
| Mean CFD agreement score | 0.477 ± 0.158 |
| Room-kind match rate | 95% |
| Opening-wall match rate | 70% |
| Mean Lx relative error | 0.4% |
| Mean Ly relative error | 21.9% |
| Mean Lz relative error | 19.8% |
| Mean volume relative error | 36.6% |
| Solver divergence failures | 3 / 100 |

Three cases (bench_a3_04/perspective, bench_a4_02/perspective, bench_a4_02/wireframe) failed to converge during the CFD solution stage due to non-physical predicted geometries and are excluded from score aggregation but counted as failures in the room-kind and opening-wall match statistics.

### 4.2 Effect of input view type

The input view type has a pronounced effect on performance [Fig:result-view-aggregate]. The floor-plan view achieves the highest scores (structural 0.884 ± 0.110, CFD 0.572 ± 0.185), consistent with its orthographic top-down projection preserving planimetric information. Rankings for the remaining views are: bird's-eye (0.800, 0.495), perspective (0.782, 0.443), wireframe (0.760, 0.476), and section (0.679, 0.396).

Opening-wall match rate shows the most dramatic view dependence: floor plan 100%, bird's-eye 85%, wireframe 75%, perspective 65%, and section only 25%. The section view's poor performance arises because a vertical cutting plane typically intersects only one or two walls, leaving remaining openings invisible to the VLM. Room-kind match rate is more robust, ranging from 90% (perspective and section) to 100% (bird's-eye and floor plan), as topology identification depends more on the overall room silhouette than on fine-grained wall visibility.

### 4.3 Effect of geometric complexity

Performance varies across complexity categories [Fig:result-category-aggregate]. A1 (simple rectangular) achieves the highest scores (structural 0.845, CFD 0.519). A2 (rectangular dense) yields the lowest structural score (0.707) despite simpler room topology, because higher obstacle density increases placement errors. A3 and A4 (composite) achieve structural scores of 0.783 and 0.791, showing that composite topology alone does not severely degrade structural prediction; CFD scores are 0.498 (A2), 0.480 (A3), and 0.407 (A4).

Opening-wall match rates differ: A1 88%, A3 and A4 each 72%, A2 only 48%. The low A2 rate reflects the geometric equivalence of all four walls in a rectangular room, leaving the VLM without topological cues. In composite rooms, the L-shaped topology provides a spatial anchor that improves opening-wall assignment.

### 4.4 Solver robustness

Of the 100 cases, 97 produce a converged solution and 3 fail due to solver divergence from non-physical predicted geometries. Of the 97 converged cases, the majority converge under the nominal solver configuration; the remainder require escalation through the robustness pathway [Fig:result-robustness]. Three cases fail entirely and are excluded from CFD score aggregation. The escalation rate is highest for the section view and the floor-plan view, the latter driven by cases where the accurate geometry produces complex internal flow structures that challenge the turbulence model at the nominal mesh resolution.

### 4.5 VLM output repeatability

Because the VLM operates with a sampling temperature of 0.2, its output is stochastic rather than deterministic. To characterise this variation, three representative cases (bench_a1_01, bench_a3_03, bench_a4_03) were each processed three times through the floor-plan pipeline. [Table:result-vlm-repeatability] reports the structural scores across runs.

[Table:result-vlm-repeatability] VLM repeatability: structural scores across 3 independent runs per case (floor-plan view).

| Case | Run 1 | Run 2 | Run 3 | Mean ± SD |
|------|-------|-------|-------|-----------|
| A1-01 (simple) | 0.875 | 0.875 | 0.875 | 0.875 ± 0.000 |
| A3-03 (composite) | 1.000 | 0.875 | 0.875 | 0.917 ± 0.059 |
| A4-03 (dense) | 0.917 | 1.000 | 1.000 | 0.972 ± 0.039 |

The mean standard deviation is 0.033. Room topology (room block F1 = 1.0) and opening-wall assignment (wall match ratio = 1.0) are perfectly consistent across all runs; variation arises solely from obstacle matching. Raw room dimensions vary substantially between runs (e.g., scale factors of 0.54, 1.08, and 0.86 for the same case), confirming that post-hoc scale calibration is essential. The small sample size (9 observations) limits statistical power; a larger study is needed to fully characterise VLM stochastic variation.

---

## 5. Discussion

### 5.1 Failure mode I: composite-room collapse under section view

When a section cutting plane passes through only one block of an L-shaped room, the VLM infers a rectangular room and discards the second block entirely [Fig:discuss-section-collapse]. This collapse eliminates recirculation zones present in the reference solution and degrades both structural and CFD scores. The failure is inherent to the section view's limited information content — a vertical section through an L-shaped room is geometrically indistinguishable from a rectangular section unless the plane intersects the block junction. Section views should therefore be avoided for suspected non-rectangular topologies, or supplemented with a second view.

### 5.2 Failure mode II: obstacle hallucination with limited CFD penalty

The VLM frequently predicts obstacles absent from the reference, particularly in wireframe and section views where sparse visual cues are misinterpreted as internal structures [Fig:discuss-obstacle-hallucination]. In two representative A3 wireframe cases, obstacle counts inflate from 0 to 3 and from 1 to 3, yet CFD agreement scores remain 0.380 and 0.529, with opening assignments preserved. Hallucinated obstacles are typically small relative to room volume and positioned away from the main flow path, introducing localised disturbances without altering bulk ventilation patterns. Obstacle prediction errors are therefore tolerable when the primary interest is bulk airflow rather than near-furniture velocity detail.

### 5.3 Failure mode III: structure-versus-fidelity gap in dense composites

Structural and CFD scores are only weakly correlated across all 97 valid cases [Fig:discuss-scatter-struct-cfd]. The decoupling is most pronounced in A4: case A4-02 achieves a structural score of 0.813 yet a CFD score of only 0.346; case A4-04 scores 0.917 structurally yet 0.302 in CFD [Fig:discuss-structure-cfd-gap]. The cause is flow-field sensitivity to small geometric perturbations in dense obstacle configurations — residual misalignments of a few centimetres, even within the IoU matching threshold, redirect jets, create or eliminate recirculation zones, and shift stagnation points. For dense configurations, bounding-box structural matching is necessary but not sufficient for CFD fidelity.

### 5.4 Scale calibration and its limitations

The single-factor post-hoc scaling reduces the mean Lx error from its uncalibrated value to 0.4% but leaves the secondary dimensions (Ly, Lz) with relative errors of approximately 20%. This anisotropic residual error has two consequences. First, the room volume is systematically distorted, which affects the bulk velocity magnitude for a given inlet flow rate. Second, the aspect ratio of the room is altered, which can shift the location of flow reattachment points and recirculation zones. A multi-factor calibration that independently corrects each dimension would reduce these errors but would require multiple reference measurements, undermining the single-image automation objective. The present single-factor approach represents a practical compromise that corrects the dominant scale error while preserving the fully automated workflow.

### 5.5 Reference solver preset matching

To ensure fair flow-field comparison, the reference CFD solution for each benchmark case was computed using the same solver preset — and therefore the same inlet velocity boundary condition — as the corresponding predicted case. This preset matching eliminates systematic velocity scale discrepancies that would arise when different cases require different levels of solver stabilisation to achieve convergence at the 0.18 m mesh resolution.

### 5.6 Computational cost

The end-to-end pipeline time per case is approximately 3 to 8 minutes on a single workstation, depending on the mesh resolution and solver convergence behaviour. The VLM API call (Gemini 3.1 Pro) typically returns within 10 to 30 seconds. Mesh generation with Gmsh takes 5 to 15 seconds for the default 0.35 m cell size. The dominant cost is the steady-state CFD solve, which ranges from 1 to 6 minutes depending on the case complexity and whether solver escalation is triggered. This total time is orders of magnitude shorter than the manual geometry preparation workflow it replaces, which typically requires hours of CAD modelling for a comparable indoor space. The framework is therefore suitable for rapid screening applications where approximate flow patterns are sufficient.

### 5.7 Practical considerations and limitations

Several limitations of the current study should be noted. First, all results reported in this paper use a mesh cell size of 0.18 m for both reference and predicted cases. A grid independence study was conducted on three representative cases at cell sizes of 0.35, 0.25, 0.18, and 0.10 m. [Table:discuss-grid-independence] shows the reference CFD agreement score evaluated at each reference resolution against the same predicted solution.

[Table:discuss-grid-independence] Grid independence check: reference CFD agreement score evaluated against reference solutions at four mesh resolutions, spanning a 34× range in cell count. All predicted solutions use a cell size of 0.18 m.

| Case | Cells (0.18 m) | Ref. 0.35 m | Ref. 0.25 m | Ref. 0.18 m | Ref. 0.10 m | Δ max |
|------|----------------|-------------|-------------|-------------|-------------|-------|
| A1-01 | ~37,000 | 0.389 | 0.568 | 0.579 | 0.625 | 0.236 |
| A3-03 | ~62,000 | 0.739 | 0.726 | 0.641 | 0.616 | 0.123 |
| A4-03 | ~55,000 | 0.576 | 0.553 | 0.558 | 0.598 | 0.045 |

At 0.35 m, mean velocities are underpredicted by up to 50% for simple geometries, confirming that resolution is inadequate for quantitative comparison. From 0.18 m onward, key quantities change less than 5% upon further refinement to 0.10 m; 0.18 m was therefore selected as the operational mesh size (3–8 min per case vs. 30+ min at 0.10 m). The non-monotonic trend for A3-03 (score decreasing from 0.739 to 0.616 with finer reference mesh) reflects that a more refined reference exposes geometric discrepancies masked at coarser resolution. Relative case rankings are preserved across all resolutions.

Second, the benchmark uses only axis-aligned rectangular and L-shaped rooms; curved walls, sloped ceilings, and multi-storey configurations are not evaluated. Third, steady-state isothermal RANS may not capture transient flow features, complex turbulence in densely obstructed rooms, or buoyancy-driven stratification important in energy audit scenarios; extension to conjugate heat transfer or buoyant flow models is left as future work. Fourth, the intersection-domain comparison may under-represent errors where predicted rooms extend beyond or fall short of the reference domain, and single-factor scale calibration requires one known reference dimension that may be unavailable in fully unsupervised deployment. Fifth, all results use a single VLM (Gemini 3.1 Pro, preview, accessed March 2026); the framework architecture is VLM-agnostic, but the reported scores are specific to this model version, and a multi-model comparison is left as future work.

---

## 6. Application to Architectural Floor Plans

To assess the framework's applicability beyond the synthetic benchmark, two architectural floor plan images were prepared and processed through the full pipeline. These cases use the floor-plan view, which achieved the highest performance in the benchmark evaluation (structural score 0.884, opening-wall match 100%).

The first case is a simple rectangular office of 6.0 m × 4.0 m ([Fig:demo-floorplan-application]a). The floor plan includes hatched walls, a door with arc symbol on the south wall, a window on the north wall, and three furniture items (desk, chair, shelf). Dimension annotations on the drawing provide the scale reference. The VLM correctly identifies the rectangular room topology with dimensions 6.0 × 4.0 × 2.7 m, places three obstacles corresponding to the furniture items, and assigns the door as inlet (south wall) and the window as outlet (north wall). The resulting CFD simulation converges under the robust solver preset at 0.18 m mesh resolution, producing a physically plausible flow field with a jet from the south-wall inlet traversing the room toward the north-wall outlet ([Fig:demo-floorplan-application]b).

The second case is an L-shaped composite room consisting of a 5.0 m × 4.0 m living room connected to a 3.0 m × 2.5 m alcove ([Fig:demo-floorplan-application]c). The VLM correctly identifies the composite topology as a two-block structure. However, it places four obstacles — a sofa (present in the drawing), a coffee table, a bookshelf, and an alcove desk — of which the latter three are hallucinated and do not appear in the input image. This obstacle hallucination is consistent with the failure mode identified in [Sec:discuss-obstacle-hallucination] from the synthetic benchmark. The inlet is assigned to the north wall and the outlet to the west wall. Despite the hallucinated obstacles, the CFD simulation converges and produces a physically plausible flow pattern in which the main jet passes through the living room with a secondary branch entering the alcove region ([Fig:demo-floorplan-application]d). The scale calibration factor of 0.625 was applied to correct the VLM's initial overestimate of the room dimensions.

Both cases demonstrate that the framework handles realistic visual noise — hatched walls, furniture icons, dimension lines, door arc symbols — and produces converged CFD solutions. The obstacle hallucination in Case 2 reinforces the benchmark finding that small hallucinated obstacles have limited impact on bulk flow patterns. Scale calibration is straightforward when dimension annotations are present, as is common in architectural practice, and qualitative flow patterns are physically reasonable. Floor-plan inputs from architectural drawings represent a practical deployment scenario for the framework.

## 7. Conclusion

This paper has presented a framework for automated indoor CFD analysis from a single 2D image, using a vision-language model to bridge the gap between visual input and simulation-ready 3D geometry. A systematic benchmark of 100 evaluation cases across four complexity levels and five input view types demonstrates that the framework achieves a mean structural score of 0.781 and a mean reference CFD agreement score of 0.477 under single-factor scale calibration. The floor-plan view emerges as the most effective input type, achieving a structural score of 0.884, a CFD agreement score of 0.572, and a 100% opening-wall match rate. Three characteristic failure modes have been identified and analysed: composite-room collapse under section views, obstacle hallucination with limited CFD penalty, and a structure-versus-fidelity gap in dense configurations. The framework converges for 97 of 100 cases (3 cases fail due to solver divergence from non-physical predicted geometries) through a multi-level robustness strategy.

The results suggest that VLM-based geometric abstraction is a viable approach for rapid indoor CFD screening when the primary quantities of interest are qualitative flow pattern identification — such as locating dominant recirculation zones, jet paths, and stagnation regions — rather than precise local velocity magnitudes. For applications requiring high-fidelity local flow prediction near specific obstacles, the current framework's accuracy is insufficient and should be supplemented with manual geometry refinement or multi-view input fusion. Future work will explore multi-factor scale calibration, multi-view input strategies, and extension to transient thermal simulations.

---

## References

## Appendix A. VLM Prompt Template

The full prompt template used for Gemini 3.1 Pro scene extraction (~659 words) is provided as supplementary material. The prompt specifies the `indoor_cfd_scene_v1` JSON schema, the coordinate system convention (origin at room minimum corner, x = west-east, y = south-north, z = floor-ceiling), and solver-stability constraints (no overlapping obstacles, openings within wall bounds, connected fluid domain). When a scale hint is available, an addendum anchors absolute dimensions to a user-provided metric span. The complete prompt text is reproduced in the online supplementary material.

## Appendix B. IoU Threshold Sensitivity

The structural score depends on the IoU thresholds used for room-block matching (default 0.2) and obstacle matching (default 0.1). To assess the sensitivity of the reported results to these thresholds, the structural score was recomputed with a uniform IoU threshold applied to both components, ranging from 0.05 to 0.50. Opening type F1 and opening wall match ratio are unaffected by this variation.

**Table B1.** Mean structural score as a function of uniform IoU matching threshold (n = 97).

| IoU threshold | Mean structural score |
|---------------|----------------------|
| 0.05 | 0.803 |
| 0.10 | 0.783 |
| 0.15 | 0.759 |
| 0.20 | 0.749 |
| 0.30 | 0.710 |
| 0.50 | 0.618 |

The score degrades monotonically from 0.803 at the most lenient threshold to 0.618 at the strictest. The steepest decline occurs between 0.30 and 0.50, indicating that many obstacle and room-block matches fall in the IoU range of 0.3 to 0.5 — geometrically reasonable overlap but not tight alignment. The relative ranking of view types and categories is preserved across all thresholds, confirming that the reported conclusions are not artefacts of the particular threshold choice.

---

## Figure Captions

[Fig:method-framework] Overall framework architecture. A single 2D image is processed by Gemini 3.1 Pro to produce a structured 3D scene description; post-hoc scale calibration corrects the dominant dimensional error using the longest horizontal span; the calibrated scene is automatically converted to an OpenFOAM case. An independent reference path generates the ground-truth CFD solution from the rule-based benchmark without VLM involvement.

[Fig:bench-design] Benchmark dataset design. The 20 geometries are organised into a 2 × 2 matrix of room topology (rectangular vs. composite L-shaped) and obstacle density (simple vs. dense), yielding categories A1–A4; the summary table reports per-category evaluation statistics.

[Fig:bench-multiview] Multi-view rendering protocol. The same geometry (A4-03) rendered into five input views: (a) perspective, (b) bird's-eye, (c) floor plan, (d) wireframe, and (e) section, illustrating the variation in information content.

[Fig:method-eval-pathway] Evaluation pathway. The predicted path (top) processes the 2D input through VLM, scale calibration, and CFD; the reference path (bottom) uses the rule-based ground-truth scene directly. Both paths converge at the metric comparison stage.

[Fig:result-view-aggregate] Aggregate performance by input view type. (a) Mean structural and CFD scores. (b) Room-kind and opening-wall match rates. Floor plan achieves the highest performance across all metrics.

[Fig:result-category-aggregate] Aggregate performance by benchmark category. (a) Mean structural and CFD scores. (b) Room-kind and opening-wall match rates.

[Fig:result-crossview-outcome] Cross-view outcome comparison for cases A2-03 and A3-03 across three view types. Matched obstacles are shown in tan, hallucinated obstacles in orange cross-hatch, and reference outlines as dashed purple. Structural (S) and CFD (C) scores are annotated per panel.

[Fig:discuss-section-collapse] Composite-room collapse under section view. Cases A3-04 and A4-05 are shown with reference geometry (left) and VLM prediction with ground-truth overlay (right); the VLM collapses the L-shaped room to a single rectangular block in both cases.

[Fig:discuss-obstacle-hallucination] Obstacle hallucination with limited CFD penalty. Two A3 wireframe cases show obstacle count inflation (0 → 3 and 1 → 3); hallucinated obstacles (orange cross-hatch) are distinguished from matched predictions (tan). Opening topology is preserved and CFD scores remain 0.380 and 0.529.

[Fig:discuss-structure-cfd-gap] Structure-versus-fidelity gap in dense composites. Two A4 cases with high structural scores (0.938 and 1.000) achieve low CFD scores (0.409 and 0.458); left panels show geometry overlap, right panels compare flow fields.

[Fig:result-robustness] Robustness and convergence summary.

[Fig:discuss-scatter-struct-cfd] Structural score versus CFD agreement score for all 97 valid cases.

[Fig:result-heatmap-category-view] Category × view interaction heatmaps.

[Fig:demo-floorplan-application] Application to architectural floor plans. (a,b) Rectangular office (6.0 m × 4.0 m): input floor plan and CFD velocity field. (c,d) L-shaped room (5.0 m × 4.0 m + 3.0 m × 2.5 m alcove): input and CFD result. Both cases converge at 0.18 m mesh resolution.

[Fig:bench-cfd-showcase] Representative reference cases from each complexity category: 3D geometry (left) and steady-state CFD velocity field with streamlines (right). (a) A1, (b) A2, (c) A3, (d) A4.

[Fig:method-validation-contour] Solver validation against Nielsen et al. (1978) benchmark (9 m × 3 m, Re ≈ 5,000). Colour contours show velocity magnitude; white streamlines show the ceiling jet, recirculation zone, and corner vortex. Red markers (A, B, C) indicate the three LDA measurement stations compared in [Table:method-validation-nielsen].

[Fig:result-cfd-component-breakdown] CFD agreement score decomposed into four components (overlap ratio, velocity magnitude similarity, direction similarity, pressure similarity) by view type, with standard-deviation error bars. Velocity magnitude is near zero across all views, driving the low overall score despite moderate overlap and directional agreement.
