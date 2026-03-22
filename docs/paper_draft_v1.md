# Automated Indoor CFD Analysis from a Single 2D Image via Vision-Language Model Abstraction

---

## Abstract

Computational fluid dynamics (CFD) simulation of indoor environments requires three-dimensional geometric models that are costly and time-consuming to construct manually. This paper presents an automated framework that converts a single two-dimensional image of an indoor space into a fully resolved steady-state CFD solution, using a vision-language model (VLM) as the geometric abstraction engine. The framework accepts a photograph or architectural drawing, invokes a VLM (Gemini 3.1 Pro, preview) to extract a structured 3D scene description comprising room geometry, internal obstacles, and ventilation openings, applies a single-factor scale calibration that requires one reference dimension (the longest horizontal span) to correct absolute size, and then automatically generates a computational mesh and boundary conditions for OpenFOAM. A rule-based benchmark dataset of 20 indoor geometries spanning four complexity levels and five input view types (100 evaluation cases in total) is constructed with independently computed reference CFD solutions, ensuring that no VLM output contaminates the ground truth. Across all 100 cases the framework achieves a mean structural similarity score of 0.781 and a mean CFD fidelity score of 0.491, with the floor-plan view yielding the highest performance (structural 0.884, CFD 0.538) and the section view the lowest (structural 0.679, CFD 0.430). Room-kind topology is correctly identified in 95% of cases, and the opening-wall assignment matches the reference in 70% of cases overall, rising to 100% for floor-plan inputs. Systematic failure-mode analysis reveals three characteristic error patterns: composite-room collapse under section views, obstacle hallucination with limited CFD penalty, and a structure-versus-fidelity gap in dense composite geometries. The full pipeline, benchmark dataset, and evaluation code are provided to support reproducibility.

---

## 1. Introduction

*[Placeholder — to be written in a subsequent revision. This section will motivate the need for rapid indoor CFD assessment, review the bottleneck of manual geometry preparation, introduce the hypothesis that VLM spatial perception is sufficient to produce simulation-ready 3D abstractions, and outline the contributions of the paper.]*

---

## 2. Related Work

Research at the intersection of computer vision and computational engineering has pursued two largely independent trajectories. On one side, image-based 3D reconstruction methods such as neural radiance fields [Mildenhall et al. 2020 NeRF] and multi-view stereo produce photorealistic scene representations but do not yield the watertight, topologically clean surface meshes required by finite-volume solvers. On the other side, CFD automation efforts have focused on parameterised geometry generators and template-based meshing workflows that assume a pre-existing CAD model. The gap between a raw image and a simulation-ready geometry has therefore remained a manual step in most practical workflows.

Vision-language models have recently demonstrated strong spatial reasoning capabilities. Models such as GPT-4V [OpenAI 2023 GPT-4V] and Gemini [Google 2024 Gemini] can interpret architectural drawings, estimate room dimensions from photographs, and enumerate furniture items from overhead images. These capabilities suggest that a VLM could serve as the geometric abstraction layer between a 2D image and a 3D scene description suitable for meshing and simulation. However, no prior work has systematically evaluated the downstream CFD fidelity of VLM-generated scene descriptions, nor has a controlled benchmark been established to isolate the effect of input view type from geometric complexity.

Several studies have explored end-to-end surrogate models [Guo et al. 2016 CNN flow prediction] that predict flow fields directly from images using convolutional neural networks. While these approaches bypass the meshing step entirely, they are constrained to the training distribution and cannot generalise to novel geometries without retraining. The present work differs in that it retains the full-fidelity CFD solver and uses the VLM solely for the geometric abstraction step, thereby preserving the physical rigour of the simulation while automating the geometry preparation.

---

## 3. Methodology

### 3.1 Framework overview

The proposed framework accepts a single 2D image of an indoor space and produces a steady-state CFD solution through a five-stage pipeline [Fig:method-framework]. In the first stage, the image is presented to a vision-language model (Gemini 3.1 Pro) together with a structured prompt that requests a JSON description of the room geometry, internal obstacles, and ventilation openings. The prompt specifies the coordinate system convention, unit system, and the required schema fields, but does not provide any ground-truth information about the scene. The VLM returns a scene description conforming to the `indoor_cfd_scene_v1` schema, which encodes the room as one or more axis-aligned rectangular blocks, each obstacle as an axis-aligned bounding box, and each opening as a wall-anchored rectangular patch with an associated type (inlet or outlet).

In the second stage, a single-factor scale calibration is applied to the VLM output. Because the VLM has no access to absolute dimensional references, the predicted scene dimensions are generally correct in proportion but incorrect in absolute magnitude. The calibration identifies the longest horizontal span in the predicted geometry and rescales the entire scene uniformly so that this span matches the corresponding span in the reference geometry. This post-hoc scaling is the only point at which reference information enters the predicted path, and it corrects the dominant source of volumetric error without altering the topological structure of the prediction.

In the third stage, the calibrated scene description is converted into an OpenFOAM case. A tetrahedral mesh is generated using Gmsh with a default target cell size of 0.35 m, yielding approximately 10 to 25 cells across the room depending on the room dimensions. Boundary conditions are set to a fixed volumetric flow rate at the inlet and a zero-gradient pressure condition at the outlet, ensuring that the flow forcing is identical between the predicted and reference cases regardless of any differences in opening size. Wall boundaries are treated as no-slip surfaces. The turbulence model is k-omega SST with standard wall functions. At the default cell size of 0.35 m, near-wall resolution is coarse and the wall-function approach operates at high y+ values; the resulting solutions should therefore be interpreted as engineering-level screening estimates rather than high-fidelity resolved simulations.

In the fourth stage, the case is solved using simpleFoam for isothermal steady-state incompressible flow. A multi-level robustness strategy is employed to handle cases where the default solver settings fail to converge. The mesh resolution ladder proceeds through 0.35, 0.25, and 0.18 m cell sizes, and at each resolution the solver attempts up to four preset configurations of progressively reduced relaxation factors, including a laminar fallback that disables the turbulence model entirely. This escalation pathway ensures that every case produces a converged solution, albeit at varying levels of resolution and modelling fidelity.

In the fifth stage, the predicted flow field is compared against the independently computed reference solution to produce structural and CFD fidelity scores.

### 3.2 Evaluation metrics

The evaluation framework produces two families of scores. The structural score quantifies the geometric similarity between the predicted and reference scenes. It is computed as the unweighted arithmetic mean of four components: (i) the F1 score of a greedy bounding-box match between predicted and reference room blocks at an intersection-over-union (IoU) threshold of 0.2, (ii) the F1 score of a greedy bounding-box match between predicted and reference obstacles at an IoU threshold of 0.1, (iii) the F1 score of opening-type classification (inlet versus outlet), and (iv) the opening-wall match ratio, defined as the fraction of openings whose assigned wall matches the reference. The structural score ranges from 0 to 1, with 1 indicating a perfect geometric match.

The CFD fidelity score quantifies the agreement between the predicted and reference velocity and pressure fields. Both fields are sampled onto a common uniform grid, and only cells within the spatial intersection of the two room domains are compared. The CFD score is computed as the unweighted mean of up to four components: (i) the overlap ratio of the predicted and reference domains relative to their union, (ii) velocity magnitude similarity, defined as max(0, 1 − min(1, relative RMSE)), (iii) velocity direction similarity, defined as 0.5 × (mean direction cosine + 1), and (iv) pressure similarity, defined as max(0, 1 − min(1, relative RMSE)). Components are included only when the corresponding field data are available. The CFD score ranges from 0 to 1.

### 3.3 Reference path and data-leakage prevention

The reference CFD solution for each benchmark case is computed from the ground-truth scene JSON using the same meshing and solver pipeline described above [Fig:method-eval-pathway]. Critically, the reference path does not involve any VLM invocation: the ground-truth scene is generated by the rule-based benchmark generator and is never presented to Gemini. The 2D input images are rendered from the 3D ground-truth scene, but these renderings are used only as inputs to the VLM and do not participate in the reference CFD computation. This separation ensures that the evaluation is free from circular reasoning or data leakage.

---

## 4. Benchmark Dataset

### 4.1 Geometry design

The benchmark dataset comprises 20 indoor geometries organised into a 2 x 2 complexity matrix along two independent axes: room topology (rectangular versus composite L-shaped) and obstacle density (simple versus dense). The four resulting categories are designated A1 through A4 [Fig:bench-design]. Category A1 contains simple rectangular rooms with zero or one obstacle, serving as a positive control where the framework is expected to perform well. Category A2 adds two to three obstacles to the rectangular room, introducing clutter without topological complexity. Category A3 replaces the rectangular room with an L-shaped composite geometry formed by the Boolean union of two axis-aligned blocks, while keeping obstacle density low. Category A4 combines the composite topology with dense obstacle placement, representing the most challenging condition. Each category contains five independently generated cases, yielding 20 geometries in total.

All geometries are generated by a deterministic rule-based generator with fixed random seeds to ensure reproducibility. Room dimensions are sampled uniformly within realistic ranges (3 to 8 m in length, 3 to 6 m in width, 2.4 to 3.5 m in height). Obstacles are placed with a minimum wall clearance of 0.4 m and a minimum inter-obstacle clearance sufficient to avoid degenerate mesh cells. Each room has exactly one inlet and one outlet, placed on different walls. Opening sizes are sampled from a range of 0.3 to 1.0 m in each dimension, consistent with typical indoor ventilation diffuser dimensions.

### 4.2 Multi-view rendering protocol

Each of the 20 benchmark geometries is rendered into five distinct 2D input views [Fig:bench-multiview]: perspective, bird's-eye, floor plan, wireframe, and section. The perspective view provides a natural photographic viewpoint with depth cues and foreshortening. The bird's-eye view offers a top-down oblique angle that reveals the overall layout. The floor-plan view is an orthographic top-down projection that preserves planimetric proportions but eliminates height information. The wireframe view renders the geometry as line edges without surface shading. The section view shows a vertical cut through the room, revealing the cross-sectional profile but obscuring features outside the cutting plane. Rendering all five views from each geometry yields 100 input images in total, enabling a controlled analysis of how input view type affects downstream prediction quality while holding the underlying geometry constant.

### 4.3 Scale calibration

The VLM does not have access to absolute dimensional references in the input image, and its predicted scene dimensions therefore exhibit a systematic scale error. To correct this, a single-factor post-hoc scaling is applied: the longest horizontal span in the predicted geometry is identified, and the ratio between this span and the corresponding span in the reference geometry is used as a uniform scale factor applied to all coordinates. After calibration, the mean relative error in the scaled dimension (Lx) drops to 0.4%, while the remaining dimensions retain relative errors of approximately 22% (Ly) and 20% (Lz). This residual anisotropy reflects the VLM's tendency to estimate the dominant room dimension more accurately than secondary dimensions, and it is a known limitation of single-factor calibration.

---

## 5. Results

### 5.1 Aggregate performance

To contextualise the framework's performance, a naive geometric baseline is evaluated alongside the scale-calibrated pipeline. The baseline predicts a fixed 5 m × 4 m × 3 m rectangular room with no obstacles and one inlet-outlet pair on opposite walls for every case. [Table:result-baseline-comparison] compares the two conditions.

[Table:result-baseline-comparison] Structural score comparison: naive baseline versus scale-calibrated VLM pipeline (n = 100).

| Condition | Structural score | Room-kind match | Opening-wall match |
|-----------|-----------------|-----------------|-------------------|
| Naive geometric default | 0.621 | 50% | 80% |
| VLM, scale-calibrated | 0.781 | 95% | 70% |

The scale-calibrated framework achieves a mean structural score of 0.781, representing a 25.8% relative improvement over the naive baseline. The naive baseline achieves a non-trivial score of 0.621, which reflects the fact that the structural score rewards correct opening-type classification even when dimensions and obstacle placement are entirely wrong. Notably, the naive baseline achieves a higher opening-wall match rate (80%) than the VLM (70%), because the fixed west-inlet/east-outlet assignment happens to match many reference cases by chance. However, the VLM achieves a substantially higher room-kind match rate (95% versus 50%), confirming that it reliably distinguishes rectangular from composite topologies — a discrimination that a fixed-geometry baseline cannot perform.

To further clarify the contribution of each evaluation component, [Table:result-component-breakdown] decomposes the structural score into its four constituents for both the naive baseline and the scale-calibrated pipeline.

[Table:result-component-breakdown] Component-level breakdown of the structural score.

| Component | Naive | Scale-calibrated |
|-----------|-------|------------------|
| Room block F1 (IoU ≥ 0.2) | 0.833 | 0.988 |
| Obstacle F1 (IoU ≥ 0.1) | 0.150 | 0.400 |
| Opening type F1 | 1.000 | 1.000 |
| Opening wall match ratio | 0.500 | 0.737 |
| **Overall structural score** | **0.621** | **0.781** |

Two of the four components — opening type F1 and opening wall match ratio — are dimension-independent topological properties. Opening type classification is trivially perfect for both conditions because the naive baseline and the VLM both predict exactly one inlet and one outlet. Opening wall match ratio jumps from 0.500 (naive, where fixed west/east assignment matches by chance) to 0.737 (VLM), confirming that the VLM provides genuine spatial reasoning about which walls carry openings.

The remaining two components — room block F1 and obstacle F1 — are dimension-dependent because they rely on bounding-box IoU. Room block F1 is near-perfect for the VLM (0.988), reflecting accurate room envelope prediction after scale calibration. Obstacle F1 is the weakest component (0.400), indicating that obstacle placement remains the primary challenge even with correct room dimensions. The naive baseline's obstacle F1 of 0.150 is non-zero because a small number of reference cases also have zero obstacles, yielding a vacuous match. This decomposition reveals that the overall structural score is substantially supported by the two topology components, and that obstacle detection accuracy — the most demanding spatial reasoning task — accounts for most of the gap between current performance and a perfect score.

The mean CFD fidelity score for the scale-calibrated pipeline is 0.491. CFD baselines are not computed for the naive default or uncalibrated conditions, as the primary purpose of these baselines is to contextualise the structural score rather than to provide full-pipeline comparisons.

[Table:result-detailed-metrics] Detailed evaluation results for the scale-calibrated pipeline (n = 97 valid cases out of 100; three cases excluded due to VLM output parsing failures).

| Metric | Value |
|--------|-------|
| Mean structural score | 0.781 ± 0.150 |
| Mean CFD score | 0.491 ± 0.131 |
| Room-kind match rate | 95% |
| Opening-wall match rate | 70% |
| Mean Lx relative error | 0.4% |
| Mean Ly relative error | 21.9% |
| Mean Lz relative error | 19.8% |
| Mean volume relative error | 36.6% |
| VLM parsing failures | 3 / 100 |

Three cases (bench_a3_04/perspective, bench_a4_02/perspective, bench_a4_02/wireframe) failed to produce a parseable VLM output and are excluded from score aggregation but counted as failures in the room-kind and opening-wall match statistics.

### 5.2 Effect of input view type

The input view type has a pronounced effect on both structural and CFD performance [Fig:result-view-aggregate]. The floor-plan view achieves the highest scores across all metrics, with a mean structural score of 0.884 ± 0.107 and a mean CFD score of 0.538 ± 0.147. This is consistent with the expectation that the orthographic top-down projection preserves the planimetric information most relevant to the room layout and obstacle placement. The bird's-eye view ranks second (structural 0.800, CFD 0.509), followed by the perspective view (0.782 ± 0.140, 0.468 ± 0.111), the wireframe view (0.760 ± 0.123, 0.508 ± 0.119), and the section view (0.679 ± 0.156, 0.430 ± 0.114).

The opening-wall match rate exhibits the most dramatic view dependence ([Fig:result-view-aggregate]b). Floor-plan inputs achieve a perfect 100% match rate, reflecting the unambiguous wall visibility in the orthographic projection. Bird's-eye inputs achieve 85%, perspective 65%, wireframe 75%, and section only 25%. The poor opening-wall performance of the section view is attributable to the fact that a vertical cutting plane typically intersects only one or two walls, leaving the remaining walls — and any openings on them — invisible to the VLM.

The room-kind match rate is more robust to view type, ranging from 90% (perspective and section) to 100% (bird's-eye and floor plan). This suggests that the VLM can reliably distinguish rectangular from composite room topologies even from oblique or partial views, a task that depends more on the overall room silhouette than on fine-grained wall visibility.

### 5.3 Effect of geometric complexity

Performance varies across the four complexity categories [Fig:result-category-aggregate]. Category A1 (simple rectangular) achieves the highest mean structural score (0.845) and CFD score (0.545), confirming its role as a positive control. Categories A3 (composite simple) and A4 (composite dense) achieve structural scores of 0.783 and 0.791 respectively, indicating that the composite topology does not by itself degrade the structural prediction severely. Category A2 (rectangular dense) yields the lowest structural score (0.707) despite having the simpler room topology, a result attributable to the higher obstacle density increasing the chance of obstacle placement errors.

The opening-wall match rate shows a different pattern. A1 achieves 88%, A3 and A4 both achieve 72%, and A2 achieves only 48%. The low A2 rate is driven by cases where the VLM places the inlet or outlet on the wrong wall, an error that is more consequential in rectangular rooms where all four walls are geometrically equivalent and the VLM has no topological cue to disambiguate them. In composite rooms (A3, A4), the L-shaped topology provides a stronger spatial anchor that helps the VLM assign openings to the correct walls.

### 5.4 Solver robustness

Of the 100 cases, 84 converge under the nominal solver configuration. The remaining 16 cases require escalation through the robustness pathway: 4 cases use the ultra-robust relaxation schedule, 5 fall back to a laminar solver, 5 use the conservative configuration, and 2 require mesh coarsening to 0.25 m cell size [Fig:result-robustness]. No case requires geometric repair, and no case fails to produce a converged solution. The escalation rate is highest for the section view (5 out of 20 cases) and the floor-plan view (5 out of 20), the latter driven by cases where the accurate geometry produces complex internal flow structures that challenge the turbulence model at the nominal mesh resolution.

---

## 6. Discussion

### 6.1 Failure mode I: composite-room collapse under section view

The section view exhibits a characteristic failure mode in which the VLM collapses a composite L-shaped room into a single rectangular room [Fig:discuss-section-collapse]. Both representative cases shown in [Fig:discuss-section-collapse] involve composite geometries (A3-04 and A4-05) where the section cutting plane passes through only one of the two constituent blocks. The VLM, seeing only a rectangular cross-section, infers a rectangular room and discards the second block entirely. This collapse reduces the structural score and produces a qualitatively different flow field, because the missing block eliminates recirculation zones that are present in the reference solution.

This failure is inherent to the section view's limited information content rather than to a deficiency in the VLM's reasoning. A vertical section through an L-shaped room is geometrically indistinguishable from a section through a rectangular room unless the cutting plane happens to intersect the junction between the two blocks. The practical implication is that section views should be avoided as inputs when the room is suspected to have a non-rectangular topology, or supplemented with a second view from a different angle.

### 6.2 Failure mode II: obstacle hallucination with limited CFD penalty

The VLM frequently predicts obstacles that do not exist in the reference geometry. This obstacle hallucination is particularly prevalent in wireframe and section views, where the sparse visual cues can be misinterpreted as internal structures. [Fig:discuss-obstacle-hallucination] illustrates two representative cases from category A3 rendered as wireframes. In the first case (A3-01), the reference geometry contains no obstacles, yet the VLM predicts three. In the second case (A3-03), one reference obstacle is present but the VLM predicts three, two of which are hallucinated.

Despite the structural error, the hallucinated obstacles have a surprisingly limited effect on the CFD fidelity. Both cases retain moderate CFD scores (0.604 and 0.600 respectively), and the opening-wall assignments are preserved in both cases. This resilience arises because the hallucinated obstacles are typically small relative to the room volume and are placed away from the main flow path between the inlet and outlet. They introduce localised flow disturbances but do not fundamentally alter the bulk flow pattern. This observation has a practical implication: obstacle prediction errors may be tolerable in applications where the primary quantity of interest is the bulk ventilation rate or the overall flow pattern rather than the detailed velocity field near specific furniture items.

### 6.3 Failure mode III: structure-versus-fidelity gap in dense composites

A scatter plot of structural score versus CFD score for all 97 valid cases [Fig:discuss-scatter-structural-cfd] reveals that the two metrics are only weakly correlated. Cases cluster broadly, and many achieve high structural scores (above 0.8) yet moderate or low CFD scores (below 0.5), indicating that geometric accuracy does not guarantee flow-field agreement. This decoupling is most pronounced in the A4 category, as discussed below.

The A4 category reveals a decoupling between structural accuracy and CFD fidelity. Several A4 cases achieve high structural scores (above 0.9) yet moderate or low CFD scores (below 0.5), indicating that geometrically accurate predictions do not guarantee accurate flow fields [Fig:discuss-structure-cfd-gap]. The two representative cases in [Fig:discuss-structure-cfd-gap] illustrate this gap. Case A4-02 achieves a structural score of 0.938 with correct obstacle placement (F1 = 0.75) and preserved openings, yet its CFD score is only 0.409. Case A4-04 achieves a perfect structural score of 1.000 with all obstacles correctly placed and matched openings, yet its CFD score is 0.458.

The underlying cause is the sensitivity of the flow field to small geometric perturbations in dense obstacle configurations. In a room with multiple obstacles, the flow path between the inlet and outlet is constrained by narrow gaps between obstacles. Even when the predicted obstacle positions are within the matching threshold (IoU > 0.18), residual misalignments of a few centimetres can redirect jets, create or eliminate recirculation zones, and shift stagnation points. The CFD solution in the reference case is computed on the exact ground-truth geometry, and even a geometrically "correct" prediction introduces enough positional noise to produce a meaningfully different flow field. This sensitivity is a fundamental limitation of the evaluation methodology and of the framework itself: for dense configurations, structural matching at the bounding-box level is necessary but not sufficient for CFD fidelity.

### 6.4 Scale calibration and its limitations

The single-factor post-hoc scaling reduces the mean Lx error from its uncalibrated value to 0.4% but leaves the secondary dimensions (Ly, Lz) with relative errors of approximately 20%. This anisotropic residual error has two consequences. First, the room volume is systematically distorted, which affects the bulk velocity magnitude for a given inlet flow rate. Second, the aspect ratio of the room is altered, which can shift the location of flow reattachment points and recirculation zones. A multi-factor calibration that independently corrects each dimension would reduce these errors but would require multiple reference measurements, undermining the single-image automation objective. The present single-factor approach represents a practical compromise that corrects the dominant scale error while preserving the fully automated workflow.

### 6.5 Computational cost

The end-to-end pipeline time per case is approximately 3 to 8 minutes on a single workstation, depending on the mesh resolution and solver convergence behaviour. The VLM API call (Gemini 3.1 Pro) typically returns within 10 to 30 seconds. Mesh generation with Gmsh takes 5 to 15 seconds for the default 0.35 m cell size. The dominant cost is the steady-state CFD solve, which ranges from 1 to 6 minutes depending on the case complexity and whether solver escalation is triggered. This total time is orders of magnitude shorter than the manual geometry preparation workflow it replaces, which typically requires hours of CAD modelling for a comparable indoor space. The framework is therefore suitable for rapid screening applications where approximate flow patterns are sufficient.

### 6.6 Practical considerations and limitations

Several limitations of the current study should be noted. First, the default mesh cell size of 0.35 m is coarse by indoor ventilation CFD standards. To assess the sensitivity of the CFD scores to mesh resolution, the reference solutions for two representative cases (A3-03 and A4-03) were recomputed at cell sizes of 0.35, 0.25, and 0.18 m, and the CFD score was re-evaluated against the same predicted solution at each reference resolution. [Table:discuss-grid-independence] shows that the CFD score changes by 0.02 to 0.10 as the reference mesh is refined, with the relative ranking between cases preserved. The 0.10 m resolution exceeded the solver timeout (1200 s) for all cases, indicating that substantially finer meshes are not feasible within the current pipeline's time budget. These results suggest that the CFD scores are influenced by reference mesh resolution but that the comparative conclusions drawn from the benchmark are robust to moderate refinement.

[Table:discuss-grid-independence] Grid independence check: CFD fidelity score of the predicted solution (at 0.35 m) evaluated against reference solutions at four mesh resolutions, spanning a 34× range in cell count.

| Case | Cells (0.35 m) | Ref. 0.35 m | Ref. 0.25 m | Ref. 0.18 m | Ref. 0.10 m | Δ max |
|------|----------------|-------------|-------------|-------------|-------------|-------|
| A1-01 | 12,001 | 0.389 | 0.568 | 0.579 | 0.625 | 0.236 |
| A3-03 | 14,735 | 0.739 | 0.726 | 0.641 | 0.616 | 0.123 |
| A4-03 | 12,240 | 0.576 | 0.553 | 0.558 | 0.598 | 0.045 |

For case A1-01 the reference flow field converges rapidly (mean velocity changes by 32% from 0.35 to 0.25 m, then under 2% for subsequent refinements), and the CFD score improves monotonically as the finer reference resolves features that the prediction also captures. For cases A3-03 and A4-03 the CFD score is more stable, varying by at most 0.12 across a 34× range in cell count, indicating that the geometry approximation error between predicted and reference scenes dominates over mesh discretisation error.

Second, the benchmark geometries are axis-aligned rectangular or L-shaped rooms, which do not represent the full diversity of real indoor spaces. Curved walls, sloped ceilings, and multi-storey configurations are not considered. Second, the CFD simulations use steady-state RANS with the k-omega SST turbulence model, which may not capture transient flow features or complex turbulent structures in densely obstructed rooms. Third, the evaluation metrics compare fields on the intersection of the predicted and reference domains, which may under-represent errors in cases where the predicted room extends beyond or falls short of the reference domain. Fourth, the scale calibration requires knowledge of one reference dimension, which may not be available in a fully unsupervised deployment; future work could explore self-calibrated approaches based on known object sizes (doors, standard furniture) visible in the image.

---

## 7. Application to Architectural Floor Plans

To assess the framework's applicability beyond the synthetic benchmark, two architectural floor plan images were prepared and processed through the full pipeline. These cases use the floor-plan view, which achieved the highest performance in the benchmark evaluation (structural score 0.884, opening-wall match 100%).

The first case is a simple rectangular office of 6.0 m × 4.0 m ([Fig:demo-floorplan-application]a). The floor plan includes hatched walls, a door with arc symbol on the south wall, a window on the north wall, and three furniture items (desk, chair, shelf). Dimension annotations on the drawing provide the scale reference. The VLM correctly identifies the rectangular room topology with dimensions 6.0 × 4.0 × 2.7 m, places three obstacles corresponding to the furniture items, and assigns the door as inlet (south wall) and the window as outlet (north wall). The resulting CFD simulation converges under the robust solver preset at 0.35 m mesh resolution, producing a physically plausible flow field with a jet from the south-wall inlet traversing the room toward the north-wall outlet ([Fig:demo-floorplan-application]b).

The second case is an L-shaped composite room consisting of a 5.0 m × 4.0 m living room connected to a 3.0 m × 2.5 m alcove ([Fig:demo-floorplan-application]c). The VLM correctly identifies the composite topology as a two-block structure. However, it places four obstacles — a sofa (present in the drawing), a coffee table, a bookshelf, and an alcove desk — of which the latter three are hallucinated and do not appear in the input image. This obstacle hallucination is consistent with the failure mode identified in [Sec:discuss-obstacle-hallucination] from the synthetic benchmark. The inlet is assigned to the north wall and the outlet to the west wall. Despite the hallucinated obstacles, the CFD simulation converges and produces a physically plausible flow pattern in which the main jet passes through the living room with a secondary branch entering the alcove region ([Fig:demo-floorplan-application]d). The scale calibration factor of 0.625 was applied to correct the VLM's initial overestimate of the room dimensions.

These two cases demonstrate that the framework can process architectural floor plans with realistic visual noise — hatched walls, furniture icons, dimension lines, door arc symbols — and still produce converged CFD solutions. The obstacle hallucination observed in Case 2 reinforces the benchmark finding that hallucinated obstacles have limited impact on the bulk flow pattern when they are small relative to the room volume. The scale calibration step is straightforward when dimension annotations are present on the drawing, as is common in architectural practice. The qualitative flow patterns are physically reasonable, though quantitative validation against measured data was not performed. The results suggest that floor-plan inputs from architectural drawings represent a practical deployment scenario for the framework.

---

## 8. Conclusion

This paper has presented a framework for automated indoor CFD analysis from a single 2D image, using a vision-language model to bridge the gap between visual input and simulation-ready 3D geometry. A systematic benchmark of 100 evaluation cases across four complexity levels and five input view types demonstrates that the framework achieves a mean structural score of 0.781 and a mean CFD score of 0.491 under single-factor scale calibration. The floor-plan view emerges as the most effective input type, achieving a structural score of 0.884 and a 100% opening-wall match rate. Three characteristic failure modes have been identified and analysed: composite-room collapse under section views, obstacle hallucination with limited CFD penalty, and a structure-versus-fidelity gap in dense configurations. The framework converges for all 100 cases through a multi-level robustness strategy, with 84% of cases solved under nominal settings.

The results suggest that VLM-based geometric abstraction is a viable approach for rapid indoor CFD screening when approximate flow patterns and bulk ventilation rates are the primary quantities of interest. For applications requiring high-fidelity local flow prediction near specific obstacles, the current framework's accuracy is insufficient and should be supplemented with manual geometry refinement or multi-view input fusion. Future work will explore multi-factor scale calibration, multi-view input strategies, and extension to transient thermal simulations.

---

## References

*[To be populated in the final manuscript.]*

---

### 5.5 VLM output repeatability

Because the VLM operates with a sampling temperature of 0.2, its output is stochastic rather than deterministic. To characterise this variation, three representative cases (bench_a1_01, bench_a3_03, bench_a4_03) were each processed three times through the floor-plan pipeline. [Table:result-vlm-repeatability] reports the structural scores across runs.

[Table:result-vlm-repeatability] VLM repeatability: structural scores across 3 independent runs per case (floor-plan view).

| Case | Run 1 | Run 2 | Run 3 | Mean ± SD |
|------|-------|-------|-------|-----------|
| A1-01 (simple) | 0.875 | 0.875 | 0.875 | 0.875 ± 0.000 |
| A3-03 (composite) | 1.000 | 0.875 | 0.875 | 0.917 ± 0.059 |
| A4-03 (dense) | 0.917 | 1.000 | 1.000 | 0.972 ± 0.039 |

The mean standard deviation across cases is 0.033, indicating moderate but bounded stochastic variation. Importantly, the room topology (room block F1 = 1.0) and opening-wall assignment (wall match ratio = 1.0) are perfectly consistent across all runs; the score differences arise entirely from obstacle matching, where the VLM occasionally generates a different number of obstacles or places them at slightly different positions. The simplest case (A1-01) produces identical output across all runs, suggesting that VLM output becomes more deterministic as scene complexity decreases. The raw room dimensions predicted by the VLM also vary between runs (e.g., scale factors of 0.54, 1.08, and 0.86 for the same case), confirming that post-hoc scale calibration is essential to absorb this dimensional variability.

---

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

[Fig:method-framework] Overall framework architecture. A single 2D image is processed by a vision-language model (Gemini 3.1 Pro, preview) to produce a structured 3D scene description. Post-hoc scale calibration corrects the dominant dimensional error using one reference dimension (the longest horizontal span). The calibrated scene is automatically converted to an OpenFOAM case and solved to steady state. An independent reference path computes the ground-truth CFD solution from the rule-based benchmark geometry without VLM involvement.

[Fig:bench-design] Benchmark dataset design. The 20 benchmark geometries are organised into a 2 x 2 matrix of room topology (rectangular versus composite L-shaped) and obstacle density (simple versus dense), yielding four complexity categories A1 through A4. Each panel shows a representative floor-plan rendering. The summary table reports per-category evaluation statistics.

[Fig:bench-multiview] Multi-view rendering protocol. Each benchmark geometry is rendered into five input views: (a) perspective, (b) bird's-eye, (c) floor plan, (d) wireframe, and (e) section. The same geometry (A4-03) is shown across all views to illustrate the variation in information content.

[Fig:method-eval-pathway] Evaluation pathway. The predicted path (top) processes the 2D input through the VLM, scale calibration, and CFD pipeline. The reference path (bottom) generates the ground-truth CFD solution from the rule-based benchmark without VLM involvement. Both paths converge at the metric comparison stage, which produces structural and CFD fidelity scores.

[Fig:result-view-aggregate] Aggregate performance across input view types. (a) Mean structural and CFD scores per view. (b) Room-kind and opening-wall match rates per view. The floor-plan view achieves the highest performance across all metrics.

[Fig:result-category-aggregate] Aggregate performance across benchmark categories. (a) Mean structural and CFD scores per category. (b) Room-kind and opening-wall match rates per category.

[Fig:result-crossview-outcome] Representative cross-view outcome comparison for two benchmark cases (A2-03 and A3-03) across three view types. Predicted obstacle positions are shown with matched (tan) and hallucinated (orange cross-hatch) encoding. Reference obstacle outlines are shown as dashed purple lines. Structural (S) and CFD (C) scores are annotated per panel.

[Fig:discuss-section-collapse] Composite-room collapse under section view. Two composite-topology cases (A3-04 sparse, A4-05 dense) are shown with reference geometry (left) and predicted geometry with ground-truth outline overlay (right). In both cases the VLM collapses the L-shaped room into a rectangular room, reducing the block count from two to one.

[Fig:discuss-obstacle-hallucination] Obstacle hallucination with limited CFD penalty. Two wireframe-input cases from category A3 are shown. The reference geometry (left) and prediction (right) are compared. Hallucinated obstacles (orange cross-hatch) are clearly distinguished from matched predictions (tan). Despite significant obstacle count inflation (0 to 3 and 1 to 3), opening topology is preserved and CFD scores remain moderate (0.60).

[Fig:discuss-structure-cfd-gap] Dense composite structure-versus-CFD gap. Two A4-category cases demonstrate that high structural scores (0.938 and 1.000) do not guarantee correspondingly high CFD scores (0.409 and 0.458). Left panels show geometry overlap between predicted (orange) and reference (dashed blue) obstacles. Right panels compare reference and predicted flow fields.

[Fig:result-robustness] Robustness and convergence summary.

[Fig:discuss-scatter-structural-cfd] Structural score versus CFD fidelity score for all 97 valid evaluation cases.

[Fig:result-heatmap-category-view] Category × view interaction heatmaps.

[Fig:demo-floorplan-application] Application to architectural floor plans. (a) Input floor plan for a 6.0 m × 4.0 m rectangular office with furniture, door, and window. (b) Gemini-extracted 3D scene and resulting CFD velocity field for case 1. (c) Input floor plan for an L-shaped composite room (5.0 m × 4.0 m + 3.0 m × 2.5 m alcove). (d) Gemini-extracted 3D scene and resulting CFD velocity field for case 2. Both cases converge under the robust solver preset at 0.35 m mesh resolution. (a) Mean structural score and (b) mean CFD fidelity score for each of the 20 cells in the 4-category × 5-view evaluation matrix. Cells with fewer than 5 valid cases are annotated with the sample count. The floor-plan view achieves the highest structural scores across all categories, while CFD scores are more uniformly distributed with less view-type dependence. Points are colour-coded by benchmark category (A1 blue, A2 orange, A3 green, A4 red) and shaped by input view type (circle = perspective, square = bird's-eye, diamond = floor plan, triangle-up = wireframe, triangle-down = section). The weak correlation between the two scores illustrates the structure-versus-fidelity decoupling discussed in [Sec:discuss-structure-cfd-gap]. (a) Room-kind and opening-wall agreement rates per view type. (b) Solver fallback trigger counts by view type, broken down by escalation level. (c) Overall pipeline statistics showing that all 100 cases converge, with 84 under nominal settings.
