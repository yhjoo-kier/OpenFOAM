# Frozen-20 CLI Evaluation Aggregate Results

> Date: 2026-03-19

## Headline

- Frozen-20 image-conditioned benchmark evaluation reached **100/100 task success**.
- Aggregate mean structural score: **0.7067**
- Aggregate mean CFD score: **0.4934**
- Room-kind match rate: **98.00%**
- Opening-wall match rate: **72.00%**
- Repair was needed in only **2** tasks; mesh-size `0.25` fallback appeared in **2** tasks.
- Non-blocking repair-sidecar warnings were preserved in **8** tasks.

## View-level results

- **perspective**: structural `0.6958`, CFD `0.4554`, room-kind match `100.00%`, opening-wall match `70.00%`
- **birdseye**: structural `0.7110`, CFD `0.5225`, room-kind match `100.00%`, opening-wall match `85.00%`
- **floorplan**: structural `0.7623`, CFD `0.5765`, room-kind match `100.00%`, opening-wall match `100.00%`
- **wireframe**: structural `0.7009`, CFD `0.4854`, room-kind match `100.00%`, opening-wall match `80.00%`
- **section**: structural `0.6633`, CFD `0.4273`, room-kind match `90.00%`, opening-wall match `25.00%`

Interpretation:

- `floorplan` was the strongest aggregate view on both structure and CFD, suggesting that plan-layout/opening placement dominates many benchmark cases more than photorealistic texture cues.
- `section` was the weakest aggregate view and the only view with room-kind collapses, indicating that single-cut geometry evidence is often insufficient for composite-room recovery.
- `perspective` remained useful but less stable than expected; it needed both of the `mesh_size=0.25` fallbacks and one repaired-scene salvage, and it also carries part of the repair-sidecar warning burden.

## Category-level results

- **A1**: structural `0.7850`, CFD `0.5567`, room-kind match `100.00%`, opening-wall match `88.00%`
- **A2**: structural `0.6727`, CFD `0.4438`, room-kind match `100.00%`, opening-wall match `48.00%`
- **A3**: structural `0.6567`, CFD `0.5466`, room-kind match `96.00%`, opening-wall match `72.00%`
- **A4**: structural `0.7124`, CFD `0.4266`, room-kind match `96.00%`, opening-wall match `80.00%`

Interpretation:

- `A1` is the cleanest regime overall and behaves as the expected positive-control category.
- `A2` is the weakest rectangular category mainly because opening-wall fidelity collapses; its low opening-wall match rate suggests blockage/layout metrics alone are not enough.
- `A3` shows the clearest structure-vs-CFD decoupling: structure scores are modest, yet CFD stays relatively strong because opening/topology fidelity can outweigh obstacle-count hallucinations in empty/light composite rooms.
- `A4` keeps room-kind fidelity mostly intact but still yields low CFD, meaning dense composite layouts remain physically hard even when coarse structure recovery looks acceptable.

## Room-kind split

- **rectangular**: structural `0.7288`, CFD `0.5003`, room-kind match `100.00%`, opening-wall match `68.00%`
- **composite**: structural `0.6845`, CFD `0.4866`, room-kind match `96.00%`, opening-wall match `76.00%`

## Derived interpretation tags

- **opening_topology_preserved_with_obstacle_hallucination**: 10 tasks
  - Composite task kept opening/topology signal strong enough for CFD despite extra hallucinated obstacles.
  - Examples: `bench_a3_01/birdseye` (struct `0.6250`, CFD `0.7315`), `bench_a3_01/floorplan` (struct `0.7500`, CFD `0.6965`), `bench_a3_01/section` (struct `0.6250`, CFD `0.5855`)
- **rectangular_blockage_layout_failure**: 17 tasks
  - Rectangular multi-obstacle task underperformed through blockage/layout distortion and/or opening-wall mismatch rather than room collapse.
  - Examples: `bench_a2_01/perspective` (struct `0.7500`, CFD `0.4543`), `bench_a2_01/section` (struct `0.7250`, CFD `0.4216`), `bench_a2_01/wireframe` (struct `0.7500`, CFD `0.3908`)
- **dense_composite_structure_physics_gap**: 10 tasks
  - Dense composite task looked structurally acceptable but still lost CFD fidelity.
  - Examples: `bench_a4_01/birdseye` (struct `0.7500`, CFD `0.4888`), `bench_a4_02/birdseye` (struct `0.7500`, CFD `0.3391`), `bench_a4_02/floorplan` (struct `0.8125`, CFD `0.4023`)
- **section_room_kind_collapse**: 2 tasks
  - Section-view task collapsed a composite room into a rectangular prediction.
  - Examples: `bench_a3_04/section` (struct `0.4167`, CFD `0.3352`), `bench_a4_05/section` (struct `0.4167`, CFD `0.3311`)
- **repair_salvaged_success**: 2 tasks
  - Task required repaired-scene salvage to finish successfully.
  - Examples: `bench_a2_04/perspective` (struct `0.5000`, CFD `0.5222`), `bench_a4_02/wireframe` (struct `0.7500`, CFD `0.2963`)
- **mesh025_fallback**: 2 tasks
  - Task needed the smaller 0.25 mesh-size fallback to succeed.
  - Examples: `bench_a2_01/perspective` (struct `0.7500`, CFD `0.4543`), `bench_a3_05/perspective` (struct `0.6250`, CFD `0.5042`)
- **ultra_robust_escalation**: 5 tasks
  - Task required ultra_robust solver escalation.
  - Examples: `bench_a1_05/birdseye` (struct `0.7500`, CFD `0.4819`), `bench_a3_01/perspective` (struct `0.6250`, CFD `0.3795`), `bench_a3_02/floorplan` (struct `0.7500`, CFD `0.4703`)
- **nonblocking_repair_sidecar_warning**: 8 tasks
  - Original scene succeeded, but a repair-sidecar attempt still left warning traces worth preserving as robustness metadata.
  - Examples: `bench_a2_04/wireframe` (struct `0.5000`, CFD `0.4580`), `bench_a2_05/floorplan` (struct `0.8333`, CFD `0.7839`), `bench_a3_03/wireframe` (struct `0.6250`, CFD `0.5680`)

Tag interpretation:

- The new tag layer makes it easier to cite benchmark behavior directly from the manifest instead of reconstructing the narrative from many per-case notes.
- `opening_topology_preserved_with_obstacle_hallucination` captures the recurring A3 pattern where extra obstacles are hallucinated but CFD remains decent because openings/topology stay usable.
- `rectangular_blockage_layout_failure` isolates the A2-style rectangular degradation mode where room-kind is fine but obstacle/opening placement still damages physics.
- `dense_composite_structure_physics_gap` flags tasks where structural recovery alone would overstate performance on dense composite rooms.
- `nonblocking_repair_sidecar_warning` keeps stabilization noise visible without misclassifying those tasks as benchmark failures.

## Hard failure modes that still matter despite 100/100 task success

- Room-kind collapse: `bench_a3_04/section` (`composite -> rectangular`), structural `0.4167`, CFD `0.3352`
- Room-kind collapse: `bench_a4_05/section` (`composite -> rectangular`), structural `0.4167`, CFD `0.3311`
- Repair salvage: `bench_a2_04/perspective` with preset `robust` at mesh `0.35`; structural `0.5000`, CFD `0.5222`
- Repair salvage: `bench_a4_02/wireframe` with preset `ultra_robust` at mesh `0.35`; structural `0.7500`, CFD `0.2963`
- Non-blocking repair-sidecar warning: `bench_a2_04/wireframe` with preset `robust`; structural `0.5000`, CFD `0.4580`
- Non-blocking repair-sidecar warning: `bench_a2_05/floorplan` with preset `robust`; structural `0.8333`, CFD `0.7839`
- Non-blocking repair-sidecar warning: `bench_a3_03/wireframe` with preset `robust`; structural `0.6250`, CFD `0.5680`
- Non-blocking repair-sidecar warning: `bench_a3_04/perspective` with preset `robust`; structural `0.3750`, CFD `0.4156`
- Non-blocking repair-sidecar warning: `bench_a3_04/wireframe` with preset `robust`; structural `0.3750`, CFD `0.4123`
- Non-blocking repair-sidecar warning: `bench_a4_01/perspective` with preset `robust`; structural `0.3750`, CFD `0.2904`
- Non-blocking repair-sidecar warning: `bench_a4_05/birdseye` with preset `robust`; structural `0.5625`, CFD `0.4280`
- Non-blocking repair-sidecar warning: `bench_a4_05/perspective` with preset `robust`; structural `0.7500`, CFD `0.3797`

Additional robustness counters:

- `ultra_robust` tasks: **5**
- `conservative` tasks: **1**
- repaired-scene successes: **2**
- non-blocking repair-sidecar warnings: **8**

## Lowest- and highest-performing cases by mean CFD

### Lowest mean CFD

- `bench_a4_02` (A4, composite): structural `0.7375`, CFD `0.3340`, opening-wall matches `5/5`, dominant tags: dense_composite_structure_physics_gap
- `bench_a2_02` (A2, rectangular): structural `0.6500`, CFD `0.3568`, opening-wall matches `1/5`, dominant tags: rectangular_blockage_layout_failure
- `bench_a4_05` (A4, composite): structural `0.6208`, CFD `0.3947`, opening-wall matches `3/5`, dominant tags: nonblocking_repair_sidecar_warning, dense_composite_structure_physics_gap
- `bench_a4_01` (A4, composite): structural `0.5917`, CFD `0.4124`, opening-wall matches `3/5`
- `bench_a2_03` (A2, rectangular): structural `0.6500`, CFD `0.4199`, opening-wall matches `1/5`, dominant tags: rectangular_blockage_layout_failure

### Highest mean CFD

- `bench_a1_03` (A1, rectangular): structural `0.7500`, CFD `0.6961`, opening-wall matches `5/5`
- `bench_a1_02` (A1, rectangular): structural `0.8250`, CFD `0.5936`, opening-wall matches `4/5`
- `bench_a3_03` (A3, composite): structural `0.7750`, CFD `0.5893`, opening-wall matches `4/5`, dominant tags: opening_topology_preserved_with_obstacle_hallucination
- `bench_a3_01` (A3, composite): structural `0.6750`, CFD `0.5882`, opening-wall matches `4/5`, dominant tags: opening_topology_preserved_with_obstacle_hallucination
- `bench_a3_02` (A3, composite): structural `0.7750`, CFD `0.5601`, opening-wall matches `4/5`, dominant tags: opening_topology_preserved_with_obstacle_hallucination

## Results/Discussion-ready takeaways

1. **Success-rate reporting alone is insufficient.** The benchmark now has 100/100 task success, yet aggregate CFD fidelity remains moderate (`~0.49`) and varies systematically by view/category.
2. **Opening/topology fidelity matters more than exact obstacle count in several composite cases.** The new interpretation tags now preserve this A3-style pattern directly in the manifest.
3. **Rectangular multi-obstacle cases need blockage-sensitive interpretation.** `A2` underperforms largely through opening-wall mistakes and obstacle-layout distortion rather than room-kind collapse.
4. **Dense composite cases expose the biggest structure-vs-physics gap.** `A4` often preserves room-kind but still loses CFD fidelity, so future discussion should not overclaim from structural scores alone.
5. **Section view should be framed as a stress input, not a generally strong modality.** Its room-kind collapses and weak opening-wall match rate make it a useful failure-analysis axis.
6. **Stabilization metadata should survive into the dataset card/manifests.** Repair-sidecar warnings and salvage/fallback traces are rare but meaningful reproducibility cues, not just log noise.

## Artifact

- JSON summary: `benchmark/manifests/evaluation_aggregate_summary.json`
- Generated from: `benchmark/evaluations/*/*/evaluation_summary.json`
