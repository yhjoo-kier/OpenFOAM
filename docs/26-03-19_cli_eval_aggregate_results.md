# Frozen-20 CLI Evaluation Aggregate Results

> Date: 2026-03-19

## Headline

- Frozen-20 image-conditioned benchmark evaluation reached **100/100 task success**.
- Aggregate mean structural score: **0.7067**
- Aggregate mean CFD score: **0.4934**
- Room-kind match rate: **98.00%**
- Opening-wall match rate: **72.00%**
- Repair was needed in only **2** tasks; mesh-size `0.25` fallback appeared in **2** tasks.

## View-level results

- **perspective**: structural `0.6958`, CFD `0.4554`, room-kind match `100.00%`, opening-wall match `70.00%`
- **birdseye**: structural `0.7110`, CFD `0.5225`, room-kind match `100.00%`, opening-wall match `85.00%`
- **floorplan**: structural `0.7623`, CFD `0.5765`, room-kind match `100.00%`, opening-wall match `100.00%`
- **wireframe**: structural `0.7009`, CFD `0.4854`, room-kind match `100.00%`, opening-wall match `80.00%`
- **section**: structural `0.6633`, CFD `0.4273`, room-kind match `90.00%`, opening-wall match `25.00%`

Interpretation:

- `floorplan` was the strongest aggregate view on both structure and CFD, suggesting that plan-layout/opening placement dominates many benchmark cases more than photorealistic texture cues.
- `section` was the weakest aggregate view and the only view with room-kind collapses, indicating that single-cut geometry evidence is often insufficient for composite-room recovery.
- `perspective` remained useful but less stable than expected; it needed both of the `mesh_size=0.25` fallbacks and one of the two repaired-scene salvages.

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

## Hard failure modes that still matter despite 100/100 task success

- Room-kind collapse: `bench_a3_04/section` (`composite -> rectangular`), structural `0.4167`, CFD `0.3352`
- Room-kind collapse: `bench_a4_05/section` (`composite -> rectangular`), structural `0.4167`, CFD `0.3311`
- Repair salvage: `bench_a2_04/perspective` with preset `robust` at mesh `0.35`; structural `0.5000`, CFD `0.5222`
- Repair salvage: `bench_a4_02/wireframe` with preset `ultra_robust` at mesh `0.35`; structural `0.7500`, CFD `0.2963`

Additional robustness counters:

- `ultra_robust` tasks: **5**
- `conservative` tasks: **1**
- repaired-scene successes: **2**

## Lowest- and highest-performing cases by mean CFD

### Lowest mean CFD

- `bench_a4_02` (A4, composite): structural `0.7375`, CFD `0.3340`, opening-wall matches `5/5`
- `bench_a2_02` (A2, rectangular): structural `0.6500`, CFD `0.3568`, opening-wall matches `1/5`
- `bench_a4_05` (A4, composite): structural `0.6208`, CFD `0.3947`, opening-wall matches `3/5`
- `bench_a4_01` (A4, composite): structural `0.5917`, CFD `0.4124`, opening-wall matches `3/5`
- `bench_a2_03` (A2, rectangular): structural `0.6500`, CFD `0.4199`, opening-wall matches `1/5`

### Highest mean CFD

- `bench_a1_03` (A1, rectangular): structural `0.7500`, CFD `0.6961`, opening-wall matches `5/5`
- `bench_a1_02` (A1, rectangular): structural `0.8250`, CFD `0.5936`, opening-wall matches `4/5`
- `bench_a3_03` (A3, composite): structural `0.7750`, CFD `0.5893`, opening-wall matches `4/5`
- `bench_a3_01` (A3, composite): structural `0.6750`, CFD `0.5882`, opening-wall matches `4/5`
- `bench_a3_02` (A3, composite): structural `0.7750`, CFD `0.5601`, opening-wall matches `4/5`

## Results/Discussion-ready takeaways

1. **Success-rate reporting alone is insufficient.** The benchmark now has 100/100 task success, yet aggregate CFD fidelity remains moderate (`~0.49`) and varies systematically by view/category.
2. **Opening/topology fidelity matters more than exact obstacle count in several composite cases.** This is especially visible in `A3`, where obstacle hallucinations often coexist with reasonable CFD alignment.
3. **Rectangular multi-obstacle cases need blockage-sensitive interpretation.** `A2` underperforms largely through opening-wall mistakes and obstacle-layout distortion rather than room-kind collapse.
4. **Dense composite cases expose the biggest structure-vs-physics gap.** `A4` often preserves room-kind but still loses CFD fidelity, so future discussion should not overclaim from structural scores alone.
5. **Section view should be framed as a stress input, not a generally strong modality.** Its room-kind collapses and weak opening-wall match rate make it a useful failure-analysis axis.

## Suggested follow-up tags / metrics

- Composite cases: add opening/topology-sensitive commentary tags when obstacle hallucination does not strongly degrade CFD.
- Rectangular multi-obstacle cases: add occupancy/blockage-sensitive tags so `A2` penalties are not hidden behind room-kind correctness.
- Empty or light composite controls: add hallucinated-obstacle burden tags separate from CFD penalties.
- Dense composite tail cases: preserve non-blocking repair-sidecar / stabilization warning tags as robustness metadata rather than benchmark failures.

## Artifact

- JSON summary: `benchmark/manifests/evaluation_aggregate_summary.json`
- Generated from: `benchmark/evaluations/*/*/evaluation_summary.json`
