# Post-hoc Scaling Full Benchmark Comparison

> Date: 2026-03-20

## Headline

- Baseline root: `/home/yhjoo/projects/OpenFOAM/benchmark/evaluations`
- Post-hoc root: `/home/yhjoo/projects/OpenFOAM/benchmark/evaluations_posthoc_scaled_longest_span`
- Post-hoc completion: **97/100 success**, **3 failed**
- Mean CFD score: **0.4934 -> 0.4909** (`-0.0025`)
- Mean structural score: **0.7067 -> 0.7813** (`+0.0746`)
- Mean room bbox relative error: `Lx 0.2487 -> 0.0039` (`-0.2447`), `Ly 0.1904 -> 0.2194` (`+0.0290`), `Lz 0.0898 -> 0.1983` (`+0.1084`)
- Opening-wall match rate: **0.72 -> 0.7216** (`+0.0016`)
- Room-kind match rate: **0.98 -> 0.9794** (`-0.0006`)
- Task-level CFD deltas among completed tasks: improved **53**, worsened **43**, unchanged **1**, severe worsened (<= -0.05) **9**

## View-level comparison

- **perspective**: CFD `0.4554 -> 0.4677` (`+0.0123`), structural `0.6958 -> 0.7824` (`+0.0866`), Lx err `0.275 -> 0.0043`
- **birdseye**: CFD `0.5225 -> 0.509` (`-0.0135`), structural `0.711 -> 0.8` (`+0.0890`), Lx err `0.3042 -> 0.0038`
- **floorplan**: CFD `0.5765 -> 0.5379` (`-0.0386`), structural `0.7623 -> 0.884` (`+0.1217`), Lx err `0.2488 -> 0.0038`
- **wireframe**: CFD `0.4854 -> 0.5082` (`+0.0227`), structural `0.7009 -> 0.7603` (`+0.0594`), Lx err `0.2816 -> 0.004`
- **section**: CFD `0.4273 -> 0.4304` (`+0.0032`), structural `0.6633 -> 0.6788` (`+0.0155`), Lx err `0.1339 -> 0.0038`

## Category-level comparison

- **A1**: success `25/25`, CFD `0.5567 -> 0.5452` (`-0.0115`), structural `0.785 -> 0.845` (`+0.0600`), opening-wall `0.88 -> 0.88`
- **A2**: success `25/25`, CFD `0.4438 -> 0.426` (`-0.0178`), structural `0.6727 -> 0.707` (`+0.0343`), opening-wall `0.48 -> 0.48`
- **A3**: success `24/25`, CFD `0.5466 -> 0.5336` (`-0.0130`), structural `0.6567 -> 0.783` (`+0.1263`), opening-wall `0.72 -> 0.75`
- **A4**: success `23/25`, CFD `0.4266 -> 0.458` (`+0.0314`), structural `0.7124 -> 0.791` (`+0.0786`), opening-wall `0.8 -> 0.7826`

## Failed post-hoc tasks

- `bench_a3_04/perspective` (A3), baseline CFD `0.4156`
- `bench_a4_02/perspective` (A4), baseline CFD `0.2875`
- `bench_a4_02/wireframe` (A4), baseline CFD `0.2963`

## Best CFD gains

- `bench_a4_01/floorplan`: `0.2935 -> 0.4802` (`+0.1867`)
- `bench_a1_02/floorplan`: `0.686 -> 0.8091` (`+0.1230`)
- `bench_a4_03/wireframe`: `0.4858 -> 0.572` (`+0.0862`)
- `bench_a1_01/birdseye`: `0.6061 -> 0.6869` (`+0.0808`)
- `bench_a1_05/birdseye`: `0.4819 -> 0.5555` (`+0.0736`)

## Worst CFD losses

- `bench_a1_05/floorplan`: `0.7405 -> 0.3252` (`-0.4152`)
- `bench_a3_01/floorplan`: `0.6965 -> 0.345` (`-0.3516`)
- `bench_a1_02/birdseye`: `0.6873 -> 0.3823` (`-0.3050`)
- `bench_a2_01/floorplan`: `0.5937 -> 0.2959` (`-0.2977`)
- `bench_a3_02/wireframe`: `0.6606 -> 0.4308` (`-0.2298`)
