# Post-hoc Scaling Smoke Comparison Note

> 작성일: 2026-03-20


## 설정
- baseline root: `benchmark/evaluations`
- post-hoc root: `benchmark/evaluations_posthoc_scaled_longest_span`
- subset: `bench_a1_01`, `bench_a2_03`, `bench_a3_03`, `bench_a4_03` × 5 views = 20 tasks
- scaling: reference GT longest horizontal span / baseline predicted longest horizontal span
- intervention point: baseline `predicted_scene.json` 이후, Gemini 재생성 없이 uniform scaling 후 downstream(mesh/solver/eval) 재실행

## 핵심 결과
- mean CFD score: **0.5181 → 0.5329** (`+0.0148`)
- mean structural score: **0.7613 → 0.8307** (`+0.0693`)
- opening-wall match rate: **0.75 → 0.75** (변화 없음)
- room-kind match rate: **1.00 → 1.00** (변화 없음)
- mean Lx relative error: **0.2737 → 0.0000**
- mean Ly relative error: **0.1492 → 0.2209**
- mean Lz relative error: **0.0714 → 0.2151**

## 해석
- uniform post-hoc scaling은 의도대로 **dominant horizontal span(Lx) 보정에는 매우 직접적**으로 작동했다.
- topology / opening-wall identity / room-kind은 subset에서 **보존**되었다.
- 다만 uniform scaling이라서 **Ly/Lz까지 함께 움직여**, baseline이 우연히 맞추고 있던 축에서는 오히려 오차가 커지는 경우가 있었다.
- 그럼에도 smoke subset 평균 CFD는 **소폭 개선(+0.0148)** 되었고, 전반적으로 “심하게 망가지는 레이어”는 아니었다.
- 개선/악화 task 수는 거의 비슷했지만(10 improved / 9 worsened / 1 unchanged), 큰 폭 개선 사례가 일부 있어 평균은 상승했다.

## 주의해서 볼 task
- 개선 상위: `bench_a4_03/wireframe` `+0.0862`, `bench_a1_01/birdseye` `+0.0808`, `bench_a1_01/wireframe` `+0.0612`, `bench_a4_03/birdseye` `+0.0538`
- 악화 상위: `bench_a3_03/birdseye` `-0.0636`, `bench_a2_03/wireframe` `-0.0191`, `bench_a2_03/floorplan` `-0.0124`, `bench_a4_03/perspective` `-0.0096`

## 1차 판단
- **prompt-side scale hint를 다시 만지지 않고도**, framework 내부 minimal post-hoc layer만으로 baseline 재사용 비교가 가능해졌다.
- smoke subset 기준으로는 **평균 개선 + topology 보존 + severe regression 부재**라서, 완전히 버릴 접근은 아니다.
- 다만 uniform single-factor scaling은 `Ly/Lz`까지 같이 흔들기 때문에, full benchmark 확대 전에 span anchor 적용 범위/축 제한 여부를 논의할 가치가 있다.
