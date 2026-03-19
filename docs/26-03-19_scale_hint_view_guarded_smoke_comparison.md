# Scale-Hint Variant Smoke Comparison (view-guarded wording refinement)

> 작성일: 2026-03-19  
> 범위: 4-case smoke subset (`bench_a1_01`, `bench_a2_03`, `bench_a3_03`, `bench_a4_03`) × 5 views = 20 tasks  
> 비교 대상:  
> - baseline: `benchmark/evaluations` (`no_scale_hint_baseline`)  
> - layout-protected: `benchmark/evaluations_scale_hinted_layout_protected` (`scale_hinted_longest_span_layout_protected_v1`)  
> - view-guarded: `benchmark/evaluations_scale_hinted_view_guarded` (`scale_hinted_longest_span_view_guarded_v1`)

## 1. 목적

`layout-protected span-only`를 유지한 채, stronger hint를 추가하지 않고도 다음 부작용 축을 줄이는 것이 이번 refinement의 목적이었다.

- `perspective`: hidden depth 과잉 추론 억제  
- `section`: unseen height / opening-wall relocation 과잉 보정 억제  
- `A4 dense composite`: opening/topology 우선, speculative obstacle detail 억제

## 2. completion status

- baseline smoke counterpart: **20/20 success**
- layout-protected smoke: **20/20 success**
- view-guarded smoke: **20/20 success**

참고: 첫 `bench_a4_03/wireframe` 시도는 dense-scene obstacle overlap으로 meshing 실패했지만, dense obstacle clearance guard를 추가한 뒤 단일 재시도로 성공했다. 따라서 최종 비교는 **completion-recovered 20/20** 기준이다.

## 3. overall mean delta

### 3.1 view-guarded - baseline
- structural score: **+0.0320**
- CFD score: **-0.0254**
- room bbox relative error Lx: **-0.2733**
- room bbox relative error Ly: **+0.0817**
- room bbox relative error Lz: **+0.0592**
- room volume relative error: **+0.0062**
- opening-wall match rate: **-0.0500**

### 3.2 view-guarded - layout-protected
- structural score: **+0.0079**
- CFD score: **-0.0091**
- room bbox relative error Lx: **+0.0000**
- room bbox relative error Ly: **+0.0363**
- room bbox relative error Lz: **+0.0001**
- room volume relative error: **+0.0245**
- opening-wall match rate: **+0.0000**

## 4. target-axis readout

### perspective
- vs baseline: structural **+0.0417**, CFD **-0.0212**, opening-wall **+0.0000**, Lx/Ly/Lz **-0.2082 / +0.1254 / -0.0016**
- vs layout-protected: structural **+0.0521**, CFD **+0.0639**, opening-wall **+0.0000**, Lx/Ly/Lz **+0.0000 / -0.0065 / -0.0183**

### section
- vs baseline: structural **-0.1479**, CFD **-0.0179**, opening-wall **-0.2500**, Lx/Ly/Lz **-0.0628 / +0.0624 / +0.0603**
- vs layout-protected: structural **-0.0021**, CFD **+0.0334**, opening-wall **+0.0000**, Lx/Ly/Lz **+0.0000 / +0.1038 / -0.1741**

### A4_all
- vs baseline: structural **+0.1448**, CFD **-0.0184**, opening-wall **+0.0000**, Lx/Ly/Lz **-0.3837 / +0.0259 / +0.0322**
- vs layout-protected: structural **+0.0400**, CFD **+0.0346**, opening-wall **+0.0000**, Lx/Ly/Lz **+0.0000 / -0.0330 / +0.0000**

### A4_perspective
- vs baseline: structural **+0.2084**, CFD **+0.0340**, opening-wall **+0.0000**, Lx/Ly/Lz **-0.3700 / +0.2215 / +0.0000**
- vs layout-protected: structural **+0.0000**, CFD **+0.1085**, opening-wall **+0.0000**, Lx/Ly/Lz **+0.0000 / -0.1108 / -0.1611**

### A4_section
- vs baseline: structural **+0.0333**, CFD **-0.0072**, opening-wall **+0.0000**, Lx/Ly/Lz **-0.0959 / +0.0062 / +0.0000**
- vs layout-protected: structural **+0.1167**, CFD **+0.1140**, opening-wall **+0.0000**, Lx/Ly/Lz **+0.0000 / +0.0062 / +0.0000**

### A4_wireframe
- vs baseline: structural **-0.0596**, CFD **+0.0359**, opening-wall **+0.0000**, Lx/Ly/Lz **-0.3014 / +0.2215 / +0.0967**
- vs layout-protected: structural **+0.0833**, CFD **+0.1388**, opening-wall **+0.0000**, Lx/Ly/Lz **+0.0000 / +0.0000 / +0.0000**

## 5. 해석

### 긍정 신호

1. **perspective / section / A4 타깃 축에서는 개선 신호가 분명하다.**  
   - perspective mean CFD는 baseline 대비 **-0.0212**로, layout-protected 대비 **+0.0639** 개선됐다.  
   - section mean CFD는 baseline 대비 **-0.0179**로, layout-protected 대비 **+0.0334** 개선됐다.  
   - A4 전체 mean CFD는 baseline 대비 **-0.0184**로, layout-protected 대비 **+0.0346** 개선됐다.

2. **특히 dense composite의 핵심 pain-point였던 A4 perspective / section / wireframe이 layout-protected 대비 모두 좋아졌다.**  
   - `A4/perspective`: CFD **+0.1085** vs layout-protected  
   - `A4/section`: CFD **+0.1140** vs layout-protected  
   - `A4/wireframe`: CFD **+0.1388** vs layout-protected

3. **Lx anchor 효과는 그대로 유지된다.**  
   - overall `Lx` relative error mean은 baseline 대비 **-0.2733**으로 layout-protected와 동일 수준이다.

### 남아 있는 문제

1. **overall mean CFD는 아직 baseline을 넘지 못한다.**  
   - view-guarded overall CFD delta는 baseline 대비 **-0.0254**다.

2. **전체 평균만 보면 layout-protected보다 오히려 약간 후퇴했다.**  
   - overall CFD는 layout-protected 대비 **-0.0091** 낮다.

3. **opening-wall match rate는 baseline보다 여전히 낮다.**  
   - overall opening-wall match rate delta는 baseline 대비 **-0.0500**이며, layout-protected 대비 추가 회복은 없었다.

4. **section의 구조 점수는 여전히 baseline보다 낮다.**  
   - section mean structural delta는 baseline 대비 **-0.1479**로, 과잉 regularization은 줄였지만 구조 회복이 완전히 해결된 것은 아니다.

## 6. 운영 결론

1. **view-guarded wording은 targeted side-effect mitigation에는 성공적**이다.  
   perspective / section / A4 dense composite 축에서는 layout-protected보다 해석상 더 맞는 방향으로 움직였다.

2. 하지만 **main-setting 승격 기준(= baseline 대비 CFD degradation 해소)** 은 아직 통과하지 못했다.

3. 따라서 현재 결론은 다음과 같다.

> **결정: frozen-20 full rerun은 계속 보류.**  
> view-guarded는 targeted follow-up 방향성은 확인했지만, overall smoke 기준으로는 baseline CFD degradation을 제거하지 못했다.

4. 다음 단계는 full rerun이 아니라, 이번 결과를 바탕으로 **targeted wording/selection refinement 또는 axis-aware acceptance rule** 을 더 좁히는 것이다. 예를 들어, overall main setting으로 바로 승격하기보다 `perspective/section/A4` 개선 신호를 별도 evidence lane으로 다루는 편이 더 타당하다.
