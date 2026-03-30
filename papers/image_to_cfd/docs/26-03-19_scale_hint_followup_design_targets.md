# Scale-Hint Follow-up Design Targets

> 작성일: 2026-03-19  
> 범위: 4-case smoke subset (`bench_a1_01`, `bench_a2_03`, `bench_a3_03`, `bench_a4_03`) 기준의 baseline / single-anchor / dual-anchor / layout-protected 비교를 바탕으로, 다음 redesign의 우선 타깃을 정리한다.

## 1. 왜 이 문서가 필요한가

scale-hint pivot 자체는 이미 타당해졌다.

- benchmark dataset는 유지된다.
- baseline vs scale-hinted bookkeeping은 분리되었다.
- scale hint가 `Lx` 절대 scale ambiguity를 줄인다는 사실도 smoke에서 반복 확인되었다.

하지만 **full frozen-20 rerun을 아직 시작하지 않는 이유**는 단순히 mean CFD가 baseline보다 낮다는 한 줄 때문만은 아니다.
다음 redesign이 어디를 겨냥해야 하는지 더 구체적으로 남겨둘 필요가 있다.

이 문서는 그 목적의 짧은 설계 메모다.

## 2. 현재 best candidate의 성격

현재 3개 hinted variant 중 가장 균형 잡힌 후보는 여전히:

- `scale_hinted_longest_span_layout_protected_v1`

이다.

이 variant는 baseline 대비:

- `Lx` relative error를 크게 줄이고,
- single-anchor / dual-anchor보다 mean CFD degradation을 줄였고,
- opening-wall fidelity도 dual-anchor보다는 덜 해친다.

하지만 **모든 뷰/카테고리에서 좋아진 것은 아니다.**
오히려 후속 redesign 타깃은 꽤 명확하다.

## 3. view-level signal (layout-protected vs baseline)

### 좋아진 축
- `floorplan`
  - structural: **+0.1250**
  - CFD: **+0.1213**
  - `Lx`: **-0.2667**
  - `Ly`: **-0.0905**
  - `Lz`: **-0.0331**
- 해석:
  - layout-protected wording은 **평면 배치 보존 + span anchor** 조합에서 가장 잘 작동한다.
  - 즉, 이 variant의 장점은 "layout-dominant view에서 absolute span을 붙였을 때" 드러난다.

### 거의 유지되는 축
- `wireframe`
  - structural: **-0.0045**
  - CFD: **-0.0039**
  - `Lx`: **-0.3543**
- 해석:
  - wireframe은 scale regularization은 잘 받지만, CFD uplift로까지 이어지지는 않는다.
  - 그래도 큰 손해는 아니라서 next redesign에서 크게 건드릴 1순위는 아니다.

### 여전히 약한 축
- `perspective`
  - structural: **-0.0104**
  - CFD: **-0.0851**
  - `Ly`: **+0.1319**
  - `Lz`: **+0.0167**
- `section`
  - structural: **-0.1459**
  - CFD: **-0.0513**
  - opening-wall match: **-0.25**
  - `Lz`: **+0.2344**
- 해석:
  - perspective에서는 span anchor가 들어가도 **depth/side-wall 해석이 안정화되지 않는다.**
  - section에서는 scale hint가 오히려 **height / opening-wall fidelity를 건드릴 가능성**이 있다.
  - 즉, 다음 redesign은 floorplan을 더 좋게 만드는 작업이 아니라, **perspective·section에서의 부작용을 줄이는 작업**이어야 한다.

## 4. category-level signal (layout-protected vs baseline)

### 개선이 확인된 카테고리
- `A2`
  - structural: **+0.0667**
  - CFD: **+0.0348**
  - `Lx`: **-0.3761**
- `A3`
  - structural: **+0.0250**
  - CFD: **+0.0252**
  - `Lx`: **-0.1938**
- 해석:
  - layout-protected wording은 **rectangular multi-obstacle(A2)** 와 **light composite(A3)** 에서는 실제로 의미 있는 개선을 낸다.
  - 즉, scale hint 자체가 잘못된 게 아니라, 특정 regime에서는 이미 유효하다.

### 여전히 약한 카테고리
- `A1`
  - structural: **-0.1000**
  - CFD: **-0.0725**
- `A4`
  - structural: **+0.1048**
  - CFD: **-0.0529**
- 해석:
  - A1은 쉬운 케이스라 scale hint 이득보다 wording side-effect가 더 크게 보인다.
  - A4는 구조 점수는 올라가도 CFD가 따라오지 않는 기존 패턴이 그대로 남는다.
  - 따라서 next redesign에서 가장 중요한 stress target은 **A4**, 그다음은 **A1 over-regularization 회피**다.

## 5. 다음 redesign이 겨냥해야 할 실패 모드

### Target 1. perspective depth regularization 실패
현재 layout-protected wording도 perspective CFD를 baseline보다 크게 깎는다.
이는 longest-span anchor가 들어가도:

- depth axis(`Ly`)가 흔들리고,
- obstacle thickness / opening placement가 함께 흔들리며,
- 결과적으로 flow path가 불안정해짐

을 시사한다.

**다음 문구 설계는 perspective에서 "front/back depth를 과도하게 재해석하지 말라"는 제약을 더 직접적으로 넣는 편이 낫다.**

예시 방향:
- visible wall-to-wall depth ordering을 유지하라
- avoid expanding hidden depth beyond what the image supports
- preserve opening wall identity before refining obstacle sizes

### Target 2. section view의 opening/height 손상
section은 원래도 stress view였지만, layout-protected variant에서도:

- opening-wall match가 회복되지 않았고,
- `Lz` error가 크게 늘었다.

즉, section에서는 span anchor가 **cut plane 밖의 3D regularization**으로 번지는 순간 오히려 손해일 수 있다.

**후속 variant에서는 section 전용으로 hint 적용을 약하게 하거나, section scenario 문구를 별도로 두는 것이 더 타당할 수 있다.**

예시 방향:
- if the image is a section/cut view, do not infer unseen ceiling height aggressively
- keep opening wall identity and cut-plane geometry before enforcing global scale

### Target 3. dense composite(A4) 구조-물리 decoupling
A4는 layout-protected에서도 structural gain은 보이지만 CFD는 baseline보다 낮다.
이는 scale hint 자체보다도:

- dense composite에서 obstacle/opening local topology 보존이 더 중요하고,
- coarse global scale regularization만으로는 CFD fidelity가 안 오른다

는 점을 다시 확인한다.

**즉, A4에 대해서는 "global span anchor"보다 "opening/topology preservation" 우선순위를 더 강하게 못 박는 wording이 필요하다.**

## 6. 권장 next experiment (small, not full rerun)

full rerun 전 다음 정도의 targeted smoke가 가장 합리적이다.

### 추천 variant name (working)
- `scale_hinted_longest_span_view_guarded_v1`

### 핵심 설계 원칙
1. 기본 anchor는 여전히 **longest horizontal span only** 유지
2. layout-protected wording 유지
3. 대신 아래 guard를 추가
   - perspective: hidden depth 과잉 추론 억제
   - section: unseen height / topology 과잉 보정 억제
   - 공통: opening wall identity를 obstacle detail보다 우선 보존

### 추천 smoke 범위
- 최소 3-case stress subset
  - `bench_a1_01` — easy rectangular over-regularization 감시
  - `bench_a4_03` — dense composite physics gap 감시
  - `bench_a3_03` 또는 `bench_a2_03` — current positive signal 유지 확인
- 또는 기존 4-case subset 재사용

## 7. 운영 결론

현재까지의 evidence를 한 줄로 요약하면 다음과 같다.

> **scale hint pivot은 맞지만, next redesign의 목적은 더 강한 scale regularization이 아니라 perspective/section/A4 부작용 억제다.**

즉,

- single-anchor → dual-anchor처럼 정보를 더 늘리는 방향은 우선순위가 낮고,
- layout-protected를 바탕으로 **view-guarded wording** 을 넣는 방향이 다음 한 수로 더 낫다.

이 결론이 유지되는 한, frozen-20 full rerun 보류 판단은 타당하다.
