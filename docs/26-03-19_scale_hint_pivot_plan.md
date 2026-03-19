# Scale-Hint Pivot Plan for Image-to-CFD Evaluation

> 작성일: 2026-03-19
> 프로젝트: OpenFOAM Image-to-CFD paper
> 배경: image-conditioned scene generation 결과를 검토하던 중, 예측 형상의 절대 길이(scale)가 GT와 체계적으로 어긋날 수 있다는 핵심 방법론 이슈가 식별되었다. 벤치마크 데이터셋 자체는 rule-based absolute metric geometry를 포함하므로 유지 가능하지만, **Gemini가 2D benchmark image로부터 형상을 예측하는 단계에는 절대 길이 기준 힌트가 필요**하다는 판단을 반영해 평가 설정을 수정한다.

---

## 1. 핵심 결론

### 유지되는 것
- frozen benchmark dataset 자체 (`scene JSON`, `reference CFD`, `5-view render bundle`)는 그대로 유지한다.
- rule-based benchmark는 이미 절대 길이(m 단위)를 갖는 정답 데이터셋이다.
- 따라서 benchmark release/artifact/reference 경로를 폐기할 필요는 없다.

### 바뀌는 것
- **image-conditioned Gemini scene generation prompt / evaluation setting** 은 수정해야 한다.
- 기존 설정은 `units = m` 출력은 강제하지만, **이미지 해석 단계에서 absolute scale anchor를 주지 않았다.**
- 앞으로는 **길이 힌트(scale hint)** 를 함께 제공하는 설정으로 재정의한다.

### 직접적 영향 범위
- 다시 해야 하는 것은 주로:
  1. Gemini prompt / generation interface
  2. benchmark evaluation task runner
  3. image-conditioned evaluation outputs / aggregate summaries
  4. 해당 결과에 의존하는 figure / manuscript interpretation
- 다시 할 필요가 없는 것은 주로:
  1. benchmark scene generation
  2. reference CFD generation
  3. benchmark rendering bundle 자체

---

## 2. 왜 이 pivot이 필요한가

현재 프레임워크는 다음을 만족한다.
- scene JSON은 m 단위 absolute geometry를 표현한다.
- validator는 `units = "m"`, room/opening/obstacle bounds를 검증한다.
- benchmark evaluation은 사후적으로 room bbox absolute/relative error를 계산한다.

하지만 현재 Gemini generation prompt는 아래를 제공하지 않는다.
- longest side가 몇 m인지
- ceiling height가 몇 m인지
- overall room extent가 어느 정도인지
- absolute scale anchor가 무엇인지

즉 현재 설정은:
- **relative structure inference는 요구하지만**
- **absolute scale inference는 모델의 암묵적 상식에 맡기는 상태**다.

이는 아래 이유로 실사용성과 논문 타당성 모두에서 약점이 된다.
- 단일 이미지로 절대 길이를 안정적으로 유추하는 것은 본질적으로 어렵다.
- CFD에서는 absolute geometry가 opening size, blockage ratio, characteristic length, flow development에 직접 영향을 준다.
- 따라서 “형상 비율은 비슷하지만 절대 길이가 틀린 scene”은 downstream CFD 관점에서 중요한 오류일 수 있다.

---

## 3. 새 평가 정의

기존 image-conditioned evaluation을 다음 두 층으로 해석한다.

### A. 기존 결과
- **no-explicit-scale-hint baseline**
- 이미지 정보만으로 구조를 추론한 결과
- 이미 완결된 100/100 benchmark 결과는 이 setting의 baseline artifact로 보존 가능
- 단, 최종 논문에서 main evidence로 그대로 둘지 여부는 재검토 필요

### B. 새 결과
- **scale-hinted image-conditioned evaluation**
- 이미지와 함께 lightweight absolute scale anchor를 제공
- 앞으로의 주력 evaluation setting은 이쪽으로 전환

현재 판단으로는 논문 본문에서는 최소한:
- no-scale-hint baseline
- scale-hinted setting
를 비교하거나,
- 적어도 scale-hinted setting을 main setting으로 제시하고 baseline은 보조/한계 분석으로 다루는 것이 타당하다.

---

## 4. 추천 scale hint 형태

### 우선 추천안: longest horizontal span hint

Gemini prompt에 아래와 같은 힌트를 추가한다.

- 예시 문구:
  - `Scale hint: the longest horizontal span of the room is approximately 6.2 m.`
  - `Use this as a global metric anchor when choosing room dimensions and obstacle sizes.`

### 이유
- 제공 정보가 과도하지 않다.
- practical deployment에서도 비교적 자연스럽다.
- benchmark에서는 reference scene으로부터 자동 계산이 쉽다.
- geometry 전체를 과도하게 정답 유도하지 않으면서도 absolute scale ambiguity를 줄인다.

### 대안 후보
1. **ceiling height + longest horizontal span**
   - 예: `ceiling height ≈ 2.8 m`, `longest horizontal span ≈ 6.2 m`
   - 더 안정적일 수 있으나 정보량이 조금 늘어난다.
2. **full bounding box prior**
   - 예: `overall room size ≈ 6.2 x 4.5 x 2.8 m`
   - 가장 강하지만, benchmark fairness 측면에서 과한 힌트가 될 수 있다.

현재 1차 구현은 **longest horizontal span only**가 가장 합리적이다.

---

## 5. 코드/파이프라인 수정 범위

### 필수 수정

1. `scripts/generate_indoor_scene_with_gemini.py`
   - prompt template에 scale hint 슬롯 추가
   - CLI/API 공통 prompt builder에서 optional scale hint 문자열을 반영
   - summary JSON에 어떤 scale hint가 사용되었는지 기록

2. `scripts/run_benchmark_evaluation_task.py`
   - reference scene에서 benchmark용 scale hint 자동 계산
   - task run 시 Gemini scenario prompt에 scale hint 포함
   - evaluation summary에 `scale_hint` 필드 기록

3. `scripts/run_benchmark_evaluation_batch.py`
   - batch rerun 시 scale-hinted mode를 명시적으로 호출할 수 있도록 옵션 추가 필요 가능

4. aggregate / manifest layer
   - 기존 no-scale-hint 결과와 새 scale-hinted 결과를 구분할 bookkeeping 필요
   - 최소한 evaluation summary와 aggregate summary에서 setting tag가 남아야 한다.

### 선택 수정

5. `benchmark/evaluations/.../task.json`
   - task scaffold에 precomputed scale hint를 넣을지 여부 검토
   - 런타임 계산으로 충분하면 필수는 아님

6. docs / manuscript
   - 논문 계획 문서와 결과 해석 문서에 scale-hint setting 반영

---

## 6. 새로 얻어야 할 데이터 범위

### 다시 생성할 필요가 없는 데이터
- benchmark scenes (`benchmark/scenes/`)
- reference CFD (`benchmark/reference_cfd/`)
- 5-view benchmark renderings (`benchmark/renderings/`)
- dataset card / integrity summary의 기본 구조

### 새로 생성해야 하는 데이터

1. **scale-hinted image-conditioned evaluation outputs**
   - per-task predicted scene JSON
   - per-task stabilization summary
   - per-task evaluation summary
   - per-task CFD metrics

2. **scale-hinted aggregate summaries**
   - 전체 100-task 집계
   - view/category split 집계
   - 필요 시 baseline 대비 비교 summary

3. **비교용 결과 문서**
   - no-scale-hint vs scale-hinted 성능 비교 note
   - scale error 개선 여부 note
   - 구조/CFD 측면에서 어느 정도 개선되었는지 정리

4. **figure regeneration input**
   - scale-hinted aggregate 숫자
   - scale-hinted representative cases
   - 기존 figures 5/6/8/9/10 중 일부의 source data 교체 가능성 검토

---

## 7. 새로 필요한 문서화 항목

반드시 남겨야 할 문서:

1. **이 문서**
   - scale issue 발견 및 pivot 배경

2. 향후 추가 문서
   - `26-03-19_scale_hint_implementation_plan.md`
   - `26-03-19_scale_hint_eval_rerun_plan.md`
   - `26-03-19_scale_hint_vs_baseline_comparison.md`

문서에 남겨야 할 핵심 포인트:
- benchmark dataset 자체는 유지되며 invalidated되지 않음
- invalidated되는 것은 기존 image-conditioned evaluation setting의 “main-setting” 지위임
- scale hint는 practical deployment에서도 자연스러운 보조 정보임
- 이는 leakage라기보다 deployment-realistic prior에 가깝다는 해석이 가능함

---

## 8. figure 작업에 대한 즉시 영향

- 현재 작동 중이던 figure 생산 cron은 중단한다. (완료)
- 기존 QC 통과 figure들도 당장은 최종본으로 보지 않는다.
- 이유는 단순 figure formatting 문제가 아니라, **source evaluation setting 자체가 재정의될 가능성**이 생겼기 때문이다.
- scale-hinted evaluation 결과가 나온 뒤 아래를 다시 결정한다.
  - 어떤 aggregate figure를 본문 main result로 둘지
  - 어떤 case figure를 representative example로 둘지
  - baseline/no-scale-hint 결과를 supplementary로 내릴지 여부

---

## 9. 현재 작업 우선순위 재정렬

### 최우선
1. scale hint 방식 확정
2. prompt / runner 수정
3. bookkeeping 방식 확정 (baseline vs scale-hinted 구분)
4. scale-hinted evaluation rerun 범위 확정

### 다음
5. 대표 subset으로 smoke test
6. full benchmark rerun
7. aggregate summary 재생성
8. figure source 재선정

### 후순위
9. paperbanana framework figure
10. 논문용 final figure full production

---

## 10. 당장 합의된 운영 결론

- benchmark dataset는 그대로 유지한다.
- Gemini가 2D benchmark image로부터 scene을 생성하는 단계에는 **길이 힌트(scale hint)** 를 넣는 방향으로 간다.
- 현재 figure 생산은 멈춘 상태에서, 먼저 방법론 수정과 재평가 범위를 식별한다.
- 이후에만 figure 작업을 재개한다.
