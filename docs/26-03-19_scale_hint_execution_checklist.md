# Scale-Hint Pivot Execution Checklist

> 작성일: 2026-03-19
> 토픽: SurfClaw / topic 2242
> 목적: image-conditioned Gemini scene generation 단계에 absolute scale hint를 도입하고, benchmark evaluation을 scale-hinted setting으로 재정의/재실행하기 위한 실행 체크리스트

## 운영 원칙

- benchmark dataset 자체(`scene JSON`, `reference CFD`, `5-view render bundle`)는 유지한다.
- 다시 해야 하는 것은 **Gemini 기반 image-conditioned generation/evaluation 경로**다.
- 기존 결과는 폐기하지 않고 **no-scale-hint baseline** 으로 보존한다.
- 새 결과는 **scale-hinted evaluation** 으로 별도 관리한다.
- 사용자 판단이 꼭 필요한 경우가 아니면 자율 진행한다.
- 상태 갱신은 이 문서를 기준으로 한다.
- 진행 보고는 한국어로 짧고 의미 있게 한다.
- 안정적 체크포인트가 되면 OpenFOAM 저장소에 커밋하되, push는 사용자 요청 시만 한다.

## 상태 표기

- [ ] 미착수
- [-] 진행중
- [x] 완료
- [~] 보류 / 추가 판단 필요
- [!] 이슈 발생 / 수정 필요

---

## A. 문제 정의와 설정 분리

- [x] benchmark dataset 자체는 유지하고 image-conditioned evaluation setting만 재정의한다는 원칙 문서화
- [x] scale issue / pivot 배경 문서 작성 (`26-03-19_scale_hint_pivot_plan.md`)
- [x] 기존 결과를 `no-scale-hint baseline`으로 명시하는 bookkeeping 전략 확정
- [x] 새 결과를 `scale-hinted` setting으로 분리 저장하는 디렉터리/manifest 전략 확정
- [-] 논문 계획 문서에 scale-hint setting 도입 사실 반영

---

## B. Scale hint specification 확정

- [x] 1차 기본안 확정: `longest horizontal span` 힌트 사용 여부 결정
- [-] 대안안 검토: `ceiling height + longest horizontal span` 필요성 판단
- [x] benchmark reference scene에서 scale hint를 자동 계산하는 규칙 정의
  - [x] rectangular case 규칙
  - [x] composite case 규칙
- [x] Gemini prompt에 들어갈 최종 문구 초안 확정
- [x] leakage / fairness 관점 메모 문서화

---

## C. 코드 수정

### C1. Gemini generation script
- [x] `scripts/generate_indoor_scene_with_gemini.py`에 scale hint 인자 추가
- [x] prompt template에 scale hint 슬롯 추가
- [x] CLI/API 공통 경로에서 scale hint가 동일하게 반영되도록 수정
- [x] generation summary에 `scale_hint` 기록

### C2. Benchmark evaluation runner
- [x] `scripts/run_benchmark_evaluation_task.py`에서 reference scene 기반 scale hint 자동 계산 추가
- [x] task run prompt에 scale hint 삽입
- [x] evaluation summary에 `setting=no_scale_hint|scale_hinted` 같은 식별자 추가
- [x] evaluation summary에 실제 사용된 `scale_hint` 값 기록

### C3. Batch / aggregate / bookkeeping
- [-] `scripts/run_benchmark_evaluation_batch.py`에 scale-hinted mode 반영 필요 여부 점검
- [x] aggregate summary 생성 스크립트가 setting 구분을 다루도록 수정 필요 여부 점검
- [x] baseline과 scale-hinted 결과가 섞이지 않도록 경로/manifest 전략 정리

---

## D. Smoke test

- [x] 대표 소규모 subset 선정
  - [x] easy rectangular 1개
  - [x] rectangular multi-obstacle 1개
  - [x] light composite 1개
  - [x] dense composite 1개
- [x] 최소 4 case × 5 views 또는 더 작은 smoke subset 범위 확정
- [-] scale-hinted smoke test 실행
- [x] geometry/CFD/room bbox error가 baseline 대비 어떤 방향으로 바뀌는지 확인
- [x] prompt wording / hint formulation의 1차 수정 필요 여부 판단

---

## E. Full rerun scope 확정 및 실행

- [~] smoke test 결과가 타당하면 frozen-20 전체 rerun 범위 확정
- [ ] full scale-hinted evaluation 재실행
- [ ] per-task outputs 정리
  - [ ] predicted scene JSON
  - [ ] stabilization summary
  - [ ] evaluation summary
  - [ ] CFD metrics
- [ ] 실패/blocked/rerun bookkeeping 정리

---

## F. 결과 집계 및 비교

- [ ] scale-hinted aggregate summary 생성
- [-] no-scale-hint baseline vs scale-hinted 비교 summary 생성
- [ ] 특히 아래 항목 비교
  - [ ] room bbox absolute / relative error
  - [ ] structural score
  - [ ] CFD score
  - [ ] room-kind match
  - [ ] opening-wall match
- [ ] scale hint가 실제로 geometry scale mismatch를 줄였는지 결론 정리
- [ ] 결과를 문서화한 비교 note 작성

---

## G. Figure / manuscript 영향 재판단

- [ ] 어떤 figures가 scale-hinted 결과로 교체되어야 하는지 식별
- [ ] 기존 figure source 중 invalid/obsolete 해진 것 표시
- [ ] scale-hinted 결과 기준 대표 case 재선정 필요 여부 판단
- [ ] paper plan / Results-Discussion 구조에 반영할 문장 초안 정리

---

## H. 현재 immediate next actions

- [x] scale hint 최종 스펙 문구 결정
- [x] generation/evaluation 코드 수정 시작
- [x] smoke subset 정의
- [-] smoke run 실행

---

## 메모

- 현재 figure 생산 cron은 중단된 상태다.
- figure 작업은 scale-hinted evaluation의 방향이 확인된 뒤 재개한다.
- 기존 100/100 evaluation은 버리지 말고 baseline artifact로 유지한다.
