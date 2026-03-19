# Post-hoc Scaling Layer 체크리스트

> 작성일: 2026-03-20
> 목적: baseline Gemini output을 재사용하면서, 프롬프트 개입 없이 framework 내부에 최소한의 post-hoc scaling layer를 추가해 성능 변화를 빠르게 검증한다.

## 0. 작업 목표

- benchmark dataset / reference CFD / 기존 no-scale-hint baseline 결과는 유지한다.
- Gemini 프롬프트 wording 수정 없이, **baseline scene output 이후 단계**에만 개입한다.
- 재구성된 room/domain geometry에서 계산한 characteristic length를 benchmark GT characteristic length와 맞추도록 **uniform scaling** 을 적용한다.
- 기존 baseline과 동일한 downstream 절차(mesh/solver/evaluation)를 다시 수행해 성능 차이를 비교한다.
- **post-hoc scaling layer 성능 비교 결과가 확보되면 본 cron/autopilot은 중단하고 사용자와 논의한다.**

---

## 1. 기준선 / 실험 정의 고정

- [x] baseline 재사용 대상 결과 경로 확인
  - baseline root: `/home/yhjoo/projects/OpenFOAM/benchmark/evaluations`
  - 재사용 입력: 각 task의 `predicted_scene.json` (Gemini 재생성 금지)
- [x] baseline과 post-hoc scaling 비교 시 사용할 evaluation subset 정의
  - smoke subset: `bench_a1_01`, `bench_a2_03`, `bench_a3_03`, `bench_a4_03` × 5 views = 20 tasks
- [x] characteristic length 정의 1안 고정
  - 채택: reconstructed room/domain의 **longest horizontal span** = `max(Lx, Ly)`
- [x] GT characteristic length 추출 경로 확인
  - `reference_scene.json`의 room overall bbox에서 동일 기준으로 추출
- [x] uniform scaling만 허용하고 topology 변경 금지 원칙 명문화
  - room/openings/obstacles에 동일 배율 적용, wall identity / room-kind / topology는 변경 금지

## 2. 구현 포인트 확정

- [x] framework에서 post-hoc scaling을 삽입할 stage 결정
  - 결정: **baseline `predicted_scene.json` 직후**, meshing 전 stage에서 별도 `scaled_scene.json` 생성
- [x] room / openings / obstacles에 동일 배율 적용되도록 대상 필드 정리
  - room blocks `origin/size`, obstacle `min/origin/size`, opening `center/size`를 uniform scaling
- [x] inlet/outlet wall identity, topology, orientation가 보존되는지 점검 항목 정의
  - 체크 포인트: `opening_wall_match`, `room_kind_match`, obstacle count delta, opening signature 유지
- [x] baseline output 재생성 없이 downstream 재실행 경로 확보
  - `run_indoor_stabilized.py --scenario <scaled_scene.json>` 경로로 Gemini generation bypass 확인

## 3. 코드 구현

- [x] post-hoc scaling helper 또는 layer 구현
  - `scripts/posthoc_scale_indoor_scene.py`
  - `scripts/run_benchmark_posthoc_scaling_task.py`
  - `scripts/scaffold_posthoc_scaled_evaluations.py`
- [x] scaling factor 계산 로직 구현
  - factor = `reference longest_horizontal_span / baseline predicted longest_horizontal_span`
- [x] scaled artifact 저장 위치 분리
  - 각 task에 `scaled_scene.json` 별도 저장
- [x] no-scale baseline과 섞이지 않도록 result namespace 분리
  - 분리 root: `benchmark/evaluations_posthoc_scaled_longest_span`
- [x] 최소 단위 smoke test 1개로 geometry validity 확인
  - `bench_a1_01/perspective` 성공, Gemini regeneration 없이 scene_json downstream 경로 확인

## 4. 비교 실험

- [x] 소규모 subset에서 no-scale vs post-hoc scaling 비교 실행
  - smoke subset 20 tasks 전부 완료
- [x] geometry 지표 비교
  - [x] room bbox / longest span
    - mean `Lx` relative error: `0.2737 -> 0.0000`
    - 대신 mean `Ly`: `0.1492 -> 0.2209`, mean `Lz`: `0.0714 -> 0.2151` 로 악화
  - [x] opening-wall fidelity
    - opening-wall match rate: `0.75 -> 0.75` (변화 없음)
    - room-kind match rate: `1.00 -> 1.00` (변화 없음)
  - [x] obstacle count / placement 영향
    - obstacle count delta / topology identity는 subset에서 유지됨
- [x] CFD 지표 비교
  - [x] mesh success
    - 20/20 성공
  - [x] solver success
    - 20/20 성공
  - [x] field error / aggregate score
    - mean CFD score: `0.5181 -> 0.5329` (`+0.0148`)
- [x] 성능 악화 여부 우선 판단
  - 전반 평균은 소폭 개선, severe global regression은 아님
  - 다만 task별 mixed result: `10 improved / 9 worsened / 1 unchanged`
  - worst case: `bench_a3_03/birdseye -0.0636`
- [x] 개선 또는 harmless 여부가 보이면 subset 확대 여부 판단
  - 결과는 확보했으므로 우선 사용자와 논의 후 full benchmark 확대 여부 결정

## 5. 해석 / 논문화 포인트 정리

- [x] 결과 해석 메모 작성
  - `docs/26-03-20_posthoc_scaling_smoke_note.md`
- [x] prompt-side scale hint 대비 framework-side calibration 장단점 정리
  - 장점: Gemini 재생성 없이 baseline 재사용 가능, topology 보존, framework-side minimal intervention
  - 한계: single-factor uniform scaling이라 `Ly/Lz`까지 같이 흔들림
- [ ] 논문 메인 방법 후보로서의 명분 정리
  - 현재 판단: 보조 calibration layer 후보는 가능하나, full benchmark 확장 전 축 제한/anchor 범위 논의 필요
- [x] post-hoc scaling layer 성능 비교 결과 확보 시 사용자와 논의 후 다음 단계 결정
  - smoke subset 결과 확보 완료. 이 cron 종료 후 사용자와 논의 필요.

---

## 운영 규칙

- 항상 기존 no-scale-hint baseline을 보존한다.
- prompt를 다시 손대는 실험보다, **baseline output 재사용 + 최소 scaling intervention** 을 우선한다.
- user decision이 꼭 필요하지 않으면 자율 진행한다.
- 중요한 상태 변화는 이 문서를 업데이트해 남긴다.
- 성능 비교 결과가 확보되면 cron을 중단하고 짧게 보고한다.
