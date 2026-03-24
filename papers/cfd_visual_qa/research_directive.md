# CFD Visual QA — Autonomous Research Directive

> Version: 1.0
> Created: 2026-03-24
> Last updated: 2026-03-24
> Authority: User may edit this file at any time. Autonomous sessions may propose edits but MUST NOT change sections marked [USER-ONLY].

---

## [USER-ONLY] Core Constraints

- **라벨링 예산**: 하루 20개, 2주 = 최대 280개
- **라벨링 방식**: 텔레그램으로 이미지 + 질문 전송, 유저가 JSON으로 응답
- **라벨링 요청 시간**: 매일 오전 9시 (KST)
- **작업 브랜치**: `paper/cfd-visual-qa`
- **커밋 규칙**: 매 작업 세션 종료 시 커밋 & 푸시
- **CFD 도메인**: 열유동 (자연대류, 강제대류, 혼합대류)
- **타겟 규모**: 이미지 ~250-300장, QA ~1,000-1,200개
- **논문 분량**: 벤치마크 + 평가 결과, 8,000 단어 이내

## Research Phases

### Phase 1: Pilot (Day 1-3)
- S4 (backward-facing step) 시나리오로 파일럿
- 올바른 케이스 2개 (laminar Re=200, turbulent Re=36000)
- 오류 케이스 4개 (E1 미수렴, E2 BC오류, E5 메시부적절, E6 인위적변형)
- 시각화: V1 velocity contour, V4 vector, V5 streamlines
- 질문 생성: L1-L5 각 레벨별 2-3개씩
- 예상 산출물: 이미지 ~30장, QA ~60개

### Phase 2: Preliminary VLM Evaluation (Day 4-5)
- 파일럿 데이터로 GPT-4o, Gemini 2.5 Pro 초기 평가
- 질문 난이도 보정, 시각화 품질 조정
- 벤치마크 설계 v2 확정

### Phase 3: Full Scenario Expansion (Day 6-10)
- S1-S10 전체 시나리오 CFD 케이스 생성
- 시각화 생성 + 질문 생성
- 매일 라벨링 배치 전송

### Phase 4: Full VLM Evaluation (Day 11-12)
- 5+ 모델 전체 평가
- 결과 테이블, 분석

### Phase 5: Paper Draft (Day 13-14)
- 논문 초안 작성

## Daily Session Protocol

매일 자율 작업 세션은 다음 순서를 따른다:

1. **상태 읽기**: `research_state.json` 읽어 현재 진행 상태 파악
2. **라벨 수거**: `labeling_queue/completed/` 에 새 라벨이 있으면 수거하여 `benchmark/` 에 통합
3. **작업 실행**: 현재 phase에 맞는 작업 수행
   - CFD 케이스 생성 (OpenFOAM via Docker)
   - 시각화 렌더링 (matplotlib + PyVista)
   - 질문 생성
   - VLM 평가 (해당 phase일 때)
4. **라벨링 배치 준비**: 다음 날 라벨링할 이미지 20개 + 질문을 `labeling_queue/pending/` 에 저장
5. **상태 업데이트**: `research_state.json` 갱신
6. **커밋 & 푸시**: 모든 산출물 커밋
7. **진행 로그**: `docs/daily_log.md` 에 당일 작업 내용 기록

## Labeling Protocol

### 라벨링 요청 형식 (텔레그램 간소화 — 방안 A)
- 이미지를 5개씩 4묶음으로 나눠 전송
- 각 이미지에 번호 + 핵심 질문 1개 (L3 물리 판단 위주)
- 답변 형식: "1:OK 2:이상 3:OK 4:이상-recirculation없음 5:OK" 식으로 한 줄
- L1/L2/L4 상세 질문은 AI가 자동 라벨링 후, 의심스러운 것만 추가 확인 요청
- 예상 소요: 10-15분

### 라벨링 응답 처리
유저 텔레그램 답변을 AI가 파싱하여 `labeling_queue/completed/day_NN/` 에 JSON으로 저장.

## Figure QC Standards

- 해상도: 1200×900 px 이상
- 컬러맵: 시나리오별 고정 (velocity=viridis, temperature=coolwarm, pressure=RdBu_r)
- 컬러바: 항상 포함, 우측 배치, 물리량+단위 라벨
- 축: x, y 라벨 포함 (SI 단위)
- 배경: 흰색
- 폰트: Arial 우선 → Liberation Sans → DejaVu Sans
- 자동 QC: Gemini CLI로 시각화 품질 확인 (가능 시)

## Direction Change Policy

자율 세션이 directive를 수정할 수 있는 경우:
- Phase 일정 미세 조정 (±1일)
- 시나리오 우선순위 변경
- 질문 난이도 보정 (예비 평가 결과 기반)
- 시각화 파라미터 조정

자율 세션이 수정할 수 없는 경우 ([USER-ONLY] 영역):
- 라벨링 예산/방식 변경
- 핵심 연구 범위 변경
- 타겟 저널 결정
- 논문 구조 대폭 변경
→ 이 경우 텔레그램으로 유저에게 제안하고 승인 대기
