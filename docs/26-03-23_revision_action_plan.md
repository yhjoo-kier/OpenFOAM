# Revision Action Plan — Codex Q1 Review Response

> Date: 2026-03-23
> Review verdict: Minor Revision
> Target: Building and Environment (Q1, IF ~7.1)

---

## Major Issues (3건)

### M1. Reference CFD 검증 — 문헌 실험 케이스 대비 검증
- **지적**: CFD-to-CFD 비교만 있고, reference CFD 자체의 타당성이 미검증
- **작업**:
  1. 문헌에서 검증된 실내 환기 실험 케이스 1개 선정 (Nielsen room / IEA Annex 20 등)
  2. 동일 솔버 설정(simpleFoam, k-ω SST, 0.18m mesh)으로 시뮬레이션
  3. 실측 데이터 대비 Hit Rate, FAC2, NMSE, FB 보고
  4. Section 2.2에 validation subsection 추가
- **난이도**: 중 (실험 데이터 확보 + 시뮬레이션 + 메트릭 산출)
- **코드 작업**: 필요 (validation 케이스 셋업 + 메트릭 계산)
- **논문 작업**: Section 2.2에 ~300 단어 + Table 1개 추가

### M2. VLM 일반화 — scope 제한 명시 또는 추가 VLM 비교
- **지적**: Gemini 3.1 Pro만 테스트, 다른 VLM 비교 없음
- **선택지**:
  - **A (추천)**: Limitations에 model-specific scope 명시 + framework가 VLM-agnostic한 이유 설명 (~200 단어)
  - **B (이상적)**: GPT-4o로 대표 케이스 4-8개 추가 실행 → 비교 표 추가 (API 비용 + 시간 필요)
- **난이도**: A는 낮 (텍스트만), B는 중-높 (API 비용 + 프롬프트 적응 + 실행)
- **논문 작업**: A는 Section 5.6에 문단 추가, B는 추가 Table + Discussion

### M3. Velocity magnitude ≈ 0 원인 분석
- **지적**: 전체 CFD agreement를 끌어내리는 핵심 컴포넌트인데 분석이 1문장
- **작업**:
  1. 차원 오차 → 단면적 변화 → inlet velocity 변화 back-of-envelope 계산
  2. 대표 케이스 2-3개에서 ref vs pred의 Umag 분포 비교 (histogram 또는 profile)
  3. Section 4.1 또는 5.3에 dedicated paragraph 추가 (~400 단어)
  4. "screening에서 속도 크기보다 유동 패턴(방향)이 중요한 경우" 실용적 의미 논의
- **난이도**: 낮-중 (계산은 간단, 추가 figure 1개 가능)
- **코드 작업**: 간단한 분석 스크립트
- **논문 작업**: ~400 단어 + 선택적 figure 1개

---

## Minor Issues (7건)

### m1. VLM 반복성 테스트 한계 명시
- **작업**: Section 4.5에 "The small sample size (3 cases × 3 runs) limits the statistical power of this analysis" 1문장 추가
- **난이도**: 매우 낮

### m2. "Second" 중복 수정 (Section 5.6)
- **작업**: 두 번째 "Second"를 "Third"로 변경, 이후 넘버링 조정
- **난이도**: 매우 낮

### m3. SD 보고 일관성 (Section 4.2)
- **작업**: 모든 뷰의 structural + CFD score에 SD 표기 통일
- **난이도**: 매우 낮

### m4. Grid independence 비단조 수렴 설명
- **작업**: A3-03의 score 감소가 "finer reference가 prediction의 오차를 더 드러내기 때문"이라는 설명 1-2문장 추가
- **난이도**: 매우 낮

### m5. Naive baseline opening match rate 오해 방지
- **작업**: Table 1에 footnote 또는 inline 설명 강화 ("the naive baseline's higher rate is due to chance alignment")
- **난이도**: 매우 낮

### m6. 등온 해석 vs 에너지 진단 동기 gap 논의
- **작업**: Section 5.6 Limitations에 "The isothermal assumption excludes thermal stratification..." 문단 추가
- **난이도**: 낮 (~100 단어)

### m7. "Photograph" 주장 vs 실제 테스트 불일치
- **작업**: Abstract의 "photograph or architectural drawing" → "architectural drawing or rendered image"로 수정, 또는 실제 사진 1장 테스트 추가
- **난이도**: 텍스트 수정만이면 매우 낮, 사진 테스트면 중

---

## 추가 개선 (리뷰에서 직접 지적은 아니나 품질 향상)

### E1. 표준 메트릭 값 본문 보고
- **작업**: 현재 supplementary에만 있는 Hit Rate, FAC2 값을 본문 Results에 1-2줄 추가
- **난이도**: 낮

### E2. "Screening-level"의 정량적 정의
- **작업**: "sufficient to identify the dominant recirculation zone location and primary jet path" 등 구체적 기술
- **난이도**: 낮

### E3. VLM API 접근 날짜 및 모델 버전 명시
- **작업**: Section 2.1에 "accessed via API on [date], model version gemini-3.1-pro-preview" 1문장
- **난이도**: 매우 낮

### E4. 고정 유량 가정 명시적 논의
- **작업**: Section 2.2에 "The fixed volume flow rate is a deliberate simplification..." 2-3문장
- **난이도**: 매우 낮

---

## 우선순위 및 의존성

```
[즉시 실행 가능 — 텍스트 수정만]
├─ m1 VLM 반복성 한계 명시
├─ m2 "Second" 중복 수정
├─ m3 SD 보고 일관성
├─ m4 Grid independence 설명
├─ m5 Naive baseline footnote
├─ m6 등온 해석 gap 논의
├─ m7 "Photograph" 수정
├─ E1 표준 메트릭 본문 보고
├─ E2 Screening-level 정의
├─ E3 API 접근 정보
└─ E4 고정 유량 논의

[분석 작업 필요]
├─ M3 Vel. magnitude 원인 분석 (back-of-envelope + optional figure)
└─ M2-A VLM scope 제한 명시 (텍스트)

[시뮬레이션 필요]
├─ M1 Reference CFD 검증 (문헌 실험 케이스)
└─ M2-B GPT-4o 비교 (선택적)
```

---

## 예상 소요

| 작업 그룹 | 항목 수 | 예상 시간 |
|----------|---------|----------|
| 텍스트 수정 (m1-m7, E1-E4) | 11 | 30분 |
| M3 분석 + 본문 | 1 | 1-2시간 |
| M2-A scope 명시 | 1 | 30분 |
| M1 검증 시뮬레이션 | 1 | 반나절 (데이터 확보 + 실행 + 분석) |
| M2-B GPT-4o 비교 (선택) | 1 | 수시간 (API 비용) |
