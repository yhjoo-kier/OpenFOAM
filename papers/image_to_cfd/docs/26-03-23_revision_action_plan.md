# Revision Action Plan — Codex Q1 Review Response

> Date: 2026-03-23
> Review verdict: Minor Revision
> Target: Building and Environment (Q1, IF ~7.1)
> Last updated: 2026-03-23 (preset matching 반영)

---

## Major Issues (3건)

### M1. Reference CFD 검증 — 문헌 실험 케이스 대비 검증
- **상태**: 🔧 **진행 중** — Nielsen (1990) 케이스 셋업 + 해석 완료, 메트릭 산출 완료
- **결과**: q=0.083, R=0.837 — profile shape는 잘 맞으나 absolute magnitude 차이 있음
- **남은 작업**: 논문 Section 2.2에 validation subsection 작성 (~300 단어 + Table)
- **참고**: `cases/validation_nielsen_1990/`, `validation_results.json`

### M2. VLM 일반화 — scope 제한 명시
- **상태**: 🔧 **M2-A 선택, 미반영** — Limitations에 model-specific scope 문단 추가 필요
- **남은 작업**: Section 5.7 Limitations에 ~200 단어 추가

### M3. Velocity magnitude ≈ 0 원인 분석
- **상태**: ✅ **해결** — inlet velocity mismatch 발견 + preset matching으로 수정
- **경과**:
  1. Scientist agent가 근본 원인 진단: solver preset별 inlet velocity 불일치 (85/97 케이스)
  2. `--force-preset` 플래그 추가 + 26개 reference 재해석
  3. vel_mag: 0.023 → 0.117, CFD agreement: 0.454 → 0.477
  4. 논문 수치 전체 업데이트 + Section 5.5 "Reference solver preset matching" 추가
- **문서**: `docs/26-03-23_inlet_velocity_mismatch_issue.md`

---

## Minor Issues (7건)

| # | 항목 | 상태 |
|---|------|------|
| m1 | VLM 반복성 테스트 한계 명시 | ✅ 완료 |
| m2 | "Second" 중복 수정 | ✅ 완료 |
| m3 | SD 보고 일관성 | ✅ 완료 |
| m4 | Grid independence 비단조 수렴 설명 | ✅ 완료 |
| m5 | Naive baseline opening match footnote | ✅ 완료 |
| m6 | 등온 해석 vs 에너지 진단 gap | ✅ 완료 |
| m7 | "Photograph" → "architectural drawing" | ✅ 완료 |

---

## 추가 개선 (E1-E4)

| # | 항목 | 상태 |
|---|------|------|
| E1 | 표준 메트릭 값 본문 보고 (q, FAC2, R) | ✅ 완료 |
| E2 | "Screening-level" 정량적 정의 | ✅ 완료 |
| E3 | VLM API 접근 날짜/버전 명시 | ✅ 완료 |
| E4 | 고정 유량 가정 논의 | ✅ 완료 |

---

## 남은 작업 요약

| 작업 | 내용 | 예상 시간 |
|------|------|----------|
| **M1 논문 반영** | Nielsen validation subsection 작성 + Table | 1시간 |
| **M2-A 논문 반영** | VLM scope 제한 문단 작성 | 30분 |
| **커밋/푸시** | 위 2건 반영 후 | 5분 |

완료율: **17/19 항목 완료 (89%)**
