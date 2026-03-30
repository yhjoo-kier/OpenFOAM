# CFD Fidelity Score 학술적 적합성 검토 — Codex GPT-5.4 (2026-03-22)

**검토 대상**: CFD Fidelity Score (4-component unweighted mean)
**검토 모델**: OpenAI Codex (GPT-5.4)

---

## 종합 판정

현재 형태의 `CFD Fidelity Score`는 합리적인 내부 벤치마크이나, **단독 출판용 CFD 검증 메트릭으로는 학술적 엄밀성이 부족**. Task-specific surrogate agreement score에 가까움.

---

## 1. 수학적 정의의 명확성

대체로 정의되어 있으나 완전히 모호하지 않음:

- **Overlap Ratio**: Jaccard-style 점유율 점수, union이 비어있지 않으면 명확
- **Velocity/Pressure**: 절단된 정규화 RMSE이나 다음을 명확히 해야 함:
  - `ref_RMS`가 겹침 셀에서만 계산되는지 여부
  - `ref_RMS = 0`이거나 매우 작을 때의 처리
  - 공간 겹침이 없을 때의 처리
- **Velocity Direction**: 한쪽 벡터만 0에 가까울 때 cosine이 미정의 → 양쪽 모두 `|U| < 1e-12`인 경우만 제외하는 것은 불충분
- **Pressure**: 비압축성 CFD에서 압력은 가산 상수까지만 정의됨 → gauge 고정 없이 raw pressure RMSE는 무의미

## 2. Relative RMSE 문제점

- 참조 RMS가 작을 때 **불안정**
- `RMSE ≥ ref_RMS`이면 모두 0점 → 나쁜 케이스 간 **구별 불가 (포화)**
- **비대칭적** (symmetric하지 않음)
- 물리적 해석 불가: 0.6점이 표준 CFD 수용 기준에 매핑되지 않음
- **분산 민감**: 큰 참조 변동을 가진 매끈한 장이 물리적으로 틀린 저분산 장보다 높은 점수

## 3. 균등 가중치의 정당성

**정당화되지 않음.**

- Overlap, 속도 크기, 방향, 압력이 동등하게 기여해야 할 물리적/통계적 근거 없음
- 별도의 Structural Score가 있으므로 Overlap을 CFD Score에 포함하면 **기하구조 이중 계산**
- 균등 가중이 실패 모드를 숨김: 압력이나 유동 방향이 수용 불가해도 총점은 괜찮을 수 있음
- **권장**: 4개 컴포넌트를 개별 보고하고, 복합 점수는 보조 요약으로만 사용

## 4. 18×18×10 샘플링 그리드

**조밀 스크리닝에는 적절, 강한 충실도 주장에는 부족.**

- 마스킹 전 3,240개 샘플 포인트, 겹침 필터 후 더 적음
- 실내 제트, 재순환 영역, 개구부 근방, 벽면 근방 경사도에 대해 조밀
- 매끈한 예측에 유리하고 좁은 특징을 놓침
- **리뷰어 예상 요구**: 36×36×20, 72×72×40 등 고밀도 민감도 분석

## 5. 기존 CFD 검증 관행과의 비교

기존 프레임워크보다 약함:

- **ASME V&V 20**: 지정 검증 포인트에서 변수별 비교 + 불확실성 기반, 단일 비가중 복합 점수 아님
- **COST 732 (건축/도시 CFD)**: Hit Rate, FAC2/FAC1.3, Fractional Bias, NMSE, 상관계수, 프로파일별 비교
- 현재 메트릭은 CFD 해를 다른 CFD 해와 비교 → ASME 의미의 정식 검증이 아닌 **참조 솔버 대비 벤치마킹**

## 6. 잠재적 편향 및 맹점

- "참조 진실"이 simpleFoam RANS k-ω SST 벽함수 + 조밀 메시일 뿐, **ground truth 아님**
- 벌크 체적 일치가 지배적, 벽면/개구부 근방 물리 과소 평가
- Direction cosine이 좋아도 크기가 심각하게 틀릴 수 있음
- RMSE 기반 항이 체계적 부호/편향 구조를 놓침
- 압력 비교가 gauge에 민감
- 정규화 박스 샘플링이 절대 스케일 효과를 제거 → 레이놀즈 수 민감 거동 숨김
- **100 케이스 = 20 기하구조 × 5 뷰 → 군집(clustered) 데이터**, 독립 표본이 아님 → CI/유의성 검정은 기하구조 군집 부트스트랩 필요

## 7. 제출 전 개선 권장사항

1. **메트릭 이름 변경**: "CFD Fidelity Score" → `Reference CFD Agreement Score` 등
2. **엣지 케이스 정의**: 빈 겹침, 0 RMS, 0속도 방향 처리 명시
3. **압력 보정**: 도메인 평균 압력 차감 또는 고정 기준점 압력 제거
4. **4개 컴포넌트 개별 보고**, 복합 점수는 보조만
5. **표준 메트릭 병행**: NRMSE/NMSE, FB, 상관계수, FAC2/FAC1.3, Hit Rate
6. **샘플링 그리드 민감도 분석** 추가
7. **참조 CFD 불확실성 근거**: 메시 독립성, 솔버 민감도, 오차 예산
8. **물리적 QoI 추가**: 개구부 유량, 평균 거주영역 속도, 압력 강하, 재순환 크기, 환기 효율
9. **클러스터 인식 통계**: 기하구조 기반 부트스트랩 CI, naive n=100 아님

## 8. 출판 가능성

- **CFD 학술지**: ❌ 1차 검증 메트릭으로 부적합. 리뷰어 거절 예상
- **ML/CV 학술지**: ⚠️ "task-specific internal benchmark for surrogate/reference agreement"로 좁게 프레이밍 + 강한 주의사항 + 표준 CFD 메트릭 병행 보고 시 수용 가능
- 정식 CFD 검증으로 제시하는 것을 중단해야 함

## 참고 자료 (Codex 제시)

- ASME V&V 20: https://www.asme.org/codes-standards/find-codes-standards/standard-for-verification-and-validation-in-computational-fluid-dynamics-and-heat-transfer
- COST 732 Best Practice: https://theairshed.com/pdf/COST%20732%20Best%20Practice%20Guideline%20May%202007.pdf
- Building and Environment (Hit Rate/FAC2 활용 사례): https://www.sciencedirect.com/science/article/pii/S0360132324011272
