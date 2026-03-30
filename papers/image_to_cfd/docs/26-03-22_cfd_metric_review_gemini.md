# CFD Fidelity Score 학술적 적합성 검토 — Gemini (2026-03-22)

**검토 대상**: CFD Fidelity Score (4-component unweighted mean)
**검토 모델**: Gemini 3.1 Pro (via Gemini CLI)

---

## 1. 수학적 정의의 명확성

메트릭은 수학적으로 정의되어 있으나 **압력(Pressure)에서 물리적 모호성**이 존재.
- 비압축성 CFD(`simpleFoam`)에서 압력은 상대값(p/ρ)이며, 절대값은 기준점에 의존
- 예측 케이스가 다른 기준점을 사용하면 gradient가 동일해도 pressure similarity가 0으로 떨어짐
- **체계적 편향(systematic bias)** 발생 가능

## 2. Relative RMSE의 문제점

- **정체(Stagnation) 편향**: 실내 환경의 넓은 사영역(dead zone)에서 `ref_RMS`가 극히 작아 (~10⁻⁴ m/s), 수치 노이즈에 과민
- **이상치 민감성**: RMSE가 국소 오차(inlet jet 위치 이동 등)를 이차적으로 과벌

## 3. 균등 가중치의 정당성

**엄밀하게 정당화되지 않음.**
- 10% overlap에서 유동장이 완벽하게 일치하면 0.775점 → 90% 틀린 기하구조에 높은 점수
- **권장**: Overlap을 곱셈 가중치로 변경 → `Score = Overlap × mean(Field Components)`

## 4. 18×18×10 샘플링 그리드

- ~3,240 포인트는 **조밀 전역 유사도 점수**로는 적절
- CFD 메시 자체가 조밀하므로(~10-25k 셀) 셀당 약 1개 샘플 포인트
- 경계층, 좁은 제트류는 놓칠 수 있음
- Epsilon 마진(5-8%)은 벽면 샘플링 아티팩트 방지에 합리적

## 5. 기존 CFD 검증 메트릭과의 비교

현재 메트릭은 **"ML 유사도"에 가깝고 "CFD 검증"에는 부족**:
- ASME V&V 20, COST Action 732 등은 다음을 선호:
  - **Hit Rate (q)**: 예측값이 특정 임계값 내에 있는 비율
  - **FAC2**: 2배 이내 예측 비율
  - **상관계수 (R)**: 구조적 유동 유사도
- CFD 리뷰어가 물리적 타당성을 신뢰하기 어려움

## 6. 잠재적 편향 및 맹점

- **체적 평균화**: 방 전체의 direction cosine 평균은 정체 공기가 고속 영역의 오차를 가림
- **스케일-충실도 결합**: 정규화 격자가 절대 치수 오차를 숨김 → 2배 큰 방도 유동 패턴 유사하면 높은 점수

## 7. 제출 전 개선 권장사항

1. **압력 정규화**: 겹침 영역에서 평균 압력 차감, 동압(0.5ρU²avg)으로 정규화
2. **가중 필드 메트릭**: direction cosine와 RMSE를 국소 속도 크기로 가중 → 고속 유동 영역 집중
3. **표준 메트릭 추가**: Hit Rate (|U_pred − U_ref| < 0.25 U_ref) 병행 보고
4. **곱셈적 Overlap**: `Overlap × mean(Field Components)`로 재정의

## 8. 종합 판정

- **공학/CFD 저널**: ❌ 현재 상태로는 출판 준비 불충분 (압력 정규화 부재, 가산적 가중치)
- **AI/CV 저널 (CVPR, ICCV 등)**: ⚠️ "reconstruction similarity score"로 프레이밍하면 수용 가능
  - Structural Score(0.781)와 CFD Score(0.491) 간 gap을 VLM의 미세 유동 세부사항 포착 불능의 발견으로 논의해야 함
