# Paper Revision Plan

> Date: 2026-03-21 (updated 2026-03-22)
> Scope: peer-review 피드백 기반 개선 방향 정리 (Introduction은 별도 후속 작업)
> Target: Building and Environment (IF ~7.1, Q1)
>
> **IMPORTANT**: Phase 2 전면 재실행 (0.18m mesh, timeout 3600s) 결과가 논문의 메인 결과를 대체한다.
> 기존 0.35m 결과는 `benchmark/evaluations_posthoc_scaled_035m_backup/`에 백업됨.
> 논문의 모든 정량 결과(Table 1-5, Figure 5-6, aggregate 통계)는 Phase 2 완료 후 재산출해야 한다.

---

## 1. 코드-본문 불일치 수정 (Critical — 제출 전 필수)

### 1.1 VLM 모델명

논문에 "Gemini 2.0 Flash"로 기술되어 있으나 실제 평가 코드는 `gemini-3-flash-preview`를 호출한다. 모든 본문, Abstract, Methodology에서 실제 모델명으로 교체하고, preview 모델 사용에 따른 재현성 한계를 Methodology 말미 또는 Limitations에 한 문장으로 명시한다.

### 1.2 Structural score 정의

현재 본문은 "weighted combination of room bounding-box relative errors, room-kind match, opening-wall assignment, obstacle F1 at IoU 0.18"로 기술하고 있으나 실제 구현은 아래 네 항목의 unweighted average이다.

- room_block_match F1
- obstacle_match F1 (IoU threshold 0.1)
- openings type F1
- openings wall_match_ratio

Section 3.2를 이 네 항목으로 정확히 재작성한다. IoU threshold를 0.1로 수정하고, 수식을 한 줄로 명시한다.

### 1.3 CFD score 정의

현재 본문의 "velocity magnitude correlation, pressure field correlation, bulk flow statistics"는 실제 구현과 다르다. 실제 구현은 아래 네 항목의 unweighted average이다.

- overlap_ratio (predicted/reference 도메인 교집합 비율)
- velocity magnitude similarity: max(0, 1 − min(1, relative RMSE))
- velocity direction similarity: 0.5 × (mean direction cosine + 1)
- pressure similarity: max(0, 1 − min(1, relative RMSE))

Section 3.2를 이에 맞게 재작성하고 수식을 명시한다.

### 1.4 메시 해상도

본문의 "target cell size of 0.05 m"은 실제 기본값 0.35 m과 7배 차이가 난다. 0.35 m로 수정하고, robustness ladder (0.35 → 0.25 → 0.18 m)를 정확히 기술한다. 이 해상도에서의 CFD 해의 한계(y+ 값, wall function 적용 조건, grid independence 미보장)를 Discussion의 Limitations 소절에서 명시적으로 논의한다.

---

## 2. 실험 설계 보강

### 2.1 Baseline 추가 (필수)

최소 두 가지 baseline을 추가한다.

**Baseline A — Naive geometric default.** 모든 케이스에 대해 5 m × 4 m × 3 m 직사각형 방, 장애물 없음, 중앙 대칭 inlet/outlet을 예측값으로 사용했을 때의 structural score 및 CFD score를 산출한다. 이 baseline이 특히 A1 카테고리에서 어느 정도의 점수를 얻는지 보여줌으로써 0.781이라는 structural score의 절대적 의미를 해석할 수 있는 맥락을 제공한다.

**Baseline B — No scale calibration ablation.** 현재 파이프라인에서 post-hoc scaling을 제거하고 VLM raw output 그대로 평가한 결과를 제시한다. 이를 통해 scale calibration의 기여도를 정량적으로 분리하고, "reference 정보를 하나 사용한다"는 설계 선택의 정당성을 데이터로 뒷받침한다.

결과 표에 Overall mean을 Naive / No-scaling / Scale-calibrated 세 조건으로 나란히 제시하면 리뷰어가 성능을 맥락적으로 해석할 수 있다.

### 2.2 통계 보강

현재 모든 결과가 평균값만 제시되어 있고 분산 정보가 없다. 다음을 추가한다.

- Figure 5, 6의 bar chart에 표준편차 whisker 추가
- Table 1에 표준편차 열 추가
- 뷰 간 성능 차이에 대해 paired Wilcoxon signed-rank test 수행 (floorplan vs section 등 주요 비교)
- n=5/셀의 한계를 Limitations에서 명시

### 2.3 Null 케이스 공개

by_room_kind에 room_kind = None인 3개 케이스가 존재한다. 이 케이스들의 발생 원인(VLM 파싱 실패 vs 집계 아티팩트)을 확인하고, Methodology 또는 Results에서 명시적으로 기술한다. 집계 통계에서 이 3개가 어떻게 처리되는지(제외 vs 0점 처리) 밝힌다.

---

## 3. Scale Calibration 프레이밍 조정

현재 논문 제목과 Abstract가 "fully automated"를 암시하나, scale calibration 단계에서 reference의 longest span 하나를 사용한다. 이 모순을 해소하기 위해 다음과 같이 프레이밍을 조정한다.

- Abstract에서 "single reference dimension을 사용한 post-hoc scale calibration"을 명시적으로 언급
- Methodology에서 이 선택의 근거를 서술: VLM이 비율은 잘 맞추지만 절대 치수를 모르므로, 하나의 기준 치수만으로 전체를 보정하는 것이 실용적 절충안
- Baseline B (no-scaling ablation)의 결과로 calibration 없이는 성능이 어느 정도 저하되는지 정량적으로 보여줌
- Discussion에서 향후 self-calibration 접근(문, 표준 가구 등 known object 기반)을 구체적으로 논의

제목 자체는 "automated"를 유지하되, Abstract 첫 단락에서 one reference dimension 의존성을 투명하게 밝히는 방식으로 정직성을 확보한다.

---

## 4. Real-Image Demonstration 섹션 추가

리뷰어의 "합성 이미지만 사용" 지적과 사용자의 구상을 반영하여, 벤치마크 평가(Section 5-6) 다음에 실제 사진 적용 사례 섹션을 추가한다.

### 4.1 구성

Section 7 (현재 Conclusion) 앞에 새로운 섹션을 삽입한다.

**Section 7. Application to Real Photographs**

벤치마크 결과에서 가장 효과적인 뷰 타입으로 확인된 floor plan 뷰를 기준으로, 실제 실내 공간의 평면도 사진 또는 건축 도면 2건을 확보하여 프레임워크를 적용한다. 이 섹션은 정성적 feasibility demonstration 성격이며, 정량 평가 (reference CFD 대비)는 수행하지 않는다.

### 4.2 케이스 선정 기준

- Case 1: 비교적 단순한 실내 공간 (사무실, 원룸 등)의 실제 평면도 — A1/A2급 난이도에 대응
- Case 2: L자형 또는 복합 구조 공간의 실제 평면도 — A3/A4급 난이도에 대응
- 소스: 인터넷 공개 평면도, 부동산 매물 도면, 또는 직접 촬영한 사진

### 4.3 제시 내용

각 케이스에 대해 다음을 제시한다.

- 입력 이미지 (실제 평면도 사진)
- VLM이 추출한 3D scene JSON의 plan-view 시각화
- 생성된 CFD 결과 (velocity magnitude contour)
- Scale calibration에 사용한 reference dimension의 출처 (예: 도면에 표기된 치수, 문 폭 0.9 m 등)

### 4.4 논의 포인트

- 합성 벤치마크 대비 실제 이미지에서 VLM 출력 품질이 어떻게 달라지는지 정성적으로 관찰
- 실제 도면의 노이즈(텍스트, 가구 아이콘, 치수선 등)가 VLM 해석에 미치는 영향
- Self-calibration 가능성: 도면에 표기된 치수 또는 문 폭 등의 known dimension 활용

### 4.5 논문 구성 변경

| 기존 | 개선 |
|------|------|
| 5. Results | 5. Results (유지) |
| 6. Discussion | 6. Discussion (유지) |
| 7. Conclusion | **7. Application to Real Photographs (신규)** |
| — | **8. Conclusion** (번호 변경) |

---

## 5. 기타 보강 사항

### 5.1 Computational cost

Methodology 또는 Results 말미에 파이프라인 소요 시간을 기술한다. VLM API 호출 시간 (케이스당 ~10-30초), 메싱 시간, 솔버 시간, 전체 end-to-end 시간을 대략적으로 보고한다. "rapid screening" 용도를 주장하려면 이 정보가 필수적이다.

### 5.2 Prompt 공개

VLM에 전달하는 structured prompt의 전체 텍스트를 Appendix 또는 Supplementary Material로 제공한다. Prompt sensitivity가 VLM 연구에서 중요한 변수이므로, 재현성을 위해 필수적이다.

### 5.3 IoU threshold sensitivity

Structural score의 obstacle matching IoU threshold를 0.1, 0.2, 0.3, 0.5로 변경했을 때 aggregate structural score가 어떻게 변하는지 간단한 sensitivity table을 Appendix에 추가한다.

### 5.4 VLM 출력 일관성

동일 이미지에 대해 VLM을 3-5회 반복 호출하여 출력의 stochastic variation을 측정한다. Temperature 설정값을 명시하고, 반복 간 structural score의 변동 범위를 보고한다. 본문에서 1-2문장으로 언급하고 상세는 Appendix로 돌린다.

### 5.5 Grid independence 언급

0.35 m 기본 해상도에서 grid independence가 보장되지 않음을 Limitations에서 명시한다. 1-2개 대표 케이스에 대해 0.35 → 0.25 → 0.18 m으로 refinement했을 때의 CFD score 변화를 보고하면 이상적이지만, 100개 케이스 전체 재실행은 비현실적이므로 대표 케이스로 한정한다.

### 5.6 Scatter plot 추가

100개 전체 케이스의 structural score vs CFD score scatter plot을 Figure로 추가한다. 카테고리별 색상, 뷰 타입별 마커로 인코딩하면 structure-fidelity decoupling을 Figure 10보다 더 직관적으로 보여줄 수 있다. 이 figure는 Discussion 도입부에 배치한다.

---

## 6. 수정 우선순위

| 순위 | 항목 | 난이도 | 영향도 |
|------|------|--------|--------|
| 1 | ~~코드-본문 불일치 수정 (1.1–1.4)~~ | ~~낮음~~ | **완료 (2026-03-21)** — 모델명 gemini-3.1-pro-preview, structural score 4항목 unweighted mean (IoU 0.2/0.1), CFD score 4항목 공식 명시, 메시 0.35m+ladder 수정, Abstract에 scale cal. 의존성 명시 |
| 2 | ~~Baseline A, B 추가 (2.1)~~ | ~~중간~~ | **완료 (2026-03-21)** — Naive default: 0.621, No-scaling: 0.707, Scaled: 0.781. Table 1에 3조건 비교 반영 완료 |
| 3 | ~~통계 보강 (2.2)~~ | ~~낮음~~ | **완료 (2026-03-21)** — 전체/뷰별/카테고리별 SD 산출, Table 2에 SD 반영, 3건 VLM 파싱 실패 공개 |
| 4 | Real-image demo 2건 (4) | 중간 (이미지 확보+실행) | 논문 완성도 |
| 5 | ~~Scale calibration 프레이밍 (3)~~ | ~~낮음~~ | **완료 (2026-03-21)** — Abstract에 one reference dimension 명시 |
| 6 | ~~Computational cost (5.1)~~ | ~~낮음~~ | **완료 (2026-03-21)** — Section 6.5에 3-8분/케이스 기술 |
| 7 | ~~Scatter plot 추가 (5.6)~~ | ~~낮음~~ | **완료 (2026-03-21)** — Fig 12 scatter plot 생성, Section 6.3에 참조 추가, caption 작성 |
| 8 | ~~Prompt 공개 (5.2)~~ | ~~낮음~~ | **완료 (2026-03-21)** — docs/appendix_vlm_prompt.md 생성, Appendix A에 참조 추가 |
| 9 | ~~IoU sensitivity (5.3)~~ | ~~중간~~ | **완료 (2026-03-21)** — IoU 0.05~0.50 범위 sensitivity table, Appendix B에 추가, 결론 강건성 확인 |
| 10 | VLM 반복 일관성 (5.4) | 중간 (API 비용) | 방법론 강건성 — 대기 (API 필요) |
| 11 | ~~Grid independence (5.5)~~ | ~~높음~~ | **완료 (2026-03-22)** — Phase 2 전면 재실행 (0.18m mesh, timeout 3600s). REF 20/20, PRED 97/100 성공. 3건 솔버 발산은 실패로 유지하고 논문에서 기술. |
| 14 | ~~Phase 2 메트릭 개선~~ | ~~중간~~ | **완료 (2026-03-22)** — 압력 gauge 보정, RMS floor, direction cosine 임계값, 고정 4-component, 표준 메트릭(Hit Rate/FAC2/NMSE/FB/R) 추가, cfd_agreement_score 병행 출력 |
| 15 | ~~Phase 2 재평가 + 논문 수치 업데이트~~ | ~~중간~~ | **완료 (2026-03-22)** — 97건 재평가, structural 0.781, CFD agreement 0.454. Section 3.2 메트릭 정의 대폭 확장 (물리적 의미, 0-1 해석 포함). 전체 논문 수치 업데이트. |
| 12 | ~~Null 케이스 공개 (2.3)~~ | ~~낮음~~ | **완료 (2026-03-21)** — Table 2에 3건 파싱 실패 명시 |
| 13 | ~~Component-level breakdown~~ | ~~중간~~ | **완료 (2026-03-21)** — Table 3 추가, dimension-dependent vs topology component 분리 해석 |

---

## 7. 예상 최종 논문 구조

1. Introduction (후속 작성)
2. Related Work
3. Methodology
4. Benchmark Dataset
5. Results
6. Discussion
7. **Application to Real Photographs (신규)**
8. Conclusion
9. References
10. Appendix: VLM Prompt, IoU Sensitivity, Repeatability

예상 분량: 6,000–7,000 words (Introduction + Real-image 섹션 포함 시)
