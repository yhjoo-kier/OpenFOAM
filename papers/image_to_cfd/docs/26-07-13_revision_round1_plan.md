# Revision Round 1 — 심사 대응 계획 & Claude Code 핸드오프 (JAAI)

- **저널**: 한국실용인공지능학회지 (JAAI)
- **원고**: Automated Indoor CFD from a Single 2D Image via VLM Abstraction (단독저자 Younghwan Joo)
- **판정 성격**: minor~moderate revision (두 심사위원 모두 게재 우호적)
- **작성일**: 2026-07-13
- **작업 방침 (사용자 확정)**:
  - **신규 시뮬레이션 없음.** 기존 벤치마크 데이터 재집계·재plot·문구 수정으로만 대응.
  - **R2-1** (opening type F1 / IoU 임계값): 지표 포함·제외 점수 병기 + 근거·참고문헌 보강.
  - **R2-4** (opening 대비 격자 조밀도): 신규 해석 없이 **한계(limitation)로 명시**.
  - 수정은 **소스 오브 트루스**(`docs/paper_draft_v1.md` + 그림 스크립트 + 데이터)에서 하고 → PDF → JAAI docx 순으로 재생성. **docx 직접 편집 금지.**

---

## 0. 핵심 리스크 요약 (먼저 볼 것)

가장 위험한 항목은 **V5 (입구 경계조건: 체적유량 vs 속도)**. `CLAUDE.md`에 "solver preset별 inlet velocity가 다름(robust=0.02, ultra_robust=0.005, laminar=0.01 m/s)"이라고 되어 있어, **실제 구현은 입구 '속도'를 preset별로 고정**한 것으로 보인다. 그런데 본문 2.2절은 "체적유량 0.05 m³/s 고정"이라고 서술 → **실제 코드와 원고 서술이 불일치할 가능성이 높다.** 심사위원 두 명이 정확히 이 지점을 지적했다. **답변서를 쓰기 전에 반드시 case의 `0/U`를 확인**하고, 확인 결과에 맞춰 2.2절·5.5절을 통일해야 한다. (잘못 답하면 신뢰도 타격)

그 외 수치 불일치(V1~V4)는 대부분 **구버전 데이터로 생성된 그림/문구**가 원인일 가능성이 높다 (`CLAUDE.md`의 preset 미매칭 `cfd_metrics.json` 경고 참조).

---

## 1. 심사평 전문 (verbatim)

### 심사위원 #1
본 논문은 비전 언어 모델을 이용하여 단일 2D 이미지로부터 실내 공간의 기하 정보를 추출하고, 이를 메시 생성과 OpenFOAM 해석으로 연결하는 자동화 프레임워크를 제안하고 있습니다. 연구 주제가 시의성 있고, 100개 사례로 구성된 벤치마크와 독립적인 기준 해석 경로를 마련한 점, 구조적 유사성과 CFD 결과를 구분하여 평가한 점, 실패 사례와 한계를 비교적 투명하게 제시한 점을 긍정적으로 평가합니다. 논문의 전체적인 구성과 연구 결과는 게재에 적합한 수준으로 판단되며, 다음 사항을 확인하고 보완해 주기 바랍니다.

1. 본문과 Table 5에서는 평균 CFD agreement score가 0.477로 제시되어 있으나, Fig. 7에는 overall mean이 0.453으로 표시되어 있습니다. 두 값의 산정 기준이 다른 경우 이를 설명하고, 표기상의 문제라면 관련 수치를 일관되게 수정해 주기 바랍니다. 또한 Table 3은 structural score를 100건 기준으로 제시하는 반면, Table 5의 제목은 97건을 기준으로 제시하고 있어 평가 대상 수가 다르게 보입니다. Structural score와 CFD agreement score의 평가 대상 수를 구분하여 명확히 표시해 주기 바랍니다.
2. Section 6의 실제 평면도 사례에서 문과 창문을 inlet과 outlet으로 지정한 것이 영상으로부터 확인된 결과인지, 사례 해석을 위해 설정한 가정인지 명확하지 않습니다. 해당 경계조건의 설정 기준을 설명하고, 실제 적용 시에는 풍량, 압력차 또는 HVAC 운전정보와 같은 추가 조건이 필요할 수 있음을 간략히 밝혀 주기 바랍니다.
3. Figure 13의 본문과 그림 내부에 제시된 수치가 캡션의 수치와 다르게 보이므로, 동일한 사례와 입력 조건을 나타내는 값인지 확인하고 필요한 경우 일관되게 수정해 주기 바랍니다. 또한 Table 7의 일부 사례에서는 0.18 m와 0.10 m 해상도 간 CFD agreement score 변화가 5%를 넘는 것으로 판단됩니다. 본문의 "less than 5%"가 다른 물리량을 기준으로 한 설명이라면 그 기준을 명확히 제시하고, 그렇지 않다면 해당 표현을 재검토해 주기 바랍니다.
4. Section 2.2에서는 동일한 체적유량을 적용하므로 입구 면적에 따라 속도가 달라지는 것으로 설명되어 있으나, Section 5.5에서는 동일한 inlet velocity boundary condition을 사용한 것으로 읽힙니다. 실제로 동일하게 적용한 조건이 체적유량인지, 입구 속도인지, 또는 경계조건의 설정 방식과 solver preset인지 명확히 설명해 주기 바랍니다.

### 심사위원 #2
이 논문은 단일 2D 이미지 → VLM 기반 구조화된 3D 형상 → 자동 격자 생성 → OpenFOAM 해석으로 이어지는 프레임워크를 제안하고 있습니다. 20개 형상을 다섯 가지 시점으로 렌더링한 100개 입력을 평가하고, 구조 유사도와 CFD 유사도를 분리해 분석했으며, 97/100개 해석이 수렴하고, 평균 structural score 0.781, CFD agreement score 0.477이 나와 충분한 가능성을 보여주었습니다. 다만, 아래와 같은 내용에 대해 수정 및 보완을 요청드립니다.

1. 성능에 대해 과장된 부분이 있다 생각됩니다. 추가 설명이나 reference를 추가해주시기 바랍니다.
   - Opening Type F1 에 대해 이미 언급하셨지만, 실제 성능보다 과장되고 있고, room IoU threshold 0.2와 obstacle threshold 0.1도 관대한 기준으로 생각됩니다.
   - Opening type F1을 종합 점수에서 제외하거나, Opening type F1을 사용해야한다면, Inlet/outlet 개수와 유형이 달라지는 benchmark를 추가하거나, 고정적으로 1.0인 항목을 포함한 점수와 제외한 점수를 모두 나타내거나, 산술평균 대신 응용 목적에 따른 가중 점수 제시하기 바랍니다. 혹은 Opening type F1을 하나의 inlet과 하나의 outlet을 갖는 논문을 reference로 제시해주시기 바랍니다.
   - 마찬가지로 room IoU threshold 0.2와 obstacle threshold 0.1에 대해서도 근거나 reference를 제시하기 바랍니다.
2. Table 5는 n = 97 valid cases라고 되어 있으나, 본문에서는 solver failure 세 건을 structural statistics에는 포함하고, unconditional structural score는 0.758이고 conditional score는 0.781이라고 서술하였습니다. 반면, 초록과 결론에서는 100개 benchmark에서 평균 structural score가 0.781라고 서술하였습니다. 독자가 혼동이 있을 수 있으므로 명확히 표시 부탁드립니다.
3. 방법론에서는 모든 predicted/reference 사례에 동일한 체적유량 0.05 m^3/s을 적용하고, inlet area에 따라 속도가 달라진다 합니다. Discussion에서는 reference case에 predicted case와 동일한 solver preset과 동일한 inlet velocity boundary condition을 사용했다고 기술합니다. Predicted geometry와 reference geometry의 opening area가 다르다면 다음 두 조건을 동시에 만족할 수 있나요? (동일한 체적유량 / 동일한 inlet velocity) 어느 조건을 일치시켰는지 명확히 설명 바랍니다.
4. opening 크기 대비 격자가 작아보입니다. 0.3m 대비 0.18m면 한개 혹은 두개 존재할 수 있는데, 이렇게 되면 정확한 해석이 불가능해 보입니다. 설명 및 수정/보완 바랍니다.

---

## 2. 코멘트별 분석 · 결정 · 조치

겹치는 코멘트는 묶어서 처리한다.

| ID | 대응 코멘트 | 유형 | 결정 & 조치 | 의존 데이터/코드 |
|----|------------|------|-------------|------------------|
| **C1. 수치 정합성 (CFD 평균)** | R1-1 | 수정 | Fig. 7 overall mean(0.453)이 본문/Table 5(0.477)와 불일치. **원인은 Fig. 7이 구버전(preset 미매칭) 데이터로 생성**되었을 가능성. 최신 매칭 집계로 Fig. 7 재생성 → 0.477로 통일. 산정 기준 각주 추가. | V1 |
| **C2. n 규약 (structural vs CFD)** | R1-1, R2-2 | 수정+명확화 | 규약 확정: **structural = 유효 예측 전체 기준, CFD = 수렴 케이스 기준**. 모든 표·그림·초록·결론에 n 명시. 초록/결론의 "100개에서 0.781" 문구를 규약에 맞게 수정. conditional/unconditional 값을 한 곳(Table 5)에 정리. | V2 |
| **C3. Fig. 13 수치 불일치** | R1-3 | 수정 | 캡션(0.938·1.000 / 0.409·0.458)이 그림 내부·본문(A4-02 0.81/0.35, A4-04 0.92/0.30)과 다름 → **캡션이 잘못된 케이스 값**으로 추정. 실제 A4-02·A4-04 값 확인 후 캡션·본문·그림 3자 통일. | V3 |
| **C4. Table 7 "<5%" 표현** | R1-3 | 명확화 | "less than 5%"는 **agreement score가 아니라 bulk velocity 기준**이었을 것. 기준 물리량 명시하거나, agreement score는 더 엄격한 차분 지표라 변동이 더 클 수 있음을 명기. | V4 |
| **C5. 입구 BC: 체적유량 vs 속도** ⚠️ | R1-4, R2-3 | 명확화(+오류 정정 가능) | **최우선.** 실제 case `0/U` 확인. `CLAUDE.md` 기준으론 preset이 **입구 속도(m/s)**를 고정 → 2.2절 "체적유량 0.05 m³/s" 서술이 코드와 불일치일 수 있음. 확인 후 2.2·5.5 통일. R2-3 답변: opening area가 다르면 둘 다 만족 불가 → **매칭한 것은 (preset을 통한) 입구 속도**임을 명확히. | V5 |
| **C6. Section 6 BC 근거** | R1-2 | 명확화 | 문=inlet/창=outlet이 VLM 판정인지 가정인지 확인 후 명시. + 실제 적용 시 풍량/차압/HVAC 운전정보 필요함을 한 문장 추가. | V8 |
| **C7. Opening Type F1 과대평가** | R2-1 | 반박+보강 | 신규 벤치마크 없이: **① opening type F1 포함/제외 두 점수 병기**(제외 시 3성분 평균 ≈ 0.708 [수치 확인]), ② 응용목적 가중 대안 언급, ③ 벤치마크가 1 inlet/1 outlet 고정임을 한계로 명기. Table 4 근처+Discussion에 반영. | V6 |
| **C8. IoU 임계값 근거** | R2-1 | 보강 | room 0.2 / obstacle 0.1 근거+참고문헌. 실내 3D 검출 표준이 mAP@0.25 IoU (VoteNet/ScanNet 등)임을 들어 room 0.2가 표준에 근접함을, obstacle 0.1은 단일영상 난이도 반영임을 서술. Appendix B 민감도로 결론 강건성 재강조. | V7 (+신규 인용 2건) |
| **C9. Opening 대비 격자** | R2-4 | 한계 명시 | 신규 해석 없이: opening 0.3–1.0 m, 셀 0.18 m → 최소 개구는 1~2셀만 걸침을 인정. **국소 미세화 부재를 명시적 한계**로 5.7절에 추가, 결과가 screening-level임을 재확인. | V9 |

---

## 3. 확인해야 할 데이터/코드 체크리스트 (Claude Code에서 수행)

> ⚠️ 데이터 규약: **preset-matched 데이터만 사용** (`CLAUDE.md`의 "Paper Data Convention" 준수). `cfd_metrics.json`(미매칭)·구버전 집계 사용 금지.

- [ ] **V1 — CFD 평균 확정 & Fig. 7 재생성**
  - 확인: `benchmark/manifests/evaluation_statistics_phase2.json`, `evaluation_aggregate_summary_phase2.json`, `evaluation_statistics_matched.json` → CFD agreement 평균(n=97)과 뷰별 평균.
  - Fig. 7 생성 스크립트(`papers/image_to_cfd/scripts/paper_figures/` 내 `fig_result_cfd_component_breakdown` 관련)가 읽는 입력 파일 확인 → 구버전이면 최신 매칭 데이터로 교체 후 재생성.
  - 산출: overall mean이 0.477과 일치하는지 확인. 그림 PDF/PNG 재생성 → `latex/figures/`.

- [ ] **V2 — structural score n 규약**
  - 확인: `evaluation_statistics_phase2.json`에서 structural 평균이 n=100인지 97인지, 그리고 conditional 0.781 / unconditional 0.758의 정의.
  - 결정 후 문구 규약표 작성(어느 표가 어느 n인지).

- [ ] **V3 — Fig. 13 케이스 값**
  - 확인: `bench_a4_02`, `bench_a4_04` (bird's-eye view) 각 케이스의 structural·CFD (eval dir 내 `cfd_metrics_matched.json` 또는 통계 JSON).
  - 캡션·본문·그림 라벨 3자 정합. 그림 재생성 필요 시 `fig_discuss_structure_cfd_gap` 스크립트 확인.

- [ ] **V4 — Table 7 "<5%" 근거**
  - 확인: grid independence 원자료(속도 등 물리량이 저장된 곳; `src/grid_study/` 산출물 또는 관련 결과 파일). "5% 미만"이 어떤 양 기준이었는지 추적.
  - 없으면 문장을 "bulk velocity magnitude 기준" 등으로 정정하고 agreement score 변동은 별도 설명.

- [ ] **V5 — 입구 경계조건 (최우선)** ⚠️
  - 확인: `cases/phase2_ref_{scene}_preset_{preset}/0/U`, `cases/phase2_pred_bench_{case}_{view}/0/U` 의 inlet 타입 (`flowRateInletVelocity` vs `fixedValue`)과 값.
  - 확인: `scripts/run_indoor_stabilized.py` 및 preset 정의(로그의 robust/ultra_robust/laminar inlet velocity).
  - 결론: 실제 고정한 물리량(속도 or 체적유량)을 확정 → 2.2·5.5 재작성. R2-3 답변 근거 확보.

- [ ] **V6 — opening type F1 포함/제외 점수**
  - 확인: 케이스별 4개 성분 값(room F1, obstacle F1, opening type F1, wall match) 저장 위치 → 제외 시 3성분 평균의 mean±SD 계산(전체 및 뷰별).
  - 산출: 병기용 수치표. (참고: 성분 평균 0.988/0.400/1.000/0.737 → 제외 평균 ≈ 0.708, 단 케이스별 재계산으로 SD 확보)

- [ ] **V7 — IoU 임계값 참고문헌**
  - 후보: VoteNet (Qi et al., 2019), ScanNet (Dai et al., 2017), SUN RGB-D (Song et al., 2015) — 실내 3D object detection이 mAP@0.25 IoU를 표준으로 사용. Zotero(`Zotero_YHJoo.bib`)에 있는지 확인, 없으면 추가.

- [ ] **V8 — Section 6 inlet/outlet 지정 근거**
  - 확인: floorplan 데모 케이스의 VLM 출력(scene JSON)에서 opening type이 VLM 판정인지, 후처리 규칙인지.

- [ ] **V9 — opening 대비 셀 수 (한계 서술용)**
  - 확인: opening 크기 분포(0.3–1.0 m) 대비 셀 0.18 m → 개구를 가로지르는 셀 수 범위 계산(약 1.7–5.6). tetra 국소 셀이 nominal보다 작을 수 있는지 여부.

---

## 4. Response-to-Reviewers 답변서 뼈대 (국문 초안)

> 수치는 `[확인: Vx]`로 표기. Claude Code에서 데이터 확정 후 채운다. 톤: 정중, 각 항목마다 "지적 요지 → 저자 응답 → 원고 반영 위치".

### 공통 서두
> 심사위원님들의 세심한 검토와 건설적인 제언에 깊이 감사드립니다. 지적해 주신 사항을 아래와 같이 반영하였으며, 수정 부분은 원고에 [색상/밑줄] 표기하였습니다.

### 심사위원 #1

**1-1 (CFD 평균 0.477 vs 0.453 / structural·CFD 평가 대상 수)**
- 응답: Fig. 7의 overall mean은 `[확인: V1 — 구버전 데이터 여부]`에 기인한 표기 오류였으며, 최신 preset-matched 집계 기준 CFD agreement 평균은 0.477 `[확인: V1]`로 통일하였습니다. 또한 structural score는 유효 예측 전체(n=`[확인: V2]`), CFD agreement score는 수렴 케이스(n=97) 기준임을 명확히 구분하여 모든 표·그림·본문에 표기하였습니다.
- 반영: Fig. 7, Table 3/5 제목, 2.3절·4.1절.

**1-2 (Section 6 inlet/outlet)**
- 응답: 해당 지정은 `[확인: V8 — VLM 판정 or 가정]`입니다. 실제 현장 적용 시에는 개구의 급기/배기 역할과 함께 풍량·차압 또는 HVAC 운전정보와 같은 경계조건이 추가로 요구됨을 6절에 명시하였습니다.
- 반영: 6절 문단 + 한 문장 추가.

**1-3 (Fig. 13 수치 / Table 7 "<5%")**
- 응답(Fig.13): 캡션 수치가 다른 케이스 값으로 잘못 기재되었습니다. 그림이 나타내는 A4-02·A4-04의 실제 값(structural `[V3]`, CFD `[V3]`)으로 캡션·본문을 일치시켰습니다.
- 응답(Table 7): "less than 5%"는 `[확인: V4 — bulk velocity magnitude]` 기준이었으며, 복합 agreement score는 더 엄격한 차분 지표로서 일부 사례에서 5%를 초과할 수 있음을 명시하였습니다. 문장을 그에 맞게 수정하였습니다.
- 반영: Fig. 13 캡션, 5.3절, 5.7절(Table 7 문단).

**1-4 (체적유량 vs 입구 속도)** ⚠️
- 응답: `[확인: V5]` 결과에 따라, 실제로 일치시킨 조건은 **`[속도 / 체적유량 중 확정]`**이며, solver preset이 이를 결정합니다. 2.2절과 5.5절의 서술을 일치시켰습니다. (opening area가 다른 경우 체적유량과 입구 속도를 동시에 일치시킬 수 없으므로, 본 연구는 `[확정 조건]`을 기준으로 매칭하였음을 명시)
- 반영: 2.2절, 5.5절.

### 심사위원 #2

**2-1 (Opening Type F1 과대평가 / IoU 임계값)**
- 응답(F1): 지적에 동의합니다. 본 벤치마크는 케이스별 1 inlet/1 outlet 구성으로 opening type F1이 항상 1.0이 되어 종합 점수를 상향시킵니다. 이를 (i) opening type F1 **포함 점수(0.781)와 제외 점수(`[확인: V6]` ≈ 0.708)를 병기**하고, (ii) 응용 목적별 가중 점수의 필요성을 논의에 추가하며, (iii) 개구 유형·개수 가변 벤치마크로의 확장을 향후 과제로 명시하였습니다. (해당 벤치마크 구성은 선행연구 `[V7 인용]`와도 부합)
- 응답(IoU): room 0.2 / obstacle 0.1 임계값의 근거를 보강하였습니다. 실내 3D 객체 검출의 표준 평가가 mAP@0.25 IoU (`[V7: VoteNet/ScanNet/SUN RGB-D 인용]`)임을 고려하면 room 0.2는 표준에 근접하며, obstacle 0.1은 단일 영상으로부터의 가구 규모 객체 추정 난이도를 반영한 완화 기준입니다. Appendix B의 임계값 민감도 분석에서 뷰·카테고리 순위가 임계값 전반에 걸쳐 보존됨을 재확인하였습니다.
- 반영: 2.3.1절, Table 4 인접 문단, Discussion, Appendix B.

**2-2 (n=97 vs 100 / 초록·결론)**
- 응답: (C2와 동일) 규약을 통일하고 초록·결론에서 "100개에서 0.781" 표현을 structural n=`[V2]` / CFD n=97로 정정하였습니다.
- 반영: Abstract, 4.1절, Conclusion, Table 5.

**2-3 (체적유량 vs 입구 속도)**
- 응답: (1-4와 동일) opening area가 다르면 두 조건 동시 만족 불가하며, 본 연구가 매칭한 것은 `[V5 확정]`입니다.
- 반영: 2.2절, 5.5절.

**2-4 (opening 대비 격자)**
- 응답: 지적에 동의합니다. nominal 셀 0.18 m 기준으로 최소 개구(0.3 m)는 1~2 셀만 걸쳐 국소 유동을 정확히 해상하지 못합니다. 본 연구는 개구 국소 격자 미세화를 적용하지 않았으며, 이는 개구 근방·근가구 국소 속도의 정확도를 제한합니다. 이를 5.7절에 명시적 한계로 추가하고, 결과가 bulk 유동 패턴 식별을 위한 screening-level 추정임을 재확인하였습니다. (개구 적응 격자는 향후 과제)
- 반영: 5.7절(한계), 결론.

---

## 5. 수정 대상 파일 & 재생성 파이프라인

**소스 오브 트루스 (여기서 수정)**
- `papers/image_to_cfd/docs/paper_draft_v1.md` — 본문(마크다운 원본). 대부분의 문구 수정.
- `papers/image_to_cfd/scripts/paper_figures/` — Fig. 7, Fig. 13 재생성 스크립트.
- (신규 인용) `C:\Vaults\Research\Zotero_YHJoo.bib` — IoU 참고문헌 추가.

**재생성 순서**
1. 데이터 재집계 / 수치 확정 (§3 체크리스트)
2. Fig. 7, Fig. 13 재생성 → `latex/figures/`
3. `paper_draft_v1.md` 수정
4. `bash papers/image_to_cfd/scripts/build_paper.sh` → `paper.pdf` (LaTeX 빌드)
5. **JAAI docx 재생성**: `papers/image_to_cfd/scripts/jaai_docx/` 의 `content.py`(본문 텍스트)·`refs.json`·`omml.json`·`build.py`·`finalize.py` 를 수정 반영 후 재실행.
   - ⚠️ `build.py`/`finalize.py`의 `BASE` 경로가 Cowork 임시폴더(`/sessions/.../outputs`)로 하드코딩됨 → 레포에서 재실행 시 경로 수정 필요. 그림 PNG는 `latex/figures/`에서, 수식은 `eq_src/`의 LaTeX에서 재생성.
   - content.py는 docx 본문의 **독립 사본**이므로, `paper_draft_v1.md` 수정분을 **수동 동기화**해야 함(자동 연결 아님).

---

## 6. Claude Code 세션 핸드오프 노트

**권장 진행 순서**
1. §3 **V5 먼저** (입구 BC) — 답변 방향 전체가 여기 걸림.
2. §3 V1~V4, V6, V8 데이터 확정 → 수치표 채우기 → §4 답변서 `[확인]` 자리 채움.
3. V7(인용), V9(한계 서술) 문구 작성.
4. Fig. 7 / Fig. 13 재생성.
5. `paper_draft_v1.md` 수정 → PDF 빌드 → 수치·캡션 최종 크로스체크(전 문서에서 0.477/0.453/0.781/n 검색).
6. JAAI docx 재생성(§5-5).

**주의**
- 데이터는 preset-matched 만 사용(`CLAUDE.md` 규약).
- 수정 후 **전 문서 grep 검증**: `0.453`, `0.477`, `0.463`, `0.758`, `0.781`, `n = 97`, `n = 100`, `0.05 m` 가 모두 일관되는지.
- 사용된 스킬 활용 권장: `/grid-study`(V4 데이터 확인용), figure 스크립트, `md2latex`, `build_paper`.

**관련 파일**
- 심사평·분석·답변서 뼈대: 본 문서
- docx 생성 파이프라인 사본: `papers/image_to_cfd/scripts/jaai_docx/`
- 현행 제출본: `papers/image_to_cfd/latex/files_for_submission_jaai/` (저자판 + blinded판 + 메타데이터 md)
