# CFD Visual QA — Benchmark Design v1

> Date: 2026-03-24
> Status: Draft for discussion

---

## 1. Design Philosophy

벤치마크의 목적은 **"VLM이 CFD 유동장 시각화에서 무엇을 이해하고 무엇을 이해하지 못하는가"**를 체계적으로 규명하는 것이다. 따라서:

- 단순한 pass/fail이 아니라 **능력의 스펙트럼**을 측정해야 함
- 올바른 유동장만으로는 부족 — **잘못된 유동장에 대한 반응**이 핵심
- 시각화 형식의 영향을 분리해야 함 (같은 데이터, 다른 렌더링)
- 재현 가능해야 함 (모든 데이터를 OpenFOAM으로 생성, 스크립트 공개)

---

## 2. Benchmark Axes (3-축 설계)

### Axis 1: 유동 시나리오 (Flow Scenarios)

각 시나리오는 독립적인 OpenFOAM 케이스로, 명확한 물리적 특성을 가짐.

| ID | 시나리오 | 유동 유형 | 핵심 물리 현상 | 난이도 |
|----|---------|----------|--------------|--------|
| S1 | Heated vertical plate | 자연대류 | 경계층 발달, 상승 플룸 | 기초 |
| S2 | Channel flow (laminar) | 강제대류 | 입구 발달 영역, 포물선 프로파일 | 기초 |
| S3 | Channel flow (turbulent) | 강제대류 | 난류 프로파일, 벽면 전단 | 기초 |
| S4 | Backward-facing step | 강제대류 | 유동 박리, recirculation, 재부착 | 중급 |
| S5 | Lid-driven cavity | 강제대류 | 1차/2차 와류, 코너 와류 | 중급 |
| S6 | Natural convection in cavity | 자연대류 | Ra수 의존 유동 구조, 열 경계층 | 중급 |
| S7 | Mixed convection in channel | 혼합대류 | Ri수에 따른 부력/관성 경쟁 | 고급 |
| S8 | Flow over heated cylinder | 혼합대류 | 와류 방출 + 열 플룸 상호작용 | 고급 |
| S9 | Rayleigh-Bénard convection | 자연대류 | 대류 셀 패턴, 다중 안정 해 | 고급 |
| S10 | Ventilated room (Nielsen-type) | 강제대류 | 제트 부착, 실내 순환 | 고급 |

**시나리오 수: 10개** — 각각 검증된 참조 해(analytical, 실험, 또는 문헌 비교)가 존재하는 것으로 선정.

### Axis 2: 시각화 유형 (Visualization Types)

동일한 CFD 결과를 다양한 방식으로 렌더링.

| ID | 시각화 유형 | 설명 | 정보 밀도 |
|----|-----------|------|----------|
| V1 | Velocity magnitude contour | 컬러맵 + 컬러바 | 중 |
| V2 | Temperature contour | 컬러맵 + 컬러바 | 중 |
| V3 | Pressure contour | 컬러맵 + 컬러바 | 중 |
| V4 | Velocity vector field | 화살표 (방향 + 크기) | 고 |
| V5 | Streamlines | 유선 (경로 + 방향) | 중 |
| V6 | Composite: contour + streamlines | V1/V2 + V5 overlay | 고 |
| V7 | Line profile plot | 특정 위치의 1D 프로파일 그래프 | 저 |

**시각화 유형 수: 7개** — 모든 시나리오에 모든 유형이 적용되지는 않음 (해당 물리량이 존재하는 경우만).

### Axis 3: 질문 난이도 (Question Levels)

| Level | 이름 | 설명 | 예시 |
|-------|------|------|------|
| L1 | Visual Reading | 시각적 정보 직접 읽기 | "최대 속도 영역의 색상은?" |
| L2 | Qualitative Pattern | 정성적 패턴 식별 | "recirculation zone이 있는가?" |
| L3 | Physical Reasoning | 물리 법칙 기반 추론 | "이 유동장에 물리적 오류가 있는가?" |
| L4 | Quantitative Estimation | 정량적 판단 | "경계층 두께가 약 얼마인가?" |
| L5 | Comparative Judgment | 두 이미지 비교 판단 | "A와 B 중 물리적으로 타당한 것은?" |

---

## 3. 오류 유동장 생성 (Error Flow Field Generation)

벤치마크의 핵심 차별점. 각 시나리오별로 올바른 유동장 1개 + 오류 유동장 N개를 생성.

### 3.1 오류 유형 분류

| ID | 오류 유형 | 생성 방법 | 물리적 위반 |
|----|----------|----------|------------|
| E1 | 미수렴 (Under-converged) | 솔버 early-stop (잔차 > 1e-2) | 질량/에너지 보존 불완전 |
| E2 | 경계조건 오류 | inlet↔outlet 교환, 벽면 조건 변경 | 유입/유출 방향 역전 |
| E3 | 물성치 오류 | 점성계수 10×, 밀도 반전 등 | 비현실적 Re/Ra |
| E4 | 난류 모델 부적합 | 층류에 난류 모델 (또는 반대) | 과도한/부족한 혼합 |
| E5 | 메시 부적절 | 경계층 해상도 극도로 부족 | 벽면 전단 왜곡, 비물리적 온도 분포 |
| E6 | 인위적 변형 | 수렴 해의 특정 영역을 코드로 조작 | 국소적 비물리적 패턴 |
| E7 | 대칭 위반 | 대칭 문제에서 비대칭 결과 강제 | 보존법칙 위반 |
| E8 | 부력 방향 반전 | 중력 방향 변경 | 자연대류 방향 역전 |

### 3.2 오류 심각도 등급

| 등급 | 설명 | 전문가 기준 |
|------|------|-----------|
| Severe | 한눈에 비물리적 | 학부생도 감지 가능 |
| Moderate | 주의 깊게 보면 감지 | 석사급 CFD 경험 필요 |
| Subtle | 전문가만 감지 가능 | 해당 유동의 물리를 깊이 이해해야 감지 |

### 3.3 시나리오 × 오류 매트릭스 (예시)

| 시나리오 | 적용 가능한 오류 유형 |
|---------|-------------------|
| S1 Heated plate | E1, E2, E3, E5, E6, E8 |
| S4 Backward step | E1, E2, E4, E5, E6 |
| S6 Cavity convection | E1, E2, E3, E5, E6, E7, E8 |
| S9 Rayleigh-Bénard | E1, E3, E5, E6, E7, E8 |

---

## 4. 라벨 구조 (Annotation Schema)

각 이미지에 대한 다층 라벨.

### 4.1 이미지 메타데이터
```json
{
  "image_id": "S4_V1_correct_001",
  "scenario": "S4",
  "scenario_name": "backward_facing_step",
  "visualization": "V1",
  "visualization_name": "velocity_magnitude_contour",
  "is_correct": true,
  "error_type": null,
  "error_severity": null,
  "solver": "simpleFoam",
  "mesh_cells": 45000,
  "turbulence_model": "kOmegaSST",
  "Re": 36000
}
```

### 4.2 오류 이미지 라벨
```json
{
  "image_id": "S4_V1_error_E2_001",
  "scenario": "S4",
  "visualization": "V1",
  "is_correct": false,
  "error_type": "E2",
  "error_type_name": "boundary_condition_error",
  "error_severity": "severe",
  "error_description": "Inlet and outlet boundary conditions swapped",
  "error_region": "entire domain — flow direction reversed",
  "physical_violation": "Mass flow direction inconsistent with geometry",
  "expert_explanation": "Flow enters from the outlet face and exits through the inlet, creating reverse recirculation pattern incompatible with the step geometry."
}
```

### 4.3 질문-답변 라벨
```json
{
  "image_id": "S4_V1_error_E2_001",
  "questions": [
    {
      "level": "L1",
      "question": "What is the approximate maximum velocity magnitude shown in this image?",
      "answer": "Approximately 1.2 m/s based on the colorbar",
      "answer_type": "free_text"
    },
    {
      "level": "L2",
      "question": "Is there a recirculation zone downstream of the step?",
      "answer": false,
      "answer_type": "boolean",
      "explanation": "No recirculation because flow direction is reversed"
    },
    {
      "level": "L3",
      "question": "Is this flow field physically plausible for a backward-facing step?",
      "answer": false,
      "answer_type": "boolean",
      "explanation": "Flow enters from the right (outlet face) which is physically impossible for the given geometry and boundary setup"
    },
    {
      "level": "L5",
      "question": "Compare with the reference image. Which aspects differ and are they physically justified?",
      "answer": "Flow direction is entirely reversed. This is not physically justified.",
      "answer_type": "free_text",
      "reference_image": "S4_V1_correct_001"
    }
  ]
}
```

---

## 5. 데이터셋 규모 추정

### 5.1 이미지 수

| 항목 | 수량 | 근거 |
|------|------|------|
| 시나리오 | 10 | Axis 1 |
| 시나리오당 올바른 유동장 | 1-2 | 파라미터 변형 (Re/Ra) |
| 시나리오당 오류 유동장 | 3-5 | 오류 유형별 1개 |
| 시나리오당 이미지 (올바른) | ~7-10 | 시각화 유형 × 올바른 유동장 수 |
| 시나리오당 이미지 (오류) | ~15-25 | 시각화 유형 × 오류 유동장 수 (일부만) |
| **총 이미지 수** | **~250-350** | |

### 5.2 질문 수

| 항목 | 수량 |
|------|------|
| 이미지당 질문 | 3-5 (level 혼합) |
| **총 질문 수** | **~1,000-1,500** |

### 5.3 CFD 케이스 수

| 항목 | 수량 |
|------|------|
| 올바른 케이스 | 15-20 (시나리오 × 파라미터) |
| 오류 케이스 | 30-50 (시나리오 × 오류 유형) |
| **총 CFD 실행** | **~50-70** |

---

## 6. 평가 프로토콜

### 6.1 평가 대상 VLM (확정)
- **Claude Opus 4.6** — 블라인드 서브에이전트로 평가 ✅ 완료
- **Gemini 3.1** — `/ask gemini` CLI로 평가
- **GPT-5.4 (Codex)** — `/ask codex` CLI로 평가
- (추가 가능: Qwen2.5-VL, LLaVA 등 오픈소스)

### 6.2 평가 메트릭

**L1 (Visual Reading):**
- 정확도 (exact match 또는 허용 오차 내)

**L2 (Qualitative Pattern):**
- 분류 정확도 (boolean/categorical)
- F1 score

**L3 (Physical Reasoning):**
- 이진 분류: correct/incorrect 판별 정확도
- Precision / Recall / F1
  - Precision: 비물리적이라 판단한 것 중 실제 비물리적인 비율
  - Recall: 실제 비물리적인 것 중 감지한 비율
- 오류 유형 식별 정확도 (multi-label)

**L4 (Quantitative Estimation):**
- 상대 오차 (estimated vs ground truth)
- 허용 오차 내 적중률

**L5 (Comparative Judgment):**
- 쌍별 비교 정확도
- 추론 근거의 물리적 타당성 (전문가 평가)

### 6.3 결과 보고 구조

```
Overall Accuracy
├── By Question Level (L1-L5)
├── By Flow Scenario (S1-S10)
├── By Visualization Type (V1-V7)
├── By Error Type (E1-E8)
├── By Error Severity (Severe/Moderate/Subtle)
└── Cross-analysis
    ├── Level × Visualization heatmap
    ├── Scenario × Level heatmap
    └── Error Detection: Precision-Recall curves
```

---

## 7. 시각화 생성 파이프라인

### 7.1 표준화 요구사항
- 모든 이미지는 동일한 렌더링 파이프라인으로 생성
- 해상도: 1200×900 px (또는 유사 비율)
- 컬러맵: 시나리오별 고정 (viridis for velocity, coolwarm for temperature 등)
- 컬러바: 항상 포함, 동일 위치/크기
- 축/라벨: 포함 (실제 CFD 후처리 도구 출력과 유사하게)
- 배경: 흰색
- 렌더링 도구: matplotlib + PyVista (일관성)

### 7.2 파이프라인 구조
```
OpenFOAM case → foamToVTK → PyVista load →
  → matplotlib contour/vector/streamline rendering → PNG export
```

---

## 8. 작업 순서 (제안)

| Phase | 작업 | 산출물 |
|-------|------|--------|
| P1 | 파일럿: S4 (backward step) 1개 시나리오 완성 | 올바른 + 오류 케이스 3-4개, 이미지 ~20장, 질문 ~60개 |
| P2 | 파일럿 평가: GPT-4o + Gemini로 예비 실험 | 초기 결과, 질문 난이도 보정 |
| P3 | 전체 시나리오 확장 (S1-S10) | 50-70 CFD 케이스 |
| P4 | 전체 시각화 + 라벨링 | ~300 이미지, ~1,200 QA |
| P5 | 전체 VLM 평가 | 결과 테이블 |
| P6 | 논문 작성 | 초안 |

---

## 9. Open Questions (논의 필요)

1. **시나리오 수**: 10개가 적절한가? 더 줄이거나 늘려야 하나?
2. **2D vs 3D**: 2D 시뮬레이션만으로 충분한가, 3D→2D slice도 포함할 것인가?
3. **컬러맵 표준화 vs 다양화**: 동일 컬러맵으로 통일 vs 다양한 컬러맵도 변수로 포함?
4. **질문 형식**: MCQ vs open-ended vs boolean 혼합 비율?
5. **비교 질문(L5)**: 두 이미지를 동시에 보여주는 것이 기술적으로 가능한 모든 VLM에서 지원되는가?
6. **라벨링 워크로드**: ~300 이미지 × 4 질문 ≈ 1,200 QA의 전문가 라벨링이 현실적인가?
7. **타겟 저널**: Building and Environment? Computers & Fluids? 또는 AI 쪽 (NeurIPS Datasets)?
