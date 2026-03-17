# Image-to-CFD 자동화 프레임워크 논문 계획

> 작성일: 2026-03-17

## 1. 핵심 메시지 (Thesis)

VLM(Vision-Language Model)의 공간 지각 능력을 활용하면, 2D 이미지 한 장으로부터
시뮬레이션 가능한 3D 실내 형상을 자동 추상화하여 CFD 해석까지 도달할 수 있다.

## 2. 논문 구성 (Draft Outline)

| 섹션 | 내용 |
|------|------|
| Introduction | 실내 환경 CFD의 필요성, 형상 모델링 병목, VLM 활용 가능성 |
| Related Work | image-to-3D 재구성, VLM 공간 지각, CFD 자동화 선행 연구 |
| Methodology | 프레임워크 아키텍처 (이미지 → scene JSON → mesh → CFD) |
| Benchmark Dataset | rule-based 형상 생성, 난이도 체계, 뷰 타입 정의 |
| Experiments | 난이도 × 뷰 타입 매트릭스 평가, 정량 지표 |
| Results & Discussion | 성공률, 물리량 오차, 실패 모드 분석 |
| Wild Image Demo | 인터넷/실제 사진 적용 사례 (정성적) |
| Conclusion | 기여 요약, 한계, 향후 과제 |

## 3. Contribution

1. **프레임워크**: 2D 이미지 → VLM → 3D 추상화 → CFD 자동화 파이프라인
2. **벤치마크 데이터셋**: 난이도별 실내 형상 + 다중 뷰 렌더링 + reference CFD (공개 가능)
3. **정량 평가**: 뷰 타입별/난이도별 성공률 및 물리량 오차 분석
4. **한계 분석**: VLM 공간 지각의 현재 한계와 실패 조건 규명

## 4. 데이터셋 설계

### 4.1 전체 흐름

```
[1] Rule-based 형상 생성 (scene JSON)
     │
[2] Gmsh 메싱 → OpenFOAM Reference CFD (정답 확보)
     │
[3] 3D → 2D 렌더링 (뷰 타입별 입력 이미지 생성)
     │
[4] 프레임워크 적용 (2D 이미지만 입력, Gemini scene generation)
     │
[5] 정답 대비 평가 (형상 + CFD 결과)
```

- Reference 경로에는 Gemini가 **개입하지 않음** → 순환 논증 없음
- 평가 대상인 Gemini는 [4]에서만 사용

### 4.2 형상 난이도 체계 (2×2 매트릭스)

두 개의 독립 변수로 구성:
- **공간 형태**: 직육면체(Rectangular) vs L자형(L-shaped)
- **장애물 복잡도**: 단순(0~1개) vs 복잡(2~3개+)

|  | Rectangular | L-shaped |
|---|---|---|
| **단순 장애물** (0~1개) | **A1** — 빈 방, 단일 가구 | **A3** — L자 빈 공간, 최소 가구 |
| **복잡 장애물** (2~3개+) | **A2** — 가구 배치 사무실 | **A4** — L자 + 밀집 가구 |

각 카테고리 3~5개 형상, 총 **12~20개 형상**

**L자형 구현**: 두 개의 직육면체 블록을 접합하는 방식 (칸막이 장애물 아님)
- scene JSON 스키마를 확장하여 composite room 지원
- 기존: `"room": { "size": { "Lx", "Ly", "Lz" } }` (단일 직육면체)
- 확장안: `"room": { "blocks": [ { "origin": {...}, "size": {...} }, ... ] }` (복수 블록 접합)
- 단일 블록이면 기존과 호환, 복수 블록이면 L자형 등 표현 가능
- Gmsh `.geo` 생성 시 블록들의 Boolean union으로 도메인 구성
- `scene_to_gmsh.py`, `validate_indoor_scene.py` 등 기존 파이프라인 코드도 확장 필요

**제외 사항**:
- T자형: inlet/outlet 배치에 따라 유동 고립 영역 발생 위험
- A4는 프레임워크 한계를 드러내는 용도로도 활용

### 4.3 뷰 타입 (독립변수)

| 뷰 타입 | 정보 특성 | 기대 난이도 |
|---------|----------|------------|
| Photorealistic rendering | 질감, 조명, 깊이감 포함 | 낮음 (정보 풍부) |
| Wireframe perspective | 구조 명확, 질감 없음 | 중간 |
| Bird's eye view (조감도) | 전체 레이아웃, 깊이 약함 | 중간 |
| Floor plan (평면도) | 높이 정보 없음 | 높음 |
| Section view (단면도) | 부분적 3D 정보 | 높음 |

- 동일 형상에서 5종 뷰를 렌더링 → **입력 이미지 유형의 영향을 분리 평가**

### 4.4 실험 매트릭스

```
(형상 난이도 4단계) × (뷰 타입 5종) × (형상 수 ~15개)
= 약 60~75 runs
```

### 4.5 경계 조건 통일

Reference CFD와 프레임워크 적용 CFD의 BC를 동일하게 고정:
- **유량 기준 통일**: inlet은 유량(volumetric flow rate)으로 지정하여 개구부 크기와 무관하게 동일 조건 보장
- inlet 온도, outlet 조건(pressure outlet) 등도 표준화
- rule-based 생성 시 BC 정보를 scene JSON에 포함하거나 별도 config로 관리

### 4.6 개구부 설계 원칙

- inlet/outlet 크기의 랜덤화 범위를 좁게 제한 (과도한 변수 다양화 방지)
- 일반적인 실내 공간 환기구 크기 기준: du, dv 약 0.3~0.6m
- diffuser 등 복잡한 급기 형태는 가정에서 제외 (단순 개구부로 통일)

## 5. Rule-Based 형상 생성기 설계

### 5.1 기반 스키마

기존 `indoor_cfd_scene_v1` JSON 포맷을 그대로 활용:

```json
{
  "schema_version": "indoor_cfd_scene_v1",
  "units": "m",
  "coordinate_system": { "origin": "room_min_corner", ... },
  "room": { "size": { "Lx": ..., "Ly": ..., "Lz": ... } },
  "obstacles": [ { "id": "...", "type": "box", "min": {...}, "size": {...} } ],
  "openings": [ { "id": "...", "type": "inlet/outlet", "wall": "...", ... } ],
  "meta": { ... }
}
```

### 5.2 생성 전략

| 파라미터 | 범위 / 규칙 |
|---------|-------------|
| 방 크기 (Lx, Ly, Lz) | 3~8m × 3~6m × 2.4~3.5m |
| 장애물 수 | L1: 0~1, L2: 2~3, L3: 1~2 + 칸막이, L4: 3+ 밀집 |
| 장애물 크기 | dx ≥ 0.6, dy ≥ 0.8, dz ≥ 1.5 (solver-friendly 최소값) |
| 벽면 이격 | ≥ 0.4m |
| 장애물 간 간격 | 겹침 불가, clearance 확보 |
| inlet/outlet | 정확히 1개씩, 서로 다른 벽면 |
| 개구부 크기 | du, dv 0.3~1.0m 범위 |
| seed | 공개하여 재현성 보장 |

### 5.3 구현 위치

```
scripts/generate_benchmark_scenes.py   # rule-based 생성기
benchmark/                             # 생성된 데이터셋
├── scenes/                            # scene JSON 파일들
├── reference_cfd/                     # reference CFD 결과
├── renderings/                        # 뷰 타입별 2D 이미지
└── evaluations/                       # 프레임워크 적용 결과 및 평가
```

## 6. 평가 지표

| 레벨 | 지표 | 설명 |
|------|------|------|
| 형상 | 장애물 수 일치, 위치/크기 오차 (RMSE) | 추론된 형상 vs 정답 형상 |
| 메시 | 메싱 성공률 | 유효한 메시가 생성되는가 |
| 솔버 | 수렴 성공률 | 발산 없이 해가 나오는가 |
| 물리량 | 온도/유속 분포 오차 (RMSE, relative error) | CFD 결과가 정답과 얼마나 가까운가 |

- 파이프라인 단계별로 **어디에서 실패하는가**를 추적하는 것이 핵심

## 7. Data Leakage 방지 체크리스트

- [ ] Reference 형상 생성에 Gemini 미사용
- [ ] 2D 렌더링에 치수 annotation, 깊이 맵 등 부가 정보 미포함
- [ ] 프레임워크 적용 시 프롬프트에 정답 형상 힌트 미포함
- [ ] 공개 건축물/도면 미사용 (Gemini 학습 데이터 오염 방지)
- [ ] VLM 모델 비교 시 정답 생성에 사용된 모델은 비교 대상에서 제외

## 8. 향후 확장 가능성

- **VLM 모델 비교**: Gemini vs GPT-4o vs Claude 등 (rule-based 정답이므로 모든 모델 평가 가능)
- **텍스트 baseline**: 이미지 없이 텍스트 프롬프트만으로 생성 → 이미지 입력의 부가가치 정량화
- **접근 A (wild images)**: 인터넷 수집 실제 사진으로 정성적 실용성 시연

## 9. 타겟 저널

- 결과 우수 시: 열유체 계통 JCR 상위 10% 저널
- 결과 보통 시: Q3 수준 저널

## 10. 미결 사항

- [ ] 2D 렌더링 도구 선정 (Blender Python API, VTK, matplotlib 3D 등)
- [ ] L자형 두 블록 접합 방식의 scene JSON 스키마 확장 설계 (composite room)
- [ ] 경계 조건 표준 세트 구체화 (inlet 유량 값, 온도, outlet 압력 조건)
- [ ] 형상 난이도 4단계 세부 정의 확정 (아래 논의 필요)
