# CFD Visual QA — VLM Benchmark for Flow Field Understanding

## Project Overview
CFD 유동장 시각화에 대한 VLM(Vision-Language Model)의 이해 능력을 체계적으로 평가하는 벤치마크 데이터셋 및 연구 프로젝트.

**핵심 가설:** 현재 프론티어 VLM은 CFD 유동장 시각화의 물리적 타당성을 체계적으로 판단하지 못하며, 이 gap을 정량적으로 규명하는 것 자체가 연구 기여가 된다.

**연구 동기:** Image-to-CFD 자동화 파이프라인(본 레포 메인 프로젝트)에서 결과 검증이 병목. VLM이 유동장의 물리적 타당성을 판단할 수 있다면 자동화 루프가 완성되지만, 현재 어디까지 가능한지 체계적으로 평가된 바 없음.

## Project Structure
```
papers/cfd_visual_qa/
├── docs/               # 연구 노트, 선행연구, 논문 초안
├── scripts/            # 벤치마크 데이터셋 생성, VLM 평가 스크립트
├── benchmark/          # 데이터셋 (유동장 이미지 + 라벨)
├── results/            # 실험 결과
├── CLAUDE.md           # 이 파일
└── README.md
```

## Research Scope

### 도메인
- 열유동 (자연대류, 강제대류, 혼합대류)
- 2D 시각화 + 3D→2D 렌더링
- OpenFOAM 기반 시뮬레이션 데이터

### 벤치마크 데이터셋 구성
1. **올바른 유동장** — 물리적으로 타당한 CFD 결과의 시각화
2. **잘못된 유동장** — 의도적으로 생성한 비물리적 결과 시각화
3. **전문가 라벨** — 올바름/오류 + 오류 유형/심각도/위치 다층 라벨

### 잘못된 유동장 생성 전략
- 수렴 실패/중간 상태 (early-stop, 발산 직전)
- 경계조건 오류 (inlet/outlet 뒤바뀜, 단열↔등온 등)
- 메시 부적절 (경계층 해상도 부족)
- 난류 모델 부적합 (층류에 난류 모델, 또는 반대)
- 인위적 조작 (recirculation zone 제거/이동, 경계층 뒤집기)

### 평가 대상 현상
| 현상 | 올바른 패턴 | 오류 패턴 예시 |
|------|-----------|--------------|
| 자연대류 (heated wall) | 벽면 상승류, 경계층 발달 | 하강류, 경계층 없음 |
| 강제대류 (channel) | 입구 발달 영역, 포물선 프로파일 | 균일 속도, 역류 |
| 혼합대류 | Richardson수에 따른 패턴 전이 | 부력 방향 반전 |
| Rayleigh-Benard | 셀 패턴, 대칭성 | 비대칭, 온도 역전 |
| 후방계단 (backward step) | 재부착 길이, recirculation | 재부착 없음, 위치 오류 |

### VLM 질문 난이도 체계
- **Level 1** (시각 읽기): "고온 벽면은 어디인가?"
- **Level 2** (패턴 인식): "recirculation zone이 존재하는가?"
- **Level 3** (물리 판단): "물리적으로 비합리적인 부분이 있는가?"
- **Level 4** (정량 추론): "경계층 두께가 이론값과 비교하여 합리적인가?"

### 시각화 유형
- Velocity contour (magnitude)
- Velocity vector field
- Streamlines
- Pressure contour
- Temperature contour
- 다중 물리량 합성 (overlay)

## Key Prior Art

| 논문 | 도메인 | CFD 시각화? | VLM 평가? | 연도 |
|------|--------|-----------|----------|------|
| MDPI Technologies (Anon.) | 실내 CFD (MR) | YES | Qwen2.5-VL 파인튜닝 | 2026 |
| Kashefi | 유체역학 생성 AI | YES (생성) | 생성 실패 분석 | 2024 |
| CFDLLMBench | CFD 지식/코드 | No | 텍스트 전용 | 2025 |
| ClimateIQA | 기상 히트맵 | 유사 | 파인튜닝 | 2024 |
| PhysBench | 물리 세계 | No CFD | 75 VLM | 2025 |
| SciFIBench | 과학 도표 | 간접 | 28 LMM | 2024 |

### 확인된 Gap
1. CFD 유동장 이미지 전용 VLM 벤치마크 없음
2. 파인튜닝 연구는 MDPI 1편뿐 (소규모, 단일 도메인)
3. 컬러맵 기반 정량적 추론 미해결
4. 다중 렌더링 스타일 간 일관성 평가 없음
5. "물리적으로 올바른가" 판단 메트릭 없음
6. CFD 결과 품질 (수치 아티팩트) 감지 평가 없음
7. 공개된 CFD 유동장 이미지 + annotation 데이터셋 없음

## Convention
- 본 프로젝트의 모든 파일은 `papers/cfd_visual_qa/` 하위에서 관리
- 메인 프로젝트(Image-to-CFD 논문)와 파일 충돌 금지
- CFD 케이스 생성 시 `../../cases/` 참조 가능하나 결과물은 본 폴더로 복사
- 코드 스타일: PEP 8
- 문서: 한국어 + 영문 혼용 허용 (논문 초안은 영문)
