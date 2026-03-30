# Literature Review — VLM Evaluation on CFD Flow Field Visualizations

> Date: 2026-03-24
> Status: Initial survey complete

---

## A. VLM Benchmarks for Scientific / Engineering Visuals

### 1. MMMU (Yue et al., CVPR 2024 Oral)
- 11,500 college-level questions, 30 subjects, 183 subfields
- Engineering/Science에서 모델 성능 최저
- **Gap:** CFD 유동장 이미지 없음. 공학 범위가 generic.
- Source: https://mmmu-benchmark.github.io/

### 2. SciFIBench (Roberts et al., NeurIPS 2024)
- 2,000 questions, 8 scientific figure categories (arXiv 논문 기반)
- 28개 LMM 평가. "과학 도표 이해 능력이 제대로 특성화되지 않았다"
- **Gap:** CFD 도메인 특화 아님. 범분야 과학 도표.
- Source: https://arxiv.org/abs/2405.08807

### 3. PhysBench (Chow et al., ICLR 2025 Oral)
- 10,002 entries, 4 domains, 19 subclasses. 75 VLM 테스트.
- GPT-4o 49.49% vs 인간 95.87%
- "physics-based dynamics"에 유체역학 포함하나 일상 영상 (물 쏟기 등)
- PhysAgent 제안 (+18.4% 향상)
- **Gap:** CFD 솔버 출력 이미지 아님. 자연 유체 영상만.
- Source: https://arxiv.org/abs/2501.16411

### 4. Visualization Literacy Test (Pandey & Ottley, EuroVis 2025)
- GPT-4, Claude, Gemini, Llama 대상 VLAT/CALVI 테스트
- 라인 차트 76-96%, 버블 차트 18-61%, 이상 탐지 25-30%
- **Insight:** CFD 컬러맵은 "다중 인코딩" 범주 → VLM 취약 영역
- **Gap:** 데이터 차트만. 컨투어/벡터/스트림라인 미포함.
- Source: https://arxiv.org/abs/2503.16632

### 5. ClimateIQA (Chen et al., 2024/2026)
- 26,280 기상 히트맵 + 762,120 instruction 샘플
- GPT-4o, Qwen-VL, LLaVA 1.6이 "정밀 색상 식별 및 공간 위치 파악에 실패"
- Climate-Zoo 파인튜닝 VLM 개발 + SPOT 알고리즘
- **Insight:** CFD 컨투어와 가장 유사한 도메인. 기상 히트맵 ≈ 온도/속도 컨투어.
- **Gap:** 기상 데이터만. 속도 벡터, 스트림라인, 압력 컨투어 미포함.
- Source: https://arxiv.org/abs/2406.09838

### 6. Multimodal ArXiv (ACL 2024)
- 6.4M 과학 figure + caption (ArXivCap)
- 유체역학 figure 부수적으로 포함
- **Gap:** CFD 특화 라벨링/큐레이션 없음.
- Source: https://aclanthology.org/2024.acl-long.775/

---

## B. AI/ML for CFD Result Validation

### 7. CFDLLMBench (Somasekharan et al., arXiv 2509.20374, NeurIPS 2025 제출)
- CFDQuery (90 MCQ) + CFDCodeBench (24 Python) + FoamBench (126 OpenFOAM)
- 개념 MCQ 92%, 코드 14%, OpenFOAM 실행 25-34%
- **Gap:** 완전히 텍스트/코드 기반. 시각적 이해 평가 제로.
- Source: https://arxiv.org/abs/2509.20374

### 8. Fine-tuning LLM for CFD (arXiv 2504.09602, 2025)
- Qwen2.5-7B를 자연어 → OpenFOAM 설정으로 파인튜닝
- 88.7% solution accuracy, 82.6% first-attempt success
- **Gap:** 설정 자동화. 결과 해석/시각화 이해 없음.
- Source: https://arxiv.org/abs/2504.09602

### 9. ★ VLM for Indoor CFD Mixed Reality (MDPI Technologies, 2026)
- **가장 직접적인 선행 연구**
- 실내 CFD 결과를 MR 이미지에 중첩, 도메인 특화 QA annotation
- Qwen2.5-VL 파인튜닝: 베이스라인 <30% → 파인튜닝 >60%
- MMBench 일반화 확인
- **Gap:** 소규모 데이터셋, 단일 도메인(실내 온도/속도), 단일 아키텍처. 압력/스트림라인/벡터 미포함. MR 프레이밍이 표준 CFD 후처리와 다름.
- Source: https://www.mdpi.com/2227-7080/14/3/157

### 10. NeurIPS ML4CFD Competition (arXiv 2506.08516, 2025)
- 240+ 팀, 2D 에어포일 서로게이트 모델링
- 정확도, 물리 충실도, 효율, OOD 일반화 평가
- **Gap:** 수치 출력 평가. 언어/시각 이해 없음.
- Source: https://arxiv.org/html/2506.08516

---

## C. Physics-Aware / Domain-Specific VLMs

### 11. Synthetic Vision (arXiv 2412.08619, 2024)
- 물리 시뮬레이터로 VLM 학습 데이터 생성
- Physics Context Builders (PCBs) — 물리 chain-of-thought
- **Gap:** 로보틱스/조작 장면. 유체 시뮬레이션 없음.
- Source: https://arxiv.org/html/2412.08619v1

### 12. Physics-Guided VLM World Models (Microsoft Research, 2024-2025)
- VLM → 물리 제약 장면 변환 (동적 제어)
- **Gap:** 로보틱스. 과학 시각화 없음.

### 13. Physics-Informed Computer Vision Survey (Banerjee et al., ACM Computing Surveys, 2024)
- 250+ 논문 리뷰. 난류 초해상도, PINN 기반 재구성 등
- **Gap:** VLM 이전 프레이밍. 언어 기반 이해/질의 없음.
- Source: https://arxiv.org/abs/2305.18035

### 14. NVIDIA Apollo (2025)
- 과학 시뮬레이션용 foundation model (neural operators, transformers, diffusion)
- **Gap:** 시뮬레이션 가속. 시각적 해석 아님.

---

## D. Flow Visualization Understanding — Classical CV/ML

### 15. Deep Learning for Flow Visualization Survey (Liu et al., Advances in Aerodynamics, 2022)
- VortexNet, Shock-Net 등 CNN 기반 유동 특징 탐지
- **Gap:** 분류/탐지 전용. 언어/의미 이해 없음.
- Source: https://link.springer.com/article/10.1186/s42774-022-00113-1

### 16. 3D Flow Field Segmentation/Classification (arXiv 2305.11884, 2023)
- 3D 와류 영역 분할, Re수 추정
- **Gap:** 이진/범주 출력만. 언어 인터페이스 없음.

### 17. VortexViz (arXiv 2404.01352, 2024)
- CNN + FCN 이중 모달리티 와류 경계 탐지
- **Gap:** 단일 태스크. 범용성 없음.

### 18. ★ "A Misleading Gallery of Fluid Motion by Generative AI" (Kashefi, 2024)
- Midjourney, DALL-E, Gemini 등에 고전 유체 현상 생성 요청
- von Karman 와류, Kelvin-Helmholtz, 초음속 충격파 등 모두 물리적 오류
- 이미지→텍스트 묘사도 실패
- **근본 원인:** 학술 저널 이미지 저작권으로 학습 데이터 부족
- **Insight:** CFD 유동장에 대한 VLM 능력의 강력한 부정적 결과(negative result)
- Source: https://arxiv.org/abs/2405.15406

---

## E. Gap Summary

| # | Gap | 관련 문헌 |
|---|-----|----------|
| 1 | CFD 유동장 이미지 전용 VLM 벤치마크 없음 | All — 아무도 안 함 |
| 2 | 파인튜닝 연구 MDPI 1편뿐 (소규모) | #9 |
| 3 | 컬러맵 기반 정량적 추론 미해결 | #4, #5 |
| 4 | 다중 렌더링 스타일 간 일관성 평가 없음 | — |
| 5 | 물리적 타당성 판단 메트릭 없음 | #18 |
| 6 | CFD 결과 품질(수치 아티팩트) 감지 평가 없음 | — |
| 7 | 공개 CFD 유동장 이미지 + annotation 데이터셋 없음 | #18 (원인 지적) |

---

## F. Positioning

본 연구는 Gap 1, 5, 6, 7을 직접 해결:
- **Gap 1 → 벤치마크 구축:** 열유동 도메인 CFD 시각화 + 다층 라벨 데이터셋
- **Gap 5 → 평가 메트릭:** 물리적 타당성 판단 정확도 체계 (Level 1-4)
- **Gap 6 → 품질 감지:** 의도적 오류 유동장으로 VLM 오류 탐지 능력 평가
- **Gap 7 → 공개 데이터셋:** 연구 커뮤니티에 공개

ClimateIQA(기상)와 CFDLLMBench(텍스트 CFD) 사이의 교차점을 채우는 위치.
