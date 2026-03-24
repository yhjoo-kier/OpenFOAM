# CFD Visual QA — Daily Research Log

---

## Day 0 (2026-03-24) — Setup + Pilot Start

### 완료
- 프로젝트 폴더 구조 생성
- CLAUDE.md 작성 (연구 가설, 벤치마크 설계, 선행연구 gap)
- 선행연구 조사 (18편 분석, 7개 gap 확인)
- 벤치마크 설계 v1 작성 (3-축: 시나리오 × 시각화 × 질문 난이도)
- 연구 directive 문서 작성
- 연구 상태 관리 시스템 구축 (research_state.json)
- 라벨링 큐 시스템 구축
- Cron 작업 등록 (메인 3 + 워치독 4 = 7개)
- 자율연구 프레임워크 문서화 (docs/26-03-24_autonomous_research_framework.md)
- **S4 파일럿 CFD 케이스 6개 완료** (correct 2 + error 4)
- **시각화 이미지 30개 생성** (6 케이스 × 5 viz types)
- **QA 질문 35개 생성** (20개 항목, L1-L5 혼합)
- **Day 1 라벨링 배치 준비 완료** (labeling_queue/pending/day_01/)

### 이슈
- Cron 작업은 세션 기반, 7일 후 자동 만료 → Day 7에 재등록 필요
- foamToVTK 타임아웃: 120s 부족 → 600s로 수정
- 이전 세션 hang된 Docker 컨테이너 정리 (checkMesh 23h, rm -rf 26h)
- pyvista 미설치 → pip install 완료

### 자체 QC 관찰
- S4_correct_turb: recirculation zone 명확, 재부착 ~6h, 물리적 타당
- S4_correct_lam: 층류 특성 (긴 recirculation), 물리적 타당
- S4_E1_underconverged: x≈10m에 비물리적 고속 영역 → severe 오류
- S4_E2_bc_swap: 유동 패턴 미묘하게 다름 → subtle~moderate 오류
- S4_E5_coarse_mesh: correct_turb와 비슷하나 BL 해상도 부족 → subtle 오류
- S4_E6_modified: 아직 인위적 변형 미적용 (수렴 해 그대로) → Day 1에서 처리

### 내일 계획 (Day 1)
- 09:02 라벨링 요청 텔레그램 발송 (20개 항목)
- S4_E6 인위적 변형 적용 (recirculation zone 제거 후 재렌더링)
- 라벨 수거 후 benchmark/labels/ 통합 시작
- S1 또는 S5 다음 시나리오 케이스 설계 시작

---

## Day 1 (2026-03-25) — Full Expansion

### 완료 (자율 탐색 cron + 메인 세션)
- S2 channel flow: 6 CFD + 30 images + 28 QA (자율 탐색 00:37)
- S6 natural convection: 6 CFD + 24 images + 28 QA (자율 탐색 02:37)
- S1 heated plate: 6 CFD + 24 images + 26 QA (자율 탐색 04:37 + 메인 06:03)
- S4_E6 인위적 변형 완료 (recirculation 1794 points 제거)
- Phase 상태: P1_pilot → P3_full_expansion 전환
- S3 turbulent channel CFD 실행 중

### 누적
- CFD 케이스: 30 (목표 50-70의 50%)
- 이미지: 138 (목표 250-300의 50%)
- QA 질문: 147
- 시나리오: 5/10 (S1, S2, S4, S5, S6)

### 이슈
- S4_E6 에이전트가 cell_data만 수정 → point_data도 직접 수정하여 해결
- S6 buoyantBoussinesqSimpleFoam: alphat BC, p_rgh 차원 문제 → 에이전트가 자체 해결

### 진행 중
- S3 turbulent channel CFD 백그라운드 실행

### 남은 작업 (Day 1 내)
- S3 렌더링 + QA
- 라벨링 배치 방안 A 형식으로 업데이트 (09:02 발송용)
- S7 mixed convection 또는 S10 ventilated room 시작
