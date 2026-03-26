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
- ~~S3 렌더링 + QA~~ ✅
- ~~라벨링 배치 방안 A 형식으로 업데이트~~ → Google Drive DOCX로 전환
- ~~S7 mixed convection 또는 S10 ventilated room 시작~~ ✅

### Day 1 최종 성과 (2026-03-25)
**MILESTONE: 10/10 시나리오 전체 완료 (원래 계획 Day 6-10)**

- CFD 케이스: 60개 (S1-S10 각 6개)
- 이미지: 258개
- QA 질문: 279개
- 전문가 라벨: 35개 (Day 1: 20, Day 2: 15)
- VLM 평가: Claude Opus blind 15개 (100% 정확도)
- 라벨링 시스템: Telegram → Google Drive DOCX 전환
- 블라인드 코드 도입 (케이스명 노출 방지)

**핵심 발견:**
1. Claude blind (100%) > 전문가 (80%) on 15-item pilot
2. Contamination paradox: 케이스명 노출 시 오히려 성능 저하 (87% vs 100%)
3. BC swap, wrong viscosity = consistently subtle errors
4. Gravity flip = consistently detectable
5. Coarse mesh = context-dependent (미세구조 유무에 따라)

### 교훈
- writeInterval > endTime → 결과 미저장. writeInterval=500 + purgeWrite 1 권장
- foamToVTK 타임아웃: 120s 부족 → 600s
- 케이스명 익명화는 평가 초기부터 적용해야 함
- VLM 평가 시 독립 서브에이전트(블라인드) 필수

---

---

## Day 2 (2026-03-26) — Paper Draft Complete

### 완료
- Gemini Day 2 블라인드 평가 (62.5% valid, 33% no-response)
- Claude Day 2 블라인드 평가 (100%, 30/30 종합)
- GPT-5.4 (Codex) 시도 → CLI 이미지 전달 실패
- 논문 아웃라인 작성 (7 sections, 7 figures, 5 tables)
- 논문 전체 초안 완료 (4,066 words, abstract + sections 1-7)
- Figure 7/7 생성 완료
- Supervisor review 반영 (directive v1.1)
- Supervisor review 확인 cron 추가

### 누적
- CFD: 60 cases, 258 images, 279 QA
- 라벨: 35/80 (Day 3-5 대기)
- VLM: Claude 30/30=100%, Gemini 16/20=80% (valid)
- 논문: 4,066 words + 7 figures

### 이슈
- Gemini CLI 무응답률 33% (image processing 불안정)
- Codex CLI 이미지 전달 불가
- Day 3-5 라벨 아직 수거 안 됨

### Day 3 계획 (2026-03-27)
- Day 3-5 라벨링 수거 (유저 편의에 따라)
- 라벨 수거 시: VLM 확대 평가 + 결과 섹션 업데이트
- 논문 polishing (figure 참조, citation 형식, table 삽입)
- Supervisor review 확인 (20:30)

---

## Day 2 (2026-03-26) — Full Evaluation + Paper Completion

### 완료
- Day 3 전문가 라벨 수거 (60%, 9/15)
- Day 4+5 전문가 라벨 수거 (각 66.7%)
- **전체 80개 전문가 라벨 완료** — 73.8% accuracy, 2 FP, 19 FN
- Claude Day 3 블라인드 평가: 93.3% (14/15) — 첫 FP (S8_correct_lam)
- Claude Day 4+5 블라인드 평가: 100% (30/30)
- **Claude 종합 75항목: 98.7%** (74/75, 1 FP, 0 FN)
- Gemini 무응답 13건 재시도: 11/11 valid, 86.2% (25/29)
- VLM 입출력 전체 문서화 (transcripts + review DOCX)
- 논문 Sections 1-7 초안 완료 (4,066 words)
- Figure 7/7 생성 완료
- 논문 아웃라인 작성
- Supervisor review 반영 (directive v1.1)

### 누적
- CFD: 60 cases, 258 images, 279 QA
- 전문가 라벨: 80/80 완료 (73.8%)
- Claude: 75 items, 98.7%
- Gemini: 29 valid items, 86.2%
- 논문: 4,066 words + 7 figures

### 이슈
- S8_correct_lam: Claude/전문가 모두 FP — 복잡 후류 패턴
- Gemini: S1_correct_lam에서 2 FP (복잡 자연대류)
- 3회 독립 시행 여부 유저 논의 중

### Day 3 계획
- 3회 시행 결정 시 추가 VLM 평가 실행
- 논문 Results 섹션 최종 수치 업데이트
- Supervisor review 확인
