# Paper 1 (Image-to-CFD) 폴더 구조 개편

**날짜**: 2026-03-30
**목적**: Paper 1 파일들을 `papers/image_to_cfd/`로 이동하여 Paper 2 (`papers/cfd_visual_qa/`)와 분리 관리

## 배경

- Paper 1이 Building and Environment에서 desk reject → 다른 저널 재투고 예정
- Paper 2 (`papers/cfd_visual_qa/`)가 별도 브랜치에서 작업 후 main에 병합됨
- 두 논문 파일이 루트에 혼재하면 경로 충돌 및 관리 혼선 우려

## 이동 대상

### papers/image_to_cfd/docs/ (논문 작성물)
- `docs/paper_draft_v1.md` — 본문 원고 (Markdown)
- `docs/cover_letter.*` — 커버레터 (md + docx)
- `docs/highlights.*` — 하이라이트 (md + docx)
- `docs/graphical_abstract/` — GA 이미지
- `docs/appendix_vlm_prompt.md` — VLM 프롬프트 부록
- `docs/introduction_drafts.md` — 서론 초안
- `docs/paper_search_results.json` — 문헌 검색 결과
- `docs/paperbanana_*.json` — PaperBanana 입력 설정
- `docs/references/` — RIS 파일들
- `docs/figure_qc/` — Figure QC 기록 + 자산
- `docs/26-03-14_*` ~ `docs/26-03-23_*` — Paper 1 작업 로그 (67개)
- `docs/26-03-17_image_to_cfd_paper_plan.md` — 논문 계획
- `docs/skills/indoor-cfd-pipeline` — 프레임워크 스킬

### papers/image_to_cfd/latex/ (LaTeX 빌드)
- `latex/` 전체 — paper.tex, figures/, files_for_submission/, 빌드 산출물, 참고 커버레터

### papers/image_to_cfd/scripts/ (빌드 도구)
- `scripts/build_paper.sh` — PDF 빌드 스크립트
- `scripts/md2latex.py` — MD→LaTeX 변환기
- `scripts/truncate_bib_authors.py` — BibTeX 저자 축약
- `scripts/paper_figures/` — Figure 생성 스크립트 (24개)

### papers/image_to_cfd/results/ (Figure 후보)
- `results/paperbanana_selected/` — 선택된 PaperBanana figures
- `results/paperbanana_fig*` — 후보 이미지 (gitignored)
- `results/paperbanana_ga*` — GA 후보 (gitignored)
- `results/paper_figures_phase2/` — Phase 2 figure 산출물

## 루트에 남기는 것

- `benchmark/` — 벤치마크 데이터셋 (공유 자산)
- `cases/` — OpenFOAM 시뮬레이션 케이스
- `src/` — 프레임워크 소스 코드 (geometry, visualization, grid_study)
- `scripts/run_indoor_stabilized.py` — 프레임워크 엔트리포인트
- `scripts/evaluate*.py`, `scripts/build_benchmark*.py` — 벤치마크 도구
- `docs/26-01-26_*` — Grid study 문서 (Paper 1 무관)

## build_paper.sh 경로 변경

Before:
```
PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
MD_SOURCE="$PROJECT_DIR/docs/paper_draft_v1.md"
LATEX_DIR="$PROJECT_DIR/latex"
SCRIPTS_DIR="$PROJECT_DIR/scripts"
FIG_SOURCE="$PROJECT_DIR/results/paper_figures_phase2"
BANANA_FIGS="$PROJECT_DIR/results/paperbanana_selected"
```

After:
```
PAPER_DIR="$(cd "$(dirname "$0")/.." && pwd)"
REPO_ROOT="$(cd "$PAPER_DIR/../.." && pwd)"
MD_SOURCE="$PAPER_DIR/docs/paper_draft_v1.md"
LATEX_DIR="$PAPER_DIR/latex"
SCRIPTS_DIR="$PAPER_DIR/scripts"
FIG_SOURCE="$PAPER_DIR/results/paper_figures_phase2"
BANANA_FIGS="$PAPER_DIR/results/paperbanana_selected"
```

BIB_SOURCE는 절대경로이므로 변경 불필요.

## CLAUDE.md 변경

- 루트 CLAUDE.md: Paper 1 전용 섹션 제거, 공통 컨벤션만 유지, `papers/` 인덱스 추가
- `papers/image_to_cfd/CLAUDE.md`: Paper Data Convention, 최종 수치, Figure 컨벤션 등 이동

## 검증

이동 완료 후 `bash papers/image_to_cfd/scripts/build_paper.sh`로 PDF 빌드 성공 확인.
