# Image-to-CFD Paper

2D 이미지로부터 VLM(Gemini)을 이용해 3D 실내 형상을 추상화하고 CFD까지 자동화하는 프레임워크의 논문.

## Status

- **제출**: Building and Environment (2026-03-27) → Desk reject
- **다음**: 다른 저널 재투고 예정
- **figshare**: https://doi.org/10.6084/m9.figshare.31866127
- **Funding**: KIER C6-2419-63

## Paper Structure
```
papers/image_to_cfd/
├── docs/                       # 원고, 커버레터, 작업 로그
│   ├── paper_draft_v1.md       # 본문 (Markdown 원본)
│   ├── cover_letter.*          # 커버레터 (md + docx)
│   ├── highlights.*            # 하이라이트 (md + docx)
│   ├── graphical_abstract/     # GA 이미지
│   └── 26-03-*_*.md            # 작업 로그
├── latex/                      # LaTeX 빌드
│   ├── paper.tex               # 생성된 LaTeX (md2latex.py 산출물)
│   ├── paper.pdf               # 최종 PDF
│   ├── figures/                # Figure 자산 (PDF + PNG)
│   └── files_for_submission/   # 제출용 패키지
├── scripts/                    # 빌드 도구
│   ├── build_paper.sh          # PDF 빌드 (bash papers/image_to_cfd/scripts/build_paper.sh)
│   ├── md2latex.py             # Markdown → LaTeX 변환기
│   ├── truncate_bib_authors.py # BibTeX 저자 축약
│   └── paper_figures/          # Figure 생성 스크립트
└── results/                    # Figure 후보
    ├── paperbanana_selected/   # 선택된 PaperBanana figures (4K)
    └── paper_figures_phase2/   # Phase 2 figure 산출물
```

## Related Root-Level Assets

프레임워크 코드와 벤치마크 데이터는 루트에 위치:
- **프레임워크**: `scripts/run_indoor_stabilized.py`
- **벤치마크**: `benchmark/` (scenes, renderings, evaluations, manifests)
- **케이스**: `cases/phase2_ref_*`, `cases/phase2_pred_*`

## Build

```bash
bash papers/image_to_cfd/scripts/build_paper.sh
```

## Paper Data Convention (Critical — 반드시 준수)

논문 결과는 **preset-matched** 데이터만 사용한다. Solver robustness preset마다 inlet velocity가 다르므로, reference와 predicted 케이스가 동일한 preset(=동일 inlet velocity)으로 수렴한 pair만 비교해야 한다.

**⚠️ 사용 금지:**
- `cfd_metrics.json` — preset 미매칭 상태의 CFD 비교 결과 (inlet velocity 불일치)
- `evaluation_aggregate_summary_posthoc_scaled_longest_span.json`의 원본 CFD score (구 버전)

**✅ 사용해야 하는 파일:**
- `cfd_metrics_matched.json` — 각 eval dir 내, preset-matched reference 기준 CFD 비교
- `benchmark/manifests/evaluation_aggregate_summary_phase2.json` — 매칭된 집계 (최신)
- `benchmark/manifests/evaluation_statistics_phase2.json` — 매칭된 통계 (최신)
- `benchmark/manifests/evaluation_statistics_matched.json` — 매칭 전용 통계

**Reference 케이스 경로:**
- Preset-matched: `cases/phase2_ref_{scene}_preset_{preset}/` (논문용)
- 원본 Phase 2: `cases/phase2_ref_{scene}/` (참고용만, 논문에 사용 금지)

**Predicted 케이스 경로:**
- `cases/phase2_pred_bench_{case}_{view}/`

**최종 수치 (2026-03-23):**
- Structural: 0.781 ± 0.151, CFD agreement: 0.477 ± 0.158 (n=97)
- Best view: floorplan (struct 0.884, CFD 0.572)

**배경:** Solver preset별 inlet velocity가 다름 (robust=0.02, ultra_robust=0.005, laminar=0.01 m/s). Phase 2 초기에는 ref/pred가 다른 preset으로 수렴하여 velocity magnitude similarity ≈ 0이었음. `--force-preset` 플래그 도입 + 26개 reference 재해석으로 해결 (2026-03-23).

## Paper Writing Convention (Markdown → LaTeX → PDF)

### Figure
- 폰트: Arial (불가 시 Liberation Sans → DejaVu Sans fallback, QC log에 기록)
- 파일 형식: PDF (벡터, 기본) + PNG (600 dpi 이상)
- 네이밍: `fig_섹션키워드_설명.pdf` (숫자 번호 사용 금지, LaTeX에서 자동 번호 할당)
- 본문 참조: `[Fig:label]` → LaTeX 변환 시 `Fig.~\ref{fig:label}`

### Table
- 캡션+레이블을 표 상단에 배치: `[Table:키워드-설명]`
- LaTeX 변환 시 `Table~\ref{tab:키워드-설명}`

### Equation
- `\tag{Eq:키워드-변수}` 형식
- LaTeX 변환 시 `Eq.~\ref{eq:키워드-변수}`

### Citation
- BibTeX: `C:\Vaults\Research\Zotero_YHJoo.bib` (Zotero BetterBibTeX 자동 내보내기)
- 초안: `[저자 연도 키워드]` → 확정 시 `[@citekey]` → LaTeX `\cite{citekey}`
- `elsarticle-num` + `natbib sort&compress`

### LaTeX 빌드
- 빌드 순서: pdflatex → bibtex → pdflatex → pdflatex
- 변환 규칙: `[Fig:label]` → `Fig.~\ref{fig:label}`, `[@citekey]` → `\cite{citekey}`
