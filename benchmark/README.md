# Image-to-CFD Benchmark Dataset

VLM 기반 2D→3D CFD 자동화 프레임워크의 정량 평가를 위한 벤치마크 데이터셋.

## 목적

- Rule-based로 생성한 실내 형상(ground truth)과 reference CFD 결과를 확보
- 동일 형상의 2D 렌더링을 프레임워크에 입력하여 결과를 정답과 비교

## 난이도 체계 (2×2)

|  | Rectangular | L-shaped |
|---|---|---|
| **단순 장애물** (0~1개) | A1 | A3 |
| **복잡 장애물** (2~3개+) | A2 | A4 |

## 디렉토리 구조 (현재)

```text
benchmark/
├── scenes/                          # rule-based 생성된 frozen scene JSON
├── manifests/
│   ├── scene_manifest.json          # category/seed/geometry manifest (현재 20 scenes)
│   ├── pilot_reference_status.json  # 8-scene pilot CFD status aggregate
│   ├── pilot_reference_status.csv
│   └── reference_batch_summary.json # 전체 frozen set aggregate reference status
├── reference_cfd/<case>             # cases/<case> 로 연결되는 symlink
├── visualizations/<case>            # results/<case> 로 연결되는 symlink
├── renderings/
│   ├── renderings_manifest.json     # benchmark input image export summary
│   └── <case>/<view>/<case>_<view>.png
├── evaluations/                     # 이미지-conditioned 평가 scaffold / 결과
└── README.md                        # 이 파일
```

현재 상태:
- 8-scene pilot set (A1/A2/A3/A4 × 2)은 reference CFD **8/8 성공**까지 회복됨.
- frozen scene subset은 **20 scenes** (A1/A2/A3/A4 × 5)까지 확장되었고 reference CFD도 **20/20 성공** 상태임.
- benchmark input view export는 `perspective`, `birdseye`, `floorplan`, `wireframe`, `section` 5종을 사용함.
- `run_benchmark_reference_batch.py`는 reference CFD 실행과 benchmark input view export를 한 흐름으로 함께 처리하며, `--collect-only`로 전체 aggregate manifest를 새로 고칠 수 있음.
- `benchmark/evaluations/`에는 frozen-20 × 5-view 기준의 **100개 evaluation task scaffold**가 준비되어 있음.
- `scripts/run_benchmark_evaluation_task.py`로 개별 scaffold task를 실제 image-conditioned pipeline에 연결할 수 있으며, task 상태와 aggregate summary/manifest를 함께 갱신함.
- 성공한 evaluation task는 geometry 비교뿐 아니라 `scripts/compute_benchmark_cfd_metrics.py`를 통해 normalized-grid 기준의 CFD result-side 비교(`cfd_metrics.json`)도 추가로 기록함.
- 현재 evaluation smoke test 결과는 **1 blocked / 99 pending**이며, blocked 상태는 benchmark failure가 아니라 Gemini backend credential 부재를 의미함.

운영상 주의:
- 부분 rerun 뒤에는 `benchmark/manifests/reference_batch_summary.json`이 최근 batch만 반영할 수 있으므로, frozen-set 전체 상태를 복구하려면 `python3 scripts/run_benchmark_reference_batch.py --collect-only`를 다시 실행하는 것이 안전함.
- 이후 `python3 scripts/scaffold_benchmark_evaluations.py`를 실행해 evaluation scaffold를 full-set 기준으로 맞춰 주는 것을 기본 절차로 삼는 것이 좋음.

## 관련 문서

- 논문 계획: `docs/26-03-17_image_to_cfd_paper_plan.md`
- frozen-20 상태 메모: `docs/26-03-18_frozen20_benchmark_status.md`
- evaluation scaffold 메모: `docs/26-03-18_benchmark_evaluation_scaffold.md`
- 프레임워크 스킬: `docs/skills/indoor-cfd-pipeline/SKILL.md`
