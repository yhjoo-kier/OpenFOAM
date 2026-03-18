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

```
benchmark/
├── scenes/                          # rule-based 생성된 frozen scene JSON
├── manifests/
│   ├── scene_manifest.json          # category/seed/geometry manifest (현재 12 scenes)
│   ├── pilot_reference_status.json  # 8-scene pilot CFD status aggregate
│   ├── pilot_reference_status.csv
│   └── reference_batch_summary.json # 최근 batch rerun 요약
├── reference_cfd/<case>             # cases/<case> 로 연결되는 symlink
├── visualizations/<case>            # results/<case> 로 연결되는 symlink
├── renderings/
│   ├── renderings_manifest.json     # benchmark input image export summary
│   └── <case>/<view>/<case>_<view>.png
├── evaluations/                     # 이미지-conditioned 평가 scaffold / 결과
└── README.md                        # 이 파일
```

현재 상태:
- 8-scene pilot set(A1/A2/A3/A4 × 2)은 reference CFD 8/8 성공까지 회복됨.
- frozen scene subset은 12 scenes(A1/A2/A3/A4 × 3)로 확장됨.
- benchmark input view export는 `perspective`, `birdseye`, `floorplan`, `wireframe`, `section` 5종을 사용함.
- `run_benchmark_reference_batch.py`는 이제 reference CFD 실행과 benchmark input view export를 한 흐름으로 함께 처리함.
- `benchmark/evaluations/`에는 frozen-12 × 5-view 기준의 60개 evaluation task scaffold가 준비되어 있음.
- `scripts/run_benchmark_evaluation_task.py`로 개별 scaffold task를 실제 image-conditioned pipeline에 연결할 수 있으며, task 상태와 aggregate summary/manifest를 함께 갱신함.

## 관련 문서

- 논문 계획: `docs/26-03-17_image_to_cfd_paper_plan.md`
- evaluation scaffold 메모: `docs/26-03-18_benchmark_evaluation_scaffold.md`
- 프레임워크 스킬: `docs/skills/indoor-cfd-pipeline/SKILL.md`
