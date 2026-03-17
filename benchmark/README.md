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

## 디렉토리 구조 (예정)

```
benchmark/
├── scenes/              # rule-based 생성된 scene JSON
├── reference_cfd/       # ground truth CFD 결과
├── renderings/          # 뷰 타입별 2D 입력 이미지
├── evaluations/         # 프레임워크 적용 결과 및 평가
├── generate_scenes.py   # rule-based 형상 생성기
└── README.md            # 이 파일
```

## 관련 문서

- 논문 계획: `docs/26-03-17_image_to_cfd_paper_plan.md`
- 프레임워크 스킬: `docs/skills/indoor-cfd-pipeline/SKILL.md`
