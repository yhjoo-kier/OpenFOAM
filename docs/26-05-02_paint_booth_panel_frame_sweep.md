# Paint Booth Panel-Frame v036 Calibration Sweep

## 목적

`paint_booth_plenum_filter_panel_frame_v036` 구조를 기준으로 supply velocity와 filter Forchheimer coefficient를 sweep하여 다음 항목을 비교했다.

- filter pressure drop proxy
- filter 아래 평균 downdraft
- work-zone 평균 downdraft
- near-car 평균 downdraft
- mass balance
- 역류 fraction

## 실행 스크립트

신규 스크립트:

```text
scripts/run_paint_booth_panel_frame_sweep.py
scripts/plot_paint_booth_panel_frame_sweep.py
```

실행 명령:

```bash
python3 scripts/run_paint_booth_panel_frame_sweep.py --force
python3 scripts/plot_paint_booth_panel_frame_sweep.py
```

주의: Docker postprocess 단계에서 PyVista/VTK cleanup이 interpreter shutdown에서 멈추는 경우가 있어 `scripts/postprocess_paint_booth_plenum_metrics.py` 마지막에 output write 이후 `os._exit(0)`로 종료하도록 보완했다.

## Sweep 범위

Supply velocity:

```text
2.30, 3.00, 3.50, 4.36 m/s
```

Filter Forchheimer coefficient:

```text
6.8e3, 1.0e4, 2.0e4, 4.0e4
```

총 case 수:

```text
16 cases
```

Case root:

```text
cases/paint_booth_panel_frame_sweep/
```

집계 결과:

```text
cases/paint_booth_panel_frame_sweep/sweep_summary.json
cases/paint_booth_panel_frame_sweep/sweep_summary.csv
```

## 검증 상태

모든 16개 case에서 다음을 확인했다.

```text
case directories: 16
post_plenum_metrics.json: 16
internal.vtu: 16
checkMesh Mesh OK: 16
simpleFoam End: 16
```

즉, 전체 sweep은 mesh/solver/VTK/postprocess 기준으로 정상 완료되었다.

## 주요 결과 요약

정렬 기준은 heuristic score이다.

- target near-car downdraft: 0.30 m/s
- target pressure drop: 50 Pa
- near-car 목표, work-zone 목표, pressure drop, near-car reverse fraction을 가중 합산

```text
case                         U_supply  f       Δp[Pa]  filter|Uz|  work|Uz|  near|Uz|  near rev
panel_frame_u4p36_f10000     4.36      10000   51.0    0.652       0.285     0.240     7.8%
panel_frame_u4p36_f6800      4.36      6800    34.9    0.652       0.291     0.247     7.8%
panel_frame_u3p50_f10000     3.50      10000   32.8    0.523       0.229     0.193     7.9%
panel_frame_u3p50_f20000     3.50      20000   65.0    0.523       0.223     0.187     8.0%
panel_frame_u3p00_f20000     3.00      20000   47.8    0.448       0.191     0.160     8.0%
panel_frame_u2p30_f40000     2.30      40000   55.7    0.343       0.144     0.121     8.1%
```

전체 heatmap:

```text
docs/figures/26-05-02_paint_booth_panel_frame_sweep_heatmaps.png
```

## 해석

### 1. Pressure drop은 supply velocity와 filter coefficient 모두에 강하게 의존

가장 작은 조건:

```text
U=2.30, f=6800 -> 약 9.7 Pa
```

가장 큰 조건:

```text
U=4.36, f=40000 -> 약 200 Pa
```

목표 50 Pa 근처 조건:

```text
U=4.36, f=10000 -> 약 51.0 Pa
U=3.00, f=20000 -> 약 47.8 Pa
U=2.30, f=40000 -> 약 55.7 Pa
```

### 2. Downdraft 속도는 supply velocity가 지배적

Filter-below 평균 속도는 같은 supply velocity 내에서 filter coefficient가 바뀌어도 큰 차이가 없었다.

예:

```text
U=2.30 -> filter_below |Uz| ≈ 0.343~0.344 m/s
U=3.00 -> filter_below |Uz| ≈ 0.447~0.449 m/s
U=3.50 -> filter_below |Uz| ≈ 0.522~0.523 m/s
U=4.36 -> filter_below |Uz| ≈ 0.650~0.652 m/s
```

이는 현재 boundary condition이 supply inlet fixed velocity이기 때문에 총 유량이 inlet velocity로 강하게 고정되기 때문이다. Filter coefficient는 주로 필요한 압력강하를 조정한다.

### 3. 목표 선택은 어떤 위치를 기준으로 하느냐에 따라 달라짐

#### Filter face 직하부를 0.35 m/s에 맞추려면

가장 적절한 후보:

```text
U_supply = 2.30 m/s
f = 40000
filter_below |Uz| ≈ 0.343 m/s
Δp ≈ 55.7 Pa
```

단점:

```text
work_zone |Uz| ≈ 0.144 m/s
near_car |Uz| ≈ 0.121 m/s
```

즉, filter 아래 목표는 맞지만 차체 주변 downdraft는 약하다.

#### Work-zone/near-car downdraft를 더 키우려면

가장 균형적인 후보:

```text
U_supply = 4.36 m/s
f = 10000
work_zone |Uz| ≈ 0.285 m/s
near_car |Uz| ≈ 0.240 m/s
Δp ≈ 51.0 Pa
```

장점:

```text
pressure drop이 50 Pa에 매우 가까움
work-zone downdraft가 0.25~0.30 m/s 수준
near-car downdraft도 sweep 내 최대권
```

단점:

```text
filter_below |Uz| ≈ 0.652 m/s로 filter 직하부는 빠름
```

## 현재 추천 후보

현재 연구 목적이 “차체 주변 도장부스 작업공간 유동”이면 우선 후보는 다음이 더 적합하다.

```text
Recommended v036 calibrated candidate:
  U_supply = 4.36 m/s
  filter Forchheimer f = 10000

Metrics:
  Δp ≈ 51.0 Pa
  work-zone |Uz| ≈ 0.285 m/s
  near-car |Uz| ≈ 0.240 m/s
  near-car reverse fraction ≈ 7.8%
```

반대로 “filter face 통과 속도 0.35 m/s”를 절대 기준으로 삼으면 다음이 적합하다.

```text
Filter-face target candidate:
  U_supply = 2.30 m/s
  filter Forchheimer f = 40000

Metrics:
  Δp ≈ 55.7 Pa
  filter-below |Uz| ≈ 0.343 m/s
  work-zone |Uz| ≈ 0.144 m/s
  near-car |Uz| ≈ 0.121 m/s
```

## 다음 작업 제안

1. `U=4.36, f=10000`을 calibrated v036 기준 case로 별도 생성
2. 해당 case의 velocity vector/streamline, pressure field, Uz field 이미지 생성
3. near-car reverse flow 감소를 위해 floor exhaust 구조를 full-floor에서 grating/slot variant로 변경
4. filter/frame을 high-resistance porous approximation에서 실제 internal baffle/wall 구조로 개선하는 v037 검토
5. pressure-drop 계산을 VTK slab proxy에서 OpenFOAM sampled surface/functionObject 기반으로 정교화
