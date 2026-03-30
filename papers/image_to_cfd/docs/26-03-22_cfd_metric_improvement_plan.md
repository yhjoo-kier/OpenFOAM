# CFD Metric 개선 계획

**작성일**: 2026-03-22
**대상 파일**: `scripts/compute_benchmark_cfd_metrics.py`
**검토 근거**: Gemini 3.1 Pro 및 Codex GPT-5.4 독립 검토 (docs/26-03-22_cfd_metric_review_*.md)
**상태**: 계획 수립 완료, 구현 대기

---

## 배경

현재 CFD Fidelity Score는 4개 컴포넌트의 단순 산술 평균으로 구성:

```
CFD Score = mean(Overlap Ratio, Velocity Mag Similarity, Velocity Dir Similarity, Pressure Similarity)
```

두 개의 독립적 AI 모델(Gemini, Codex) 검토 결과, 리프레이밍(서술/주장 범위 축소)과 별개로 메트릭 정의 자체에 4가지 기술적 개선이 필요한 것으로 확인됨.

---

## 개선 항목

### 1. 압력 Gauge 보정 [Critical]

**문제**: 비압축성 솔버(simpleFoam)에서 압력은 p/rho 형태의 상대값. 참조 해와 예측 해의 gauge 기준점이 다르면, gradient가 동일해도 raw RMSE가 크게 나옴. 현재 코드(`scalar_metrics`, line 154-166)는 raw 압력값을 직접 비교.

**현재 코드** (`compare_cases`, line 260-262):
```python
if "p" in ref_sample and "p" in pred_sample:
    field_metrics["pressure"] = scalar_metrics(ref_sample["p"][overlap], pred_sample["p"][overlap])
```

**개선안**: 겹침 영역에서 각각의 평균 압력을 차감한 후 비교

```python
if "p" in ref_sample and "p" in pred_sample:
    ref_p = ref_sample["p"][overlap]
    pred_p = pred_sample["p"][overlap]
    # Gauge correction: subtract domain-mean pressure from each field
    ref_p_corrected = ref_p - np.mean(ref_p)
    pred_p_corrected = pred_p - np.mean(pred_p)
    field_metrics["pressure"] = scalar_metrics(ref_p_corrected, pred_p_corrected)
    field_metrics["pressure"]["gauge_corrected"] = True
    field_metrics["pressure"]["reference_mean_removed"] = round_or_none(float(np.mean(ref_p)))
    field_metrics["pressure"]["predicted_mean_removed"] = round_or_none(float(np.mean(pred_p)))
```

**영향 범위**:
- `compare_cases()` 함수 내 압력 비교 블록 (line 260-262)
- 모든 기존 evaluation_summary.json의 pressure_similarity 값이 변경됨
- Phase 2 재실행 시 자동 반영되므로, Phase 2 전에 코드 수정 권장

**검증 방법**: 기존 결과와 보정 후 결과를 나란히 출력하여 pressure_similarity 변화 확인

---

### 2. 엣지 케이스 정의 보완 [Important]

**문제**: 다음 경계 조건에서 메트릭이 미정의이거나 불안정:

| 엣지 케이스 | 현재 처리 | 문제점 |
|------------|----------|--------|
| `ref_RMS = 0` | `relative_rmse = None` → 컴포넌트 스킵 | 컴포넌트 수가 케이스마다 달라져 점수 비교 불가 |
| 한쪽 벡터만 `|U| ≈ 0` | 비교에 포함됨 | cosine 미정의, 수치 불안정 |
| 겹침 없음 (`overlap = 0`) | `cfd_score = 0.0` | 적절 (유지) |
| `ref_RMS` 극소 (정체영역 지배) | 정상 계산 | relative_rmse 폭발 → 항상 0점 |

**개선안**:

#### 2a. ref_RMS = 0 처리
```python
# scalar_metrics / vector_metrics 내부
if ref_rms is not None and ref_rms < RMS_FLOOR:
    ref_rms = RMS_FLOOR  # e.g., 1e-6
    flagged_rms_floor = True
```
- `RMS_FLOOR` 상수 도입 (기본값 `1e-6`)
- 플래그 기록: `"rms_floored": true`

#### 2b. 단측 영벡터 direction cosine 처리
현재 코드 (line 139-141):
```python
denom = ref_mag * pred_mag
direction_mask = denom > 1e-12  # 양쪽 곱이 0이면 제외
```

개선: 한쪽이라도 0에 가까우면 제외
```python
direction_mask = (ref_mag > 1e-8) & (pred_mag > 1e-8)
```
- 임계값을 `1e-12`(곱) → `1e-8`(각각)으로 변경하여 수치 안정성 확보

#### 2c. 컴포넌트 수 고정
현재는 `None`인 컴포넌트를 스킵하여 분모가 가변적:
```python
# 현재: components 리스트에 조건부 append
components.append(max(0.0, 1.0 - min(1.0, umag_rel_rmse)))
```

개선: 항상 4개 컴포넌트 유지, 계산 불가 시 0.0 + 플래그
```python
# 컴포넌트가 None이면 0.0으로 처리하되 플래그 기록
umag_sim = max(0.0, 1.0 - min(1.0, umag_rel_rmse)) if umag_rel_rmse is not None else 0.0
components.append(umag_sim)
if umag_rel_rmse is None:
    notes.append("velocity_magnitude: ref_RMS unavailable, scored as 0.0")
```

**영향 범위**:
- `vector_metrics()`, `scalar_metrics()`, `compare_cases()` 함수
- 점수 비교의 일관성 향상 (케이스 간 동일 분모)

---

### 3. 컴포넌트 개별 보고 + 복합 점수 보조화 [Important]

**문제**: 현재 `cfd_score` 단일 값이 주 결과로 보고됨. 4개 컴포넌트의 균등 평균은 물리적 정당성이 없으며, 실패 모드를 숨김.

**개선안 — 논문 서술 + 코드 출력 구조 변경**:

#### 3a. 출력 JSON 구조 확장
```json
{
  "aggregate_score": {
    "cfd_agreement_score": 0.491,
    "scoring_note": "Unweighted mean of 4 components. Report components individually for detailed analysis.",
    "components": {
      "overlap_ratio_vs_union": 0.938,
      "velocity_magnitude_similarity": 0.000,
      "velocity_direction_similarity": 0.615,
      "pressure_similarity": 0.675
    },
    "component_count": 4,
    "components_floored": []
  }
}
```

#### 3b. 메트릭 이름 변경
코드 내 키 이름 변경:
- `cfd_score` → `cfd_agreement_score`
- 논문 내: "CFD Fidelity Score" → "Reference CFD Agreement Score"

#### 3c. 논문 결과 표 형식 제안
| Case | Overlap | Vel. Mag. | Vel. Dir. | Pressure | Agreement |
|------|---------|-----------|-----------|----------|-----------|
| a1_01 | 0.938 | 0.000 | 0.615 | 0.675 | 0.557 |
| ... | ... | ... | ... | ... | ... |
| **Mean** | **0.93** | **0.12** | **0.58** | **0.65** | **0.491** |

- 4개 컴포넌트 열이 주 결과, Agreement 열은 보조 요약

**영향 범위**:
- `compare_cases()` 출력 구조
- `build_benchmark_evaluation_aggregate.py` (집계 스크립트)
- 논문 draft Section 5 (Results)
- 하위 호환성: 기존 키 `cfd_score` 유지하되 `cfd_agreement_score`도 병행 출력

---

### 4. 표준 CFD 메트릭 병행 보고 [Recommended]

**문제**: 현재 메트릭이 CFD 커뮤니티의 표준 검증 메트릭(ASME V&V 20, COST 732)과 괴리.

**추가할 표준 메트릭**:

#### 4a. Hit Rate (q)
COST 732에서 권장하는 CFD 검증 메트릭:
```
q = (1/N) * sum( |pred_i - ref_i| < max(W, D * |ref_i|) )
```
- W: 절대 허용 오차 (속도의 경우 통상 inlet velocity의 비율)
- D: 상대 허용 오차 (통상 0.25)
- 실내 CFD에서 D=0.25가 일반적

```python
def hit_rate(ref: np.ndarray, pred: np.ndarray, D: float = 0.25, W: float = 0.0) -> float:
    """COST 732 hit rate metric."""
    diff = np.abs(pred - ref)
    threshold = np.maximum(W, D * np.abs(ref))
    return float(np.mean(diff < threshold))
```

#### 4b. FAC2 (Factor of 2)
```
FAC2 = fraction of points where 0.5 <= pred_i/ref_i <= 2.0
```

```python
def fac2(ref: np.ndarray, pred: np.ndarray) -> float:
    """Fraction of predictions within factor of 2."""
    valid = np.abs(ref) > 1e-10  # avoid division by near-zero
    if not np.any(valid):
        return 1.0
    ratio = pred[valid] / ref[valid]
    return float(np.mean((ratio >= 0.5) & (ratio <= 2.0)))
```

#### 4c. NMSE (Normalized Mean Square Error)
```
NMSE = mean((pred - ref)^2) / (mean(pred) * mean(ref))
```

#### 4d. FB (Fractional Bias)
```
FB = 2 * (mean(ref) - mean(pred)) / (mean(ref) + mean(pred))
```

#### 4e. 상관계수 R
```
R = corr(ref, pred)
```

**구현 위치**: `compute_benchmark_cfd_metrics.py`의 `compare_cases()` 함수 내, 기존 `field_metrics` 딕셔너리에 추가

**출력 구조 확장**:
```json
{
  "field_metrics": {
    "velocity_magnitude": {
      "mae": 0.012,
      "rmse": 0.018,
      "relative_rmse": 0.85,
      "hit_rate_D025": 0.42,
      "fac2": 0.55,
      "nmse": 0.72,
      "fb": -0.05,
      "correlation_R": 0.68
    },
    "pressure": { ... }
  }
}
```

**영향 범위**:
- `compute_benchmark_cfd_metrics.py` (메트릭 함수 추가 + compare_cases 확장)
- 논문 Results 표에 Hit Rate, FAC2 열 추가
- 기존 결과에 영향 없음 (추가 메트릭만 병행)

---

## 구현 우선순위 및 의존성

```
[Phase A: 코드 수정] ─────────────────────────────────
│
├─ 1. 압력 gauge 보정 (Critical, 독립적)
├─ 2. 엣지 케이스 보완 (Important, 독립적)
├─ 3. 출력 구조 확장 + 이름 변경 (Important, 1·2와 병행 가능)
└─ 4. 표준 메트릭 함수 추가 (Recommended, 독립적)

[Phase B: 검증] ──────────────────────────────────────
│
├─ 기존 케이스 1-2개로 before/after 비교 실행
├─ 압력 보정 효과 확인 (pressure_similarity 변화)
└─ 표준 메트릭 값 sanity check

[Phase C: 전체 재평가] ───────────────────────────────
│
├─ Phase 2 (0.18m 메시) 재실행 시 자동 반영
│  또는 기존 0.35m 결과에 대해 메트릭만 재계산
└─ 집계 스크립트 업데이트 (build_benchmark_evaluation_aggregate.py)

[Phase D: 논문 반영] ─────────────────────────────────
│
├─ Results 표 구조 변경 (컴포넌트 개별 + 표준 메트릭)
├─ Methodology 섹션 메트릭 정의 업데이트
├─ 이름 변경: CFD Fidelity Score → Reference CFD Agreement Score
└─ Limitation 섹션에 CFD-to-CFD 비교 한계 명시
```

---

## 하위 호환성 고려

- 기존 `evaluation_summary.json`의 `cfd_score` 키는 유지 (deprecated 마킹)
- 새 키 `cfd_agreement_score` 병행 출력
- 기존 집계/시각화 스크립트가 `cfd_score`를 참조하므로 즉시 삭제하지 않음
- Phase 2 완료 후 일괄 마이그레이션

---

## 참고 문서

- Gemini 검토: `docs/26-03-22_cfd_metric_review_gemini.md`
- Codex 검토: `docs/26-03-22_cfd_metric_review_codex.md`
- ASME V&V 20: Standard for Verification and Validation in CFD and Heat Transfer
- COST 732: Best Practice Guideline for CFD Simulation of Flows in Urban Environment
