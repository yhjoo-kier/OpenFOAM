# Paint Booth L2 Diagnostic Visualization QC

## 배경

사용자가 이전 diagnostic visualization이 엉망이라고 지적했다. 시각 QC를 수행한 결과, 기존 그림들은 OpenFOAM VTK unstructured/polyhedral mesh를 regular grid로 직접 point-sampling하면서 강한 row-wise artifact가 발생한 것으로 판단했다.

QC 대상 기존 그림:

- `26-05-02_paint_booth_l2_diagnostic_uz_deviation_horizontal_stack.png`
- `26-05-02_paint_booth_l2_diagnostic_reverse_flow_mask_midplane.png`
- `26-05-02_paint_booth_l2_diagnostic_clipped_uz_midplane.png`
- `26-05-02_paint_booth_l2_diagnostic_vorticity_midplane.png`

## 기존 그림 QC 판정

기존 diagnostic set은 **QC fail**로 판단한다.

주요 문제:

1. **수평 stripe/banding artifact**
   - 모든 패널에서 y 또는 z 방향 row-wise 줄무늬가 강하게 나타났다.
   - 물리적 coherent structure라기보다 sampling/interpolation artifact 가능성이 높다.

2. **solid car 내부 값 표시 문제**
   - reverse-flow/vorticity 그림에서 차체 내부에도 field/mask가 보이는 듯했다.
   - solid mask가 제대로 적용되지 않아 해석 혼란을 유발했다.

3. **reverse-flow mask 가독성 문제**
   - `Uz > 0` 영역이 면이 아니라 얇은 검은 선분처럼 보였다.
   - 배경 color, geometry edge, grid line과 구분되지 않았다.
   - threshold noise를 역류로 과장했을 가능성이 있다.

4. **vorticity plot QC fail**
   - raw VTK derivative 기반 `|curl(U)|`는 수평 streak와 흰색 고값/NaN 유사 artifact가 지배적이었다.
   - 차체 주변 shear/wake를 안정적으로 해석하기 어려웠다.

5. **색상/스케일/annotation 부족**
   - 상대편차 평균 기준, clipping, mask 처리, solid 영역 제외 여부가 충분히 명확하지 않았다.

따라서 기존 diagnostic 이미지는 정성적 참고용으로도 조심해야 하며, 물리 해석용으로 사용하지 않는 것이 좋다.

## 수정 방향

새 시각화 스크립트:

- `scripts/render_paint_booth_qc_safe_visualizations.py`

수정 원칙:

1. VTK point sampling 대신 **cell-centered U 원자료** 사용
2. 무한히 얇은 plane 대신 **finite-thickness slab** 선택
3. regular grid point sampling 대신 **cell-center 2D binning** 사용
4. 차체 side/top footprint를 명시적으로 **solid mask 처리**
5. reverse-flow는 `Uz > 0` 대신 유의미한 threshold `Uz > +0.02 m/s` 사용
6. raw vorticity는 QC fail로 보고, 대신 `|∇Uz|` shear proxy를 참고용으로 생성
7. 도메인 경계/차체 마스크 근처는 해석 주의 대상으로 명시

## 개선 이미지

생성된 QC-safe images:

- `docs/figures/26-05-02_paint_booth_l2_qc_safe_montage.png`
- `docs/figures/26-05-02_paint_booth_l2_qc_safe_uz_deviation_horizontal_stack.png`
- `docs/figures/26-05-02_paint_booth_l2_qc_safe_reverse_flow_mask_midplane.png`
- `docs/figures/26-05-02_paint_booth_l2_qc_safe_clipped_uz_midplane.png`
- `docs/figures/26-05-02_paint_booth_l2_qc_safe_shear_proxy_midplane.png`

## 개선 이미지 QC 판정

개선된 montage는 기존보다 artifact가 줄었고, 정성적 diagnostic 목적으로는 사용 가능하다.

개선점:

- 차체 solid 영역이 명확히 마스킹됨
- row-wise sampling artifact가 크게 감소
- reverse-flow는 thresholded magenta overlay로 표시되어 이전보다 명확함
- work-zone clipped `Uz`는 차체 주변 downwash 구조를 상대적으로 잘 보여줌
- raw vorticity 대신 `|∇Uz|` proxy를 사용해 artifact를 줄임

남은 한계:

1. 도메인 외곽과 경계부는 여전히 artifact 가능성이 높음
2. 차체 마스크 경계 바로 근처는 binning/부분 점유 cell 영향이 있음
3. reverse-flow patch는 매우 국소적이며, 지면/차체 경계 근처라 물리적 역류로 확정하기 어렵다
4. shear proxy는 vorticity가 아니므로 정량 해석에는 부적합
5. binning 기반이므로 blocky/pixelated appearance는 남아 있음

## 현재 해석 가능 범위

안전하게 해석 가능한 것:

- 전체적으로 downdraft가 지배적이라는 점
- 필터 아래 평균성분을 제거하면 residual nonuniformity가 보인다는 점
- 차체 높이로 내려오면서 차체 footprint 주변에서 `Uz` 편차가 커진다는 점
- work-zone 중간 단면에서 큰 규모의 연속 reverse-flow bubble은 명확하지 않다는 점
- 차체 주변 전단/gradient가 클 가능성이 있는 위치를 `|∇Uz|` proxy로 확인할 수 있다는 점

조심해야 할 것:

- 지면/차체 경계 근처 작은 reverse-flow patch의 물리적 의미
- 도메인 외곽 stripe/band 구조
- shear proxy의 정량적 크기
- raw vorticity 기반 물리 해석

## 다음 개선 제안

1. ParaView에서 직접 cell data slice와 solid clipping을 검증
2. y-slab 두께 민감도 비교: ±0.02, ±0.04, ±0.08 m
3. bin size 민감도 비교
4. reverse-flow threshold 민감도: `Uz > 0.01, 0.02, 0.05 m/s`
5. 차체 주변 local cutout view 생성
6. floor exhaust/side-wall 근처는 별도 mask 또는 crop 적용
