# Indoor CFD Auto-Stabilization Strategy

## 목적
Gemini가 생성한 실내 형상으로부터 자동 생성된 OpenFOAM 케이스는 자주 발산한다. 따라서 파이프라인은 단순 실행이 아니라, **실패 감지 -> 완화 -> 재시도** 루프를 내장해야 한다.

## 현재 파이프라인
1. Gemini scene JSON 생성
2. scene validator 통과
3. Gmsh geometry/mesh 생성
4. OpenFOAM case scaffold 생성
5. `gmshToFoam` / `checkMesh` / `simpleFoam`
6. 결과 시각화

## 왜 안정화 루프가 필요한가
자동 생성 형상은 아래 문제를 자주 가진다.
- opening 주변의 급격한 형상 변화
- obstacle edge 주변 품질 저하
- 높은 non-orthogonality
- wall/patch 분류 누락
- turbulence 초기조건 과도
- inlet velocity 대비 너무 공격적인 relaxation

따라서 발산은 예외가 아니라 정상적인 failure mode로 취급해야 한다.

## failure 신호
다음은 자동 감지해야 한다.
- `FOAM FATAL ERROR`
- `Floating point exception`
- residual이 반복적으로 증가
- continuity error 급증
- `bounding omega`, `bounding k`, `bounding U`가 반복적으로 심각
- `defaultFaces`가 과도하게 많음
- `checkMesh`에서 severe non-orth faces가 많음

## stabilization ladder
자동 재시도는 아래 순서를 기본으로 한다.

### Level 0 — aggressive / baseline
- inlet velocity: user/default
- internal field U: inlet과 동일
- relaxation: p=0.3, U/k/omega=0.7
- nonOrthogonalCorrectors=1

### Level 1 — conservative
- inlet velocity scale: 0.5x 또는 0.2x
- internal field U: `(0 0 0)`
- relaxation: p=0.2, U/k/omega=0.3
- nonOrthogonalCorrectors=3
- defaultFaces는 `empty` 대신 보수적으로 `wall/noSlip` 또는 `calculated` 처리

### Level 2 — robust
- inlet velocity scale: 0.1x
- 더 작은 turbulence intensity
- `k`, `omega` 초기값 추가 완화
- solver relTol 강화
- 필요 시 더 많은 pressure correction

### Level 3 — mesh-aware retry
- opening 주변 local refinement
- obstacle edge refinement
- mesh size 감소
- 필요 시 geometry 단순화

## 권장 자동 정책
1. baseline 시도
2. 실패 시 conservative preset 적용
3. 여전히 실패 시 robust preset 적용
4. 여전히 실패 시 mesh-aware retry
5. 최종 실패 시 로그/원인 리포트 생성

## 현재 검증된 사실
- baseline 성격의 설정은 발산 가능
- conservative starter 설정은 동일 geometry에서 수렴 완료 사례 확인

## framework 산출물
최종 산출물은 아래를 포함해야 한다.
- scene JSON
- `.geo_unrolled`
- `.msh`
- OpenFOAM case
- solver log
- stabilization attempt summary
- 1x2 comparison image
  - 좌: Gemini scene geometry
  - 우: OpenFOAM flow field

## 다음 확장
- retry 전략을 JSON/YAML policy로 외부화
- flow quality score 기반 자동 승급/중단
- geometry validity score와 solver stability score의 연결
- case archive zip 자동 생성
