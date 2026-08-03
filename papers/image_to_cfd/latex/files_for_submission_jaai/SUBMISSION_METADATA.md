# 투고용 메타데이터 — Journal of Korean Applied Artificial Intelligence (JAAI)

- **저널**: 한국실용인공지능학회지 / *Journal of Korean Applied Artificial Intelligence* (약어: *J. Korean Appl. Artif. Intell.*)
- **E-ISSN**: 3091-5902
- **원고 종류**: 정규 연구논문 (Research Article)
- **원고 파일**: `Image-to-CFD_JAAI.docx`
- **본문 언어**: 영문 (국문 제목·저자·소속 병기)
- **분량**: 22쪽 / 그림 17, 표 8, 수식 14, 참고문헌 33
- **작성일**: 2026-07-13

---

## 제목

**국문**
비전-언어 모델 추상화를 이용한 단일 2D 이미지 기반 실내 전산유체역학 자동 해석

**영문**
Automated Indoor Computational Fluid Dynamics Analysis from a Single 2D Image via Vision-Language Model Abstraction

**러닝 타이틀 (short title)**
Image-to-CFD via VLM geometric abstraction

---

## 저자

단독저자 (single author)

| 항목 | 국문 | 영문 |
|------|------|------|
| 성명 | 주영환 | Younghwan Joo |
| 소속 1 | 한국에너지기술연구원 에너지효율연구본부 | Energy Efficiency Research Division, Korea Institute of Energy Research |
| 소속 2 | 과학기술연합대학원대학교 에너지공학 | Energy Engineering, University of Science and Technology |

- **교신저자**: 주영환 (Younghwan Joo)
- **E-mail**: yhjoo@kier.re.kr
- **주소**: 대전광역시 유성구 가정로 152, 한국에너지기술연구원 (34129)

---

## 초록 (Abstract)

### 영문 (본문 수록본, 150 words)

Computational fluid dynamics (CFD) simulation of indoor environments requires three-dimensional geometric models that take hours to days of manual effort to build. This paper presents a framework that converts a single 2D image of an indoor space into a steady-state CFD solution using a vision-language model (VLM) as the geometric abstraction engine. The VLM (Gemini 3.1 Pro) extracts a structured 3D scene description, which is scale-calibrated from one reference dimension and automatically meshed and solved with OpenFOAM. A rule-based benchmark of 20 indoor geometries and five input view types (100 cases) provides independently computed reference solutions with no VLM involvement in the ground truth. Of 100 cases, 97 converge, achieving a mean structural score of 0.781 and a mean CFD agreement score of 0.477. Floor-plan inputs perform best (structural 0.884, CFD 0.572). Three failure modes are identified: composite-room collapse, obstacle hallucination, and a structure-versus-fidelity gap. The pipeline, benchmark, and evaluation code are publicly available.

### 국문 (투고 시스템 입력용)

실내 환경의 전산유체역학(CFD) 해석은 3차원 형상 모델을 필요로 하며, 이를 수작업으로 구축하는 데 수 시간에서 수일이 소요된다. 본 연구는 비전-언어 모델(VLM)을 기하 추상화 엔진으로 활용하여 실내 공간의 단일 2D 이미지로부터 정상상태 CFD 해를 자동으로 도출하는 프레임워크를 제시한다. VLM(Gemini 3.1 Pro)이 구조화된 3차원 장면 기술을 추출하면, 이를 하나의 기준 치수로 스케일 보정한 뒤 자동으로 격자를 생성하여 OpenFOAM으로 해석한다. 규칙 기반으로 생성한 20개의 실내 형상과 5종의 입력 뷰(총 100 케이스)로 구성된 벤치마크가 VLM이 개입하지 않은 독립적인 기준 해를 제공한다. 100개 케이스 중 97개가 수렴하였으며, 평균 구조 점수 0.781, 평균 CFD 일치도 점수 0.477을 달성하였다. 평면도 입력이 가장 우수한 성능을 보였다(구조 0.884, CFD 0.572). 복합 형상 붕괴, 장애물 환각, 구조-정확도 괴리의 세 가지 실패 모드를 규명하였다. 파이프라인, 벤치마크 및 평가 코드는 공개되어 있다.

---

## 키워드 (Key Words)

**영문 (본문 수록본)**
Indoor computational fluid dynamics, Vision-language model, Geometric abstraction, Automated simulation, Benchmark

**국문 (투고 시스템 입력용)**
실내 전산유체역학, 비전-언어 모델, 기하 추상화, 자동화 시뮬레이션, 벤치마크

---

## 후기 / 사사 (Acknowledgement)

> 본 문안은 현재 docx 본문에 미포함 상태입니다. 블라인드 심사가 아닌 경우 참고문헌 앞에 삽입하시면 됩니다.

본 연구는 한국에너지기술연구원의 주요사업의 지원으로 수행되었음. (No. C6-2419-63)

---

## 연구윤리 및 공개 사항

- **이해상충 (Declaration of competing interest)**: 저자는 본 논문의 내용에 영향을 줄 수 있는 경쟁적 재정적 이해관계나 개인적 관계가 없음을 선언한다. *(docx 본문 수록)*
- **데이터 가용성 (Data availability)**: 벤치마크 데이터셋, 평가 스크립트, 프레임워크 코드는 논문 게재 승인 시 공개 예정. *(docx 본문 수록)*
- **데이터 저장소 (figshare)**: https://doi.org/10.6084/m9.figshare.31866127
- **생성형 AI 사용 고지**: 연구 대상 자체가 VLM(Gemini 3.1 Pro)을 형상 추출 엔진으로 사용하는 프레임워크이며, 해당 사용은 본문 2장(Methodology)에 명시되어 있음.

---

## 원고 구성

| 장 | 제목 |
|----|------|
| 1 | Introduction (1.1 Motivation / 1.2 Related work / 1.3 Contributions) |
| 2 | Methodology (2.1 Framework overview / 2.2 Governing equations and CFD setup / 2.3 Evaluation metrics / 2.4 Solver validation / 2.5 Reference path and data-leakage prevention) |
| 3 | Benchmark Dataset (3.1 Geometry design / 3.2 Multi-view rendering protocol / 3.3 Scale calibration) |
| 4 | Results (4.1 Aggregate performance / 4.2 Input view type / 4.3 Geometric complexity / 4.4 Solver robustness / 4.5 VLM repeatability) |
| 5 | Discussion (5.1–5.3 Failure modes / 5.4 Scale calibration / 5.5 Preset matching / 5.6 Computational cost / 5.7 Limitations) |
| 6 | Application to Architectural Floor Plans |
| 7 | Conclusion |
| — | Appendix A. VLM Prompt Template / Appendix B. IoU Threshold Sensitivity |

**핵심 정량 결과**: 구조 점수 0.781 ± 0.150, CFD 일치도 0.477 ± 0.158 (n = 97), 수렴 97/100, 최우수 뷰 = 평면도(floorplan)

---

## 서식 준수 체크리스트

- [x] JAAI 제출용 워드 템플릿 양식 적용 (A4, 상단 마스트헤드, 2단 조판)
- [x] 국문·영문 제목 병기, 저자·소속·교신저자 이메일 표기
- [x] 영문 초록 150단어 이내 + Key Words
- [x] 그림·표 캡션 영문 작성 (그림 캡션 하단, 표 캡션 상단)
- [x] 그림·표 번호를 본문에서 모두 인용
- [x] 참고문헌 영문 작성, 본문 내 상첨자 인용, JAAI 서식 `(n) 저자, 연도, "제목," 저널, Vol., No., pp.`
- [x] 수식은 Word 네이티브 수식(편집 가능), 우측 정렬 번호 (1)–(14)
- [x] 6쪽 이상 (22쪽)
- [ ] 투고 시스템에 국문 초록·국문 키워드 입력 (위 내용 사용)
- [ ] 후기(사사) 문안 본문 삽입 여부 결정

---

## 폴더 파일 목록

| 파일 | 설명 |
|------|------|
| `Image-to-CFD_JAAI.docx` | 투고 원고 (본 서식 적용 완료) |
| `SUBMISSION_METADATA.md` | 본 문서 |
| `제출용 논문 템플릿.docx` | 학회 제공 템플릿 (저자 정보 미포함) |
| `참고용 논문 템플릿.docx` | 학회 제공 템플릿 (저자 정보 포함) |
