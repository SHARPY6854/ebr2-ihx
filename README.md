# ebr2-ihx

SMR 안전해석 팀프로젝트. EBR-II 1차계통 1D 열수력 코드 중 **IHX(중간 열교환기)** 파트 작업 저장소.

- 대상 시험: SHRT-17 (보호 유량상실, PLOF) → SHRT-45R (무보호 유량상실, ULOF)
- 기준 자료: ANL-ARC-226 Rev.1, *Benchmark Specifications and Data Requirements for EBR-II Shutdown Heat Removal Tests SHRT-17 and SHRT-45R* (2012)

---

## 진행 현황

| 단계 | 내용 | 상태 | 결과물 |
|---|---|---|---|
| **1단계** | IHX 1D 모델 설계서 및 노달라이제이션 작성 | 완료 | `03_IHX_모델_설계서.md`, `도면/ihx_nodalization.svg` |
| **2단계** | 파이썬 모듈 `ihx.py` 구현 및 단위 검산 | 완료 | `ihx.py`, `04_IHX_코드_사용법_및_검산.md` |
| **3단계** | SHRT-17 정상상태 점검 (약 57 MW, 온도상승 약 145 K) | 진행 예정 | - |
| **4단계** | SHRT-17 단독 과도 계산 (0~900 s) | 진행 예정 | - |
| **5단계** | SHRT-45R 2차측 초기 유량 불일치 판단 | 진행 예정 | - |

---

## 파일 및 폴더 구성

| 경로 | 내용 |
|---|---|
| `01_프로젝트_컨텍스트.md` | 프로젝트 결정 사항, 물리 배경, 검증된 수치, 사양서 오류 정리. **먼저 읽을 것** |
| `02_IHX_작업지시_프롬프트.md` | 단계별 작업 지시서 |
| `03_IHX_모델_설계서.md` | **IHX 20노드 대향류 모델 설계서** (노달라이제이션 그림, 지배방정식, 물성식, 상관식, 수치해법) |
| `04_IHX_코드_사용법_및_검산.md` | **`ihx.py` 입출력 설명, 타 파트 연결 예시, 단위 검산 결과** |
| `ihx.py` | **IHX 1D 정상상태 및 과도 해석 파이썬 모듈** (외부 라이브러리 불필요) |
| `도면/ihx_nodalization.svg` | IHX 20노드 대향류 노달라이제이션 도식 이미지 |
| `경계조건_CSV/` | 사양서 표를 SI 단위로 변환한 시간별 경계조건 8개 (원본 표 대조 완료) |
| `스크립트/` | IHX 형상(`ihx_geometry.json`), 경계조건 추출(`ihx_bc.py`), 변환 및 검산(`calc.py`) |
| `보고서/` | 노심, 펌프, 루프 운동량 자료조사 보고서 (docx, pdf) |

### 경계조건 CSV 목록

| 파일 | 원본 출처 (사양서) |
|---|---|
| `bc_ihx_secondary_flow_SHRT17.csv`, `bc_ihx_secondary_temp_SHRT17.csv` | 표 3.1.4.2.1~2, p.21 |
| `bc_ihx_secondary_flow_SHRT45R.csv`, `bc_ihx_secondary_temp_SHRT45R.csv` | 표 3.2.4.2.1 p.36, 3.2.4.2.2 p.37 |
| `bc_pump_speed_SHRT17.csv`, `bc_pump_speed_SHRT45R.csv` | 표 3.1.4.1.1 p.14, 3.2.4.1.1 p.29 |
| `bc_power_SHRT17.csv` | 표 3.1.4.2.1, p.19 |
| `bc_em_pump_current_SHRT45R.csv` | 표 3.2.4.1.2, p.31 |

유량은 gpm 기준온도(1차측 800°F, 중간 루프 582°F)의 나트륨 밀도(SAS4A/SASSYS-1 식 12.13-5)로 질량유량 변환.

---

## 바로 실행해 보기

```bash
python3 ihx.py                          # IHX 모듈 자체 단위 확인 및 물리 검산 실행
python3 스크립트/ihx_geometry.py         # IHX 형상 파생값 계산
python3 스크립트/calc.py                 # 전체 단위 변환 및 펌프 수두 검산
```

---

## 주의 사항 및 참고 출처

- 공개 저장소라서 원본 사양서 PDF(`ANL-ARC-226 Rev.1`), 벤치마크 논문(`NURETH-16`, `IAEA-TECDOC-1819`), 팀 내부 문서는 제외되어 있다.
- SHRT-45R IHX 2차측 초기 유량이 사양서 안에서 두 값으로 나온다 (유량표 5079.5 gpm, 초기조건 표 5474.5 gpm). 5단계에서 열수지로 판정 예정.
- 펌프 정격 수두는 실측 토출압과 일치하는 358,530 Pa를 쓴다 (37.7 m 표기값과 불일치).
