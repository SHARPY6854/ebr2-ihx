# ebr2-ihx

SMR 안전해석 팀프로젝트. EBR-II 1차계통 1D 열수력 코드 중 **IHX(중간 열교환기)** 파트 작업 저장소.

- 대상 시험: SHRT-17 (보호 유량상실, PLOF) → SHRT-45R (무보호 유량상실, ULOF)
- 기준 자료: ANL-ARC-226 Rev.1, *Benchmark Specifications and Data Requirements for EBR-II Shutdown Heat Removal Tests SHRT-17 and SHRT-45R* (Sumner and Wei, 2012)

## 폴더 구성

| 경로 | 내용 |
|---|---|
| `01_프로젝트_컨텍스트.md` | 지금까지의 결정 사항, 물리 배경, 검증된 수치, 사양서 오류, IHX 상세. **먼저 읽을 것** |
| `02_IHX_작업지시_프롬프트.md` | AI(Gemini 등)에 작업을 시킬 때 쓰는 프롬프트 |
| `경계조건_CSV/` | 사양서 표를 SI로 변환한 시간별 경계조건 8개. 원본 표 이미지와 한 줄씩 대조 완료 |
| `스크립트/` | IHX 형상(`ihx_geometry.json`), 경계조건 추출(`ihx_bc.py`), 변환 및 검산(`calc.py`) |
| `보고서/` | 노심, 펌프, 루프 운동량 자료조사 보고서 (docx, pdf) |

### 경계조건 CSV

| 파일 | 원본 (사양서) |
|---|---|
| `bc_ihx_secondary_flow_SHRT17.csv`, `bc_ihx_secondary_temp_SHRT17.csv` | 표 3.1.4.2.1~2, p.21 |
| `bc_ihx_secondary_flow_SHRT45R.csv`, `bc_ihx_secondary_temp_SHRT45R.csv` | 표 3.2.4.2.1 p.36, 3.2.4.2.2 p.37 |
| `bc_pump_speed_SHRT17.csv`, `bc_pump_speed_SHRT45R.csv` | 표 3.1.4.1.1 p.14, 3.2.4.1.1 p.29 |
| `bc_power_SHRT17.csv` | 표 3.1.4.2.1, p.19 |
| `bc_em_pump_current_SHRT45R.csv` | 표 3.2.4.1.2, p.31 |

유량은 gpm 기준온도(1차측 800°F, 중간 루프 582°F)의 나트륨 밀도(SAS4A/SASSYS-1 식 12.13-5)로 질량유량 변환.

## 이 저장소에 없는 자료

공개 저장소라서 아래는 올리지 않았다. 팀 공유 드라이브나 원 출처에서 받을 것.

- 사양서 ANL-ARC-226 Rev.1 원본 PDF: https://publications.anl.gov/anlpubs/2012/06/73647.pdf
- NURETH-16 Phase 2 벤치마크 논문 (2015)
- 팀 내부 자료 (보고서 계획, 회의 결과, 담당별 자료조사)

## 스크립트 실행

```bash
python 스크립트/ihx_geometry.py                 # IHX 형상 파생값 → ihx_geometry.json
python 스크립트/ihx_bc.py <사양서 PDF 경로>        # IHX 2차측 경계조건 CSV 재생성 (pdftotext 필요)
python 스크립트/calc.py                          # 단위 변환, 검산
```

## 주의 사항

- SHRT-45R IHX 2차측 초기 유량이 사양서 안에서 두 값으로 나온다 (유량표 5079.5 gpm, 초기조건 표 5474.5 gpm). 유량표 헤더의 기준온도도 "XX°F"로 빠져 있어 582°F로 가정했다
- 펌프 정격 수두는 358,530 Pa를 쓴다 (37.7 m와 불일치, 실측 토출압과는 358,530 Pa가 일치)
- 실측 데이터는 사양서에 없다. Phase 2 논문 그림 디지타이징 또는 IAEA-TECDOC-1819 필요
