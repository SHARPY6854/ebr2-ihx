# `ihx.py` 코드 설명 및 단위 검산 결과 (2단계)

`ihx.py`는 [03_IHX_모델_설계서.md](03_IHX_모델_설계서.md)를 파이썬 코드로 구현한 모듈이다. `numpy` 같은 외부 패키지 설치 없이 **파이썬 기본 내장 기능만으로 바로 실행**된다.

---

## 1. 실행 방법

```bash
python3 ihx.py
```

위 명령어를 치면 기하 형상, 물성치, 원통 벽 열저항, 등온 테스트, 과도 에너지 보존, CSV 경계조건 보간까지 6가지 자동 검산이 실행된다.

---

## 2. 다른 파트와 연결하는 법 (사용 예시)

Z-pipe, 저온풀, 루프 운동량 파트 코드를 합칠 때 아래 3개 함수만 부르면 된다. **모든 단위는 SI 단위(m, kg, s, K, W, Pa)** 다.

```python
from ihx import create_ihx, load_secondary_bc, solve_ihx_steady, step_ihx

# 1) IHX 생성 (기본 20노드) 및 2차측 CSV 경계조건 불러오기
ihx = create_ihx(n_nodes=20)
bc_sec = load_secondary_bc("SHRT17")  # 시각 t(초)를 넣으면 (T2_in_K, mdot2_kg_s) 반환

# 2) 0초 정상상태 초기화
T2_in_0, mdot2_0 = bc_sec(0.0)
init_res = solve_ihx_steady(ihx, T1_in_K=723.5, mdot1_kg_s=452.0, T2_in_K=T2_in_0, mdot2_kg_s=mdot2_0)

# 3) 시간 루프에서 매 dt초마다 1스텝씩 진행
dt = 0.1
T2_in, mdot2 = bc_sec(0.1)
res = step_ihx(ihx, dt_s=dt, T1_in_K=723.0, mdot1_kg_s=440.0, T2_in_K=T2_in, mdot2_kg_s=mdot2)

print("저온풀로 나가는 1차 출구 온도 [K]:", res["T1_out_K"])
print("2차측 출구 온도 [K]:", res["T2_out_K"])
print("열교환량 [MW]:", res["Q_W"] / 1e6)
```

---

## 3. 입출력 변수 정리

### 입력값
| 변수 | 단위 | 설명 | 어디서 받나 |
|---|---|---|---|
| `T1_in_K` | $\text{K}$ | 1차측 IHX 입구 온도 | Z-pipe 출구 (단독 검증 시 실측/초기값) |
| `mdot1_kg_s` | $\text{kg/s}$ | 1차측 질량유량 | 루프 운동량 방정식 |
| `T2_in_K` | $\text{K}$ | 2차측 IHX 입구 온도 | `load_secondary_bc` (사양서 CSV 자동 보간) |
| `mdot2_kg_s` | $\text{kg/s}$ | 2차측 질량유량 | `load_secondary_bc` (사양서 CSV 자동 보간) |

### 출력값 (`res` 딕셔너리)
| 키 이름 | 단위 | 설명 | 어디로 넘기나 |
|---|---|---|---|
| `T1_out_K` | $\text{K}$ | 1차측 출구 온도 | 저온풀 (Cold Pool) 입구 |
| `T2_out_K` | $\text{K}$ | 2차측 출구 온도 | 벤치마크 비교 결과 |
| `Q_W` | $\text{W}$ | 튜브 벽에서 2차측으로 전달된 총 열량 | 열수지 확인 |
| `T1_K`, `Tw_K`, `T2_K` | $\text{K}$ (길이 20 리스트) | 위(`0`)에서 아래(`19`)까지 층별 온도 분포 | 내부 상태 확인 및 그래프 |
| `rho1_kg_m3`, `z_nodes_m`, `dz_m` | $\text{kg/m}^3$, $\text{m}$, $\text{m}$ | 층별 1차측 나트륨 밀도와 높이 | 루프 운동량 부력항 계산 |
| `L_over_A_m_inv` | $\text{m}^{-1}$ | IHX 1차측 관성 계수 ($4.9373\text{ m}^{-1}$) | 루프 운동량 $\sum (L/A)$ |
| `dp_loss_Pa`, `dp_grav_Pa` | $\text{Pa}$ | 1차측 마찰+형상 압력손실, 중력 수두 | 루프 운동량 압력수지 |

*(참고: 만약 팀 전체 코드가 `scipy` 등 외부 미분방정식 솔버를 쓴다면 `ihx.derivatives(T1_in, mdot1, T2_in, mdot2)`로 미분값 `dT/dt`만 바로 꺼낼 수도 있다.)*

---

## 4. 단위 및 정합성 자동 검산 결과 (`python3 ihx.py` 실측치)

| 검산 항목 | 확인 내용 | 결과 수치 | 판정 |
|---|---|---|---|
| **기하 치수** | `ihx_geometry.json` 로드 및 면적, 수력직경, $L/A$ 확인 | $A_1=0.63536\text{ m}^2$, $A_2=0.42674\text{ m}^2$, $L/A=4.9373\text{ m}^{-1}$ | 통과 |
| **기준 밀도 대조** | `calc_results.json`의 $800^\circ\text{F}$, $582^\circ\text{F}$ 나트륨 밀도와 비교 | $849.13\text{ kg/m}^3$, $876.75\text{ kg/m}^3$ (오차 $10^{-12}$ 미만) | 통과 |
| **물성 단위 확인** | $400^\circ\text{C}$($673.15\text{ K}$) 나트륨 및 304SS 물성치 범위 점검 | $c_{p,\text{Na}}=1273.08\text{ J/(kg K)}$, $k_{\text{Na}}=72.86\text{ W/(m K)}$, $k_{304}=19.01\text{ W/(m K)}$ | 통과 |
| **원통 벽 열저항** | 반벽 분할 저항 합($G_{1w}^{-1} + G_{w2}^{-1}$)과 원통 벽 이론 엄밀해 비교 | 상대 오차 $1.58 \times 10^{-16}$ | 통과 |
| **등온 영점 테스트** | 양쪽 입구를 똑같이 $650\text{ K}$로 넣었을 때 온도 변화와 열교환량이 0인지 확인 | 최대 온도 편차 $4.21 \times 10^{-12}\text{ K}$, 열량 $1.95 \times 10^{-7}\text{ W}$ | 통과 |
| **과도 에너지 보존** | 과도 1스텝($\Delta t = 0.1\text{ s}$)에서 `[저장 에너지 증가율 - 순 유입 열량]` 오차 | 잔차 $-9.02 \times 10^{-7}\text{ W}$ (57 MW 대비 사실상 0) | 통과 |
| **CSV 보간 확인** | SHRT-17 2차측 CSV $0\text{ s}$, $900\text{ s}$ 로드 및 단위 변환 확인 | $0\text{ s}$: $301.12^\circ\text{C}$, $310.84\text{ kg/s}$ / $900\text{ s}$: $292.13^\circ\text{C}$, $17.20\text{ kg/s}$ | 통과 |
