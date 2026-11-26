"""EBR-II 중간열교환기(IHX) 1D 대향류 과도 열수력 모듈.

기준 자료:
- 형상 및 경계조건: ANL-ARC-226 Rev.1, 5.2.7절 표 5.2.7.1 (p.128), 표 3.1.4.2.1~2 (p.21), 표 3.2.4.2.1~2 (p.36~37)
- 나트륨 물성: SAS4A/SASSYS-1 매뉴얼 (ANL/NE-16/19, Ch.12, 식 12.13-5 외)
- 304 스테인리스강 물성: ANL-75-55 (SI 환산)
- 열전달 상관식: 쉘측 Mikityuk (2009) / Ushakov (1977) / Schad-Modified (1976), 튜브측 Seban-Shimazaki (1951)

모든 내부 계산과 입출력은 SI 단위(m, kg, s, K, W, Pa)를 사용한다.
외부 라이브러리 없이 파이썬 표준 라이브러리만으로 동작한다.
"""

from __future__ import annotations

import csv
import json
import math
import os
from typing import Callable, Dict, List, Optional, Tuple

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_GEOM_PATH = os.path.join(BASE_DIR, "스크립트", "ihx_geometry.json")
DEFAULT_BC_DIR = os.path.join(BASE_DIR, "경계조건_CSV")
GRAVITY = 9.80665  # 표준 중력가속도 [m/s^2]

# 격자 기본값. 1차 상류 차분은 NTU(약 9)가 큰 이 IHX에서 수치 확산이 커서
# N=20은 수렴값보다 열량을 약 2.6~2.9% 작게 준다. N=80이면 수렴값과 약 0.7% 이내.
DEFAULT_N_NODES: int = 80

# 1차측 실제 유동 구간 (ANL-ARC-226 Rev.1 5.2.7절 p.124 본문, 그림 5.2.7.3 p.127)
# - 1차측 입구관(Z-pipe 출구)은 노심 하단 기준 4.939 m 높이로 들어온다 (본문 명시).
#   그 위(약 5.369 m 상부 튜브시트까지)는 튜브시트로 막혀 1차측 하향 유동이 없다.
# - 1차측 출구는 하부 튜브시트 바로 위 환형 개구부 (그림 판독값 약 2.372 m, 확인 필요).
PRIMARY_INLET_ELEV_M: float = 4.939
PRIMARY_OUTLET_ELEV_M: float = 2.372


# =============================================================================
# 1. 물성치 함수 (입력 온도 T 단위: K, 반환값: 모두 SI 단위)
# =============================================================================

def rho_na(T: float) -> float:
    """액체 나트륨 밀도 [kg/m^3]. 출처: SAS4A/SASSYS-1 식 12.13-5 (371 < T < 2509 K)."""
    return 1004.23 - 0.21390 * T - 1.1046e-5 * T * T


def cp_na(T: float) -> float:
    """액체 나트륨 정압비열 [J/(kg K)]. 출처: SAS4A/SASSYS-1 Ch.12 (400 < T < 2200 K)."""
    T2 = T * T
    return 1972.0 - 2.4560 * T + 3.0921e-3 * T2 - 1.7191e-6 * T2 * T + 3.7782e-10 * T2 * T2


def k_na(T: float) -> float:
    """액체 나트륨 열전도도 [W/(m K)]. 출처: SAS4A/SASSYS-1 Ch.12 (371 < T < 1500 K)."""
    T2 = T * T
    return 110.45 - 6.5112e-2 * T + 1.5430e-5 * T2 - 2.4617e-9 * T2 * T


def mu_na(T: float) -> float:
    """액체 나트륨 점성계수 [Pa s]. 출처: SAS4A/SASSYS-1 Ch.12 (371 < T < 1200 K)."""
    invT = 1.0 / T
    return 3.6522e-5 + 0.16626 * invT - 45.6877 * invT * invT + 2.8733e4 * invT * invT * invT


def pr_na(T: float) -> float:
    """액체 나트륨 프란틀 수 [-]. Pr = cp * mu / k."""
    return cp_na(T) * mu_na(T) / k_na(T)


def rho_ss304(T: float) -> float:
    """Type 304 스테인리스강 밀도 [kg/m^3]. 출처: ANL-75-55."""
    return 7984.1 - 0.26506 * T - 1.1580e-4 * T * T


def cp_ss304(T: float) -> float:
    """Type 304 스테인리스강 비열 [J/(kg K)]. 출처: ANL-75-55."""
    return 469.44 + 0.13481 * T


def k_ss304(T: float) -> float:
    """Type 304 스테인리스강 열전도도 [W/(m K)]. 출처: ANL-75-55."""
    return 8.116 + 0.01618 * T


# =============================================================================
# 2. 열전달 및 마찰계수 상관식
# =============================================================================

def nusselt_shell(Pe: float, p_over_d: float, corr: str = "mikityuk") -> float:
    """쉘측(삼각 배열 봉다발 평행류) 누셀트 수 Nu [-].

    corr 옵션:
    - 'mikityuk': Mikityuk (2009), 1.1 <= P/D <= 1.95, 30 <= Pe <= 5000 (저유량 수렴값 포함)
    - 'ushakov' : Ushakov et al. (1977), 1.3 <= P/D <= 2.0, 1 <= Pe <= 4000
    - 'schad'   : Schad-Modified (1976), 1.1 <= P/D <= 1.5, Pe <= 1000 (Pe <= 150 상수식 포함)
    """
    Pe_pos = max(Pe, 0.0)
    x = p_over_d
    key = corr.lower()
    if key == "mikityuk":
        return 0.047 * (1.0 - math.exp(-3.8 * (x - 1.0))) * (Pe_pos ** 0.77 + 250.0)
    if key == "ushakov":
        return 7.55 * x - 20.0 * (x ** -13.0) + (0.041 / (x * x)) * (Pe_pos ** (0.56 + 0.19 * x))
    if key == "schad":
        geom = -16.15 + 24.96 * x - 8.55 * x * x
        return geom * (Pe_pos ** 0.3 if Pe_pos > 150.0 else 4.496)
    raise ValueError(f"지원하지 않는 쉘측 상관식: {corr}")


def nusselt_tube(Pe: float, corr: str = "seban") -> float:
    """튜브측(원형관 내부 액체금속) 누셀트 수 Nu [-].

    corr 옵션:
    - 'seban': Seban-Shimazaki (1951), Nu = 5.0 + 0.025 * Pe^0.8 (등벽온 경계조건)
    - 'lyon' : Lyon-Martinelli (1947), Nu = 7.0 + 0.025 * Pe^0.8 (등열유속 경계조건)
    """
    Pe_pos = max(Pe, 0.0)
    key = corr.lower()
    if key == "seban":
        return 5.0 + 0.025 * (Pe_pos ** 0.8)
    if key == "lyon":
        return 7.0 + 0.025 * (Pe_pos ** 0.8)
    raise ValueError(f"지원하지 않는 튜브측 상관식: {corr}")


def friction_factor_darcy(Re: float) -> float:
    """매끈한 유로 기준 Darcy 마찰계수 f [-]. 층류(Re<=2000), 난류(Re>=4000), 천이구간 선형보간."""
    if Re <= 1.0e-12:
        return 0.0
    if Re <= 2000.0:
        return 64.0 / Re
    if Re >= 4000.0:
        return 0.3164 / (Re ** 0.25)
    f_lam = 64.0 / 2000.0
    f_tur = 0.3164 / (4000.0 ** 0.25)
    w = (Re - 2000.0) / 2000.0
    return (1.0 - w) * f_lam + w * f_tur


# =============================================================================
# 3. 선형 연립방정식 해법 (띠행렬 가우스 소거법, 외부 라이브러리 불필요)
# =============================================================================

def _solve_banded(mat: List[List[float]], rhs: List[float], bw: int = 3) -> List[float]:
    """대각선 주변 대역폭 bw 이내인 정방행렬 A x = b 풀이 (부분 피벗팅 포함)."""
    n = len(rhs)
    a = [row[:] for row in mat]
    b = rhs[:]

    for k in range(n):
        row_end = min(n, k + bw + 1)
        col_end = min(n, k + 2 * bw + 1)

        # 부분 피벗팅
        piv = k
        max_val = abs(a[k][k])
        for i in range(k + 1, row_end):
            if abs(a[i][k]) > max_val:
                max_val = abs(a[i][k])
                piv = i
        if piv != k:
            a[k], a[piv] = a[piv], a[k]
            b[k], b[piv] = b[piv], b[k]

        diag = a[k][k]
        if abs(diag) < 1.0e-30:
            raise ZeroDivisionError(f"특이 행렬 발생 (행 {k})")

        for i in range(k + 1, row_end):
            factor = a[i][k] / diag
            if factor == 0.0:
                continue
            a[i][k] = 0.0
            for j in range(k + 1, col_end):
                a[i][j] -= factor * a[k][j]
            b[i] -= factor * b[k]

    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        col_end = min(n, i + 2 * bw + 1)
        s = b[i]
        for j in range(i + 1, col_end):
            s -= a[i][j] * x[j]
        x[i] = s / a[i][i]
    return x


# =============================================================================
# 4. 경계조건 CSV 로더 및 시간 보간기
# =============================================================================

def interp_linear(t: float, times: List[float], values: List[float]) -> float:
    """시간 배열 times에 대한 선형 보간 (범위 밖은 양 끝값 유지)."""
    if t <= times[0]:
        return values[0]
    if t >= times[-1]:
        return values[-1]
    lo, hi = 0, len(times) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if times[mid] <= t:
            lo = mid
        else:
            hi = mid
    dt_seg = times[hi] - times[lo]
    if dt_seg <= 0.0:
        return values[lo]
    w = (t - times[lo]) / dt_seg
    return values[lo] + w * (values[hi] - values[lo])


# SHRT-45R 2차측 초기 유량 보정 비율 (초기조건 표 5474.5 gpm / 과도 유량표 0초 5079.5 gpm)
SHRT45R_FLOW_SCALE: float = 5474.5 / 5079.5


def load_secondary_bc(
    test_name: str = "SHRT17",
    bc_dir: str = DEFAULT_BC_DIR,
    flow_scale: float = 1.0,
) -> Callable[[float], Tuple[float, float]]:
    """사양서 변환 CSV에서 2차측 입구 경계조건 보간 함수 bc_func(t_s) -> (T2_in_K, mdot2_kg_s) 생성.

    test_name : 'SHRT17' 또는 'SHRT45R'
    flow_scale: 2차측 질량유량 배율 (기본값 1.0, SHRT-45R 초기유량 정합 시 SHRT45R_FLOW_SCALE 사용)
    """
    tag = test_name.upper().replace("-", "")
    flow_path = os.path.join(bc_dir, f"bc_ihx_secondary_flow_{tag}.csv")
    temp_path = os.path.join(bc_dir, f"bc_ihx_secondary_temp_{tag}.csv")

    t_flow, mdot_flow = [], []
    with open(flow_path, encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            t_flow.append(float(row["time_s"]))
            mdot_flow.append(float(row["mdot_kg_s"]) * float(flow_scale))

    t_temp, T_temp = [], []
    with open(temp_path, encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            t_temp.append(float(row["time_s"]))
            T_temp.append(float(row["T_K"]))

    def bc_func(t_s: float) -> Tuple[float, float]:
        T2_in_K = interp_linear(t_s, t_temp, T_temp)
        mdot2_kg_s = interp_linear(t_s, t_flow, mdot_flow)
        return T2_in_K, mdot2_kg_s

    return bc_func


# =============================================================================
# 5. IHX 1D 대향류 과도/정상상태 모델
# =============================================================================

class IHXModel:
    """EBR-II IHX 축방향 N노드 대향류 모델.

    노드 인덱스 j = 0, 1, ..., N-1:
    - j = 0   : 최상단 (1차측 Z-pipe 유입단 / 2차측 상부 출구단)
    - j = N-1 : 최하단 (1차측 저온풀 배출단 / 2차측 하부 유입단)
    """

    def __init__(
        self,
        n_nodes: int = DEFAULT_N_NODES,
        geom_path: str = DEFAULT_GEOM_PATH,
        shell_corr: str = "mikityuk",
        tube_corr: str = "seban",
        use_secondary_plena: bool = False,
        k_form_shell: float = 0.0,
        T_init_K: float = 625.0,
        length_mode: str = "active",
    ) -> None:
        """length_mode:
        - 'active' (기본): 1차측이 실제로 튜브를 따라 흐르는 구간만 열전달 길이로 사용
          (1차 입구 4.939 m ~ 1차 출구 약 2.372 m, 길이 2.567 m)
        - 'full'         : 튜브 전체 길이 3.137 m 사용 (민감도 비교용)
        """
        with open(geom_path, encoding="utf-8") as f:
            raw = json.load(f)
        g = raw["given"]
        d = raw["derived"]

        self.N = int(n_nodes)
        self.shell_corr = shell_corr
        self.tube_corr = tube_corr
        self.use_secondary_plena = use_secondary_plena
        self.k_form_shell = float(k_form_shell)

        # 기하 치수 (모두 SI 단위: m, m^2, m^3)
        self.N_tubes: int = int(g["N_tubes"])
        self.Do: float = float(g["tube_OD_m"])
        self.Di: float = float(d["tube_ID_m"])
        self.wall_thk: float = float(g["tube_wall_m"])
        self.pitch: float = float(g["pitch_m"])
        self.P_over_D: float = float(d["P_over_D"])
        self.length_mode = length_mode.lower()
        if self.length_mode == "active":
            self.z_top: float = PRIMARY_INLET_ELEV_M
            self.z_bot: float = PRIMARY_OUTLET_ELEV_M
        elif self.length_mode == "full":
            self.z_top = float(g["fig_upper_tubesheet_elev_m"])
            self.z_bot = float(g["fig_lower_tubesheet_bottom_elev_m"])
        else:
            raise ValueError(f"지원하지 않는 length_mode: {length_mode}")
        self.L_tube: float = float(g["tube_length_m"])
        self.L: float = self.z_top - self.z_bot
        self.dz: float = self.L / self.N

        self.A1: float = float(d["shell_flow_area_m2_assume_center_pipe"])
        self.A2: float = float(d["tube_flow_area_m2"])
        self.Aw: float = self.N_tubes * 0.25 * math.pi * (self.Do ** 2 - self.Di ** 2)
        self.Dh1: float = float(d["shell_Dh_m"])
        self.Dh2: float = float(d["tube_Dh_m"])

        # 단일 노드 전열면적 및 체적
        self.Ao_node: float = self.N_tubes * math.pi * self.Do * self.dz
        self.Ai_node: float = self.N_tubes * math.pi * self.Di * self.dz
        self.V1_node: float = self.A1 * self.dz
        self.V2_node: float = self.A2 * self.dz
        self.Vw_node: float = self.Aw * self.dz

        # 원통형 반벽 열전도 저항 계산용 반경 (기하평균 r_m = sqrt(r_o * r_i))
        self.ro: float = 0.5 * self.Do
        self.ri: float = 0.5 * self.Di
        self.rm: float = math.sqrt(self.ro * self.ri)
        # R_w_half = ln(r_ratio) / (2 * pi * k_w * dz * N_tubes) = geom_factor / k_w
        self._rw_outer_geom: float = math.log(self.ro / self.rm) / (2.0 * math.pi * self.dz * self.N_tubes)
        self._rw_inner_geom: float = math.log(self.rm / self.ri) / (2.0 * math.pi * self.dz * self.N_tubes)

        # 높이 정보 (노심 유효하단 z = 0 기준 [m], 루프 부력 계산용). z_top, z_bot은 위에서 설정
        self.z_nodes: List[float] = [
            self.z_top - (j + 0.5) * self.dz for j in range(self.N)
        ]
        self.L_over_A: float = self.L / self.A1  # [m^-1]

        # 2차측 상하부 플레넘 체적 [m^3]
        self.V_pl_in: float = float(g["sec_inlet_plenum_m3"])
        self.V_pl_out: float = float(g["sec_outlet_plenum_m3"])

        # 상태변수 초기화 [K]
        self.T1: List[float] = [float(T_init_K)] * self.N
        self.Tw: List[float] = [float(T_init_K)] * self.N
        self.T2: List[float] = [float(T_init_K)] * self.N
        self.T_pl_in: float = float(T_init_K)
        self.T_pl_out: float = float(T_init_K)

    def _node_conductances(
        self, T1_j: float, Tw_j: float, T2_j: float, mdot1: float, mdot2: float
    ) -> Tuple[float, float, float, float]:
        """노드별 열전달 컨덕턴스 (G_1w, G_w2 [W/K]) 및 대류 열전달계수 (h1, h2 [W/(m^2 K)]) 계산."""
        # 1차측(쉘측)
        mu1, cp1, k1 = mu_na(T1_j), cp_na(T1_j), k_na(T1_j)
        Re1 = abs(mdot1) * self.Dh1 / (self.A1 * mu1)
        Pe1 = Re1 * (cp1 * mu1 / k1)
        Nu1 = nusselt_shell(Pe1, self.P_over_D, self.shell_corr)
        h1 = Nu1 * k1 / self.Dh1

        # 2차측(튜브측)
        mu2, cp2, k2 = mu_na(T2_j), cp_na(T2_j), k_na(T2_j)
        Re2 = abs(mdot2) * self.Dh2 / (self.A2 * mu2)
        Pe2 = Re2 * (cp2 * mu2 / k2)
        Nu2 = nusselt_tube(Pe2, self.tube_corr)
        h2 = Nu2 * k2 / self.Dh2

        # 304SS 튜브 벽 반경방향 전도 저항 결합
        kw = k_ss304(Tw_j)
        R_1w = (1.0 / (h1 * self.Ao_node)) + (self._rw_outer_geom / kw)
        R_w2 = (self._rw_inner_geom / kw) + (1.0 / (h2 * self.Ai_node))
        return 1.0 / R_1w, 1.0 / R_w2, h1, h2

    def _build_and_solve_linear_system(
        self,
        T1_in: float,
        mdot1: float,
        T2_bundle_in: float,
        mdot2: float,
        inv_dt: float,
        T1_ref: List[float],
        Tw_ref: List[float],
        T2_ref: List[float],
    ) -> Tuple[List[float], List[float], List[float], List[float], List[float]]:
        """반음해 후방 오일러(inv_dt = 1/dt) 또는 정상상태(inv_dt = 0) 3N 선형계 조립 및 풀이."""
        n_vars = 3 * self.N
        mat = [[0.0] * n_vars for _ in range(n_vars)]
        rhs = [0.0] * n_vars

        G1w_list = [0.0] * self.N
        Gw2_list = [0.0] * self.N

        m1_abs = abs(mdot1)
        m2_abs = abs(mdot2)

        for j in range(self.N):
            T1_j, Tw_j, T2_j = T1_ref[j], Tw_ref[j], T2_ref[j]
            G1w, Gw2, _, _ = self._node_conductances(T1_j, Tw_j, T2_j, mdot1, mdot2)
            G1w_list[j] = G1w
            Gw2_list[j] = Gw2

            cp1 = cp_na(T1_j)
            cpw = cp_ss304(Tw_j)
            cp2 = cp_na(T2_j)

            C1_dt = rho_na(T1_j) * self.V1_node * cp1 * inv_dt
            Cw_dt = rho_ss304(Tw_j) * self.Vw_node * cpw * inv_dt
            C2_dt = rho_na(T2_j) * self.V2_node * cp2 * inv_dt

            W1 = m1_abs * cp1
            W2 = m2_abs * cp2

            r1 = 3 * j
            rw = 3 * j + 1
            r2 = 3 * j + 2

            # (1) 1차측 에너지식 (정방향 mdot1 >= 0: 상단 j=0 -> 하단 j=N-1)
            mat[r1][r1] = C1_dt + W1 + G1w
            mat[r1][rw] = -G1w
            rhs[r1] = C1_dt * self.T1[j]
            if mdot1 >= 0.0:
                if j == 0:
                    rhs[r1] += W1 * T1_in
                else:
                    mat[r1][3 * (j - 1)] = -W1
            else:
                if j == self.N - 1:
                    rhs[r1] += W1 * T1_in
                else:
                    mat[r1][3 * (j + 1)] = -W1

            # (2) 304SS 튜브 벽 에너지식
            mat[rw][r1] = -G1w
            mat[rw][rw] = Cw_dt + G1w + Gw2
            mat[rw][r2] = -Gw2
            rhs[rw] = Cw_dt * self.Tw[j]

            # (3) 2차측 에너지식 (정방향 mdot2 >= 0: 하단 j=N-1 -> 상단 j=0 상향류)
            mat[r2][rw] = -Gw2
            mat[r2][r2] = C2_dt + W2 + Gw2
            rhs[r2] = C2_dt * self.T2[j]
            if mdot2 >= 0.0:
                if j == self.N - 1:
                    rhs[r2] += W2 * T2_bundle_in
                else:
                    mat[r2][3 * (j + 1) + 2] = -W2
            else:
                if j == 0:
                    rhs[r2] += W2 * T2_bundle_in
                else:
                    mat[r2][3 * (j - 1) + 2] = -W2

        sol = _solve_banded(mat, rhs, bw=3)
        T1_new = [sol[3 * j] for j in range(self.N)]
        Tw_new = [sol[3 * j + 1] for j in range(self.N)]
        T2_new = [sol[3 * j + 2] for j in range(self.N)]
        return T1_new, Tw_new, T2_new, G1w_list, Gw2_list

    def hydraulics(self, mdot1: float) -> Dict[str, float]:
        """루프 운동량 방정식 전달용 1차측 압력손실, 중력수두, 관성계수 계산."""
        rho_nodes = [rho_na(T) for T in self.T1]
        rho_mean = sum(rho_nodes) / self.N
        dp_grav = sum(r * GRAVITY * self.dz for r in rho_nodes)

        f_sum = 0.0
        for j in range(self.N):
            Re_j = abs(mdot1) * self.Dh1 / (self.A1 * mu_na(self.T1[j]))
            f_sum += friction_factor_darcy(Re_j) * (self.dz / self.Dh1)

        dyn_head_signed = mdot1 * abs(mdot1) / (2.0 * rho_mean * self.A1 * self.A1)
        dp_fric = f_sum * dyn_head_signed
        dp_form = self.k_form_shell * dyn_head_signed
        return {
            "L_over_A_m_inv": self.L_over_A,
            "dp_grav_Pa": dp_grav,
            "dp_fric_Pa": dp_fric,
            "dp_form_Pa": dp_form,
            "dp_loss_Pa": dp_fric + dp_form,
        }

    def _package_outputs(
        self,
        T1_in: float,
        mdot1: float,
        T2_in: float,
        mdot2: float,
        G1w_list: List[float],
        Gw2_list: List[float],
        energy_balance_residual_W: float = 0.0,
    ) -> Dict[str, object]:
        """타 파트 연동 및 검산용 표준 출력 딕셔너리 생성."""
        T1_out = self.T1[-1] if mdot1 >= 0.0 else self.T1[0]
        T2_bundle_out = self.T2[0] if mdot2 >= 0.0 else self.T2[-1]
        T2_out = self.T_pl_out if self.use_secondary_plena else T2_bundle_out

        # 열전달량 계산 [W]
        Q_1toW = sum(G1w_list[j] * (self.T1[j] - self.Tw[j]) for j in range(self.N))
        Q_Wto2 = sum(Gw2_list[j] * (self.Tw[j] - self.T2[j]) for j in range(self.N))

        # 노드별 대류 엔탈피 수송량 합산 [W]
        T2_bundle_in = self.T_pl_in if self.use_secondary_plena else T2_in
        Q_1_adv = 0.0
        Q_2_adv = 0.0
        for j in range(self.N):
            cp1_j = cp_na(self.T1[j])
            cp2_j = cp_na(self.T2[j])
            T1_up = T1_in if (mdot1 >= 0.0 and j == 0) else (
                self.T1[j - 1] if mdot1 >= 0.0 else (T1_in if j == self.N - 1 else self.T1[j + 1])
            )
            T2_up = T2_bundle_in if (mdot2 >= 0.0 and j == self.N - 1) else (
                self.T2[j + 1] if mdot2 >= 0.0 else (T2_bundle_in if j == 0 else self.T2[j - 1])
            )
            Q_1_adv += abs(mdot1) * cp1_j * (T1_up - self.T1[j])
            Q_2_adv += abs(mdot2) * cp2_j * (self.T2[j] - T2_up)

        hyd = self.hydraulics(mdot1)
        return {
            "T1_in_K": T1_in,
            "T1_out_K": T1_out,
            "T2_in_K": T2_in,
            "T2_out_K": T2_out,
            "T1_K": self.T1[:],
            "Tw_K": self.Tw[:],
            "T2_K": self.T2[:],
            "Q_W": Q_Wto2,
            "Q_1toW_W": Q_1toW,
            "Q_Wto2_W": Q_Wto2,
            "Q_1_adv_W": Q_1_adv,
            "Q_2_adv_W": Q_2_adv,
            "energy_residual_W": energy_balance_residual_W,
            "rho1_kg_m3": [rho_na(T) for T in self.T1],
            "z_nodes_m": self.z_nodes[:],
            "dz_m": self.dz,
            **hyd,
        }

    def solve_steady_state(
        self,
        T1_in: float,
        mdot1: float,
        T2_in: float,
        mdot2: float,
        tol: float = 1.0e-9,
        max_iter: int = 50,
    ) -> Dict[str, object]:
        """정상상태(d/dt = 0) 온도 분포를 고정점 반복으로 직접 계산."""
        T1_guess = [0.5 * (T1_in + T2_in)] * self.N
        Tw_guess = [0.5 * (T1_in + T2_in)] * self.N
        T2_guess = [0.5 * (T1_in + T2_in)] * self.N

        G1w_list = [0.0] * self.N
        Gw2_list = [0.0] * self.N

        for _ in range(max_iter):
            T1_new, Tw_new, T2_new, G1w_list, Gw2_list = self._build_and_solve_linear_system(
                T1_in=T1_in,
                mdot1=mdot1,
                T2_bundle_in=T2_in,
                mdot2=mdot2,
                inv_dt=0.0,
                T1_ref=T1_guess,
                Tw_ref=Tw_guess,
                T2_ref=T2_guess,
            )
            diff = max(
                max(abs(a - b) for a, b in zip(T1_new, T1_guess)),
                max(abs(a - b) for a, b in zip(Tw_new, Tw_guess)),
                max(abs(a - b) for a, b in zip(T2_new, T2_guess)),
            )
            T1_guess, Tw_guess, T2_guess = T1_new, Tw_new, T2_new
            if diff < tol:
                break

        self.T1 = T1_guess
        self.Tw = Tw_guess
        self.T2 = T2_guess
        self.T_pl_in = T2_in
        self.T_pl_out = self.T2[0] if mdot2 >= 0.0 else self.T2[-1]

        out = self._package_outputs(T1_in, mdot1, T2_in, mdot2, G1w_list, Gw2_list)
        # 정상상태에서 1차측 방출열과 2차측 흡수열 차이 [W]
        out["energy_residual_W"] = float(out["Q_1_adv_W"]) - float(out["Q_2_adv_W"])
        return out

    def step(
        self,
        dt: float,
        T1_in: float,
        mdot1: float,
        T2_in: float,
        mdot2: float,
    ) -> Dict[str, object]:
        """과도 상태 1 시간스텝(dt [s]) 적분 (반음해 후방 오일러법)."""
        if dt <= 0.0:
            raise ValueError("시간간격 dt는 양수여야 합니다.")

        T1_old = self.T1[:]
        Tw_old = self.Tw[:]
        T2_old = self.T2[:]

        # 2차측 입구 플레넘 혼합 지연 (옵션)
        if self.use_secondary_plena:
            M_in = rho_na(self.T_pl_in) * self.V_pl_in
            alpha_in = abs(mdot2) * dt / M_in
            self.T_pl_in = (self.T_pl_in + alpha_in * T2_in) / (1.0 + alpha_in)
            T2_bundle_in = self.T_pl_in
        else:
            self.T_pl_in = T2_in
            T2_bundle_in = T2_in

        inv_dt = 1.0 / dt
        T1_new, Tw_new, T2_new, G1w_list, Gw2_list = self._build_and_solve_linear_system(
            T1_in=T1_in,
            mdot1=mdot1,
            T2_bundle_in=T2_bundle_in,
            mdot2=mdot2,
            inv_dt=inv_dt,
            T1_ref=T1_old,
            Tw_ref=Tw_old,
            T2_ref=T2_old,
        )

        # 이산 에너지 보존 잔차 검산 [W]: dE_stored/dt - (Q_1_adv - Q_2_adv)
        dE_dt = 0.0
        Q1_adv = 0.0
        Q2_adv = 0.0
        for j in range(self.N):
            cp1 = cp_na(T1_old[j])
            cpw = cp_ss304(Tw_old[j])
            cp2 = cp_na(T2_old[j])
            C1 = rho_na(T1_old[j]) * self.V1_node * cp1
            Cw = rho_ss304(Tw_old[j]) * self.Vw_node * cpw
            C2 = rho_na(T2_old[j]) * self.V2_node * cp2
            dE_dt += (
                C1 * (T1_new[j] - T1_old[j])
                + Cw * (Tw_new[j] - Tw_old[j])
                + C2 * (T2_new[j] - T2_old[j])
            ) * inv_dt

            T1_up = T1_in if (mdot1 >= 0.0 and j == 0) else (
                T1_new[j - 1] if mdot1 >= 0.0 else (T1_in if j == self.N - 1 else T1_new[j + 1])
            )
            T2_up = T2_bundle_in if (mdot2 >= 0.0 and j == self.N - 1) else (
                T2_new[j + 1] if mdot2 >= 0.0 else (T2_bundle_in if j == 0 else T2_new[j - 1])
            )
            Q1_adv += abs(mdot1) * cp1 * (T1_up - T1_new[j])
            Q2_adv += abs(mdot2) * cp2 * (T2_new[j] - T2_up)

        residual_W = dE_dt - (Q1_adv - Q2_adv)

        self.T1 = T1_new
        self.Tw = Tw_new
        self.T2 = T2_new

        # 2차측 출구 플레넘 혼합 지연 (옵션)
        T2_bundle_out = self.T2[0] if mdot2 >= 0.0 else self.T2[-1]
        if self.use_secondary_plena:
            M_out = rho_na(self.T_pl_out) * self.V_pl_out
            alpha_out = abs(mdot2) * dt / M_out
            self.T_pl_out = (self.T_pl_out + alpha_out * T2_bundle_out) / (1.0 + alpha_out)
        else:
            self.T_pl_out = T2_bundle_out

        return self._package_outputs(
            T1_in, mdot1, T2_in, mdot2, G1w_list, Gw2_list, energy_balance_residual_W=residual_W
        )

    def derivatives(
        self, T1_in: float, mdot1: float, T2_in: float, mdot2: float
    ) -> Tuple[List[float], List[float], List[float]]:
        """외부 상미분방정식(ODE) 솔버 결합용 시간 미분항 (dT1/dt, dTw/dt, dT2/dt) [K/s] 반환."""
        dT1_dt = [0.0] * self.N
        dTw_dt = [0.0] * self.N
        dT2_dt = [0.0] * self.N
        T2_bundle_in = self.T_pl_in if self.use_secondary_plena else T2_in

        for j in range(self.N):
            T1_j, Tw_j, T2_j = self.T1[j], self.Tw[j], self.T2[j]
            G1w, Gw2, _, _ = self._node_conductances(T1_j, Tw_j, T2_j, mdot1, mdot2)
            cp1 = cp_na(T1_j)
            cpw = cp_ss304(Tw_j)
            cp2 = cp_na(T2_j)
            C1 = rho_na(T1_j) * self.V1_node * cp1
            Cw = rho_ss304(Tw_j) * self.Vw_node * cpw
            C2 = rho_na(T2_j) * self.V2_node * cp2

            T1_up = T1_in if (mdot1 >= 0.0 and j == 0) else (
                self.T1[j - 1] if mdot1 >= 0.0 else (T1_in if j == self.N - 1 else self.T1[j + 1])
            )
            T2_up = T2_bundle_in if (mdot2 >= 0.0 and j == self.N - 1) else (
                self.T2[j + 1] if mdot2 >= 0.0 else (T2_bundle_in if j == 0 else self.T2[j - 1])
            )

            q1w = G1w * (T1_j - Tw_j)
            qw2 = Gw2 * (Tw_j - T2_j)
            dT1_dt[j] = (abs(mdot1) * cp1 * (T1_up - T1_j) - q1w) / C1
            dTw_dt[j] = (q1w - qw2) / Cw
            dT2_dt[j] = (abs(mdot2) * cp2 * (T2_up - T2_j) + qw2) / C2

        return dT1_dt, dTw_dt, dT2_dt


# =============================================================================
# 6. 타 파트 연동용 간결 함수 인터페이스
# =============================================================================

def create_ihx(
    n_nodes: int = DEFAULT_N_NODES,
    geom_path: str = DEFAULT_GEOM_PATH,
    shell_corr: str = "mikityuk",
    tube_corr: str = "seban",
    use_secondary_plena: bool = False,
    k_form_shell: float = 0.0,
    length_mode: str = "active",
) -> IHXModel:
    """IHX 모델 인스턴스 생성."""
    return IHXModel(
        n_nodes=n_nodes,
        geom_path=geom_path,
        shell_corr=shell_corr,
        tube_corr=tube_corr,
        use_secondary_plena=use_secondary_plena,
        k_form_shell=k_form_shell,
        length_mode=length_mode,
    )


def effectiveness_ntu_check(ihx: IHXModel, T1_in: float, mdot1: float, T2_in: float, mdot2: float) -> Dict[str, float]:
    """현재 정상상태 해의 노드별 컨덕턴스로 UA를 구해 대향류 epsilon-NTU 해석해와 비교.

    격자가 충분히 촘촘하면 수치해 Q와 해석해 Q_eNTU가 일치해야 한다 (구현 검증용).
    """
    UA = 0.0
    for j in range(ihx.N):
        G1w, Gw2, _, _ = ihx._node_conductances(ihx.T1[j], ihx.Tw[j], ihx.T2[j], mdot1, mdot2)
        UA += 1.0 / (1.0 / G1w + 1.0 / Gw2)
    T1_out = ihx.T1[-1]
    T2_out = ihx.T2[0]
    C1 = abs(mdot1) * cp_na(0.5 * (T1_in + T1_out))
    C2 = abs(mdot2) * cp_na(0.5 * (T2_in + T2_out))
    Cmin, Cmax = min(C1, C2), max(C1, C2)
    Cr = Cmin / Cmax
    NTU = UA / Cmin
    if abs(1.0 - Cr) < 1.0e-9:
        eps = NTU / (1.0 + NTU)
    else:
        e = math.exp(-NTU * (1.0 - Cr))
        eps = (1.0 - e) / (1.0 - Cr * e)
    return {"UA_W_K": UA, "NTU": NTU, "Cr": Cr, "eps": eps, "Q_eNTU_W": eps * Cmin * (T1_in - T2_in)}


def solve_ihx_steady(
    ihx: IHXModel, T1_in_K: float, mdot1_kg_s: float, T2_in_K: float, mdot2_kg_s: float
) -> Dict[str, object]:
    """IHX 초기 정상상태 계산 후 결과 딕셔너리 반환."""
    return ihx.solve_steady_state(T1_in_K, mdot1_kg_s, T2_in_K, mdot2_kg_s)


def step_ihx(
    ihx: IHXModel, dt_s: float, T1_in_K: float, mdot1_kg_s: float, T2_in_K: float, mdot2_kg_s: float
) -> Dict[str, object]:
    """IHX 과도 1스텝(dt_s) 진행 후 결과 딕셔너리 반환."""
    return ihx.step(dt_s, T1_in_K, mdot1_kg_s, T2_in_K, mdot2_kg_s)


# =============================================================================
# 7. 모듈 자체 단위 확인 및 물리 정합성 검산 (직접 실행 시 동작)
# =============================================================================

def _run_self_checks() -> None:
    print("=" * 68)
    print("EBR-II IHX 모듈 (ihx.py) 단위 일관성 및 물리 정합성 검산")
    print("=" * 68)

    ihx = create_ihx()
    print(f"[1] 기하 형상 로드 완료 (N = {ihx.N}, length_mode = {ihx.length_mode}, L = {ihx.L:.3f} m, dz = {ihx.dz:.5f} m)")
    print(f"    1차측 유동면적 A1 = {ihx.A1:.5f} m^2, 수력직경 Dh1 = {ihx.Dh1*1e3:.3f} mm")
    print(f"    2차측 유동면적 A2 = {ihx.A2:.5f} m^2, 수력직경 Dh2 = {ihx.Dh2*1e3:.3f} mm")
    print(f"    전열면적 합계: 외측 = {ihx.Ao_node*ihx.N:.2f} m^2, 내측 = {ihx.Ai_node*ihx.N:.2f} m^2")
    print(f"    1차측 관성계수 L/A = {ihx.L_over_A:.4f} m^-1")
    print(f"    노드 중심 높이 범위: z_1 = {ihx.z_nodes[0]:.3f} m ~ z_{ihx.N} = {ihx.z_nodes[-1]:.3f} m")

    # 기준온도 밀도 검산 (팀 공통 calc_results.json 대조)
    T_800F = (800.0 - 32.0) * 5.0 / 9.0 + 273.15
    T_582F = (582.0 - 32.0) * 5.0 / 9.0 + 273.15
    rho_800 = rho_na(T_800F)
    rho_582 = rho_na(T_582F)
    assert abs(rho_800 - 849.1295097687316) < 1.0e-9, "800F 밀도 불일치"
    assert abs(rho_582 - 876.7455749408058) < 1.0e-9, "582F 밀도 불일치"
    T_ref = 673.15  # 400 C
    print(f"\n[2] 물성식 단위 검산 (기준온도 대조 통과)")
    print(f"    rho_Na(800F={T_800F:.2f} K) = {rho_800:.2f} kg/m^3, rho_Na(582F={T_582F:.2f} K) = {rho_582:.2f} kg/m^3")
    print(
        f"    400C(673.15 K) 나트륨: rho={rho_na(T_ref):.2f} kg/m^3, cp={cp_na(T_ref):.2f} J/(kg K), "
        f"k={k_na(T_ref):.2f} W/(m K), mu={mu_na(T_ref):.4e} Pa s, Pr={pr_na(T_ref):.5f}"
    )
    print(
        f"    400C(673.15 K) 304SS : rho={rho_ss304(T_ref):.2f} kg/m^3, cp={cp_ss304(T_ref):.2f} J/(kg K), "
        f"k={k_ss304(T_ref):.2f} W/(m K)"
    )

    # 원통 벽 전도 저항 엄밀해 대조
    G1w, Gw2, h1, h2 = ihx._node_conductances(680.0, 660.0, 640.0, 452.0, 310.6)
    R_sum = (1.0 / G1w) + (1.0 / Gw2)
    R_exact = (
        1.0 / (h1 * ihx.Ao_node)
        + math.log(ihx.Do / ihx.Di) / (2.0 * math.pi * k_ss304(660.0) * ihx.dz * ihx.N_tubes)
        + 1.0 / (h2 * ihx.Ai_node)
    )
    rel_err_R = abs(R_sum - R_exact) / R_exact
    assert rel_err_R < 1.0e-15
    print(f"\n[3] 열전달 계수 및 원통 벽 열저항 엄밀해 검산")
    print(f"    단일 노드 h1 = {h1:.1f} W/(m^2 K), h2 = {h2:.1f} W/(m^2 K)")
    print(f"    G_1w = {G1w/1e3:.2f} kW/K, G_w2 = {Gw2/1e3:.2f} kW/K, 원통 저항 오차 = {rel_err_R:.2e}")

    # 등온 영점 테스트 (Null test)
    res_iso = solve_ihx_steady(ihx, T1_in_K=650.0, mdot1_kg_s=450.0, T2_in_K=650.0, mdot2_kg_s=310.0)
    max_iso_dev = max(
        max(abs(t - 650.0) for t in res_iso["T1_K"]),
        max(abs(t - 650.0) for t in res_iso["Tw_K"]),
        max(abs(t - 650.0) for t in res_iso["T2_K"]),
    )
    assert max_iso_dev < 1.0e-10 and abs(res_iso["Q_W"]) < 1.0e-5
    print(f"\n[4] 등온 영점 테스트(650 K) 통과: 최대 온도 편차 = {max_iso_dev:.2e} K, Q = {res_iso['Q_W']:.2e} W")

    # 임의 과도 스텝 이산 에너지 보존 검산
    solve_ihx_steady(ihx, T1_in_K=720.0, mdot1_kg_s=450.0, T2_in_K=575.0, mdot2_kg_s=310.0)
    res_step = step_ihx(ihx, dt_s=0.1, T1_in_K=710.0, mdot1_kg_s=300.0, T2_in_K=575.0, mdot2_kg_s=220.0)
    assert abs(res_step["energy_residual_W"]) < 1.0e-4
    print(f"[5] 과도 1스텝(dt=0.1 s) 이산 에너지 보존 잔차 = {res_step['energy_residual_W']:.3e} W (기계 정밀도 수준)")

    # epsilon-NTU 해석해 대조 (촘촘한 격자에서 수치해가 해석해로 수렴하는지)
    fine = create_ihx(n_nodes=400)
    r_f = solve_ihx_steady(fine, T1_in_K=723.5, mdot1_kg_s=452.0, T2_in_K=574.2, mdot2_kg_s=310.6)
    chk = effectiveness_ntu_check(fine, 723.5, 452.0, 574.2, 310.6)
    rel_entu = abs(r_f["Q_W"] - chk["Q_eNTU_W"]) / chk["Q_eNTU_W"]
    assert rel_entu < 2.0e-3, f"epsilon-NTU 해석해와 불일치: {rel_entu:.3e}"
    print(f"[6] epsilon-NTU 해석해 대조 (N=400): Q = {r_f['Q_W']/1e6:.3f} MW, Q_eNTU = {chk['Q_eNTU_W']/1e6:.3f} MW, "
          f"상대차 {rel_entu:.2e} (NTU = {chk['NTU']:.2f})")

    # 경계조건 CSV 보간기 검산
    bc17 = load_secondary_bc("SHRT17")
    T2_0, m2_0 = bc17(0.0)
    T2_900, m2_900 = bc17(900.0)
    print(f"[7] SHRT-17 경계조건 CSV 보간 검산:")
    print(f"    t =   0 s -> T2_in = {T2_0:.3f} K ({T2_0-273.15:.2f} C), mdot2 = {m2_0:.3f} kg/s")
    print(f"    t = 900 s -> T2_in = {T2_900:.3f} K ({T2_900-273.15:.2f} C), mdot2 = {m2_900:.3f} kg/s")
    print("=" * 68)
    print("모든 검산 항목 정상 통과!")


if __name__ == "__main__":
    _run_self_checks()
