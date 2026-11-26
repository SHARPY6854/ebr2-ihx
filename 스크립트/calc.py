"""EBR-II 노심 + 펌프/루프 입력값 SI 변환 및 검산.
출처: ANL-ARC-226 Rev.1 (Sumner & Wei, 2012). 물성: SAS4A/SASSYS-1 Eq. 12.13-5.
"""
import csv, json, math, os

OUT = os.path.dirname(os.path.abspath(__file__))
GPM = 6.30901964e-5  # m3/s per US gpm

def F2K(f): return (f - 32.0) * 5.0 / 9.0 + 273.15
def rho_na(T):  # SAS4A/SASSYS-1 Eq. 12.13-5, T in K
    return 1.00423e3 - 0.21390 * T - 1.1046e-5 * T**2

R = {}
T800, T582 = F2K(800), F2K(582)
R["rho_800F"] = rho_na(T800); R["rho_582F"] = rho_na(T582)

# --- 펌프 정격 검산 (표 3.1.4.1.2, 식 3.7)
QR, WR, HR_m, HR_Pa = 0.2946, 250.2, 37.7, 358530.0
R["rho_implied_by_wr_over_QR"] = WR / QR
R["HR_Pa_from_37.7m_x_rho800"] = HR_m * 9.80665 * R["rho_800F"]
R["rho_implied_by_HRPa_over_g_HRm"] = HR_Pa / (9.80665 * HR_m)

# --- 초기 유량 (표 3.1.1.1, 3.2.1.1, 5.3.x)
def gpm800(q): return q * GPM * R["rho_800F"]
def gpm582(q): return q * GPM * R["rho_582F"]
init = {
 "SHRT-17": dict(P_MW=57.3, core_gpm=8500.0, pumps_gpm=8699.1, p1=4352.6, p2=4346.5,
                 Tin_F=664.8, ihx2_gpm=5614.6, ihx2_T_F=573.8, pdis1_psi=42.9, pdis2_psi=41.9, pUP_psi=6.36),
 "SHRT-45R": dict(P_MW=60.0, core_gpm=8766.6, pumps_gpm=8971.9, p1=4489.1, p2=4482.8,
                  Tin_F=650.8, ihx2_gpm=5474.5, ihx2_T_F=552.1, pdis1_psi=41.84, pdis2_psi=41.84, pUP_psi=6.51),
}
PSI = 6894.757
for k, d in init.items():
    d["core_kgs"] = gpm800(d["core_gpm"]); d["pumps_kgs"] = gpm800(d["pumps_gpm"])
    d["p1_kgs"] = gpm800(d["p1"]); d["p2_kgs"] = gpm800(d["p2"])
    d["Tin_C"] = F2K(d["Tin_F"]) - 273.15
    d["ihx2_kgs"] = gpm582(d["ihx2_gpm"]); d["ihx2_T_C"] = F2K(d["ihx2_T_F"]) - 273.15
    d["pdis1_kPa"] = d["pdis1_psi"] * PSI / 1e3; d["pdis2_kPa"] = d["pdis2_psi"] * PSI / 1e3
    d["pUP_kPa"] = d["pUP_psi"] * PSI / 1e3
    # 노심 평균 온도상승 (cp ~ 1270 J/kg/K 근방, 검산용 근사)
    d["dT_core_est_K"] = d["P_MW"] * 1e6 / (d["core_kgs"] * 1270.0)
R["init"] = init

# --- 노심 형상 (표 4.1.1, 4.1.2, MARK-II AI)
N, Dc, tc, Df, Dw, Hw, Lf = 91, 4.41e-3, 0.315e-3, 3.36e-3, 1.24e-3, 0.1524, 0.3429
Fin = 56.1e-3
Dci = Dc - 2 * tc
P = Dc + Dw
A_hex = math.sqrt(3) / 2 * Fin**2
A_pin = N * math.pi / 4 * Dc**2
A_wire = N * math.pi / 4 * Dw**2  # 와이어 경사 무시
A_flow = A_hex - A_pin - A_wire
Pw = N * math.pi * Dc + N * math.pi * Dw + 6 * Fin / math.sqrt(3)
geo = dict(clad_ID_mm=Dci * 1e3, bond_gap_mm=(Dci - Df) / 2 * 1e3, pitch_mm=P * 1e3,
           P_over_D=P / Dc, H_over_D=Hw / Dc, A_hex_mm2=A_hex * 1e6, A_flow_mm2=A_flow * 1e6,
           wetted_perim_mm=Pw * 1e3, Dh_mm=4 * A_flow / Pw * 1e3,
           bundle_fit_flat_mm=(math.sqrt(3) * 5 * P + Dc + 2 * Dw) * 1e3)
R["geo"] = geo

# 적용범위 점검
R["range_UCTS_ok"] = (19 <= N <= 217) and (1.025 <= geo["P_over_D"] <= 1.42) and (8 <= geo["H_over_D"] <= 50)
R["range_Mikityuk_PD_ok"] = 1.1 <= geo["P_over_D"] <= 1.95

# --- 펌프 수두 모델 (식 3.1~3.7)
b1, b2, b3m, b3t, b3l, b4t, b4l, wt, em = 1.1740, 0.0818, -0.2558, -0.5923, -0.0471, 1.9, 1.0, 0.06, 0.55
def head(s, w):
    """s, w: 무차원(정격 대비). 반환: 무차원 수두 H/Hr."""
    if s > em * w:
        b3, b4 = b3m, b4t
    elif w >= wt:
        b3, b4 = b3t, b4t
    else:
        b3, b4 = b3l, b4l
    sgn = 0.0 if w == 0 else (1.0 if w > 0 else -1.0)
    return b1 * s**2 + b2 * s * w + b3 * sgn * abs(w)**b4
R["head_rated"] = head(1.0, 1.0)
R["head_shutoff"] = head(1.0, 0.0)
# 초기점: SHRT-45R 펌프1
s0 = 792.05 / 870.0; w0 = init["SHRT-45R"]["p1_kgs"] / WR
R["head_45R_p1_init_ratio"] = head(s0, w0); R["head_45R_p1_init_kPa"] = head(s0, w0) * HR_Pa / 1e3
s0 = 799.03 / 870.0; w0 = init["SHRT-17"]["p1_kgs"] / WR
R["head_17_p1_init_kPa"] = head(s0, w0) * HR_Pa / 1e3

# --- 경계조건 CSV
pump17 = [(0,799.03,764.92),(2,727.77,695.42),(4,566.79,532.61),(6.5,407.49,379.64),(9,323.18,301.59),
 (12.5,244.06,230.18),(17.5,172.50,167.39),(24.5,109.25,113.98),(33,56.46,74.62),(42,0.69,46.54),
 (51,0.76,2.29),(90,0.76,2.44),(120,0.61,2.44),(180,0.76,2.44),(240,0.84,2.59),(300,0.61,2.44),
 (360,0.61,2.44),(420,0.61,2.44),(480,0.61,2.44),(540,0.76,2.44),(600,0.76,2.75),(660,0.69,2.44),
 (720,0.69,2.44),(780,0.76,2.59),(840,0.76,2.44),(900,0.92,2.44)]
pump45 = [(0,792.05,762.27),(1,769.10,739.07),(1.5,731.19,701.21),(3,606.46,575.65),(4.5,511.19,482.21),
 (7,397.55,373.96),(11,283.90,266.69),(17,213.66,198.62),(25.5,154.38,144.98),(36.5,105.14,100.40),
 (48.5,72.67,68.91),(64.5,46.45,42.42),(82,25.45,19.69),(90,17.13,1.38),(95,10.39,-0.67),(98.5,-0.67,-0.67),
 (120,-1.06,-0.99),(180,-0.72,-0.67),(240,-0.67,-0.53),(300,-1.26,-1.23),(360,-1.53,-1.38),(420,-0.87,-0.77),
 (480,-0.70,-0.50),(540,-1.60,-1.50),(600,-1.55,-1.57),(660,-1.55,-1.50),(720,-1.33,-1.31),(780,-0.96,-0.87),
 (840,-1.72,-1.62),(900,-1.28,-1.14)]
pow17 = [(0,53.93,3.36,57.29),(1,8.24,3.12,11.36),(2,6.86,2.95,9.81),(3.5,5.61,2.77,8.38),(5,4.76,2.63,7.39),
 (7.5,3.82,2.46,6.29),(10,3.23,2.33,5.56),(15,2.47,2.15,4.62),(22.5,1.85,1.95,3.80),(30,1.46,1.82,3.28),
 (40,1.10,1.70,2.79),(50,0.84,1.60,2.44),(60,0.66,1.52,2.19),(90,0.38,1.35,1.73),(120,0.25,1.24,1.49),
 (180,0.13,1.09,1.22),(240,0.10,1.00,1.10),(300,0.10,0.93,1.02),(360,0.08,0.88,0.96),(420,0.08,0.83,0.92),
 (480,0.08,0.80,0.88),(540,0.08,0.77,0.84),(600,0.08,0.74,0.82),(660,0.08,0.71,0.79),(720,0.07,0.69,0.76),
 (780,0.08,0.67,0.75),(840,0.07,0.65,0.72),(900,0.07,0.63,0.70)]
em45 = [(0.0,5797.14),(13.5,6121.30),(47.0,6200.54),(79.5,6186.53),(121.5,6167.36),(123.2,6172.51),
 (147.0,6130.74),(254.5,6075.52),(372.5,6026.88),(532.0,5999.42),(597.5,5989.97),(616.0,6002.85),
 (623.0,6064.94),(628.0,6579.34),(633.0,6597.94),(640.5,7322.64),(644.5,7445.66),(651.0,8936.83),
 (659.5,8945.13),(667.0,9119.65),(671.5,9438.08),(680.0,9442.66),(681.5,9335.09),(685.0,9359.97),
 (686.5,9203.76),(701.5,9210.34),(702.5,9111.64),(741.0,9118.79),(743.0,9203.19),(792.0,9244.39),
 (842.0,9261.55),(899.5,9238.09)]

# 표 자체 검산: 출력표 합계
R["pow17_sum_maxerr_MW"] = max(abs(f + d - t) for _, f, d, t in pow17)
R["n_rows"] = dict(pump17=len(pump17), pump45=len(pump45), pow17=len(pow17), em45=len(em45))

def wcsv(name, header, rows):
    with open(os.path.join(OUT, name), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f); w.writerow(header); w.writerows(rows)
wcsv("bc_pump_speed_SHRT17.csv", ["time_s", "pump1_rpm", "pump2_rpm", "pump1_rad_s", "pump2_rad_s", "pump1_s_norm", "pump2_s_norm"],
     [(t, a, b, round(a*2*math.pi/60, 4), round(b*2*math.pi/60, 4), round(a/870, 5), round(b/870, 5)) for t, a, b in pump17])
wcsv("bc_pump_speed_SHRT45R.csv", ["time_s", "pump1_rpm", "pump2_rpm", "pump1_rad_s", "pump2_rad_s", "pump1_s_norm", "pump2_s_norm"],
     [(t, a, b, round(a*2*math.pi/60, 4), round(b*2*math.pi/60, 4), round(a/870, 5), round(b/870, 5)) for t, a, b in pump45])
wcsv("bc_power_SHRT17.csv", ["time_s", "fission_W", "decay_W", "total_W"],
     [(t, f*1e6, d*1e6, tt*1e6) for t, f, d, tt in pow17])
wcsv("bc_em_pump_current_SHRT45R.csv", ["time_s", "current_A"], em45)

json.dump(R, open(os.path.join(OUT, "calc_results.json"), "w"), indent=1, ensure_ascii=False)
for k, v in R.items():
    if k != "init": print(k, v)
for k, d in init.items():
    print(k, {kk: (round(vv, 2) if isinstance(vv, float) else vv) for kk, vv in d.items()})
