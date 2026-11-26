"""4단계: SHRT-17 IHX 모듈 과도 동작 시험 스크립트 (0 ~ 900 s).

주의: 기본 실행의 두 케이스는 1차측 입구 온도를 초기값으로 고정한 "모듈 동작 시험"이다.
SHRT-17은 스크램 시험이라 실제 IHX 1차 입구 온도는 시간에 따라 떨어지므로, 기본 결과를
SHRT-17 예측값으로 쓰면 안 된다. 실측 검증은 --primary-csv 로 실측 IHX 1차 입구 온도와 유량을 넣어 수행한다.

기능:
- 2차측 입구 유량과 온도는 사양서 변환 CSV (표 3.1.4.2.1~2, p.21)를 선형 보간해 입력
- 1차측 입력 모드:
  1) 기본 실행 시 두 가지 기준 케이스를 0~900 s 계산하여 비교:
     - 케이스 A (1차 유량 고정): T1_in = 450.39 C 고정, mdot1 = 452.07 kg/s 고정 (2차측 CSV 단독 영향 확인)
     - 케이스 B (1차 펌프 감속 + 자연순환 가정): T1_in = 450.39 C 고정, mdot1은 사양서 펌프 속도 CSV
       (표 3.1.4.1.1, p.14) 비례 감속 후 정격의 4.5% (20.34 kg/s) 자연순환 하한 적용
  2) --primary-csv <파일경로> 지정 시: 디지타이징 실측 또는 타 파트 계산 CSV(time_s, T1_in_C/T1_in_K, mdot1_kg_s) 사용
- 결과물:
  - 스크립트/transient_shrt17_results.csv (시간별 계산 결과)
  - 도면/shrt17_transient_0_900s.svg (유량, 열교환량, 출구 온도 3단 그래프)
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from typing import Callable, Dict, List, Tuple

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import ihx

GPM_TO_M3S = 6.30901964e-5


def f_to_k(temp_f: float) -> float:
    return (temp_f - 32.0) * 5.0 / 9.0 + 273.15


def load_pump_coastdown_flow(mdot1_init: float, nc_floor_frac: float = 0.045) -> Callable[[float], float]:
    """사양서 표 3.1.4.1.1 (p.14) 펌프 회전수 CSV 기반 1차 유량 감속 + 자연순환 하한(가정) 함수."""
    pump_csv = os.path.join(ROOT, "경계조건_CSV", "bc_pump_speed_SHRT17.csv")
    times, avg_rpm = [], []
    with open(pump_csv, encoding="utf-8") as f:
        for row in csv.DictReader(f):
            times.append(float(row["time_s"]))
            avg_rpm.append(0.5 * (float(row["pump1_rpm"]) + float(row["pump2_rpm"])))
    rpm0 = avg_rpm[0]
    rel_speed = [r / rpm0 for r in avg_rpm]

    def mdot1_func(t_s: float) -> float:
        sr = ihx.interp_linear(t_s, times, rel_speed)
        return mdot1_init * max(sr, nc_floor_frac)

    return mdot1_func


def load_custom_primary_csv(
    csv_path: str, default_T1_K: float, default_mdot1: float
) -> Callable[[float], Tuple[float, float]]:
    """외부 디지타이징 또는 타 파트 연동용 1차측 CSV 로더."""
    times, t1_k_list, m1_list = [], [], []
    with open(csv_path, encoding="utf-8") as f:
        reader = csv.DictReader(f)
        fields = reader.fieldnames or []
        for row in reader:
            times.append(float(row["time_s"]))
            if "T1_in_K" in fields:
                t1_k_list.append(float(row["T1_in_K"]))
            elif "T1_in_C" in fields:
                t1_k_list.append(float(row["T1_in_C"]) + 273.15)
            else:
                t1_k_list.append(default_T1_K)

            if "mdot1_kg_s" in fields:
                m1_list.append(float(row["mdot1_kg_s"]))
            else:
                m1_list.append(default_mdot1)

    def prim_func(t_s: float) -> Tuple[float, float]:
        return ihx.interp_linear(t_s, times, t1_k_list), ihx.interp_linear(t_s, times, m1_list)

    return prim_func


def run_single_case(
    prim_bc: Callable[[float], Tuple[float, float]],
    sec_bc: Callable[[float], Tuple[float, float]],
    t_end: float = 900.0,
) -> Tuple[List[Dict[str, float]], float]:
    """0 ~ t_end 초까지 IHX 과도 적분 실행 후 1초 간격 기록 반환."""
    model = ihx.create_ihx(shell_corr="mikityuk", tube_corr="seban")
    T1_0, m1_0 = prim_bc(0.0)
    T2_0, m2_0 = sec_bc(0.0)
    r0 = ihx.solve_ihx_steady(model, T1_0, m1_0, T2_0, m2_0)

    history: List[Dict[str, float]] = [
        {
            "time_s": 0.0,
            "mdot1_kg_s": m1_0,
            "mdot2_kg_s": m2_0,
            "T1_in_C": T1_0 - 273.15,
            "T1_out_C": float(r0["T1_out_K"]) - 273.15,
            "T2_in_C": T2_0 - 273.15,
            "T2_out_C": float(r0["T2_out_K"]) - 273.15,
            "Q_MW": float(r0["Q_W"]) / 1.0e6,
            "energy_residual_W": abs(float(r0["energy_residual_W"])),
        }
    ]

    t = 0.0
    max_res_W = abs(float(r0["energy_residual_W"]))
    next_save_t = 1.0

    while t < t_end - 1.0e-9:
        dt = 0.1 if t < 100.0 else 0.5
        t = round(t + dt, 6)
        T1_in, m1 = prim_bc(t)
        T2_in, m2 = sec_bc(t)
        res = ihx.step_ihx(model, dt, T1_in, m1, T2_in, m2)
        err_w = abs(float(res["energy_residual_W"]))
        if err_w > max_res_W:
            max_res_W = err_w

        if t >= next_save_t - 1.0e-6:
            history.append(
                {
                    "time_s": round(t, 2),
                    "mdot1_kg_s": m1,
                    "mdot2_kg_s": m2,
                    "T1_in_C": T1_in - 273.15,
                    "T1_out_C": float(res["T1_out_K"]) - 273.15,
                    "T2_in_C": T2_in - 273.15,
                    "T2_out_C": float(res["T2_out_K"]) - 273.15,
                    "Q_MW": float(res["Q_W"]) / 1.0e6,
                    "energy_residual_W": err_w,
                }
            )
            next_save_t += 1.0

    return history, max_res_W


def write_svg_plots(hist_const: List[Dict[str, float]], hist_coast: List[Dict[str, float]], svg_path: str) -> None:
    """시간(0~900 s)에 따른 유량, 열교환량, 출구 온도 3단 비교 SVG 그래프 생성."""
    W, H = 820, 740
    x0, x1 = 88, 785
    t_max = 900.0

    def sx(t: float) -> float:
        return x0 + (t / t_max) * (x1 - x0)

    def sy(val: float, vmin: float, vmax: float, py0: float, py1: float) -> float:
        return py1 - (val - vmin) / (vmax - vmin) * (py1 - py0)

    def poly(hist: List[Dict[str, float]], key: str, vmin: float, vmax: float, py0: float, py1: float) -> str:
        step = 2
        pts = hist[::step]
        if pts[-1]["time_s"] != hist[-1]["time_s"]:
            pts = pts + [hist[-1]]
        return " ".join(f"{sx(r['time_s']):.1f},{sy(r[key], vmin, vmax, py0, py1):.1f}" for r in pts)

    out = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" font-family="sans-serif">',
        f'<rect width="{W}" height="{H}" fill="#ffffff"/>',
        '<text x="435" y="28" text-anchor="middle" font-size="17" font-weight="bold" fill="#111827">IHX 모듈 과도 동작 시험 (SHRT-17 2차측 경계조건, 1차 입구 온도 고정, 0 ~ 900 s)</text>',
        '<text x="435" y="47" text-anchor="middle" font-size="11" fill="#4b5563">실선: 케이스 B (1차 유량 펌프 속도 비례 감소, 하한 4.5% 가정) | 점선: 케이스 A (1차 유량 고정) | 1차 입구 온도는 두 케이스 모두 고정</text>',
    ]

    panels = [
        ("1) 질량유량 변화 [kg/s]", 72.0, 242.0, 0.0, 500.0, [0, 100, 200, 300, 400, 500]),
        ("2) IHX 열교환량 Q 변화 [MW]", 292.0, 462.0, 0.0, 60.0, [0, 15, 30, 45, 60]),
        ("3) 1차 및 2차 출구 온도 변화 [°C]", 512.0, 692.0, 290.0, 460.0, [300, 340, 380, 420, 460]),
    ]

    for title, py0, py1, vmin, vmax, yticks in panels:
        out.append(f'<text x="{x0}" y="{py0-8:.1f}" font-size="12" font-weight="bold" fill="#111827">{title}</text>')
        out.append(f'<rect x="{x0}" y="{py0}" width="{x1-x0}" height="{py1-py0}" fill="#f9fafb" stroke="#9ca3af" stroke-width="1"/>')
        for yt in yticks:
            yy = sy(float(yt), vmin, vmax, py0, py1)
            out.append(f'<line x1="{x0}" y1="{yy:.1f}" x2="{x1}" y2="{yy:.1f}" stroke="#e5e7eb" stroke-width="1"/>')
            out.append(f'<text x="{x0-8}" y="{yy+4:.1f}" text-anchor="end" font-size="10" fill="#374151">{yt}</text>')
        for xt in range(0, 901, 100):
            xx = sx(float(xt))
            out.append(f'<line x1="{xx:.1f}" y1="{py0}" x2="{xx:.1f}" y2="{py1}" stroke="#e5e7eb" stroke-width="1"/>')
            out.append(f'<text x="{xx:.1f}" y="{py1+15:.1f}" text-anchor="middle" font-size="10" fill="#374151">{xt}</text>')

    # Panel 1 curves: mdot1 (const dashed, coast solid), mdot2 (blue solid)
    _, py0, py1, vmin, vmax, _ = panels[0]
    out.append(f'<polyline fill="none" stroke="#dc2626" stroke-width="1.8" stroke-dasharray="6,4" points="{poly(hist_const, "mdot1_kg_s", vmin, vmax, py0, py1)}"/>')
    out.append(f'<polyline fill="none" stroke="#dc2626" stroke-width="2.4" points="{poly(hist_coast, "mdot1_kg_s", vmin, vmax, py0, py1)}"/>')
    out.append(f'<polyline fill="none" stroke="#2563eb" stroke-width="2.4" points="{poly(hist_coast, "mdot2_kg_s", vmin, vmax, py0, py1)}"/>')
    out.append(f'<text x="535" y="{py0+40:.1f}" font-size="10" font-weight="bold" fill="#dc2626">빨간 실선: 1차 유량(감속) / 점선: 1차 유량(고정 {hist_const[0]["mdot1_kg_s"]:.0f} kg/s)</text>')
    out.append(f'<text x="535" y="{py0+56:.1f}" font-size="10" font-weight="bold" fill="#2563eb">파란 실선: 2차 유량 (사양서 CSV, 310.8 → 17.2 kg/s)</text>')

    # Panel 2 curves: Q_MW
    _, py0, py1, vmin, vmax, _ = panels[1]
    out.append(f'<polyline fill="none" stroke="#059669" stroke-width="1.8" stroke-dasharray="6,4" points="{poly(hist_const, "Q_MW", vmin, vmax, py0, py1)}"/>')
    out.append(f'<polyline fill="none" stroke="#059669" stroke-width="2.5" points="{poly(hist_coast, "Q_MW", vmin, vmax, py0, py1)}"/>')
    out.append(f'<text x="500" y="{py0+22:.1f}" font-size="10" font-weight="bold" fill="#047857">초록 실선: 케이스 B 열교환량 ({hist_coast[0]["Q_MW"]:.2f} → {hist_coast[-1]["Q_MW"]:.2f} MW)</text>')
    out.append(f'<text x="500" y="{py0+38:.1f}" font-size="10" fill="#065f46">초록 점선: 케이스 A 열교환량 ({hist_const[0]["Q_MW"]:.2f} → {hist_const[-1]["Q_MW"]:.2f} MW)</text>')

    # Panel 3 curves: T1_out, T2_out, T2_in
    _, py0, py1, vmin, vmax, _ = panels[2]
    out.append(f'<polyline fill="none" stroke="#9ca3af" stroke-width="1.5" stroke-dasharray="3,3" points="{poly(hist_coast, "T2_in_C", vmin, vmax, py0, py1)}"/>')
    out.append(f'<polyline fill="none" stroke="#dc2626" stroke-width="1.8" stroke-dasharray="6,4" points="{poly(hist_const, "T1_out_C", vmin, vmax, py0, py1)}"/>')
    out.append(f'<polyline fill="none" stroke="#dc2626" stroke-width="2.4" points="{poly(hist_coast, "T1_out_C", vmin, vmax, py0, py1)}"/>')
    out.append(f'<polyline fill="none" stroke="#2563eb" stroke-width="1.8" stroke-dasharray="6,4" points="{poly(hist_const, "T2_out_C", vmin, vmax, py0, py1)}"/>')
    out.append(f'<polyline fill="none" stroke="#2563eb" stroke-width="2.4" points="{poly(hist_coast, "T2_out_C", vmin, vmax, py0, py1)}"/>')
    out.append(f'<text x="445" y="{py0+42:.1f}" font-size="10" font-weight="bold" fill="#2563eb">파란선: 2차 출구 T2_out (실선 케이스 B, 점선 케이스 A)</text>')
    out.append(f'<text x="445" y="{py0+58:.1f}" font-size="10" font-weight="bold" fill="#dc2626">빨간선: 1차 출구 T1_out (실선 케이스 B, 점선 케이스 A)</text>')
    out.append(f'<text x="445" y="{py0+74:.1f}" font-size="10" fill="#4b5563">회색 점선: 2차 입구 T2_in (CSV 경계조건, 301.1 → 292.1 °C)</text>')

    out.append(f'<text x="{0.5*(x0+x1):.1f}" y="{H-12}" text-anchor="middle" font-size="12" font-weight="bold" fill="#111827">시간 [s]</text>')
    out.append("</svg>")

    with open(svg_path, "w", encoding="utf-8") as f:
        f.write("\n".join(out))


def main() -> None:
    parser = argparse.ArgumentParser(description="SHRT-17 IHX 단독 과도 계산 (0~900 s)")
    parser.add_argument("--primary-csv", type=str, default="", help="1차측 실측/연동 CSV 파일 경로 (선택)")
    args = parser.parse_args()

    rho_800f = ihx.rho_na(f_to_k(800.0))
    mdot_core0 = 8500.0 * GPM_TO_M3S * rho_800f
    # 표 5.3.1.1 (p.129): IHX 번들 유량 = IHX 출구 8421.0 gpm (입구 누설 17.5 gpm 제외, 번들 우회로 가정)
    mdot1_init = 8421.0 * GPM_TO_M3S * rho_800f  # 451.13 kg/s
    T_core_in_K = f_to_k(664.8)

    T1_init_K = T_core_in_K + 57.29e6 / (mdot_core0 * ihx.cp_na(T_core_in_K + 50.0))
    for _ in range(15):
        T1_init_K = T_core_in_K + 57.29e6 / (mdot_core0 * ihx.cp_na(0.5 * (T_core_in_K + T1_init_K)))

    sec_bc = ihx.load_secondary_bc("SHRT17")

    if args.primary_csv:
        prim_custom = load_custom_primary_csv(args.primary_csv, T1_init_K, mdot1_init)
        hist_custom, max_err = run_single_case(prim_custom, sec_bc, 900.0)
        print(f"사용자 지정 1차측 CSV 실행 완료 (최대 에너지 오차 = {max_err:.3e} W)")
        print("시각[s] | Q[MW] | T1_out[C] | T2_out[C]")
        for r in hist_custom:
            if int(r["time_s"]) in (0, 10, 30, 60, 120, 300, 600, 900) and r["time_s"] == int(r["time_s"]):
                print(f"{int(r['time_s']):6d} | {r['Q_MW']:5.2f} | {r['T1_out_C']:9.2f} | {r['T2_out_C']:9.2f}")
        return

    # 케이스 A: 1차측 입구 온도, 유량 모두 초기값 고정
    hist_const, err_const = run_single_case(lambda _t: (T1_init_K, mdot1_init), sec_bc, 900.0)

    # 케이스 B: 1차측 입구 온도 초기값 고정 + 펌프 회전수 CSV 감속 및 4.5% 자연순환 하한(가정)
    coast_flow = load_pump_coastdown_flow(mdot1_init, nc_floor_frac=0.045)
    hist_coast, err_coast = run_single_case(lambda t: (T1_init_K, coast_flow(t)), sec_bc, 900.0)

    # 결과 CSV 저장
    csv_out = os.path.join(ROOT, "스크립트", "transient_shrt17_results.csv")
    with open(csv_out, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow([
            "time_s",
            "mdot2_kg_s",
            "T2_in_C",
            "caseA_mdot1_kg_s",
            "caseA_Q_MW",
            "caseA_T1_out_C",
            "caseA_T2_out_C",
            "caseB_mdot1_kg_s",
            "caseB_Q_MW",
            "caseB_T1_out_C",
            "caseB_T2_out_C",
        ])
        for ra, rb in zip(hist_const, hist_coast):
            w.writerow([
                ra["time_s"],
                round(ra["mdot2_kg_s"], 4),
                round(ra["T2_in_C"], 3),
                round(ra["mdot1_kg_s"], 4),
                round(ra["Q_MW"], 4),
                round(ra["T1_out_C"], 3),
                round(ra["T2_out_C"], 3),
                round(rb["mdot1_kg_s"], 4),
                round(rb["Q_MW"], 4),
                round(rb["T1_out_C"], 3),
                round(rb["T2_out_C"], 3),
            ])

    svg_out = os.path.join(ROOT, "도면", "shrt17_transient_0_900s.svg")
    write_svg_plots(hist_const, hist_coast, svg_out)

    print("=" * 78)
    print("4단계: SHRT-17 IHX 단독 과도 계산 완료 (0 ~ 900 s)")
    print("=" * 78)
    print(f"최대 이산 에너지 보존 잔차: 케이스 A = {err_const:.3e} W, 케이스 B = {err_coast:.3e} W")
    print("\n[주요 시각별 계산 결과 비교]")
    print("시각[s] | 2차유량[kg/s] | 케이스A Q[MW] T1_out[C] T2_out[C] | 케이스B 1차유량 Q[MW] T1_out[C] T2_out[C]")
    print("-" * 78)
    idx_map = {int(r["time_s"]): i for i, r in enumerate(hist_const) if r["time_s"] == int(r["time_s"])}
    for ts in [0, 5, 15, 30, 50, 100, 300, 600, 900]:
        ra = hist_const[idx_map[ts]]
        rb = hist_coast[idx_map[ts]]
        print(
            f"{ts:6d}  |    {ra['mdot2_kg_s']:7.2f}    |"
            f"   {ra['Q_MW']:6.2f}     {ra['T1_out_C']:7.2f}   {ra['T2_out_C']:7.2f}  |"
            f"  {rb['mdot1_kg_s']:6.2f}  {rb['Q_MW']:5.2f}   {rb['T1_out_C']:7.2f}   {rb['T2_out_C']:7.2f}"
        )
    print("=" * 78)
    print(f"CSV 저장 완료 : {csv_out}")
    print(f"그래프 저장 완료: {svg_out}")


if __name__ == "__main__":
    main()
