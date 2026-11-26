"""5단계: SHRT-45R 2차측 초기 유량 불일치 판정 스크립트.

사양서(ANL-ARC-226 Rev.1) 내 불일치:
- 초기조건 표 3.2.1.1 (p.24) : 5474.5 gpm (582 F 기준 -> 302.82 kg/s)
- 과도 유량표 3.2.4.2.1 (p.36): 0.0 s 값 5079.5 gpm ("gpm at XX F" -> 582 F 가정 시 280.97 kg/s)

정상상태 열수지(60.0 MW)와 ihx.py 대향류 모델(기본 N=80, active 길이)로 어느 값이 열수지와 맞는지 판정한다.
전제: (1) IHX 열량 = 노심 출력 60.0 MW (탱크 열손실, 펌프 발열 무시), (2) IHX 1차 입구 온도 = 노심 혼합 출구 온도 (Z-pipe 단열).
"""

from __future__ import annotations

import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import ihx

GPM_TO_M3S = 6.30901964e-5


def f_to_k(temp_f: float) -> float:
    return (temp_f - 32.0) * 5.0 / 9.0 + 273.15


def required_t2_out_for_power(T2_in_K: float, mdot2: float, Q_target_W: float) -> float:
    """주어진 유량 mdot2로 Q_target_W를 모두 흡수할 때 필요한 2차측 출구 온도 [K]."""
    T_out = T2_in_K + Q_target_W / (mdot2 * ihx.cp_na(T2_in_K + 75.0))
    for _ in range(15):
        T_out = T2_in_K + Q_target_W / (mdot2 * ihx.cp_na(0.5 * (T2_in_K + T_out)))
    return T_out


def main() -> None:
    rho_800f = ihx.rho_na(f_to_k(800.0))
    rho_582f = ihx.rho_na(f_to_k(582.0))

    # SHRT-45R 사양서 기준값 (표 3.2.1.1 p.24, p.26 및 표 5.3.2.1 p.130)
    Q_target_W = 60.0e6
    T_core_in_K = f_to_k(650.8)                  # 616.93 K (343.78 C)
    mdot_core = 8766.6 * GPM_TO_M3S * rho_800f   # 469.64 kg/s
    # 표 5.3.2.1 (p.130): IHX 입구 8703.2 gpm 중 18.1 gpm은 입구 누설(번들 우회로 가정) -> 번들 유량은 IHX 출구 8685.1 gpm
    mdot1_ihx = 8685.1 * GPM_TO_M3S * rho_800f   # 465.28 kg/s
    T2_in_K = f_to_k(552.1)                      # 562.09 K (288.94 C)

    # 노심 60.0 MW 열수지 기준 혼합 출구 온도 (= 단열 Z-pipe 가정 시 IHX 1차측 입구 온도)
    T1_in_K = T_core_in_K + Q_target_W / (mdot_core * ihx.cp_na(T_core_in_K + 50.0))
    for _ in range(15):
        T1_in_K = T_core_in_K + Q_target_W / (mdot_core * ihx.cp_na(0.5 * (T_core_in_K + T1_in_K)))

    candidates = [
        ("표 3.2.4.2.1 0초 값 (5079.5 gpm)", 5079.5),
        ("초기조건 표 3.2.1.1 값 (5474.5 gpm)", 5474.5),
    ]

    print("=" * 76)
    print("5단계: SHRT-45R 2차측 초기 유량 불일치 판정 (60.0 MW 정상상태 열수지)")
    print("=" * 76)
    print(f"[기준 조건]")
    print(f"  초기 열출력              : {Q_target_W/1e6:.2f} MW (표 3.2.1.1, p.24)")
    print(f"  노심 입구(=저온풀) 온도  : {T_core_in_K-273.15:.2f} C (650.8 F)")
    print(f"  노심 혼합 출구 온도      : {T1_in_K-273.15:.2f} C (IHX 1차 입구 기준 온도)")
    print(f"  IHX 1차측 유량           : {mdot1_ihx:.2f} kg/s (번들 유량 8685.1 gpm, 표 5.3.2.1 p.130)")
    print(f"  IHX 2차측 입구 온도      : {T2_in_K-273.15:.2f} C (552.1 F, 표 3.2.1.1 p.26)")

    results_summary = {}

    for label, gpm in candidates:
        mdot2 = gpm * GPM_TO_M3S * rho_582f
        # 1) 단순 열수지(60 MW 전량 흡수 가정) 시 요구되는 2차측 출구 온도
        T2_req_K = required_t2_out_for_power(T2_in_K, mdot2, Q_target_W)
        pinch_K = T1_in_K - T2_req_K  # 음수면 전제 하에서 T2_out > T1_in 이 필요 (열수지 성립 불가)

        # 2) 이론적 최대 가능 열전달량 Q_max = mdot2 * cp2 * (T1_in - T2_in)
        cp2_mid = ihx.cp_na(0.5 * (T2_in_K + T1_in_K))
        Q_max_W = mdot2 * cp2_mid * (T1_in_K - T2_in_K)

        # 3) ihx.py 모델 실제 정상상태 계산 (기본 격자)
        model = ihx.create_ihx()
        res = ihx.solve_ihx_steady(model, T1_in_K, mdot1_ihx, T2_in_K, mdot2)

        # 4) 정확히 60.0 MW를 전달하려면 1차측 입구 온도가 몇 도여야 하는지 역산
        T1_needed_K = T1_in_K
        m_inv = ihx.create_ihx()
        for _ in range(25):
            r_inv = ihx.solve_ihx_steady(m_inv, T1_needed_K, mdot1_ihx, T2_in_K, mdot2)
            T1_needed_K += (Q_target_W - float(r_inv["Q_W"])) / 3.7e5

        print("\n" + "-" * 76)
        print(f"[{label}] -> 질량유량 = {mdot2:.2f} kg/s (582 F 기준)")
        print("-" * 76)
        print(f"  1) 60.0 MW 흡수에 필요한 2차 출구 온도 : {T2_req_K-273.15:.2f} C (dT2 = +{T2_req_K-T2_in_K:.2f} K)")
        print(f"     고온단 온도차 (T1_in - T2_req)      : {pinch_K:+.2f} K  "
              f"{'-> [전제 하에서 성립 불가: 온도 역전 필요]' if pinch_K < 0 else '-> [성립 가능]'}")
        print(f"  2) 무한 면적 가정 최대 가능 열량 Q_max : {Q_max_W/1e6:.2f} MW")
        print(f"  3) ihx.py(N={model.N}) 실제 계산 열교환량 Q   : {float(res['Q_W'])/1e6:.2f} MW (60 MW 대비 {float(res['Q_W'])/1e6 - 60.0:+.2f} MW)")
        print(f"     이때 1차 출구 온도 T1_out           : {float(res['T1_out_K'])-273.15:.2f} C "
              f"(저온풀 343.78 C 대비 {float(res['T1_out_K'])-T_core_in_K:+.2f} K)")
        print(f"     이때 2차 출구 온도 T2_out           : {float(res['T2_out_K'])-273.15:.2f} C")
        print(f"  4) 60 MW 달성에 필요한 1차 입구 온도   : {T1_needed_K-273.15:.2f} C "
              f"(노심 출구 {T1_in_K-273.15:.2f} C 대비 {T1_needed_K-T1_in_K:+.2f} K)")

        results_summary[str(gpm)] = {
            "mdot2_kg_s": round(mdot2, 3),
            "T2_req_for_60MW_C": round(T2_req_K - 273.15, 2),
            "pinch_T1in_minus_T2req_K": round(pinch_K, 2),
            "second_law_valid": pinch_K > 0.0,
            "Q_max_theoretical_MW": round(Q_max_W / 1e6, 2),
            "ihx_Q_MW": round(float(res["Q_W"]) / 1e6, 2),
            "ihx_T1_out_C": round(float(res["T1_out_K"]) - 273.15, 2),
            "ihx_T2_out_C": round(float(res["T2_out_K"]) - 273.15, 2),
            "required_T1_in_for_60MW_C": round(T1_needed_K - 273.15, 2),
        }

    # 추가 확인: 만약 "gpm at XX F"에서 기준온도가 582 F가 아니라 최저 온도(예: 응고점 근처 250 F = 394 K)였다면?
    rho_250f = ihx.rho_na(f_to_k(250.0))
    mdot2_extreme = 5079.5 * GPM_TO_M3S * rho_250f
    T2_req_ext_K = required_t2_out_for_power(T2_in_K, mdot2_extreme, Q_target_W)
    print("\n" + "-" * 76)
    print("[보조 확인] 만약 표 3.2.4.2.1의 빠진 기준온도(XX F)가 극단적으로 낮았다면?")
    print(f"  250 F(121 C) 밀도({rho_250f:.1f} kg/m^3)를 써도 5079.5 gpm = {mdot2_extreme:.2f} kg/s")
    print(f"  이때도 60 MW 흡수 시 T2_out = {T2_req_ext_K-273.15:.2f} C > T1_in({T1_in_K-273.15:.2f} C) 이므로 여전히 성립 불가")
    print(f"  비율 확인: 5474.5 / 5079.5 = {5474.5/5079.5:.5f} (+7.78%)")
    # CSV 보간기 배율 옵션(SHRT45R_FLOW_SCALE) 검산
    bc_raw = ihx.load_secondary_bc("SHRT45R")
    bc_scaled = ihx.load_secondary_bc("SHRT45R", flow_scale=ihx.SHRT45R_FLOW_SCALE)
    _, mdot_raw_0 = bc_raw(0.0)
    _, mdot_scaled_0 = bc_scaled(0.0)
    print(f"  CSV 0초 유량 (원본 / 보정배율 적용): {mdot_raw_0:.3f} kg/s -> {mdot_scaled_0:.3f} kg/s")
    # 보조 근거: SHRT-17은 두 표 값이 거의 같다 (45R만 튀는 값)
    print(f"  [보조 근거] SHRT-17: 초기조건 표 5614.6 gpm vs 유량표 0초 5619.55 gpm (차이 {(5619.55/5614.6-1)*100:+.2f}%)")
    print(f"              SHRT-45R: 초기조건 표 5474.5 gpm vs 유량표 0초 5079.5 gpm (차이 {(5079.5/5474.5-1)*100:+.2f}%)")
    print("=" * 76)

    out_json = os.path.join(ROOT, "스크립트", "shrt45r_flow_check.json")
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(results_summary, f, indent=2, ensure_ascii=False)
    print(f"판정 수치 JSON 저장 완료: {out_json}")

    # 시각 비교 SVG 생성 (좌: 60 MW 흡수 요구 출구온도 및 온도역전 여부, 우: ihx.py 실제 열교환량)
    svg_path = os.path.join(ROOT, "도면", "shrt45r_flow_comparison.svg")
    c5079 = results_summary["5079.5"]
    c5474 = results_summary["5474.5"]
    T1_in_C = T1_in_K - 273.15

    w_svg, h_svg = 760, 410
    # 좌측 패널 (온도 범위 425 ~ 465 C)
    lx0, lx1, ly0, ly1 = 75, 355, 80, 315
    tmin, tmax = 425.0, 465.0

    def sy_t(tc: float) -> float:
        return ly1 - (tc - tmin) / (tmax - tmin) * (ly1 - ly0)

    # 우측 패널 (열교환량 범위 50 ~ 62 MW)
    rx0, rx1, ry0, ry1 = 445, 725, 80, 315
    qmin, qmax = 50.0, 62.0

    def sy_q(qmw: float) -> float:
        return ry1 - (qmw - qmin) / (qmax - qmin) * (ry1 - ry0)

    y_t1in = sy_t(T1_in_C)
    y_req_5079 = sy_t(c5079["T2_req_for_60MW_C"])
    y_req_5474 = sy_t(c5474["T2_req_for_60MW_C"])

    y_q60 = sy_q(60.0)
    y_q_5079 = sy_q(c5079["ihx_Q_MW"])
    y_q_5474 = sy_q(c5474["ihx_Q_MW"])

    svg_lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w_svg} {h_svg}" width="{w_svg}" height="{h_svg}" font-family="sans-serif">',
        f'<rect width="{w_svg}" height="{h_svg}" fill="#ffffff"/>',
        '<text x="380" y="28" text-anchor="middle" font-size="16" font-weight="bold" fill="#111827">SHRT-45R 2차측 초기 유량 후보 비교 (5079.5 gpm vs 5474.5 gpm)</text>',
        '<text x="380" y="47" text-anchor="middle" font-size="11" fill="#4b5563">사양서 초기 출력 60.0 MW 및 노심 혼합 출구 온도 444.05 °C 기준 정합성 판정</text>',
        # 좌측 패널 배경 및 눈금
        '<text x="215" y="68" text-anchor="middle" font-size="12" font-weight="bold" fill="#1f2937">[1] 60.0 MW 흡수에 필요한 2차 출구 온도 [°C]</text>',
        f'<rect x="{lx0}" y="{ly0}" width="{lx1-lx0}" height="{ly1-ly0}" fill="#f9fafb" stroke="#9ca3af" stroke-width="1"/>',
        # 온도 역전 금지 구역 (T > T1_in) 연분홍 음영
        f'<rect x="{lx0+1}" y="{ly0+1}" width="{lx1-lx0-2}" height="{y_t1in-ly0-1:.1f}" fill="#fee2e2" opacity="0.7"/>',
        f'<text x="{0.5*(lx0+lx1):.1f}" y="{ly0+18}" text-anchor="middle" font-size="10" font-weight="bold" fill="#991b1b">성립 불가 구역: 1차 입구(444.05 °C)보다 뜨거움</text>',
    ]
    for tc_tick in range(425, 466, 5):
        yy = sy_t(float(tc_tick))
        svg_lines.append(f'<line x1="{lx0}" y1="{yy:.1f}" x2="{lx1}" y2="{yy:.1f}" stroke="#e5e7eb" stroke-width="1"/>')
        svg_lines.append(f'<text x="{lx0-8}" y="{yy+4:.1f}" text-anchor="end" font-size="10" fill="#374151">{tc_tick}</text>')

    # 좌측 막대 2개
    svg_lines.extend([
        f'<rect x="115" y="{y_req_5079:.1f}" width="76" height="{ly1-y_req_5079:.1f}" fill="#ef4444" rx="3"/>',
        f'<text x="153" y="{y_req_5079-7:.1f}" text-anchor="middle" font-size="11" font-weight="bold" fill="#b91c1c">{c5079["T2_req_for_60MW_C"]:.2f} °C</text>',
        f'<text x="153" y="{y_req_5079+16:.1f}" text-anchor="middle" font-size="10" font-weight="bold" fill="#ffffff">역전 ({-c5079["pinch_T1in_minus_T2req_K"]:+.1f} K)</text>',
        f'<rect x="235" y="{y_req_5474:.1f}" width="76" height="{ly1-y_req_5474:.1f}" fill="#2563eb" rx="3"/>',
        f'<text x="273" y="{y_req_5474+18:.1f}" text-anchor="middle" font-size="11" font-weight="bold" fill="#ffffff">{c5474["T2_req_for_60MW_C"]:.2f} °C</text>',
        f'<text x="273" y="{y_req_5474+33:.1f}" text-anchor="middle" font-size="10" fill="#dbeafe">성립 ({-c5474["pinch_T1in_minus_T2req_K"]:+.2f} K)</text>',
        # 1차 입구 한계선
        f'<line x1="{lx0}" y1="{y_t1in:.1f}" x2="{lx1}" y2="{y_t1in:.1f}" stroke="#991b1b" stroke-width="2" stroke-dasharray="6,4"/>',
        f'<text x="{lx1-6}" y="{y_t1in-5:.1f}" text-anchor="end" font-size="10" font-weight="bold" fill="#991b1b">T1_in = {T1_in_C:.2f} °C</text>',
        '<text x="153" y="335" text-anchor="middle" font-size="11" font-weight="bold" fill="#111827">5079.5 gpm</text>',
        '<text x="153" y="350" text-anchor="middle" font-size="10" fill="#6b7280">(표 3.2.4.2.1 0초)</text>',
        '<text x="273" y="335" text-anchor="middle" font-size="11" font-weight="bold" fill="#1e40af">5474.5 gpm (채택)</text>',
        '<text x="273" y="350" text-anchor="middle" font-size="10" fill="#1e40af">(초기조건 표 3.2.1.1)</text>',
    ])

    # 우측 패널 배경 및 눈금
    svg_lines.extend([
        f'<text x="585" y="68" text-anchor="middle" font-size="12" font-weight="bold" fill="#1f2937">[2] ihx.py (N = {ihx.DEFAULT_N_NODES}) 정상상태 열교환량 [MW]</text>',
        f'<rect x="{rx0}" y="{ry0}" width="{rx1-rx0}" height="{ry1-ry0}" fill="#f9fafb" stroke="#9ca3af" stroke-width="1"/>',
    ])
    for q_tick in range(50, 63, 2):
        yy = sy_q(float(q_tick))
        svg_lines.append(f'<line x1="{rx0}" y1="{yy:.1f}" x2="{rx1}" y2="{yy:.1f}" stroke="#e5e7eb" stroke-width="1"/>')
        svg_lines.append(f'<text x="{rx0-8}" y="{yy+4:.1f}" text-anchor="end" font-size="10" fill="#374151">{q_tick}</text>')

    svg_lines.extend([
        f'<rect x="485" y="{y_q_5079:.1f}" width="76" height="{ry1-y_q_5079:.1f}" fill="#f87171" rx="3"/>',
        f'<text x="523" y="{y_q_5079-7:.1f}" text-anchor="middle" font-size="11" font-weight="bold" fill="#b91c1c">{c5079["ihx_Q_MW"]:.2f} MW</text>',
        f'<text x="523" y="{y_q_5079+18:.1f}" text-anchor="middle" font-size="10" font-weight="bold" fill="#ffffff">{c5079["ihx_Q_MW"]-60.0:+.2f} MW</text>',
        f'<rect x="605" y="{y_q_5474:.1f}" width="76" height="{ry1-y_q_5474:.1f}" fill="#2563eb" rx="3"/>',
        f'<text x="643" y="{y_q_5474-7:.1f}" text-anchor="middle" font-size="11" font-weight="bold" fill="#1e40af">{c5474["ihx_Q_MW"]:.2f} MW</text>',
        f'<text x="643" y="{y_q_5474+18:.1f}" text-anchor="middle" font-size="10" font-weight="bold" fill="#ffffff">{c5474["ihx_Q_MW"]-60.0:+.2f} MW</text>',
        # 60 MW 목표선
        f'<line x1="{rx0}" y1="{y_q60:.1f}" x2="{rx1}" y2="{y_q60:.1f}" stroke="#111827" stroke-width="2" stroke-dasharray="6,4"/>',
        f'<text x="{rx0+8}" y="{y_q60-6:.1f}" font-size="10" font-weight="bold" fill="#111827">목표 초기 출력 = 60.00 MW</text>',
        '<text x="523" y="335" text-anchor="middle" font-size="11" font-weight="bold" fill="#111827">5079.5 gpm</text>',
        '<text x="523" y="350" text-anchor="middle" font-size="10" fill="#6b7280">(280.97 kg/s)</text>',
        '<text x="643" y="335" text-anchor="middle" font-size="11" font-weight="bold" fill="#1e40af">5474.5 gpm (채택)</text>',
        '<text x="643" y="350" text-anchor="middle" font-size="10" fill="#1e40af">(302.82 kg/s)</text>',
        '<text x="380" y="378" text-anchor="middle" font-size="11" fill="#374151">결론: 5079.5 gpm으로는 60 MW 열수지가 성립하지 않아 초기조건 표의 5474.5 gpm(302.82 kg/s)을 채택</text>',
        '<text x="380" y="394" text-anchor="middle" font-size="10" fill="#6b7280">(전제: IHX 열량 = 노심 출력 60 MW, IHX 1차 입구 온도 = 노심 혼합 출구 온도)</text>',
        "</svg>",
    ])
    with open(svg_path, "w", encoding="utf-8") as f:
        f.write("\n".join(svg_lines))
    print(f"비교 그래프 SVG 저장 완료: {svg_path}")


if __name__ == "__main__":
    main()
