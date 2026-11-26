"""3단계: SHRT-17 초기조건 정상상태 점검 스크립트.

점검 목표:
1. SHRT-17 초기 운전 조건에서 열교환량이 약 57 MW 근처로 재현되는지 확인
2. 2차측 나트륨 온도상승(T2_out - T2_in)이 약 145 K 근처인지 확인
3. 1차측 방출열, 벽 통과열, 2차측 흡수열 간 에너지 보존 오차 정량 출력
4. 격자 수렴(N = 10~320), epsilon-NTU 해석해 대조, 열전달 길이(active/full) 및 상관식 민감도 확인

기준 출처: ANL-ARC-226 Rev.1 표 3.1.1.1 (p.9, p.11), 표 3.1.4.2.1~2 (p.19, p.21), 표 5.3.1.1 (p.129)
"""

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


def main() -> None:
    rho_800f = ihx.rho_na(f_to_k(800.0))
    rho_582f = ihx.rho_na(f_to_k(582.0))

    # 1. 사양서 SHRT-17 초기조건 (표 3.1.1.1 p.9, p.11 및 표 3.1.4.2.1 p.19)
    Q_spec_W = 57.29e6  # 초기 출력: 핵분열 53.93 MW + 붕괴열 3.36 MW = 57.29 MW (표 3.1.1.1은 57.3 MW)
    mdot_core = 8500.0 * GPM_TO_M3S * rho_800f  # 455.36 kg/s
    # 표 5.3.1.1 (p.129, SHRT-17): IHX 입구 8438.5 gpm, IHX 출구 8421.0 gpm.
    # 두 값의 차이 17.5 gpm은 IHX 입구 누설(p.129~130 누설 목록)이다. 누설이 번들을 거치지 않는다고 가정하고 출구 값 8421.0 gpm을 번들 유량으로 쓴다.
    mdot1_ihx = 8421.0 * GPM_TO_M3S * rho_800f   # 451.13 kg/s
    T_core_in_K = f_to_k(664.8)                 # 624.71 K (351.56 C)

    # 노심 혼합 평균 출구 온도 계산 (Z-pipe 단열 가정 시 IHX 1차측 입구 온도 추정값)
    T_core_out_K = T_core_in_K + Q_spec_W / (mdot_core * ihx.cp_na(T_core_in_K + 50.0))
    for _ in range(15):
        T_mid = 0.5 * (T_core_in_K + T_core_out_K)
        T_core_out_K = T_core_in_K + Q_spec_W / (mdot_core * ihx.cp_na(T_mid))

    # 2차측 초기조건: 표 3.1.1.1 (5614.6 gpm, 573.8 F) vs CSV 0초 (5619.55 gpm, 574.02 F)
    T2_in_tab_K = f_to_k(573.8)
    mdot2_tab = 5614.6 * GPM_TO_M3S * rho_582f
    bc17 = ihx.load_secondary_bc("SHRT17")
    T2_in_csv_K, mdot2_csv = bc17(0.0)

    print("=" * 72)
    print("3단계: SHRT-17 초기조건 정상상태 점검 결과")
    print("=" * 72)
    print(f"[입력 조건 요약]")
    print(f"  사양서 초기 열출력         : {Q_spec_W/1e6:.2f} MW (표 3.1.4.2.1 p.19 / 표 3.1.1.1: 57.3 MW)")
    print(f"  노심 입구(=저온풀) 온도    : {T_core_in_K-273.15:.2f} C ({T_core_in_K:.2f} K, 664.8 F)")
    print(f"  노심 혼합 출구 온도(계산)  : {T_core_out_K-273.15:.2f} C ({T_core_out_K:.2f} K)")
    print(f"  1차측 유량 (분배율 반영)   : {mdot1_ihx:.2f} kg/s (노심 전량 기준 시 {mdot_core:.2f} kg/s)")
    print(f"  2차측 입구 (표 3.1.1.1)    : {T2_in_tab_K-273.15:.2f} C, {mdot2_tab:.2f} kg/s (5614.6 gpm)")
    print(f"  2차측 입구 (CSV 0.0 s)     : {T2_in_csv_K-273.15:.2f} C, {mdot2_csv:.2f} kg/s (5619.55 gpm)")

    # [검산 1] 노심 혼합 출구 온도를 그대로 IHX 1차측 입구에 넣었을 때
    model_a = ihx.create_ihx()
    res_a = ihx.solve_ihx_steady(model_a, T_core_out_K, mdot1_ihx, T2_in_tab_K, mdot2_tab)
    dT2_a = res_a["T2_out_K"] - res_a["T2_in_K"]
    dT1_a = res_a["T1_in_K"] - res_a["T1_out_K"]

    print("\n" + "-" * 72)
    print(f"[검산 1] 노심 출구 온도({T_core_out_K-273.15:.2f} C) 직결 입력 시 IHX 정상상태 계산 결과 (N = {model_a.N}, {model_a.length_mode}, L = {model_a.L:.3f} m)")
    print("-" * 72)
    print(f"  열교환량 Q (벽 -> 2차측)   : {res_a['Q_W']/1e6:.3f} MW  (사양서 57.29 MW 대비 {(res_a['Q_W']-Q_spec_W)/Q_spec_W*100:+.2f}%)")
    print(f"  1차측 대류 방출열 Q_1_adv  : {res_a['Q_1_adv_W']/1e6:.6f} MW")
    print(f"  튜브벽 전달열 Q_1toW       : {res_a['Q_1toW_W']/1e6:.6f} MW")
    print(f"  2차측 대류 흡수열 Q_2_adv  : {res_a['Q_2_adv_W']/1e6:.6f} MW")
    print(f"  에너지 보존 오차 |Q1 - Q2| : {abs(res_a['energy_residual_W']):.3e} W (상대오차 {abs(res_a['energy_residual_W'])/res_a['Q_W']*100:.2e}%)")
    print(f"  1차측 온도 변화            : {res_a['T1_in_K']-273.15:.2f} C -> {res_a['T1_out_K']-273.15:.2f} C (dT1 = -{dT1_a:.2f} K, 저온풀 351.56 C 대비 {res_a['T1_out_K']-T_core_in_K:+.2f} K)")
    print(f"  2차측 온도 변화            : {res_a['T2_in_K']-273.15:.2f} C -> {res_a['T2_out_K']-273.15:.2f} C (온도상승 dT2 = +{dT2_a:.2f} K)")

    # [검산 2] CSV 0초 경계조건 기준 동일 계산
    model_b = ihx.create_ihx()
    res_b = ihx.solve_ihx_steady(model_b, T_core_out_K, mdot1_ihx, T2_in_csv_K, mdot2_csv)
    dT2_b = res_b["T2_out_K"] - res_b["T2_in_K"]
    print(f"\n[검산 2] CSV 0.0초 경계조건 입력 시:")
    print(f"  열교환량 Q = {res_b['Q_W']/1e6:.3f} MW, 2차측 출구 = {res_b['T2_out_K']-273.15:.2f} C (dT2 = +{dT2_b:.2f} K), 에너지 오차 = {abs(res_b['energy_residual_W']):.3e} W")

    # [검산 3] 정확히 Q = 57.29 MW가 되도록 하는 1차측 입구 온도 역산
    T1_match_K = T_core_out_K
    model_c = ihx.create_ihx()
    for _ in range(20):
        res_c = ihx.solve_ihx_steady(model_c, T1_match_K, mdot1_ihx, T2_in_tab_K, mdot2_tab)
        T1_match_K += (Q_spec_W - res_c["Q_W"]) / 3.8e5
    dT2_c = res_c["T2_out_K"] - res_c["T2_in_K"]
    print(f"\n[검산 3] 정확히 Q = 57.29 MW 전달 시 역산 결과:")
    print(f"  필요한 1차측 입구 온도     : {T1_match_K-273.15:.2f} C (노심 혼합출구 {T_core_out_K-273.15:.2f} C 대비 {T1_match_K-T_core_out_K:+.2f} K)")
    print(f"  이때 1차측 출구 온도       : {res_c['T1_out_K']-273.15:.2f} C (저온풀 351.56 C 대비 {res_c['T1_out_K']-T_core_in_K:+.2f} K 차이)")
    print(f"  이때 2차측 출구 온도       : {res_c['T2_out_K']-273.15:.2f} C (2차측 온도상승 dT2 = +{dT2_c:.2f} K)")

    # [검산 4] 격자 수렴 + epsilon-NTU 해석해 대조 (열전달 길이 active / full)
    print("\n" + "-" * 72)
    print("[검산 4] 격자 수렴 및 epsilon-NTU 해석해 대조 (1차 상류 차분)")
    print("-" * 72)
    grid_rows = []
    for mode in ("active", "full"):
        q_prev = None
        for n_nodes in (10, 20, 40, 80, 160, 320):
            m_g = ihx.create_ihx(n_nodes=n_nodes, length_mode=mode)
            r_g = ihx.solve_ihx_steady(m_g, T_core_out_K, mdot1_ihx, T2_in_tab_K, mdot2_tab)
            row = dict(length_mode=mode, L_m=round(m_g.L, 4), N=n_nodes, Q_MW=round(r_g["Q_W"] / 1e6, 4),
                       T1_out_C=round(r_g["T1_out_K"] - 273.15, 3), T2_out_C=round(r_g["T2_out_K"] - 273.15, 3))
            if n_nodes == 320:
                chk = ihx.effectiveness_ntu_check(m_g, T_core_out_K, mdot1_ihx, T2_in_tab_K, mdot2_tab)
                row["Q_eNTU_MW"] = round(chk["Q_eNTU_W"] / 1e6, 4)
                row["NTU"] = round(chk["NTU"], 3)
            grid_rows.append(row)
            print(f"  {mode:6s} L={m_g.L:.3f} m  N={n_nodes:3d} -> Q={r_g['Q_W']/1e6:.3f} MW "
                  f"({(r_g['Q_W']-Q_spec_W)/Q_spec_W*100:+.2f}% vs 57.29 MW), T1_out={r_g['T1_out_K']-273.15:.2f} C"
                  + (f", Q_eNTU={row['Q_eNTU_MW']:.3f} MW (NTU {row['NTU']:.2f})" if n_nodes == 320 else ""))
        q160 = [r for r in grid_rows if r["length_mode"] == mode and r["N"] == 160][0]["Q_MW"]
        q320 = [r for r in grid_rows if r["length_mode"] == mode and r["N"] == 320][0]["Q_MW"]
        print(f"  {mode:6s} Richardson 외삽(1차 정확도, N=160/320): Q_inf = {2*q320-q160:.3f} MW")

    print("\n" + "-" * 72)
    print(f"[검산 4-2] 상관식 민감도 (N = {ihx.DEFAULT_N_NODES}, active)")
    print("-" * 72)
    sens_rows = []
    for s_corr, t_corr, mode in [
        ("mikityuk", "seban", "active"),
        ("ushakov", "seban", "active"),
        ("schad", "seban", "active"),
        ("mikityuk", "lyon", "active"),
        ("mikityuk", "seban", "full"),
    ]:
        m_s = ihx.create_ihx(shell_corr=s_corr, tube_corr=t_corr, length_mode=mode)
        r_s = ihx.solve_ihx_steady(m_s, T_core_out_K, mdot1_ihx, T2_in_tab_K, mdot2_tab)
        dt2_s = r_s["T2_out_K"] - T2_in_tab_K
        print(
            f"  쉘={s_corr:8s}, 튜브={t_corr:5s}, 길이={mode:6s} -> "
            f"Q={r_s['Q_W']/1e6:.3f} MW, T1_out={r_s['T1_out_K']-273.15:.2f} C, "
            f"T2_out={r_s['T2_out_K']-273.15:.2f} C (dT2={dt2_s:.2f} K)"
        )
        sens_rows.append(
            dict(
                N=m_s.N,
                length_mode=mode,
                shell_corr=s_corr,
                tube_corr=t_corr,
                Q_MW=round(r_s["Q_W"] / 1e6, 4),
                T1_out_C=round(r_s["T1_out_K"] - 273.15, 3),
                T2_out_C=round(r_s["T2_out_K"] - 273.15, 3),
                dT2_K=round(dt2_s, 3),
                energy_err_W=abs(r_s["energy_residual_W"]),
            )
        )

    # [검산 5] 층별 온도 분포 요약 출력 (상단 j=1, 5, 10, 15, 20)
    print("\n" + "-" * 72)
    print(f"[검산 5] {model_a.N}층 내부 온도 분포 요약 (검산 1 기준, 일부 층만 표시)")
    print("-" * 72)
    print("  층(j) | 높이 z[m] | 1차측 T1[C] | 튜브벽 Tw[C] | 2차측 T2[C] | T1-T2[K]")
    nA = model_a.N
    for idx in sorted(set([0, nA // 4, nA // 2, 3 * nA // 4, nA - 1])):
        t1c = res_a["T1_K"][idx] - 273.15
        twc = res_a["Tw_K"][idx] - 273.15
        t2c = res_a["T2_K"][idx] - 273.15
        z_m = res_a["z_nodes_m"][idx]
        print(f"  {idx+1:4d}  |   {z_m:5.3f}   |    {t1c:6.2f}   |    {twc:6.2f}    |    {t2c:6.2f}   |  {t1c-t2c:6.2f}")

    out_json = {
        "direct_core_outlet_case": {
            "T1_in_C": round(T_core_out_K - 273.15, 3),
            "mdot1_kg_s": round(mdot1_ihx, 3),
            "T2_in_C": round(T2_in_tab_K - 273.15, 3),
            "mdot2_kg_s": round(mdot2_tab, 3),
            "Q_MW": round(res_a["Q_W"] / 1e6, 4),
            "T1_out_C": round(res_a["T1_out_K"] - 273.15, 3),
            "T2_out_C": round(res_a["T2_out_K"] - 273.15, 3),
            "dT2_K": round(dT2_a, 3),
            "energy_residual_W": abs(res_a["energy_residual_W"]),
        },
        "exact_57p29MW_match_case": {
            "required_T1_in_C": round(T1_match_K - 273.15, 3),
            "T1_out_C": round(res_c["T1_out_K"] - 273.15, 3),
            "T2_out_C": round(res_c["T2_out_K"] - 273.15, 3),
            "dT2_K": round(dT2_c, 3),
            "energy_residual_W": abs(res_c["energy_residual_W"]),
        },
        "grid_convergence": grid_rows,
        "sensitivity": sens_rows,
    }
    out_path = os.path.join(ROOT, "스크립트", "steady_shrt17_results.json")
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(out_json, f, indent=2, ensure_ascii=False)

    # SVG 그래프 생성 (층별 온도 분포: 가로축 온도 [C], 세로축 높이 z [m])
    svg_path = os.path.join(ROOT, "도면", "shrt17_steady_profile.svg")
    w_svg, h_svg = 720, 440
    px0, px1 = 95, 660
    py0, py1 = 65, 375
    t_min, t_max = 295.0, 460.0
    z_min, z_max = model_a.z_bot, model_a.z_top

    def sx(tc: float) -> float:
        return px0 + (tc - t_min) / (t_max - t_min) * (px1 - px0)

    def sy(zm: float) -> float:
        return py1 - (zm - z_min) / (z_max - z_min) * (py1 - py0)

    z_arr = res_a["z_nodes_m"]
    t1_arr = [t - 273.15 for t in res_a["T1_K"]]
    tw_arr = [t - 273.15 for t in res_a["Tw_K"]]
    t2_arr = [t - 273.15 for t in res_a["T2_K"]]

    def poly(ts: list[float]) -> str:
        return " ".join(f"{sx(tc):.1f},{sy(zm):.1f}" for tc, zm in zip(ts, z_arr))

    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w_svg} {h_svg}" width="{w_svg}" height="{h_svg}" font-family="sans-serif">',
        f'<rect width="{w_svg}" height="{h_svg}" fill="#ffffff"/>',
        f'<text x="375" y="28" text-anchor="middle" font-size="16" font-weight="bold" fill="#111827">SHRT-17 초기 정상상태 IHX 높이별 온도 분포 (N = {model_a.N}, Q = {res_a["Q_W"]/1e6:.2f} MW)</text>',
        '<text x="375" y="46" text-anchor="middle" font-size="11" fill="#4b5563">빨간선: 1차측(쉘, 하향류) | 회색 점선: 304SS 튜브 벽 | 파란선: 2차측(튜브, 상향류)</text>',
        f'<rect x="{px0}" y="{py0}" width="{px1-px0}" height="{py1-py0}" fill="#f9fafb" stroke="#9ca3af" stroke-width="1"/>',
    ]
    for tc_tick in range(300, 461, 20):
        xx = sx(float(tc_tick))
        lines.append(f'<line x1="{xx:.1f}" y1="{py0}" x2="{xx:.1f}" y2="{py1}" stroke="#e5e7eb" stroke-width="1"/>')
        lines.append(f'<text x="{xx:.1f}" y="{py1+18}" text-anchor="middle" font-size="11" fill="#374151">{tc_tick}</text>')
    for zm_tick in [2.5, 3.0, 3.5, 4.0, 4.5]:
        yy = sy(zm_tick)
        lines.append(f'<line x1="{px0}" y1="{yy:.1f}" x2="{px1}" y2="{yy:.1f}" stroke="#e5e7eb" stroke-width="1"/>')
        lines.append(f'<text x="{px0-10}" y="{yy+4:.1f}" text-anchor="end" font-size="11" fill="#374151">{zm_tick:.1f}</text>')

    lines.append(f'<polyline fill="none" stroke="#dc2626" stroke-width="2.5" points="{poly(t1_arr)}"/>')
    lines.append(f'<polyline fill="none" stroke="#4b5563" stroke-width="2.0" stroke-dasharray="6,4" points="{poly(tw_arr)}"/>')
    lines.append(f'<polyline fill="none" stroke="#2563eb" stroke-width="2.5" points="{poly(t2_arr)}"/>')

    step_c = max(1, len(z_arr) // 20)
    for tc, zm in list(zip(t1_arr, z_arr))[::step_c]:
        lines.append(f'<circle cx="{sx(tc):.1f}" cy="{sy(zm):.1f}" r="3" fill="#dc2626"/>')
    for tc, zm in list(zip(t2_arr, z_arr))[::step_c]:
        lines.append(f'<circle cx="{sx(tc):.1f}" cy="{sy(zm):.1f}" r="3" fill="#2563eb"/>')

    # 주석 라벨
    lines.append(f'<text x="{sx(t2_arr[0])-14:.1f}" y="{sy(z_arr[0])+34:.1f}" text-anchor="end" font-size="10" font-weight="bold" fill="#991b1b">1차 입구단 ({t1_arr[0]:.1f} °C)</text>')
    lines.append(f'<text x="{sx(t2_arr[0])-8:.1f}" y="{sy(z_arr[0])+14:.1f}" text-anchor="end" font-size="10" font-weight="bold" fill="#1e40af">2차 출구 ({t2_arr[0]:.1f} °C)</text>')
    lines.append(f'<text x="{sx(t1_arr[-1])+8:.1f}" y="{sy(z_arr[-1])-4:.1f}" font-size="10" font-weight="bold" fill="#991b1b">1차 출구 ({t1_arr[-1]:.1f} °C)</text>')
    lines.append(f'<text x="{sx(t2_arr[-1])+8:.1f}" y="{sy(z_arr[-1])+4:.1f}" font-size="10" font-weight="bold" fill="#1e40af">2차 입구단 ({t2_arr[-1]:.1f} °C)</text>')

    lines.append(f'<text x="{0.5*(px0+px1):.1f}" y="{h_svg-14}" text-anchor="middle" font-size="12" font-weight="bold" fill="#111827">온도 [°C]</text>')
    lines.append(f'<text x="24" y="{0.5*(py0+py1):.1f}" text-anchor="middle" font-size="12" font-weight="bold" fill="#111827" transform="rotate(-90 24 {0.5*(py0+py1):.1f})">높이 z [m] (노심 하단 기준)</text>')
    lines.append("</svg>")
    with open(svg_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))

    print("=" * 72)
    print(f"결과 JSON 저장 완료: {out_path}")
    print(f"온도분포 그래프 저장 완료: {svg_path}")


if __name__ == "__main__":
    main()
