"""EBR-II IHX 형상 및 파생값 계산.
출처: ANL-ARC-226 Rev.1 5.2.7절, 표 5.2.7.1 (p.128), 그림 5.2.7.3 (p.127).
'fig_' 접두 값은 그림 판독값이라 원문 그림으로 재확인 필요.
"""
import json, math
G = dict(N_tubes=3026, tube_OD_m=1.588e-2, tube_wall_m=0.124e-2, pitch_m=2.064e-2, tube_length_m=3.137,
         casing_length_m=5.639, casing_ID_m=1.295, sec_pipe_OD_m=0.324, sec_pipe_wall_m=0.635e-2,
         sec_inlet_plenum_m3=0.774, sec_outlet_plenum_m3=0.774, tubesheet_thk_m=0.0762, material="304SS",
         primary_inlet_elev_m=4.939, casing_bottom_elev_m=1.559,
         fig_lower_tubesheet_bottom_elev_m=2.232, fig_upper_tubesheet_elev_m=5.369,
         fig_primary_outlet_elev_m_approx="2.302~2.372")
Do, t, P, L, N = G['tube_OD_m'], G['tube_wall_m'], G['pitch_m'], G['tube_length_m'], G['N_tubes']
Di = Do - 2 * t
D = dict(tube_ID_m=Di, P_over_D=P / Do,
         tube_flow_area_m2=N * math.pi / 4 * Di**2,
         HT_area_outer_m2=N * math.pi * Do * L, HT_area_inner_m2=N * math.pi * Di * L,
         casing_area_m2=math.pi / 4 * G['casing_ID_m']**2,
         tubes_outer_area_m2=N * math.pi / 4 * Do**2,
         center_pipe_area_m2=math.pi / 4 * G['sec_pipe_OD_m']**2)
D['shell_flow_area_m2_assume_center_pipe'] = D['casing_area_m2'] - D['tubes_outer_area_m2'] - D['center_pipe_area_m2']
Pw = N * math.pi * Do + math.pi * G['casing_ID_m'] + math.pi * G['sec_pipe_OD_m']
D['shell_Dh_m'] = 4 * D['shell_flow_area_m2_assume_center_pipe'] / Pw
D['tube_Dh_m'] = Di
D['check_tube_span_from_fig_m'] = G['fig_upper_tubesheet_elev_m'] - G['fig_lower_tubesheet_bottom_elev_m']
json.dump(dict(given=G, derived=D), open('ihx_geometry.json', 'w'), indent=1, ensure_ascii=False)
for k, v in D.items(): print(k, round(v, 5))
