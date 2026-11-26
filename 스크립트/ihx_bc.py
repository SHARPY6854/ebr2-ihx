"""IHX 2차측(중간 루프) 입구 경계조건을 사양서 PDF 표에서 직접 읽어 SI CSV로 변환.
출처: ANL-ARC-226 Rev.1 표 3.1.4.2.1~2 (p.21), 3.2.4.2.1 (p.36), 3.2.4.2.2 (p.37).
SHRT-45R 유량표 헤더는 'gpm at XX°F'로 기준온도가 빠져 있음 → SHRT-17과 같은 582°F로 가정.
"""
import csv, re, subprocess, sys
P = sys.argv[1]
GPM = 6.30901964e-5
def F2K(f): return (f - 32) * 5 / 9 + 273.15
def rho(T): return 1004.23 - 0.21390 * T - 1.1046e-5 * T**2   # SAS4A Eq.12.13-5
R582 = rho(F2K(582))
def rows(page, title):
    t = subprocess.run(['pdftotext', '-layout', '-f', str(page), '-l', str(page), P, '-'], capture_output=True, text=True).stdout
    seg = t.split(title, 1)[1]
    out = []
    for line in seg.splitlines()[2:]:
        n = re.findall(r'-?\d+\.?\d*', line)
        if len(n) == 2 and not re.search(r'[A-Za-z]', line): out.append([float(x) for x in n])
        elif out and line.strip(): break
    return out
spec = {
 'SHRT17_flow': (37, 'Table 3.1.4.2.1.'), 'SHRT17_temp': (37, 'Table 3.1.4.2.2.'),
 'SHRT45R_flow': (52, 'Table 3.2.4.2.1.'), 'SHRT45R_temp': (53, 'Table 3.2.4.2.2.'),
}
for case in ('SHRT17', 'SHRT45R'):
    fl = rows(*spec[case + '_flow']); tp = rows(*spec[case + '_temp'])
    with open(f'bc_ihx_secondary_flow_{case}.csv', 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['time_s', 'flow_gpm_at_582F', 'flow_m3_s', 'mdot_kg_s'])
        for t, q in fl: w.writerow([t, q, round(q * GPM, 6), round(q * GPM * R582, 3)])
    with open(f'bc_ihx_secondary_temp_{case}.csv', 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['time_s', 'T_F', 'T_K', 'T_C'])
        for t, T in tp: w.writerow([t, T, round(F2K(T), 3), round(F2K(T) - 273.15, 3)])
    print(case, 'flow rows', len(fl), 'temp rows', len(tp), 'first', fl[0], tp[0], 'last', fl[-1], tp[-1])
print('rho_582F', round(R582, 2))
