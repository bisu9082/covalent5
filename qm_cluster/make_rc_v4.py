#!/usr/bin/env python3
"""
covalent5 · RC opt v4 — 촉매 H결합 복원 (His NE2 - Ser Ogamma ~2.9A)
v3 수렴구조에서 His447이 회전해 NE2가 Ser Ogamma에서 4.59A로 멀어짐 -> 구속으로 복원.
사용: python make_rc_v4.py  ->  rc_vx_p3s_v4.inp
"""
import math
GEOM = "rc_vx_p3s_v3.xyz"
L = [l.split() for l in open(GEOM).read().splitlines()[2:] if l.strip()]
xyz = [(a[0], list(map(float, a[1:4]))) for a in L]
P = xyz[30][1]
anchors = [0, 5, 6, 16, 17, 26]
far = [i for i in range(28, len(xyz)) if math.dist(xyz[i][1], P) > 4.0]
frozen = sorted(set(anchors + far))
lines = ["# covalent5 VX(P3S) RC opt v4 -- 촉매 H결합 복원",
         "! wB97X-D3 def2-SVP def2/J RIJCOSX LooseOpt CPCM(water) TightSCF defgrid2",
         "%pal nprocs 8 end", "%maxcore 3500", "%geom Constraints"]
for i in frozen:
    lines.append(f"  {{C {i} C}}")
lines.append("  {B 30 1 3.30 C}")
lines.append("  {B 11 1 2.90 C}")
lines.append(" end end")
lines.append(f"* xyzfile -1 1 {GEOM}")
open("rc_vx_p3s_v4.inp", "w").write("\n".join(lines) + "\n")
print(f"고정 {len(frozen)} + 구속 2개(P-Ogamma 3.3, His NE2-Ogamma 2.9)")
print("-> rc_vx_p3s_v4.inp")
