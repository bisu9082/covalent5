#!/usr/bin/env python3
"""
covalent5 · P...Ogamma relaxed scan 입력 생성 — VX(P3S) 알콕사이드
rc_alkoxide.xyz(수렴 RC)에서 P(29)-Ser O-(1)를 3.3->1.6A 스캔(18점). 각 점 제약최적화.
frozen=삼중쌍앵커[0,4,5,15,16,25]+먼VX(>4A). 순전하 -1. 출력: scan_alkoxide.inp
"""
import math
L = [l.split() for l in open("rc_alkoxide.xyz").read().splitlines()[2:] if l.strip()]
xyz = [(a[0], list(map(float, a[1:4]))) for a in L]
P, OG = 29, 1
anchors = [0, 4, 5, 15, 16, 25]
far = [i for i in range(27, 69) if math.dist(xyz[i][1], xyz[P][1]) > 4.0]
frozen = sorted(set(anchors + far))
lines = ["# covalent5 VX(P3S) alkoxide P...Ogamma scan 3.3->1.6A. net charge -1.",
         "! wB97X-D3 def2-SVP def2/J RIJCOSX LooseOpt CPCM(water) TightSCF defgrid2",
         "%pal nprocs 8 end", "%maxcore 3500",
         "%geom",
         f"  Scan B {P} {OG} = 3.30, 1.60, 18 end",
         "  Constraints"]
for i in frozen:
    lines.append(f"    {{C {i} C}}")
lines.append("   end")
lines.append(" end")
lines.append("* xyzfile -1 1 rc_alkoxide.xyz")
open("scan_alkoxide.inp", "w").write("\n".join(lines) + "\n")
print(f"scan B {P} {OG} = 3.3->1.6 (18 pts), frozen {len(frozen)} (anchor6 + farVX {len(far)})")
print("-> scan_alkoxide.inp")
