#!/usr/bin/env python3
"""
covalent5 · RC 최적화 v3 입력 생성 — floppy 모드 제거
삼중쌍 앵커(0,5,6,16,17,26) + 반응중심(P=30)에서 먼(>4A) VX 원자 고정 -> 로터 정지.
P-Ogamma(30-1) 3.3A 구속, LooseOpt, 순전하 -1. 입력 구조=rc_vx_p3s_v2.xyz(이전 53사이클 결과).
사용: python make_rc_v3.py
출력: rc_vx_p3s_v3.inp
"""
import math

GEOM = "rc_vx_p3s_v2.xyz"
L = [l.split() for l in open(GEOM).read().splitlines()[2:] if l.strip()]
xyz = [(a[0], list(map(float, a[1:4]))) for a in L]
P = xyz[30][1]
anchors = [0, 5, 6, 16, 17, 26]
far = [i for i in range(28, len(xyz)) if math.dist(xyz[i][1], P) > 4.0]
frozen = sorted(set(anchors + far))

lines = ["# covalent5 VX(P3S) RC opt v3 -- floppy 로터 고정 + LooseOpt",
         f"# 고정 {len(frozen)}원자. P-Ogamma 3.3A 구속. 순전하 -1.",
         "! wB97X-D3 def2-SVP def2/J RIJCOSX LooseOpt CPCM(water) TightSCF defgrid2",
         "%pal nprocs 8 end",
         "%maxcore 3500",
         "%geom Constraints"]
for i in frozen:
    lines.append(f"  {{C {i} C}}")
lines.append("  {B 30 1 3.30 C}")
lines.append(" end end")
lines.append(f"* xyzfile -1 1 {GEOM}")
open("rc_vx_p3s_v3.inp", "w").write("\n".join(lines) + "\n")
print(f"고정 원자수: {len(frozen)} (앵커6 + 먼VX {len(far)})")
print("-> rc_vx_p3s_v3.inp")
