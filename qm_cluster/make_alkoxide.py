#!/usr/bin/env python3
"""
covalent5 · Ser-알콕사이드 모델 클러스터 + RC 입력 생성 (A안)
v3 수렴구조에서 Ser Ogamma-H(HG, orig idx4) 제거 -> Ser-O- ; His NE2에 HE2 추가 -> HIP.
순전하 -1. frozen=삼중쌍앵커+먼VX, P-Ogamma 3.3 구속.
사용: python make_alkoxide.py
"""
import math

def unit(v):
    n = math.sqrt(sum(c * c for c in v))
    return [c / n for c in v]

orig = []
for ln in open("rc_vx_p3s_v3.xyz").read().splitlines()[2:]:
    if ln.strip():
        p = ln.split()
        orig.append((p[0], list(map(float, p[1:4]))))

keep = [i for i in range(len(orig)) if i != 4]
new = [orig[i] for i in keep]
o2n = {o: n for n, o in enumerate(keep)}
OG = o2n[1]
P = o2n[30]
anchors = [o2n[i] for i in [0, 5, 6, 16, 17, 26]]
vx_new = [o2n[i] for i in range(28, len(orig))]

ne2 = orig[11][1]
cd2 = orig[9][1]
ce1 = orig[10][1]
mid = [(cd2[k] + ce1[k]) / 2 for k in range(3)]
d = unit([ne2[k] - mid[k] for k in range(3)])
he2 = [ne2[k] + 1.0 * d[k] for k in range(3)]
new.append(("H", he2))

far = [i for i in vx_new if math.dist(new[i][1], new[P][1]) > 4.0]
frozen = sorted(set(anchors + far))

with open("cluster_alkoxide.xyz", "w") as f:
    f.write(f"{len(new)}\nSer-alkoxide: Ser-O- + His(HIP) + Glu- + VX. net charge -1\n")
    for el, c in new:
        f.write(f"{el:2s} {c[0]:12.6f} {c[1]:12.6f} {c[2]:12.6f}\n")

lines = ["# covalent5 VX(P3S) Ser-alkoxide RC opt. net charge -1.",
         "! wB97X-D3 def2-SVP def2/J RIJCOSX LooseOpt CPCM(water) TightSCF defgrid2",
         "%pal nprocs 8 end", "%maxcore 3500", "%geom Constraints"]
for i in frozen:
    lines.append(f"  {{C {i} C}}")
lines.append(f"  {{B {P} {OG} 3.30 C}}")
lines.append(" end end")
lines.append("* xyzfile -1 1 cluster_alkoxide.xyz")
open("rc_alkoxide.inp", "w").write("\n".join(lines) + "\n")

print(f"atoms {len(new)} (Ser HG removed, His HE2 added)")
print(f"new index: P={P}, Ser Ogamma={OG}")
print(f"frozen {len(frozen)} (anchor6 + farVX {len(far)})")
print("-> cluster_alkoxide.xyz, rc_alkoxide.inp")
