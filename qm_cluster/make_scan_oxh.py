#!/usr/bin/env python3
"""
covalent5 · His447 full-QM 1D relaxed scan 입력 생성
수렴 RC(rc_oxh.xyz)에서 P...Ser-Ogamma를 3.30->1.70 A 스캔(18점).
※ 양성자(Ser Ogamma-H)는 **자유** — 협주적 PT가 스캔 중 자발 발생하도록 구속하지 않음.
freeze = caps/anchors/farVX (oxh_indices.json). His Ne2...Ogamma 구속은 제거(PT 허용).
사용: python make_scan_oxh.py        (rc_oxh.xyz, oxh_indices.json 필요)
출력: scan_oxh.inp
"""
import json
m = json.load(open("oxh_indices.json"))
P, OG = m["P"], m["Ogamma"]
frozen = m["frozen"]
L = ["# covalent5 His447 full-QM P...Ogamma scan 3.3->1.7A (proton free). net -1.",
     "! wB97X-D3 def2-SVP def2/J RIJCOSX LooseOpt CPCM(water) TightSCF defgrid2",
     "%pal nprocs 8 end", "%maxcore 3500", "%geom",
     f"  Scan B {P} {OG} = 3.30, 1.70, 18 end", "  Constraints"]
for i in frozen:
    L.append(f"    {{C {i} C}}")
L += ["   end", " end", "* xyzfile -1 1 rc_oxh.xyz"]
open("scan_oxh.inp", "w").write("\n".join(L) + "\n")
print(f"scan B {P} {OG} = 3.3->1.7 (18 pts), proton free, frozen {len(frozen)} -> scan_oxh.inp")
