#!/usr/bin/env python3
"""
covalent5 · 도킹 best pose 정렬 렌더 (PyMOL)
모든 비공유 best pose를 hAChE 활성부위(Ser203) 기준 동일 시점으로 렌더 → 비교 가능한 패널.
요구: pymol-open-source
실행 (반드시 docking/ 에서):
  export COV5="/mnt/c/Users/kkaan/Downloads/claude_research/covalent5"
  cd "$COV5/docking"
  pymol -cq render_poses.py
출력: ../figure/poses/<sid>.png
"""
import os, glob, gzip, tempfile
from pymol import cmd

# PyMOL에서 __file__이 불안정 → 현재 작업 디렉토리(docking/) 기준
HERE = os.getcwd()
REC = os.path.join(HERE, "receptor_4EY7_chainA.pdb")
POSES = sorted(glob.glob(os.path.join(HERE, "out_noncovalent", "*_nc.sdf.gz")))
OUT = os.path.join(os.path.dirname(HERE), "figure", "poses")
os.makedirs(OUT, exist_ok=True)

if not os.path.exists(REC):
    raise SystemExit(f"수용체 없음: {REC}\n→ docking/ 폴더에서 실행하세요 (cd \"$COV5/docking\")")
if not POSES:
    raise SystemExit("out_noncovalent/*_nc.sdf.gz 없음 → 도킹 먼저")

cmd.set("ray_opaque_background", 0)
cmd.set("ray_shadows", 0)
cmd.bg_color("white")

cmd.load(REC, "rec")
cmd.hide("everything")
cmd.show("cartoon", "rec")
cmd.color("gray85", "rec")
cmd.set("cartoon_transparency", 0.25, "rec")
cmd.select("triad", "rec and resi 203+334+447")
cmd.show("sticks", "triad")
cmd.color("salmon", "triad")
cmd.set("stick_radius", 0.18, "triad")
cmd.orient("rec and resi 203+334+447+86+286")
cmd.turn("y", 15)
view = cmd.get_view()

n = 0
for gz in POSES:
    sid = os.path.basename(gz).replace("_nc.sdf.gz", "")
    tmp = tempfile.NamedTemporaryFile(suffix=".sdf", delete=False)
    with gzip.open(gz, "rt") as f:
        tmp.write(f.read().encode()); tmp.close()
    cmd.load(tmp.name, "lig", state=1)
    cmd.hide("everything", "lig")
    cmd.show("sticks", "lig")
    cmd.color("teal", "lig and elem C")
    cmd.set("stick_radius", 0.20, "lig")
    cmd.set_view(view)
    png = os.path.join(OUT, sid + ".png")
    cmd.ray(1400, 1050)
    cmd.png(png, dpi=300)
    cmd.delete("lig")
    os.unlink(tmp.name)
    n += 1
    print(f"[{n}/{len(POSES)}] {sid} -> figure/poses/{sid}.png")

print(f"완료: {n}개 -> figure/poses/")
