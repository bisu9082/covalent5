#!/usr/bin/env python3
"""
covalent5 · DFT 최적화 구조 3D 렌더 (PyMOL, ball-and-stick, 수소 제외, 원소색 고정) — Fig 2 (a)용
실행 ($COV5 에서):
  cd "$COV5"
  export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"
  pymol -cq render_molecules.py
출력: figure/mol3d/<name>.png
원소색(범례와 일치): C #999999, N blue, O red, P orange, F cyan
"""
import os
from pymol import cmd

HERE = os.getcwd()
OUT = os.path.join(HERE, "figure", "mol3d")
os.makedirs(OUT, exist_ok=True)

reps = [("A-230", "A230", "A230_01_P8S_C5N7E"),
        ("A-232", "A232", "A232_00_P8R_C5N7E"),
        ("A-234", "A234", "A234_00_P3S_N6C7E"),
        ("A-242", "A242", "A242_00_P7R")]

cmd.bg_color("white")
cmd.set("ray_opaque_background", 0)
cmd.set("ray_shadows", 0)
cmd.set("valence", 1)

for name, code, sid in reps:
    xyz = os.path.join(HERE, "dft", code, sid, "opt.xyz")
    if not os.path.exists(xyz):
        print(f"[없음] {xyz}"); continue
    cmd.delete("all")
    cmd.load(xyz, "mol")
    cmd.remove("hydro")
    cmd.hide("everything")
    cmd.show("spheres")
    cmd.show("sticks")
    cmd.set("sphere_scale", 0.25)
    cmd.set("stick_radius", 0.14)
    # 원소색 고정 (범례와 일치)
    cmd.color("grey60", "elem C")
    cmd.color("blue", "elem N")
    cmd.color("red", "elem O")
    cmd.color("orange", "elem P")
    cmd.color("cyan", "elem F")
    cmd.orient("mol")
    cmd.zoom("mol", 1.0)
    png = os.path.join(OUT, name + ".png")
    cmd.ray(1200, 900)
    cmd.png(png, dpi=300)
    print(f"{name} ({sid}) -> figure/mol3d/{name}.png")

print("완료 -> figure/mol3d/")
