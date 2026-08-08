#!/usr/bin/env python3
"""
covalent5 · ML 특성행렬 조립
data/descriptors.csv (DFT) + docking/docking_results.csv (비공유) -> data/features.csv
키: stereo_id. AL/GPR 입력. 타깃(pKi)은 anchor_ki.csv 확보 후 결합.
사용: cd "$COV5" && python build_features.py
"""
import csv, os
HERE = os.path.dirname(os.path.abspath(__file__))

def load(path, key="stereo_id"):
    with open(path, encoding="utf-8") as f:
        return {r[key]: r for r in csv.DictReader(f)}

desc = load(os.path.join(HERE, "data", "descriptors.csv"))
dock_rows = [r for r in csv.DictReader(open(os.path.join(HERE, "docking", "docking_results.csv"), encoding="utf-8"))
             if r["mode"] == "noncovalent"]
dock = {r["stereo_id"]: r for r in dock_rows}

cols = ["stereo_id", "code", "cls", "role",
        "HOMO_eV", "LUMO_eV", "gap_eV", "dipole_D", "qP_mulliken", "qP_loewdin",
        "dock_affinity", "dock_CNNpose", "dock_CNNaffinity"]
out = []
for sid, d in desc.items():
    k = dock.get(sid, {})
    out.append({
        "stereo_id": sid, "code": d["code"], "cls": d["cls"], "role": d["role"],
        "HOMO_eV": d["HOMO_eV"], "LUMO_eV": d["LUMO_eV"], "gap_eV": d["gap_eV"],
        "dipole_D": d["dipole_D"], "qP_mulliken": d["qP_mulliken"], "qP_loewdin": d["qP_loewdin"],
        "dock_affinity": k.get("affinity_kcal_mol", ""),
        "dock_CNNpose": k.get("CNNpose_score", ""),
        "dock_CNNaffinity": k.get("CNNaffinity", ""),
    })
with open(os.path.join(HERE, "data", "features.csv"), "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(out)
print(f"features.csv: {len(out)} rows x {len(cols)} cols")
miss = [o["stereo_id"] for o in out if not o["dock_affinity"]]
if miss:
    print("docking missing:", ", ".join(miss))
