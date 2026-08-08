#!/usr/bin/env python3
"""
covalent5 · G3-2(b) 컨포머 안정성 판정
각 분자의 conf*/sp.out 을 extract_descriptors.py 와 동일 파서로 읽어
분자별 HOMO/LUMO/gap/dipole/qP 의 컨포머간 범위(range)·표준편차 산출.
범위가 디스크립터 척도 대비 작으면 '단일 컨포머 디스크립터가 robust' 지지.
사용: cd $COV5/conf_stability && python parse_conf_stability.py
"""
import os, re, csv, statistics as st

HERE = os.path.dirname(os.path.abspath(__file__))


def parse_sp(path):
    d = dict(HOMO_eV=None, LUMO_eV=None, gap_eV=None, dipole_D=None, qP=None)
    if not os.path.exists(path):
        return d
    txt = open(path, errors="ignore").read()
    mo = re.search(r"ORBITAL ENERGIES\s*\n.*?\n\s*NO\s+OCC.*?\n(.*?)\n\s*\n", txt, re.S)
    if mo:
        homo = lumo = None
        for ln in mo.group(1).splitlines():
            p = ln.split()
            if len(p) >= 4:
                try:
                    occ = float(p[1]); eev = float(p[3])
                except ValueError:
                    continue
                if occ > 0.5:
                    homo = eev
                elif lumo is None:
                    lumo = eev
        d["HOMO_eV"], d["LUMO_eV"] = homo, lumo
        if homo is not None and lumo is not None:
            d["gap_eV"] = round(lumo - homo, 4)
    md = re.findall(r"Magnitude \(Debye\)\s*:\s*(-?\d+\.\d+)", txt)
    if md:
        d["dipole_D"] = float(md[-1])
    ml = re.search(r"LOEWDIN ATOMIC CHARGES.*?\n(.*?)\n\s*\n", txt, re.S)
    if ml:
        pc = re.findall(r"^\s*\d+\s+P\s*:\s*(-?\d+\.\d+)", ml.group(1), re.M)
        if pc:
            d["qP"] = float(pc[0])
    return d


def main():
    mols = sorted([x for x in os.listdir(HERE)
                   if os.path.isdir(os.path.join(HERE, x)) and "_" in x])
    keys = ["HOMO_eV", "LUMO_eV", "gap_eV", "dipole_D", "qP"]
    rows = []
    print(f"{'molecule':16s} {'n':>2s}  " + "  ".join(f"{k}(range)" for k in keys))
    for m in mols:
        vals = {k: [] for k in keys}
        confs = sorted(x for x in os.listdir(os.path.join(HERE, m)) if x.startswith("conf"))
        for c in confs:
            d = parse_sp(os.path.join(HERE, m, c, "sp.out"))
            for k in keys:
                if d[k] is not None:
                    vals[k].append(d[k])
        n = max(len(vals[k]) for k in keys) if keys else 0
        if n < 2:
            print(f"{m:16s}  <2 컨포머 파싱됨 — DFT 완료 확인"); continue
        rng = {k: (max(vals[k]) - min(vals[k])) if len(vals[k]) >= 2 else float("nan") for k in keys}
        sd = {k: (st.pstdev(vals[k]) if len(vals[k]) >= 2 else float("nan")) for k in keys}
        print(f"{m:16s} {n:>2d}  " +
              "  ".join(f"{rng[k]:+.3f}" for k in keys))
        rows.append([m, n] + [round(rng[k], 4) for k in keys] + [round(sd[k], 4) for k in keys])
    out = os.path.join(HERE, "conf_stability_summary.csv")
    with open(out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["molecule", "n_conf"] + [f"{k}_range" for k in keys] + [f"{k}_sd" for k in keys])
        w.writerows(rows)
    print(f"\n-> {os.path.relpath(out, HERE)}")
    print("판정: gap range < ~0.2 eV, qP range < ~0.03 e 이면 컨포머-robust.")


if __name__ == "__main__":
    main()
