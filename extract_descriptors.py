#!/usr/bin/env python3
"""
covalent5 · DFT 디스크립터 추출
dft/job_manifest.csv의 각 잡 sp.out에서 전자구조 디스크립터 추출 -> data/descriptors.csv
추출: FINAL SP energy, HOMO/LUMO/gap(eV), dipole(Debye), P의 Mulliken/Loewdin 전하.
사용 (30개 DFT 완료 후):
  export COV5="/mnt/c/Users/kkaan/Downloads/claude_research/covalent5"
  cd "$COV5" && python extract_descriptors.py
"""
import os, re, csv

HERE = os.path.dirname(os.path.abspath(__file__))
JOBMAN = os.path.join(HERE, "dft", "job_manifest.csv")
OUT = os.path.join(HERE, "data", "descriptors.csv")


def parse_sp(path):
    d = dict(E_sp_Eh=None, HOMO_eV=None, LUMO_eV=None, gap_eV=None,
             dipole_D=None, qP_mulliken=None, qP_loewdin=None)
    if not os.path.exists(path):
        return d
    txt = open(path, errors="ignore").read()

    m = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", txt)
    if m:
        d["E_sp_Eh"] = float(m[-1])

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

    mm = re.search(r"MULLIKEN ATOMIC CHARGES.*?\n(.*?)\n\s*Sum of atomic charges", txt, re.S)
    if mm:
        pc = re.findall(r"^\s*\d+\s+P\s*:\s*(-?\d+\.\d+)", mm.group(1), re.M)
        if pc:
            d["qP_mulliken"] = float(pc[0])
    ml = re.search(r"LOEWDIN ATOMIC CHARGES.*?\n(.*?)\n\s*\n", txt, re.S)
    if ml:
        pc = re.findall(r"^\s*\d+\s+P\s*:\s*(-?\d+\.\d+)", ml.group(1), re.M)
        if pc:
            d["qP_loewdin"] = float(pc[0])
    return d


def main():
    rows = list(csv.DictReader(open(JOBMAN, encoding="utf-8")))
    out, miss = [], []
    for r in rows:
        sp = os.path.join(HERE, r["jobdir"], "sp.out")
        d = parse_sp(sp)
        if d["E_sp_Eh"] is None:
            miss.append(r["stereo_id"])
        out.append(dict(stereo_id=r["stereo_id"], code=r["code"], cls=r["class"],
                        role=r["role"], **d))
    fields = ["stereo_id", "code", "cls", "role", "E_sp_Eh", "HOMO_eV", "LUMO_eV",
              "gap_eV", "dipole_D", "qP_mulliken", "qP_loewdin"]
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=fields); wr.writeheader(); wr.writerows(out)
    ok = sum(1 for o in out if o["E_sp_Eh"] is not None)
    print(f"추출 {len(out)}건 (에너지 성공 {ok}) -> {os.path.relpath(OUT, HERE)}")
    if miss:
        print("sp.out 미완/누락:", ", ".join(miss))


if __name__ == "__main__":
    main()
