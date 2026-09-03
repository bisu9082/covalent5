#!/usr/bin/env python3
"""
covalent5 · 도킹 sensitivity 결과 → 변형별 Spearman rho (vs 실측 pKi)
out_sens_<label>/*_nc.log 를 파싱(mode-1: affinity/intramol/CNNpose/CNNaffinity),
../data/anchor_ki.csv 의 pKi 와 병합, 변형별 Vina affinity·GNINA CNN 의 Spearman rho 산출.
사용:  cd $COV5/docking && python parse_sensitivity.py
출력:  sensitivity_summary.csv  (label, n, rho_vina, rho_cnn)
"""
import os, re, csv, glob
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
ROW = re.compile(r"^\s*1\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)")


def parse_log(p):
    for ln in open(p):
        m = ROW.match(ln)
        if m:
            aff, intramol, cnnpose, cnnaff = map(float, m.groups())
            return aff, cnnaff
    return None, None


def load_ki():
    p = os.path.join(HERE, "..", "data", "anchor_ki.csv")
    return {r["stereo_id"]: float(r["pKi"]) for r in csv.DictReader(open(p)) if r.get("pKi")}


def main():
    ki = load_ki()
    rows = []
    for d in sorted(glob.glob(os.path.join(HERE, "out_sens_*"))):
        label = os.path.basename(d).replace("out_sens_", "")
        aff, cnn, y = [], [], []
        for log in glob.glob(os.path.join(d, "*_nc.log")):
            sid = os.path.basename(log).replace("_nc.log", "")
            if sid not in ki:
                continue
            a, c = parse_log(log)
            if a is None:
                continue
            aff.append(a); cnn.append(c); y.append(ki[sid])
        if len(y) >= 4:
            rv = spearmanr(aff, y)[0]; rc = spearmanr(cnn, y)[0]
            rows.append((label, len(y), rv, rc))
            print(f"{label:10s} n={len(y):2d}  rho_Vina={rv:+.2f}  rho_CNN={rc:+.2f}")
    out = os.path.join(HERE, "sensitivity_summary.csv")
    with open(out, "w", newline="") as f:
        w = csv.writer(f); w.writerow(["variant", "n", "rho_vina", "rho_cnn"]); w.writerows(rows)
    print(f"\n-> {os.path.relpath(out, HERE)}")
    if rows:
        rv = [r[2] for r in rows]; rc = [r[3] for r in rows]
        print(f"Vina rho range [{min(rv):+.2f}, {max(rv):+.2f}]  |  CNN rho range [{min(rc):+.2f}, {max(rc):+.2f}]")
        print("=> 모든 변형에서 |rho| 이 약하게 유지되면 '단일 셋업 특이적 아님' 결론 지지.")


if __name__ == "__main__":
    main()
