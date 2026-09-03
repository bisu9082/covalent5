#!/usr/bin/env python3
"""
covalent5 · G3-1 결정수 재도킹 결과 → ρ 비교 (vs 물 제거 baseline)
out_waters/*_wet.log 파싱(mode-1: affinity/intramol/CNNpose/CNNaffinity),
../data/anchor_ki.csv pKi 와 병합 → Vina·GNINA CNN Spearman ρ.
사용: cd $COV5/docking && python parse_waters.py
"""
import os, re, csv, glob, math

HERE = os.path.dirname(os.path.abspath(__file__))
ROW = re.compile(r"^\s*1\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)")


def parse_log(p):
    for ln in open(p):
        m = ROW.match(ln)
        if m:
            g = list(map(float, m.groups()))
            return g[0], g[3]           # Vina affinity, CNNaffinity
    return None, None


def rank(v):
    idx = sorted(range(len(v)), key=lambda i: v[i]); r = [0] * len(v); i = 0
    while i < len(v):
        j = i
        while j + 1 < len(v) and v[idx[j + 1]] == v[idx[i]]: j += 1
        for k in range(i, j + 1): r[idx[k]] = (i + j) / 2 + 1
        i = j + 1
    return r


def pear(a, b):
    n = len(a); ma = sum(a) / n; mb = sum(b) / n
    num = sum((a[i] - ma) * (b[i] - mb) for i in range(n))
    da = math.sqrt(sum((x - ma) ** 2 for x in a)); db = math.sqrt(sum((x - mb) ** 2 for x in b))
    return num / (da * db) if da * db else float("nan")


def spearman(a, b):
    return pear(rank(a), rank(b))


def main():
    ki = {r["stereo_id"]: float(r["pKi"]) for r in
          csv.DictReader(open(os.path.join(HERE, "..", "data", "anchor_ki.csv"))) if r.get("pKi")}
    vina, cnn, y, ids = [], [], [], []
    for sid in ki:
        log = os.path.join(HERE, "out_waters", f"{sid}_wet.log")
        if not os.path.exists(log): continue
        a, c = parse_log(log)
        if a is None: continue
        vina.append(a); cnn.append(c); y.append(ki[sid]); ids.append(sid)
    print(f"n={len(y)} (물 포함 재도킹)")
    print(f"  Vina  rho = {spearman(vina, y):+.2f}   (baseline +0.12)")
    print(f"  GNINA CNN rho = {spearman(cnn, y):+.2f}   (baseline -0.40)")
    print("=> 물 포함해도 ρ가 약하게 유지되면 '결론이 물 처리에 robust' 지지.")
    with open(os.path.join(HERE, "waters_summary.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["stereo_id", "pKi", "vina_wet", "cnn_wet"])
        for i in range(len(y)): w.writerow([ids[i], y[i], vina[i], cnn[i]])
    print("-> waters_summary.csv")


if __name__ == "__main__":
    main()
