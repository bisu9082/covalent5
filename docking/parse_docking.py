#!/usr/bin/env python3
"""
covalent5 · GNINA 도킹 결과 파싱 (v1.3.2 포맷)
GNINA 결과 테이블(4컬럼: affinity / intramol / CNN pose score / CNN affinity)에서
best pose(mode 1) 추출 -> docking_results.csv. intramol 음수 허용.
사용: python parse_docking.py
"""
import os, re, csv, glob

HERE = os.path.dirname(os.path.abspath(__file__))
# mode-1 행: "    1   -5.78   -0.03   0.9426   4.250"  (음수 허용, 4개 float)
ROW = re.compile(r"^\s*1\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)")


def parse_log(path):
    with open(path) as f:
        for ln in f:
            m = ROW.match(ln)
            if m:
                aff, intramol, cnnpose, cnnaff = map(float, m.groups())
                return aff, intramol, cnnpose, cnnaff
    return None, None, None, None


def collect(outdir, tag):
    rows = []
    for log in sorted(glob.glob(os.path.join(HERE, outdir, "*.log"))):
        sid = re.sub(r"_(nc|cov)$", "", os.path.basename(log).rsplit(".", 1)[0])
        aff, intramol, cnnpose, cnnaff = parse_log(log)
        rows.append({"stereo_id": sid, "mode": tag,
                     "affinity_kcal_mol": aff, "intramol_kcal_mol": intramol,
                     "CNNpose_score": cnnpose, "CNNaffinity": cnnaff})
    return rows


def main():
    rows = collect("out_noncovalent", "noncovalent") + collect("out_covalent", "covalent")
    out = os.path.join(HERE, "docking_results.csv")
    with open(out, "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=["stereo_id", "mode", "affinity_kcal_mol",
                                           "intramol_kcal_mol", "CNNpose_score", "CNNaffinity"])
        wr.writeheader(); wr.writerows(rows)
    ok = sum(1 for r in rows if r["affinity_kcal_mol"] is not None)
    print(f"파싱 {len(rows)}건 (점수 추출 성공 {ok}건) -> {os.path.relpath(out, HERE)}")
    miss = [r["stereo_id"] for r in rows if r["affinity_kcal_mol"] is None]
    if miss:
        print("점수 없음:", ", ".join(miss))


if __name__ == "__main__":
    main()
