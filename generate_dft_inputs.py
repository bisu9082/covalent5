#!/usr/bin/env python3
"""
covalent5 · Step 4-② ORCA DFT 입력 생성
==================================================
molecule_manifest.csv의 각 입체분리 구조에 대해 ORCA 입력 2종 생성:
  1) opt.inp  : 지오메트리 최적화 + 진동수  (wB97X-D3 / def2-SVP, CPCM water)
  2) sp.inp   : 단일점 에너지 (opt 결과 좌표 사용)  (wB97X-D3 / def2-TZVP, CPCM water)
좌표는 인라인 삽입(self-contained). 모든 분자 중성·닫힌껍질(0,1).

※ ORCA 병렬(nprocs>1)은 풀 경로 호출 필요 → run.sh가 ${ORCA_BIN:-orca} 사용.
  실행 전: export ORCA_BIN=/home/k9/orca/orca

사용:
  python generate_dft_inputs.py --nprocs 8 --maxcore 3500
"""
import os, csv, argparse, stat

HERE = os.path.dirname(os.path.abspath(__file__))
MANIFEST = os.path.join(HERE, "data", "molecule_manifest.csv")
STRUCT = os.path.join(HERE, "data", "structures")
DFT = os.path.join(HERE, "dft")

OPT_KEYS = "! wB97X-D3 def2-SVP def2/J RIJCOSX OPT FREQ CPCM(water) TightSCF defgrid2"
SP_KEYS  = "! wB97X-D3 def2-TZVP def2/J RIJCOSX SP CPCM(water) TightSCF defgrid3"


def read_xyz_body(path):
    with open(path) as f:
        lines = f.read().splitlines()
    n = int(lines[0].strip())
    return "\n".join(lines[2:2 + n])


def orca_input(keywords, nprocs, maxcore, coord_block=None, xyzfile=None,
               charge=0, mult=1):
    head = [keywords,
            f"%pal nprocs {nprocs} end",
            f"%maxcore {maxcore}"]
    if coord_block is not None:
        head.append(f"* xyz {charge} {mult}")
        head.append(coord_block)
        head.append("*")
    else:
        head.append(f"* xyzfile {charge} {mult} {xyzfile}")
    return "\n".join(head) + "\n"


def run_sh(stereo_id):
    return f"""#!/bin/bash
# {stereo_id} — opt+freq 후 single-point. ORCA는 병렬 시 풀 경로 필요.
# 실행 전: export ORCA_BIN=/home/k9/orca/orca
set -e
ORCA="${{ORCA_BIN:-orca}}"
"$ORCA" opt.inp > opt.out 2>&1
"$ORCA" sp.inp  > sp.out  2>&1
echo "DONE {stereo_id}"
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nprocs", type=int, default=8)
    ap.add_argument("--maxcore", type=int, default=3500, help="MB per core")
    args = ap.parse_args()

    os.makedirs(DFT, exist_ok=True)
    rows = list(csv.DictReader(open(MANIFEST, encoding="utf-8")))
    jobs = []
    for r in rows:
        code, sid = r["code"], r["stereo_id"]
        q, mult = int(r["formal_charge"]), 1
        xyz_path = os.path.join(HERE, r["xyz"])
        if not os.path.exists(xyz_path):
            print(f"[MISS] {sid}: xyz 없음 {r['xyz']}")
            continue
        body = read_xyz_body(xyz_path)
        jdir = os.path.join(DFT, code, sid)
        os.makedirs(jdir, exist_ok=True)
        with open(os.path.join(jdir, "opt.inp"), "w") as f:
            f.write(orca_input(OPT_KEYS, args.nprocs, args.maxcore,
                               coord_block=body, charge=q, mult=mult))
        with open(os.path.join(jdir, "sp.inp"), "w") as f:
            f.write(orca_input(SP_KEYS, args.nprocs, args.maxcore,
                               xyzfile="opt.xyz", charge=q, mult=mult))
        rsh = os.path.join(jdir, "run.sh")
        with open(rsh, "w") as f:
            f.write(run_sh(sid))
        os.chmod(rsh, os.stat(rsh).st_mode | stat.S_IEXEC)
        jobs.append({"code": code, "stereo_id": sid, "class": r["class"],
                     "role": r["role"], "charge": q, "mult": mult,
                     "n_atoms": r["n_atoms"],
                     "jobdir": os.path.relpath(jdir, HERE)})

    with open(os.path.join(DFT, "run_all.sh"), "w") as f:
        f.write("#!/bin/bash\n# covalent5 전체 DFT 잡 순차 실행\n")
        f.write("# 실행 전: export ORCA_BIN=/home/k9/orca/orca\nset -e\n")
        f.write('ROOT="$(cd "$(dirname "$0")" && pwd)"\n')
        for j in jobs:
            rel = os.path.relpath(os.path.join(HERE, j["jobdir"]), DFT)
            f.write(f'echo ">>> {j["stereo_id"]}"; (cd "$ROOT/{rel}" && bash run.sh)\n')
    os.chmod(os.path.join(DFT, "run_all.sh"),
             os.stat(os.path.join(DFT, "run_all.sh")).st_mode | stat.S_IEXEC)

    with open(os.path.join(DFT, "job_manifest.csv"), "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=["code", "stereo_id", "class", "role",
                                           "charge", "mult", "n_atoms", "jobdir"])
        wr.writeheader(); wr.writerows(jobs)

    print(f"=== DFT 입력 생성 완료 ===")
    print(f"잡 수: {len(jobs)}  (각 opt+freq & single-point)")
    print(f"러너: dft/run_all.sh (export ORCA_BIN 먼저) | 매니페스트: dft/job_manifest.csv")


if __name__ == "__main__":
    main()
