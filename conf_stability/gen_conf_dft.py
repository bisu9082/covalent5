#!/usr/bin/env python3
"""
covalent5 · G3-2(b) 컨포머 민감도 — 유연 분자 4종 각 3 컨포머 DFT 입력 생성
리뷰어 R1-2 우려("단일 컨포머 DFT 디스크립터")에 답: 유연한 분자에서
컨포머를 바꿔도 HOMO/LUMO/gap/dipole/qP 가 안정함을 보인다.

대상(가장 회전결합 많은 4종): VX, VR, DETAB(diethyl-tabun), FEN(fenamiphos).
각 분자: RDKit ETKDGv3 다중 임베드 -> MMFF 최적화 -> 에너지순 -> RMSD로
구별되는 저에너지 컨포머 3개 -> ORCA opt(def2-SVP)+sp(def2-TZVP) 입력.

실행(WSL, covalent5 env; rdkit 필요):
  export COV5=/mnt/c/Users/kkaan/Downloads/claude_research/covalent5
  cd "$COV5/conf_stability" && python gen_conf_dft.py
그다음:  bash run_confs.sh   (ORCA_BIN 설정 필요)
"""
import os, csv, math
from rdkit import Chem
from rdkit.Chem import AllChem

HERE = os.path.dirname(os.path.abspath(__file__))
MAN = os.path.join(HERE, "..", "data", "molecule_manifest.csv")
TARGETS = ["VX_00_P3S", "VR_00_P8S", "DETAB_00_P3R", "FEN_00_P3R"]
NCONF = 3            # 분자당 유지할 구별 컨포머 수
NEMBED = 40         # 초기 임베드 후보 수
RMS_MIN = 0.5       # 구별 기준(Å)

OPT_HDR = "! wB97X-D3 def2-SVP def2/J RIJCOSX OPT CPCM(water) TightSCF defgrid2\n%pal nprocs 8 end\n%maxcore 3500\n"
SP_HDR = "! wB97X-D3 def2-TZVP def2/J RIJCOSX SP CPCM(water) TightSCF defgrid3\n%pal nprocs 8 end\n%maxcore 3500\n"


def load_smiles():
    m = {}
    for r in csv.DictReader(open(MAN)):
        if r["stereo_id"] in TARGETS:
            m[r["stereo_id"]] = r["canonical_smiles"]
    return m


def distinct_low_conformers(mol):
    """MMFF 최적화 후 에너지순으로 RMSD 구별되는 저에너지 컨포머 반환."""
    params = AllChem.ETKDGv3(); params.randomSeed = 0xC0FFEE; params.pruneRmsThresh = 0.3
    cids = list(AllChem.EmbedMultipleConfs(mol, numConfs=NEMBED, params=params))
    res = AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=2000)
    energies = [(cid, res[i][1]) for i, cid in enumerate(cids) if res[i][0] == 0]
    energies.sort(key=lambda x: x[1])
    kept = []
    for cid, e in energies:
        ok = True
        for kcid in kept:
            if AllChem.GetBestRMS(mol, mol, prbId=cid, refId=kcid) < RMS_MIN:
                ok = False; break
        if ok:
            kept.append(cid)
        if len(kept) >= NCONF:
            break
    return kept, {cid: e for cid, e in energies}


def xyz_block(mol, cid):
    conf = mol.GetConformer(cid); lines = []
    for at in mol.GetAtoms():
        p = conf.GetAtomPosition(at.GetIdx())
        lines.append(f"{at.GetSymbol():2s} {p.x:14.8f} {p.y:14.8f} {p.z:14.8f}")
    return "\n".join(lines)


def main():
    sm = load_smiles(); run_lines = ["#!/bin/bash", "set -e", 'ORCA="${ORCA_BIN:-orca}"']
    summary = []
    for sid in TARGETS:
        mol = Chem.AddHs(Chem.MolFromSmiles(sm[sid]))
        kept, emap = distinct_low_conformers(mol)
        print(f"{sid}: {len(kept)} 컨포머 (ΔE_MMFF max {max(emap[c] for c in kept)-min(emap[c] for c in kept):.2f} kcal/mol)")
        for k, cid in enumerate(kept):
            d = os.path.join(HERE, sid, f"conf{k}"); os.makedirs(d, exist_ok=True)
            body = f"* xyz 0 1\n{xyz_block(mol, cid)}\n*\n"
            open(os.path.join(d, "opt.inp"), "w").write(OPT_HDR + body)
            open(os.path.join(d, "sp.inp"), "w").write(SP_HDR + "* xyzfile 0 1 opt.xyz\n")
            rel = os.path.relpath(d, HERE)
            run_lines += [f'cd "{rel}"',
                          '"$ORCA" opt.inp > opt.out 2>&1',
                          '"$ORCA" sp.inp  > sp.out  2>&1',
                          f'echo "DONE {sid}/conf{k}"', 'cd "$OLDPWD" 2>/dev/null || cd ..; cd ..']
            summary.append((sid, k, emap[cid]))
    # 안전한 run 스크립트(절대경로 기반)
    run = ["#!/bin/bash", "set -e", 'ORCA="${ORCA_BIN:-orca}"', f'ROOT="{HERE}"']
    for sid, k, _ in summary:
        run += [f'cd "$ROOT/{sid}/conf{k}"',
                '"$ORCA" opt.inp > opt.out 2>&1',
                '"$ORCA" sp.inp  > sp.out  2>&1',
                f'echo "DONE {sid}/conf{k}"']
    run.append('echo "=== 전체 완료. python parse_conf_stability.py ==="')
    open(os.path.join(HERE, "run_confs.sh"), "w").write("\n".join(run) + "\n")
    with open(os.path.join(HERE, "conf_manifest.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["stereo_id", "conf", "E_mmff_kcal"])
        w.writerows(summary)
    print(f"\n총 {len(summary)} DFT 잡 -> run_confs.sh, conf_manifest.csv")
    print("실행: export ORCA_BIN=/home/k9/orca/orca; bash run_confs.sh")


if __name__ == "__main__":
    main()
