#!/usr/bin/env python3
"""
covalent5 · Step 4-③ 리간드 준비 (도킹용)
==================================================
molecule_manifest.csv 기반 두 세트 생성:
  1) noncovalent/<sid>.sdf : intact 작용제 (RQ1 비공유 도킹 — '잘못된 프록시' 축)
  2) covalent/<sid>.sdf    : leaving-group 제거 warhead. 이탈기 자리를 더미원자[*]로
                             치환하여 Ser203 Oγ 부착점 표시 (P–[*] = 형성될 공유결합)
이탈기 규칙:
  G/Novichok = P–F의 F, Tabun = C#N, V = P–S의 S(+thiocholine 가지),
  OP-oxon = 방향족/비닐 O(아릴옥시).  티온(P=S)은 직접 인산화 X → 공유 warhead 미생성.
또한 ligand_manifest.csv 출력(스크립트가 구버전 파일 무시하고 이것만 사용).
※ 공유결합 도킹 거동은 GNINA 버전 의존 → README의 6CQT/6CQX 검증 절차 필수.

사용: python prepare_ligands.py
"""
import os, csv
from rdkit import Chem
from rdkit.Chem import AllChem

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
MANIFEST = os.path.join(ROOT, "data", "molecule_manifest.csv")
LIGDIR = os.path.join(HERE, "ligands")
NONCOV = os.path.join(LIGDIR, "noncovalent")
COV = os.path.join(LIGDIR, "covalent")
LIGMAN = os.path.join(LIGDIR, "ligand_manifest.csv")


def find_leaving(mol, cls):
    P = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == "P"]
    if not P:
        return None, None, "no_P"
    P = P[0]
    nbrs = mol.GetAtomWithIdx(P).GetNeighbors()
    if cls in ("G", "Novichok"):
        for n in nbrs:
            if n.GetSymbol() == "F":
                return P, n.GetIdx(), "F"
        for n in nbrs:  # Tabun형 시안화물 이탈
            if n.GetSymbol() == "C" and any(
                    b.GetBondType() == Chem.BondType.TRIPLE for b in n.GetBonds()):
                return P, n.GetIdx(), "CN"
    if cls == "V":
        for n in nbrs:
            if n.GetSymbol() == "S":
                return P, n.GetIdx(), "S"
    if cls == "OP-oxon":
        for n in nbrs:
            if n.GetSymbol() == "O" and n.GetDegree() == 2:
                other = [a for a in n.GetNeighbors() if a.GetIdx() != P][0]
                if other.GetIsAromatic() or any(
                        b.GetBondType() == Chem.BondType.DOUBLE for b in other.GetBonds()):
                    return P, n.GetIdx(), "OAr/Ovinyl"
    return P, None, "leaving_not_found"


def make_warhead(smiles, cls):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None, "parse_fail"
    P, lv, tag = find_leaving(m, cls)
    if lv is None:
        return None, tag
    rw = Chem.RWMol(m)
    to_del = set(); stack = [lv]; seen = {P}
    while stack:
        a = stack.pop()
        if a in seen:
            continue
        seen.add(a); to_del.add(a)
        for nb in rw.GetAtomWithIdx(a).GetNeighbors():
            if nb.GetIdx() not in seen:
                stack.append(nb.GetIdx())
    dummy = rw.AddAtom(Chem.Atom(0))
    rw.AddBond(P, dummy, Chem.BondType.SINGLE)
    for idx in sorted(to_del, reverse=True):
        rw.RemoveAtom(idx)
    mol = rw.GetMol()
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return None, f"sanitize_fail:{e}"
    return mol, tag


def embed_write(mol, path, name):
    mol = Chem.AddHs(mol)
    p = AllChem.ETKDGv3(); p.randomSeed = 0xC0FFEE
    if AllChem.EmbedMolecule(mol, p) != 0:
        p.useRandomCoords = True
        AllChem.EmbedMolecule(mol, p)
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=1000)
    except Exception:
        pass
    mol.SetProp("_Name", name)
    w = Chem.SDWriter(path); w.write(mol); w.close()


def main():
    os.makedirs(NONCOV, exist_ok=True)
    os.makedirs(COV, exist_ok=True)
    rows = list(csv.DictReader(open(MANIFEST, encoding="utf-8")))
    man = []
    for r in rows:
        sid, cls, smi = r["stereo_id"], r["class"], r["canonical_smiles"]
        m = Chem.MolFromSmiles(smi)
        embed_write(m, os.path.join(NONCOV, sid + ".sdf"), sid)
        wh, tag = make_warhead(smi, cls)
        has_cov = 0
        if wh is not None:
            embed_write(wh, os.path.join(COV, sid + ".sdf"), sid + "_warhead")
            has_cov = 1
        man.append({"stereo_id": sid, "code": r["code"], "class": cls,
                    "role": r["role"], "noncov": f"noncovalent/{sid}.sdf",
                    "cov": f"covalent/{sid}.sdf" if has_cov else "",
                    "leaving": tag, "has_covalent": has_cov})
    with open(LIGMAN, "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=["stereo_id", "code", "class", "role",
                                           "noncov", "cov", "leaving", "has_covalent"])
        wr.writeheader(); wr.writerows(man)
    n_cov = sum(x["has_covalent"] for x in man)
    print(f"비공유 intact: {len(man)}  |  공유 warhead: {n_cov}")
    print(f"ligand_manifest: {os.path.relpath(LIGMAN, HERE)}")


if __name__ == "__main__":
    main()
