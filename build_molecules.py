#!/usr/bin/env python3
"""
covalent5 · Step 4-① 분자 집합 구축
==================================================
compound_set.csv를 읽어 (1) SMILES 검증, (2) 입체이성질체 열거(R/S + C=N E/Z 라벨),
(3) 3D 좌표 생성(ETKDGv3 + MMFF94 최적화), (4) XYZ/SDF 출력, (5) manifest 작성.

데이터 무결성:
 - verify_flag=1 행은 기본 SKIP. 확정 후 --include-provisional 로만 포함.
 - SMILES 공란 행은 자동 SKIP(구조 미확정).
 - 모든 출처(source 열)를 manifest에 보존.

사용:
  python build_molecules.py
  python build_molecules.py --include-provisional
  python build_molecules.py --max-stereo 8
출력: data/structures/<code>/<stereo_id>.{xyz,sdf}, data/molecule_manifest.csv
"""
import os, csv, argparse
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions

HERE = os.path.dirname(os.path.abspath(__file__))
CSV_IN = os.path.join(HERE, "data", "compound_set.csv")
OUT_DIR = os.path.join(HERE, "data", "structures")
MANIFEST = os.path.join(HERE, "data", "molecule_manifest.csv")


def assign_label(mol):
    """입체중심 R/S + 이중결합 E/Z 라벨 (C=N 기하 포함)."""
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    labels = []
    for a in mol.GetAtoms():
        if a.HasProp("_CIPCode"):
            labels.append(f"{a.GetSymbol()}{a.GetIdx()}{a.GetProp('_CIPCode')}")
    ez = {Chem.BondStereo.STEREOE: "E", Chem.BondStereo.STEREOZ: "Z",
          Chem.BondStereo.STEREOCIS: "Z", Chem.BondStereo.STEREOTRANS: "E"}
    for b in mol.GetBonds():
        st = b.GetStereo()
        if st in ez:
            i, j = b.GetBeginAtom(), b.GetEndAtom()
            labels.append(f"{i.GetSymbol()}{i.GetIdx()}{j.GetSymbol()}{j.GetIdx()}{ez[st]}")
    return "_".join(labels) if labels else "achiral"


def embed3d(mol, seed=0xC0FFEE):
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    if AllChem.EmbedMolecule(mol, params) != 0:
        params.useRandomCoords = True
        if AllChem.EmbedMolecule(mol, params) != 0:
            return None, "embed_failed"
    try:
        if AllChem.MMFFHasAllMoleculeParams(mol):
            AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)
            engine = "MMFF94"
        else:
            AllChem.UFFOptimizeMolecule(mol, maxIters=2000)
            engine = "UFF(fallback)"
    except Exception as e:
        return None, f"opt_failed:{e}"
    return mol, engine


def write_xyz(mol, path, comment=""):
    conf = mol.GetConformer()
    with open(path, "w") as f:
        f.write(f"{mol.GetNumAtoms()}\n{comment}\n")
        for atom in mol.GetAtoms():
            p = conf.GetAtomPosition(atom.GetIdx())
            f.write(f"{atom.GetSymbol():2s} {p.x:14.8f} {p.y:14.8f} {p.z:14.8f}\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--include-provisional", action="store_true")
    ap.add_argument("--max-stereo", type=int, default=8)
    args = ap.parse_args()

    os.makedirs(OUT_DIR, exist_ok=True)
    rows = list(csv.DictReader(open(CSV_IN, encoding="utf-8")))
    manifest = []
    n_ok = n_skip = n_fail = 0

    for r in rows:
        name, code, smi = r["name"], r["code"], (r["smiles"] or "").strip()
        flag = r["verify_flag"].strip() == "1"
        if not smi:
            print(f"[SKIP] {name:16} 구조 미확정 (SMILES 공란) — {r['source']}")
            n_skip += 1; continue
        if flag and not args.include_provisional:
            print(f"[SKIP] {name:16} verify_flag=1 (검증 필요) — {r['source']}")
            n_skip += 1; continue
        m = Chem.MolFromSmiles(smi)
        if m is None:
            print(f"[FAIL] {name:16} SMILES 파싱 실패: {smi}")
            n_fail += 1; continue

        cdir = os.path.join(OUT_DIR, code)
        os.makedirs(cdir, exist_ok=True)
        opts = StereoEnumerationOptions(maxIsomers=args.max_stereo, onlyUnassigned=True)
        isomers = list(EnumerateStereoisomers(m, options=opts))
        print(f"[OK]   {name:16} class={r['class']:9} stereoisomers={len(isomers)}"
              f"{'  *PROVISIONAL*' if flag else ''}")

        for i, iso in enumerate(isomers):
            iso = Chem.MolFromSmiles(Chem.MolToSmiles(iso))
            label = assign_label(iso)
            mol3d, engine = embed3d(iso)
            if mol3d is None:
                print(f"        isomer {i} 3D 실패: {engine}")
                n_fail += 1; continue
            charge = Chem.GetFormalCharge(mol3d)
            sid = f"{code}_{i:02d}_{label}"
            xyz = os.path.join(cdir, sid + ".xyz")
            sdf = os.path.join(cdir, sid + ".sdf")
            write_xyz(mol3d, xyz, comment=f"{name} {label} q={charge}")
            w = Chem.SDWriter(sdf); w.write(mol3d); w.close()
            manifest.append({
                "name": name, "code": code, "class": r["class"], "role": r["role"],
                "stereo_id": sid, "cip_label": label,
                "canonical_smiles": Chem.MolToSmiles(iso),
                "n_atoms": mol3d.GetNumAtoms(), "formal_charge": charge,
                "ff_engine": engine, "source": r["source"],
                "provisional": int(flag), "xyz": os.path.relpath(xyz, HERE),
            })
            n_ok += 1

    default_fields = ["name", "code", "class", "role", "stereo_id", "cip_label",
                      "canonical_smiles", "n_atoms", "formal_charge", "ff_engine",
                      "source", "provisional", "xyz"]
    fields = list(manifest[0].keys()) if manifest else default_fields
    with open(MANIFEST, "w", newline="", encoding="utf-8") as f:
        wr = csv.DictWriter(f, fieldnames=fields)
        wr.writeheader()
        wr.writerows(manifest)

    print(f"\n=== 요약 === 생성 구조 {n_ok} | 화합물 SKIP {n_skip} | 실패 {n_fail}")
    print(f"manifest: {os.path.relpath(MANIFEST, HERE)}  ({len(manifest)} 구조)")


if __name__ == "__main__":
    main()
