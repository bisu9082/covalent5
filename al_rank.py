#!/usr/bin/env python3
"""
covalent5 · ML 랭킹 재정의 + 구조특성 보강
전자(DFT) + RDKit 구조특성(이탈기/MW/TPSA) 결합 -> logki 랭킹(Spearman, LOO).
사용: cd "$COV5" && python al_rank.py
"""
import csv, os
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from scipy.stats import spearmanr, pearsonr
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestRegressor
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, WhiteKernel, ConstantKernel
from sklearn.model_selection import LeaveOneOut

HERE = os.path.dirname(os.path.abspath(__file__))
def rd(p): return list(csv.DictReader(open(p, encoding="utf-8")))

manifest = {r["stereo_id"]: r for r in rd(os.path.join(HERE, "data", "molecule_manifest.csv"))}
desc = {r["stereo_id"]: r for r in rd(os.path.join(HERE, "data", "descriptors.csv"))}
ki = {r["stereo_id"]: float(r["pKi"]) for r in rd(os.path.join(HERE, "data", "anchor_ki.csv")) if r.get("pKi")}

ELEC = ["HOMO_eV", "LUMO_eV", "gap_eV", "dipole_D", "qP_mulliken", "qP_loewdin"]
SF = ["MW", "TPSA", "nRot", "hasF", "hasPS", "hasCN", "hasSR", "hasArOxy"]

def struct_feats(smi):
    m = Chem.MolFromSmiles(smi)
    hasF = 1.0 if any(a.GetSymbol() == "F" for a in m.GetAtoms()) else 0.0
    hasPS = 1.0 if any(b.GetBondType() == Chem.BondType.DOUBLE and
                       {b.GetBeginAtom().GetSymbol(), b.GetEndAtom().GetSymbol()} == {"P", "S"}
                       for b in m.GetBonds()) else 0.0
    hasCN = 1.0 if m.HasSubstructMatch(Chem.MolFromSmarts("P-C#N")) else 0.0
    hasSR = 1.0 if m.HasSubstructMatch(Chem.MolFromSmarts("P-S-C")) else 0.0
    hasAr = 1.0 if m.HasSubstructMatch(Chem.MolFromSmarts("P-O-c")) else 0.0
    return [Descriptors.MolWt(m), Descriptors.TPSA(m),
            rdMolDescriptors.CalcNumRotatableBonds(m), hasF, hasPS, hasCN, hasSR, hasAr]

ids = [s for s in ki if s in desc and s in manifest]
y = np.array([ki[s] for s in ids])
Xe = np.array([[float(desc[s][c]) for c in ELEC] for s in ids])
Xs = np.array([struct_feats(manifest[s]["canonical_smiles"]) for s in ids])
X = np.hstack([Xe, Xs]); names = ELEC + SF
print(f"n={len(y)}  features {X.shape[1]} (elec{len(ELEC)}+struct{len(SF)})")

print("\nsingle feature vs logki:")
for i, nm in enumerate(names):
    print(f"  {nm:10} r={pearsonr(X[:, i], y)[0]:+.2f} rho={spearmanr(X[:, i], y)[0]:+.2f}")

def loo(model_fn):
    pred = np.zeros(len(y)); sc = StandardScaler()
    for tr, te in LeaveOneOut().split(X):
        Xt = sc.fit_transform(X[tr]); model = model_fn().fit(Xt, y[tr])
        pred[te] = model.predict(sc.transform(X[te]))
    return pred

rf = loo(lambda: RandomForestRegressor(n_estimators=400, random_state=0))
k = ConstantKernel() * RBF(np.ones(X.shape[1])) + WhiteKernel()
gp = loo(lambda: GaussianProcessRegressor(kernel=k, normalize_y=True, n_restarts_optimizer=2))
print(f"\nLOO ranking Spearman rho: RF={spearmanr(y, rf)[0]:+.2f}  GPR={spearmanr(y, gp)[0]:+.2f}")
print(f"LOO Pearson r:            RF={pearsonr(y, rf)[0]:+.2f}  GPR={pearsonr(y, gp)[0]:+.2f}")
