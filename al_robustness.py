#!/usr/bin/env python3
"""
covalent5 · 물리정보 모델 통계 강건성 (리뷰어 대응 재현 스크립트)
al_rank.py와 동일한 14-피처(DFT 6 + RDKit 구조 8) GPR-LOO 파이프라인을 재사용하여:
  (1) GPR/RF LOO Spearman rho 재현
  (2) GPR rho 부트스트랩 95% CI  (pair-resample)
  (3) GPR rho 순열검정 p          (라벨셔플 후 LOO 재수행)
  (4) VX leave-BOTH-enantiomers-out (입체선택성 순환성 제거 테스트)
결과를 표준출력 + data/ml_robustness_out.md 로 기록.
사용:  cd "$COV5" && python al_robustness.py [--nperm 1000] [--nboot 10000] [--seed 0]
요구:  numpy, scipy, scikit-learn, rdkit
"""
import csv, os, argparse, warnings
import numpy as np
warnings.filterwarnings("ignore")
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from scipy.stats import spearmanr, pearsonr
from sklearn.preprocessing import StandardScaler
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, WhiteKernel, ConstantKernel
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import LeaveOneOut

HERE = os.path.dirname(os.path.abspath(__file__))
ELEC = ["HOMO_eV", "LUMO_eV", "gap_eV", "dipole_D", "qP_mulliken", "qP_loewdin"]
SF = ["MW", "TPSA", "nRot", "hasF", "hasPS", "hasCN", "hasSR", "hasArOxy"]


def rd(p):
    return list(csv.DictReader(open(os.path.join(HERE, p), encoding="utf-8")))


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


def make_kernel(d):
    return ConstantKernel() * RBF(np.ones(d)) + WhiteKernel()


def gpr_loo(X, y, restarts=2):
    """Leave-one-out GPR predictions."""
    pred = np.zeros(len(y)); sc = StandardScaler()
    k = make_kernel(X.shape[1])
    for tr, te in LeaveOneOut().split(X):
        g = GaussianProcessRegressor(kernel=k, normalize_y=True,
                                     n_restarts_optimizer=restarts).fit(sc.fit_transform(X[tr]), y[tr])
        pred[te] = g.predict(sc.transform(X[te]))
    return pred


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nperm", type=int, default=1000)
    ap.add_argument("--nboot", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)

    manifest = {r["stereo_id"]: r for r in rd("data/molecule_manifest.csv")}
    desc = {r["stereo_id"]: r for r in rd("data/descriptors.csv")}
    ki = {r["stereo_id"]: float(r["pKi"]) for r in rd("data/anchor_ki.csv") if r.get("pKi")}
    ids = [s for s in ki if s in desc and s in manifest]
    y = np.array([ki[s] for s in ids])
    Xe = np.array([[float(desc[s][c]) for c in ELEC] for s in ids])
    Xs = np.array([struct_feats(manifest[s]["canonical_smiles"]) for s in ids])
    X = np.hstack([Xe, Xs]); n = len(y)

    # (1) reproduce LOO rho
    gp = gpr_loo(X, y)
    sc = StandardScaler(); rf = np.zeros(n)
    for tr, te in LeaveOneOut().split(X):
        rf[te] = RandomForestRegressor(400, random_state=0).fit(
            sc.fit_transform(X[tr]), y[tr]).predict(sc.transform(X[te]))
    rho_gp = spearmanr(y, gp)[0]; rho_rf = spearmanr(y, rf)[0]
    r_gp = pearsonr(y, gp)[0]; r_rf = pearsonr(y, rf)[0]

    # (2) bootstrap CI (GPR rho, pair resample)
    bs = []
    for _ in range(args.nboot):
        idx = rng.integers(0, n, n)
        if len(set(y[idx])) < 3:
            continue
        r = spearmanr(y[idx], gp[idx])[0]
        if not np.isnan(r):
            bs.append(r)
    lo, hi = np.percentile(bs, [2.5, 97.5])

    # (3) permutation p (shuffle y, redo LOO GPR; restarts=0 for speed)
    null = np.empty(args.nperm)
    for i in range(args.nperm):
        null[i] = spearmanr(*(lambda yp: (yp, gpr_loo(X, yp, restarts=0)))(rng.permutation(y)))[0]
    pperm = (np.sum(null >= rho_gp) + 1) / (args.nperm + 1)

    # (4) leave-BOTH-VX-out
    vx = [i for i, s in enumerate(ids) if s.startswith("VX_")]
    mask = np.ones(n, bool); mask[vx] = False
    sc2 = StandardScaler(); Xt = sc2.fit_transform(X[mask])
    g = GaussianProcessRegressor(kernel=make_kernel(X.shape[1]), normalize_y=True,
                                 n_restarts_optimizer=2).fit(Xt, y[mask])
    pv, sv = g.predict(sc2.transform(X[vx]), return_std=True)
    vx_rows = [(ids[j], y[j], pv[i], sv[i]) for i, j in enumerate(vx)]

    lines = []
    lines.append(f"n = {n} anchors, {X.shape[1]} features ({len(ELEC)} DFT + {len(SF)} structural)")
    lines.append(f"(1) LOO Spearman rho:  GPR = {rho_gp:+.3f}   RF = {rho_rf:+.3f}")
    lines.append(f"    LOO Pearson  r  :  GPR = {r_gp:+.3f}   RF = {r_rf:+.3f}")
    lines.append(f"(2) GPR rho bootstrap 95% CI = [{lo:+.2f}, {hi:+.2f}]  "
                 f"(N={len(bs)}; excludes zero = {not (lo < 0 < hi)})")
    lines.append(f"(3) GPR rho permutation p = {pperm:.4f}  "
                 f"(N={args.nperm}; null mean {null.mean():+.3f}, 95th pct {np.percentile(null,95):+.3f})")
    lines.append("(4) leave-BOTH-VX-out (circularity control):")
    for sid, ya, pa, sa in vx_rows:
        lines.append(f"    {sid:14s} actual={ya:.2f}  pred={pa:.2f} +/- {sa:.2f}")
    sp = next(p for s, _, p, _ in vx_rows if s.endswith("P3S"))
    rp = next(p for s, _, p, _ in vx_rows if s.endswith("P3R"))
    lines.append(f"    -> predicts S_P > R_P out-of-sample: {sp > rp}  ({sp:.2f} vs {rp:.2f})")
    out = "\n".join(lines)
    print(out)
    with open(os.path.join(HERE, "data", "ml_robustness_out.md"), "w", encoding="utf-8") as f:
        f.write("# ml_robustness_out (auto-generated by al_robustness.py)\n\n```\n" + out + "\n```\n")
    print("\n-> data/ml_robustness_out.md")


if __name__ == "__main__":
    main()
