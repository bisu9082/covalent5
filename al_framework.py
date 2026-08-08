#!/usr/bin/env python3
"""
covalent5 · 물리정보 능동학습(GPR) 프레임워크
입력: data/features.csv (DFT+도킹) + data/anchor_ki.csv (실험 pKi, 라벨 부분집합)
기능: GPR(불확실도)+LOO-CV+RF baseline+미측정 예측(평균,sd)+acquisition+입체선택성.
요구: scikit-learn
사용: cd "$COV5" && python al_framework.py
anchor_ki.csv: stereo_id,ki,unit,pKi,source  (pKi 사용)
"""
import csv, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FEATS = ["HOMO_eV","LUMO_eV","gap_eV","dipole_D","qP_mulliken","qP_loewdin"]


def read_csv(p):
    with open(p, encoding="utf-8") as f:
        return list(csv.DictReader(f))


def main():
    feats = read_csv(os.path.join(HERE, "data", "features.csv"))
    ki_path = os.path.join(HERE, "data", "anchor_ki.csv")
    if not os.path.exists(ki_path):
        print("anchor_ki.csv 없음 — 문헌 ki 추출 후 생성 필요 (stereo_id,ki,unit,pKi,source).")
        print("프레임워크는 준비됨; 라벨 확보시 바로 실행.")
        return
    ki = {r["stereo_id"]: float(r["pKi"]) for r in read_csv(ki_path) if r.get("pKi")}

    X, y, ids, Xu, idu = [], [], [], [], []
    for r in feats:
        try:
            x = [float(r[c]) for c in FEATS]
        except ValueError:
            continue
        if r["stereo_id"] in ki:
            X.append(x); y.append(ki[r["stereo_id"]]); ids.append(r["stereo_id"])
        else:
            Xu.append(x); idu.append(r["stereo_id"])
    X, y = np.array(X), np.array(y)
    print(f"labeled {len(y)} / unlabeled {len(idu)}")
    if len(y) < 5:
        print("labeled <5 — need more anchor data.")
        return

    from sklearn.preprocessing import StandardScaler
    from sklearn.gaussian_process import GaussianProcessRegressor
    from sklearn.gaussian_process.kernels import RBF, WhiteKernel, ConstantKernel
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.model_selection import LeaveOneOut
    from sklearn.metrics import r2_score, mean_squared_error

    sc = StandardScaler().fit(X)
    Xs = sc.transform(X)
    kernel = ConstantKernel() * RBF(length_scale=np.ones(X.shape[1])) + WhiteKernel()

    loo = LeaveOneOut()
    gp_pred, rf_pred = [], []
    for tr, te in loo.split(Xs):
        g = GaussianProcessRegressor(kernel=kernel, normalize_y=True, n_restarts_optimizer=2).fit(Xs[tr], y[tr])
        gp_pred.append(g.predict(Xs[te])[0])
        rf = RandomForestRegressor(n_estimators=300, random_state=0).fit(Xs[tr], y[tr])
        rf_pred.append(rf.predict(Xs[te])[0])
    gp_pred, rf_pred = np.array(gp_pred), np.array(rf_pred)
    print(f"GPR LOO R2={r2_score(y,gp_pred):.3f} RMSE={mean_squared_error(y,gp_pred)**0.5:.3f}")
    print(f"RF  LOO R2={r2_score(y,rf_pred):.3f} RMSE={mean_squared_error(y,rf_pred)**0.5:.3f}")

    gp = GaussianProcessRegressor(kernel=kernel, normalize_y=True, n_restarts_optimizer=4).fit(Xs, y)
    if idu:
        Xus = sc.transform(np.array(Xu))
        mu, sd = gp.predict(Xus, return_std=True)
        out = [{"stereo_id": s, "pKi_pred": round(float(m), 3), "uncertainty_sd": round(float(u), 3)}
               for s, m, u in zip(idu, mu, sd)]
        with open(os.path.join(HERE, "data", "al_predictions.csv"), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["stereo_id", "pKi_pred", "uncertainty_sd"])
            w.writeheader(); w.writerows(out)
        nxt = sorted(out, key=lambda d: -d["uncertainty_sd"])[:3]
        print("acquisition(max uncertainty) next:", [d["stereo_id"] for d in nxt])
        print(f"-> data/al_predictions.csv ({len(out)})")


if __name__ == "__main__":
    main()
