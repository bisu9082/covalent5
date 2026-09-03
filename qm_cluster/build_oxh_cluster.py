#!/usr/bin/env python3
"""
covalent5 · His447 full-QM 확장 클러스터 빌더 (triad + oxyanion hole)
================================================================
중성 Michaelis 반응물(RC)을 oxyanion hole 포함으로 구축한다.
- triad 사이드체인: Ser203(-OH 중성), His447(HID 중성), Glu334(카복실레이트 -1)  [build_reactant_cluster2.py와 동일]
- oxyanion hole: Gly121 / Gly122 / Ala204 백본 아미드 N-H를 **formamide 단편**으로 모델
  (선행잔기 C=O + 공여 N-H, 기하학적 H캡 → 중성 amide, 결정 H결합 방향 보존)
- VX: vx_pose_H.sdf (수소 포함 포즈)
순전하 = -1 (Glu334만 음전하; formamide 전부 중성).

입력:
  rec_A_H.pdb   : reduce/obabel로 양성자화한 receptor_4EY7_chainA.pdb (표준번호, chain A 유지)
                  예) reduce receptor_4EY7_chainA.pdb > rec_A_H.pdb
                  ※ heavy 좌표 불변 → oxyanion 백본도 같은 프레임에서 추출됨(정합)
  vx_pose_H.sdf : H 추가된 VX 도킹 포즈 (기존 파일)
사용:
  python build_oxh_cluster.py rec_A_H.pdb vx_pose_H.sdf [chain=A]
출력:
  cluster_reactant_oxh.xyz   확장 중성 RC 클러스터
  rc_oxh.inp                 ORCA opt 입력 (caps freeze + P-Ogamma 3.3 + His Ne2-Ogamma 2.9)
  oxh_indices.json           ORCA 제약용 0-based 인덱스 (P, Ogamma, Ne2, donors, frozen)
"""
import sys, json, math

KEEP = {("SER", "203"): ["CB", "OG", "HB2", "HB3", "HG"],
        ("HIS", "447"): ["CB", "CG", "ND1", "CD2", "CE1", "NE2", "HB2", "HB3", "HD2", "HE1"],
        ("GLU", "334"): ["CB", "CG", "CD", "OE1", "OE2", "HB2", "HB3", "HG2", "HG3"]}
# oxyanion hole donor: (donor_res, donor_resi, preceding_res, preceding_resi)
OXH = [("GLY", "121", "GLY", "120"), ("GLY", "122", "GLY", "121"), ("ALA", "204", "SER", "203")]


def u(v):
    n = math.sqrt(sum(c * c for c in v)) or 1.0
    return [c / n for c in v]


def sub(a, b):
    return [a[i] - b[i] for i in range(3)]


def addp(a, d, s):
    return [a[i] + s * d[i] for i in range(3)]


def read_pdb(pdb, chain):
    A = {}
    for ln in open(pdb):
        if ln[:4] == "ATOM" and ln[21] == chain:
            nm = ln[12:16].strip(); res = ln[17:20].strip(); ri = ln[22:26].strip()
            el = (ln[76:78].strip() or nm[0])
            A[(res, ri, nm)] = (el, float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
    return A


def main():
    pdb, sdf = sys.argv[1], sys.argv[2]
    chain = sys.argv[3] if len(sys.argv) > 3 else "A"
    A = read_pdb(pdb, chain)
    out = []          # (element, x, y, z)
    cap_idx = []      # 좌표고정할 cap/anchor 인덱스
    meta = {}

    # ---- 1) triad 사이드체인 + Cbeta cap + His HID ----
    for (res, resi), keep in KEEP.items():
        for nm in keep:
            k = (res, resi, nm)
            if k in A:
                el, x, y, z = A[k]
                out.append((el, x, y, z))
        # Cbeta cap (CA 자리에 H) — 이 H를 고정해 골격 강성 모사
        if (res, resi, "CB") in A and (res, resi, "CA") in A:
            cb = list(A[(res, resi, "CB")][1:]); ca = list(A[(res, resi, "CA")][1:])
            d = u(sub(ca, cb))
            out.append(("H", *addp(cb, d, 1.09)))
            cap_idx.append(len(out) - 1)
    # His447 HID: ND1에 HD1 (고리 평면 외부)
    nd1 = list(A[("HIS", resi_of := "447", "ND1")][1:])
    cg = list(A[("HIS", "447", "CG")][1:]); ce1 = list(A[("HIS", "447", "CE1")][1:])
    mid = [(cg[i] + ce1[i]) / 2 for i in range(3)]
    out.append(("H", *addp(nd1, u(sub(nd1, mid)), 1.0)))
    # key triad indices (0-based, in build order)
    # Ser OG, His NE2 위치 재탐색
    meta["Ogamma"] = next(i for i, (el, x, y, z) in enumerate(out)
                          if el == "O" and abs(x - A[("SER", "203", "OG")][1]) < 1e-3)
    meta["Ne2"] = next(i for i, (el, x, y, z) in enumerate(out)
                       if el == "N" and abs(x - A[("HIS", "447", "NE2")][1]) < 1e-3)
    n_triad = len(out)

    # ---- 2) oxyanion hole formamide 단편 3개 ----
    donor_H = []
    for dr, di, pr, pi in OXH:
        N = list(A[(dr, di, "N")][1:]); CA = list(A[(dr, di, "CA")][1:])
        C = list(A[(pr, pi, "C")][1:]); O = list(A[(pr, pi, "O")][1:]); CAj = list(A[(pr, pi, "CA")][1:])
        # 결정 방향 아미드 H (외부 이등분선), 1.01
        dHn = u([-(u(sub(C, N))[k] + u(sub(CA, N))[k]) for k in range(3)])
        Hamide = addp(N, dHn, 1.01)
        Hcap_N = addp(N, u(sub(CA, N)), 1.01)       # N-CA(donor) 자리 cap
        Hformyl = addp(C, u(sub(CAj, C)), 1.09)     # C-CA(preceding) 자리 cap
        base = len(out)
        out += [("N", *N), ("H", *Hamide), ("H", *Hcap_N), ("C", *C), ("O", *O), ("H", *Hformyl)]
        # 고정: cap 2개(Hcap_N, Hformyl) + 카보닐 C,O (백본 강성)
        cap_idx += [base + 2, base + 3, base + 4, base + 5]
        donor_H.append(base + 1)                    # Hamide = oxyanion 공여
    meta["oxh_donor_H"] = donor_H

    # ---- 3) VX (H 포함 sdf) ----
    lines = open(sdf).read().splitlines()
    na = int(lines[3][:3])
    p_idx = None
    vx_start = len(out)
    for j, ln in enumerate(lines[4:4 + na]):
        s = ln.split()
        el = s[3]
        out.append((el, float(s[0]), float(s[1]), float(s[2])))
        if el == "P":
            p_idx = len(out) - 1
    meta["P"] = p_idx
    # 먼 VX 고정 (P에서 4 A 초과)
    P = out[p_idx][1:]
    far = [i for i in range(vx_start, len(out))
           if math.dist(out[i][1:], P) > 4.0]
    frozen = sorted(set(cap_idx + far))
    meta["frozen"] = frozen
    meta["net_charge"] = -1
    meta["n_atoms"] = len(out)

    # ---- write xyz ----
    with open("cluster_reactant_oxh.xyz", "w") as f:
        f.write(f"{len(out)}\nHis447 full-QM RC: triad(HID,Glu-,Ser-OH) + oxyanion(Gly121/122,Ala204 formamide) + VX. net -1\n")
        for el, x, y, z in out:
            f.write(f"{el:2s} {x:12.6f} {y:12.6f} {z:12.6f}\n")

    # ---- write ORCA opt input ----
    L = ["# covalent5 His447 full-QM RC opt (oxyanion hole). net -1.",
         "! wB97X-D3 def2-SVP def2/J RIJCOSX LooseOpt CPCM(water) TightSCF defgrid2",
         "%pal nprocs 8 end", "%maxcore 3500", "%geom Constraints"]
    for i in frozen:
        L.append(f"  {{C {i} C}}")
    L.append(f"  {{B {meta['P']} {meta['Ogamma']} 3.30 C}}")     # 친핵 접근 거리 약구속
    L.append(f"  {{B {meta['Ne2']} {meta['Ogamma']} 2.90 C}}")   # His Ne2...Ser Ogamma H결합
    L += [" end end", f"* xyzfile -1 1 cluster_reactant_oxh.xyz"]
    open("rc_oxh.inp", "w").write("\n".join(L) + "\n")

    json.dump(meta, open("oxh_indices.json", "w"), indent=2)
    print(f"원자 {len(out)} = triad부 {n_triad} + oxyanion 18 + VX {na}")
    print(f"P={meta['P']}  Ogamma={meta['Ogamma']}  His Ne2={meta['Ne2']}  oxyanion 공여H={donor_H}")
    print(f"frozen {len(frozen)}개 (caps+anchors+farVX)")
    print("-> cluster_reactant_oxh.xyz, rc_oxh.inp, oxh_indices.json")


if __name__ == "__main__":
    main()
