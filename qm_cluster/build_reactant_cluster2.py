#!/usr/bin/env python3
"""
covalent5 · QM 반응물 클러스터 v2 — 양성자화된 수용체 사용 + His447 HID 수동 부여
입력: reduce 처리된 receptor_H.pdb + H 추가된 VX 포즈(vx_pose_H.sdf)
삼중쌍 사이드체인(+reduce H) , Cβ 절단·H캡, His447 ND1에 HD1 추가(HID), Glu 탈양성자.
순전하 = -1 (Glu334 카복실레이트).
사용: python build_reactant_cluster2.py receptor_H.pdb vx_pose_H.sdf [chain]
출력: cluster_reactant.xyz
"""
import sys, math

# 사이드체인 유지 원자 (backbone N/CA/C/O/H/HA 제외)
KEEP = {("SER", "203"): ["CB", "OG", "HB2", "HB3", "HG"],
        ("HIS", "447"): ["CB", "CG", "ND1", "CD2", "CE1", "NE2", "HB2", "HB3", "HD2", "HE1"],
        ("GLU", "334"): ["CB", "CG", "CD", "OE1", "OE2", "HB2", "HB3", "HG2", "HG3"]}


def atoms(pdb, chain):
    A = {}
    for ln in open(pdb):
        if ln[:4] == "ATOM" and ln[21] == chain:
            A[(ln[17:20].strip(), ln[22:26].strip(), ln[12:16].strip())] = (
                (ln[76:78].strip() or ln[12:16].strip()[0]),
                float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
    return A


def unit(v):
    n = math.sqrt(sum(c * c for c in v))
    return [c / n for c in v]


def main():
    pdb, sdf = sys.argv[1], sys.argv[2]
    chain = sys.argv[3] if len(sys.argv) > 3 else "A"
    A = atoms(pdb, chain)
    out = []
    for (res, resi), keep in KEEP.items():
        for nm in keep:
            k = (res, resi, nm)
            if k in A:
                el, x, y, z = A[k]
                out.append((el, x, y, z))
        # Cβ 캡 (CA 자리에 H)
        if (res, resi, "CB") in A and (res, resi, "CA") in A:
            cb = list(A[(res, resi, "CB")][1:])
            ca = list(A[(res, resi, "CA")][1:])
            d = unit([ca[i] - cb[i] for i in range(3)])
            out.append(("H", *[cb[i] + 1.09 * d[i] for i in range(3)]))
    # His447 HID: ND1에 H 추가 (고리 외부 방향, 평면내 ≈ ND1에서 CG·CE1 중점 반대쪽)
    nd1 = A[("HIS", "447", "ND1")][1:]
    cg = A[("HIS", "447", "CG")][1:]
    ce1 = A[("HIS", "447", "CE1")][1:]
    mid = [(cg[i] + ce1[i]) / 2 for i in range(3)]
    d = unit([nd1[i] - mid[i] for i in range(3)])
    out.append(("H", *[nd1[i] + 1.0 * d[i] for i in range(3)]))  # HD1
    n_triad = len(out)
    # VX (H 포함 sdf)
    lines = open(sdf).read().splitlines()
    na = int(lines[3][:3])
    for ln in lines[4:4 + na]:
        p = ln.split()
        out.append((p[3], float(p[0]), float(p[1]), float(p[2])))
    with open("cluster_reactant.xyz", "w") as f:
        f.write(f"{len(out)}\nVX reactant cluster: triad(HID,Glu-,Ser-OH,capped) + VX. net charge -1\n")
        for el, x, y, z in out:
            f.write(f"{el:2s} {x:12.6f} {y:12.6f} {z:12.6f}\n")
    print(f"삼중쌍부 {n_triad}원자 (HID HD1 포함) + VX {na} = {len(out)}원자")
    print("순전하 -1 (Glu334 카복실레이트). -> cluster_reactant.xyz")


if __name__ == "__main__":
    main()
