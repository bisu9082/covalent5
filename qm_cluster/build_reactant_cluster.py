#!/usr/bin/env python3
"""
covalent5 · QM 반응물(Michaelis) 클러스터 조립 — VX 파일럿
수용체 촉매 삼중쌍 사이드체인(Cβ 절단·H캡) + 도킹 VX 포즈(intact, H포함) 결합.
출력: triad_heavy.xyz, cluster_reactant_raw.xyz
사용: python build_reactant_cluster.py <receptor.pdb> <vx_pose.sdf> [chain]
"""
import sys, math

KEEP = {("SER", "203"): ["CB", "OG"],
        ("HIS", "447"): ["CB", "CG", "ND1", "CD2", "CE1", "NE2"],
        ("GLU", "334"): ["CB", "CG", "CD", "OE1", "OE2"]}


def receptor_atoms(pdb, chain):
    A = {}
    for ln in open(pdb):
        if ln[:4] == "ATOM" and ln[21] == chain:
            A[(ln[17:20].strip(), ln[22:26].strip(), ln[12:16].strip())] = (
                (ln[76:78].strip() or ln[12:16].strip()[0]),
                float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
    return A


def cap_h(cb, ca):
    d = [ca[i] - cb[i] for i in range(3)]
    n = math.sqrt(sum(v * v for v in d))
    return [cb[i] + 1.09 * d[i] / n for i in range(3)]


def read_sdf(sdf):
    out = []
    lines = open(sdf).read().splitlines()
    na = int(lines[3][:3])
    for ln in lines[4:4 + na]:
        p = ln.split()
        out.append((p[3], float(p[0]), float(p[1]), float(p[2])))
    return out


def main():
    pdb = sys.argv[1]
    sdf = sys.argv[2]
    chain = sys.argv[3] if len(sys.argv) > 3 else "A"
    A = receptor_atoms(pdb, chain)
    triad = []
    for (res, resi), keep in KEEP.items():
        for nm in keep:
            if (res, resi, nm) in A:
                el, x, y, z = A[(res, resi, nm)]
                triad.append((el, x, y, z))
        if (res, resi, "CB") in A and (res, resi, "CA") in A:
            cb = A[(res, resi, "CB")][1:]
            ca = A[(res, resi, "CA")][1:]
            h = cap_h(list(cb), list(ca))
            triad.append(("H", h[0], h[1], h[2]))
    with open("triad_heavy.xyz", "w") as f:
        f.write(f"{len(triad)}\ntriad sidechains (heavy+capH)\n")
        for el, x, y, z in triad:
            f.write(f"{el:2s} {x:12.6f} {y:12.6f} {z:12.6f}\n")
    lig = read_sdf(sdf)
    allc = triad + lig
    with open("cluster_reactant_raw.xyz", "w") as f:
        f.write(f"{len(allc)}\ntriad(heavy+cap)+VX(heavy+H), no protein H\n")
        for el, x, y, z in allc:
            f.write(f"{el:2s} {x:12.6f} {y:12.6f} {z:12.6f}\n")
    print(f"triad {len(triad)} (heavy+cap) + VX {len(lig)} = {len(allc)} atoms")
    print("-> triad_heavy.xyz, cluster_reactant_raw.xyz")


if __name__ == "__main__":
    main()
