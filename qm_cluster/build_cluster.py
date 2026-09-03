#!/usr/bin/env python3
"""
covalent5 · QM 클러스터 추출 (최소 모델: 촉매 삼중쌍 + 부가물/기질)
PDB(체인 A)에서 Ser203/His447/Glu334 사이드체인(Cα-Cβ 절단·H캡) + 지정 기질 잔기 추출.
출력: cluster_heavy.xyz, frozen_atoms.txt(캡 원자 인덱스, 1-based)
※ 결정구조엔 H 없음 → 출력은 heavy+capH. 잔기/기질 H 추가는 다음(인터랙티브) 단계.
사용: python build_cluster.py <pdb> <substrate_resname> [chain]
예:   python build_cluster.py 6cqx.pdb VX A
"""
import sys, math

KEEP = {  # 사이드체인 유지 원자 (Cβ부터)
 ("SER","203"): ["CB","OG"],
 ("HIS","447"): ["CB","CG","ND1","CD2","CE1","NE2"],
 ("GLU","334"): ["CB","CG","CD","OE1","OE2"],
}

def atoms(pdb, chain):
    out=[]
    for ln in open(pdb):
        if ln[:6] in ("ATOM  ","HETATM") and ln[21]==chain:
            out.append(dict(rec=ln[:6].strip(), name=ln[12:16].strip(),
                res=ln[17:20].strip(), resi=ln[22:26].strip(),
                x=float(ln[30:38]), y=float(ln[38:46]), z=float(ln[46:54]), el=ln[76:78].strip()))
    return out

def cap_h(cb, ca):  # CB 위치에서 CA 방향 1.09Å H
    d=[ca[i]-cb[i] for i in range(3)]; n=math.sqrt(sum(v*v for v in d))
    return [cb[i]+1.09*d[i]/n for i in range(3)]

def main():
    pdb=sys.argv[1]; sub=sys.argv[2]; chain=sys.argv[3] if len(sys.argv)>3 else "A"
    A=atoms(pdb,chain)
    cluster=[]  # (el,x,y,z,tag)
    frozen=[]
    # 삼중쌍 사이드체인 + 캡
    for (res,resi),keep in KEEP.items():
        sub_at={a["name"]:a for a in A if a["res"]==res and a["resi"]==resi}
        if not sub_at:
            print(f"  [경고] {res}{resi} 없음"); continue
        for nm in keep:
            if nm in sub_at:
                a=sub_at[nm]; cluster.append((a["el"] or nm[0], a["x"],a["y"],a["z"], f"{res}{resi}:{nm}"))
        if "CB" in sub_at and "CA" in sub_at:
            h=cap_h([sub_at["CB"][k] for k in "xyz"],[sub_at["CA"][k] for k in "xyz"])
            cluster.append(("H",h[0],h[1],h[2], f"{res}{resi}:capH")); frozen.append(len(cluster))
    # 기질/부가물
    ns=sum(1 for a in A if a["res"]==sub)
    for a in A:
        if a["res"]==sub:
            cluster.append((a["el"] or a["name"][0], a["x"],a["y"],a["z"], f"{sub}:{a['name']}"))
    # 출력
    with open("cluster_heavy.xyz","w") as f:
        f.write(f"{len(cluster)}\ncovalent5 cluster: triad+{sub} (heavy+capH, no substrate H yet)\n")
        for el,x,y,z,tag in cluster:
            f.write(f"{el:2s} {x:12.6f} {y:12.6f} {z:12.6f}\n")
    open("frozen_atoms.txt","w").write(" ".join(map(str,frozen))+"\n")
    print(f"클러스터 원자: {len(cluster)} (삼중쌍 사이드체인+캡 + {sub} {ns}원자)")
    print(f"고정(캡H) 인덱스: {frozen}")
    print("→ cluster_heavy.xyz, frozen_atoms.txt  (※ 잔기/기질 H 추가는 다음 단계)")

main()
