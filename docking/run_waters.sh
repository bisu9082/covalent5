#!/bin/bash
# covalent5 · G3-1: 결정수(crystallographic waters) 포함 재도킹 (리뷰어 R1-3)
# 4EY7 원본의 활성부위 물을 보존한 wet 수용체로 14 앵커 재도킹 → ρ 비교(vs 물 제거 baseline).
# 요구: wget, obabel, gnina (covalent5 env). 실행:
#   cd $COV5/docking && conda activate covalent5
#   export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
#   bash run_waters.sh 2>&1 | tee waters_run.log
set -e
cd "$(dirname "$0")"
CX=-5.337; CY=-41.735; CZ=29.543; SZ=22.5   # Ser203 Ogamma 중심(baseline과 동일)

# 1) 4EY7 원본 다운로드 (없으면)
[ -f 4EY7.pdb ] || wget -q https://files.rcsb.org/download/4EY7.pdb -O 4EY7.pdb
echo "4EY7.pdb HOH: $(grep -c 'HOH' 4EY7.pdb)"

# 2) chain A 단백질 + Ser Ogamma 6A 내 결정수 보존 (좌표 프레임 동일)
python3 - <<PYEOF
import math
OG=($CX,$CY,$CZ); keep=[]; nwat=0
for ln in open("4EY7.pdb"):
    rec=ln[:6]
    if rec=="ATOM  " and ln[21]=="A":
        keep.append(ln)                       # chain A protein
    elif rec=="HETATM" and ln[17:20]=="HOH":
        try: x,y,z=float(ln[30:38]),float(ln[38:46]),float(ln[46:54])
        except: continue
        if math.dist((x,y,z),OG)<6.0:         # 활성부위 물만
            keep.append(ln); nwat+=1
open("receptor_wet_raw.pdb","w").write("".join(keep)+"END\n")
print(f"보존한 활성부위 물: {nwat}개 -> receptor_wet_raw.pdb")
PYEOF

# 3) 양성자화 + pdbqt (baseline과 동일 절차)
obabel receptor_wet_raw.pdb -O receptor_wet.pdb -p 7.4
obabel receptor_wet.pdb -O receptor_wet.pdbqt -xr

# 4) 14 앵커 재도킹 (앵커만: anchor_ki.csv 기준)
OUT=out_waters; mkdir -p "$OUT"
ANCHORS=$(tail -n +2 ../data/anchor_ki.csv | cut -d, -f1)
MAN=ligands/ligand_manifest.csv
for sid in $ANCHORS; do
  noncov=$(awk -F, -v s="$sid" '$1==s{print $5}' "$MAN")
  lig="ligands/$noncov"
  [ -f "$lig" ] || { echo "  [skip] $sid"; continue; }
  gnina -r receptor_wet.pdbqt -l "$lig" \
    --center_x $CX --center_y $CY --center_z $CZ \
    --size_x $SZ --size_y $SZ --size_z $SZ \
    --exhaustiveness 16 --num_modes 9 --cnn_scoring rescore --seed 42 \
    -o "$OUT/${sid}_wet.sdf.gz" --log "$OUT/${sid}_wet.log" >/dev/null 2>&1 \
    && echo "  ok $sid" || echo "  FAIL $sid"
done
echo "=== 완료. python parse_waters.py 로 ρ 비교 ==="
