#!/bin/bash
# covalent5 · 비공유 도킹 (RQ1 '잘못된 프록시' 축)
# ligand_manifest.csv 기준 (구버전 파일 무시). hAChE gorge에 도킹, CNN affinity 산출.
set -e
cd "$(dirname "$0")"
REC=receptor.pdbqt
CX=-5.337; CY=-41.735; CZ=29.543; SZ=22.5
OUT=out_noncovalent; mkdir -p "$OUT"
MAN=ligands/ligand_manifest.csv

tail -n +2 "$MAN" | while IFS=, read sid code cls role noncov cov leaving hascov; do
  lig="ligands/$noncov"
  [ -f "$lig" ] || { echo "  [건너뜀] $sid 파일없음"; continue; }
  echo ">>> [noncov] $sid"
  gnina -r "$REC" -l "$lig" \
    --center_x $CX --center_y $CY --center_z $CZ \
    --size_x $SZ --size_y $SZ --size_z $SZ \
    --exhaustiveness 16 --num_modes 9 \
    --cnn_scoring rescore \
    -o "$OUT/${sid}_nc.sdf.gz" --log "$OUT/${sid}_nc.log"
done
echo "비공유 도킹 완료 -> $OUT/  (parse_docking.py로 점수 추출)"
