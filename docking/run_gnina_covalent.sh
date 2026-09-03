#!/bin/bash
# covalent5 · 공유 도킹 (P -> Ser203 Oγ 인산화 부가물)
# ligand_manifest.csv에서 has_covalent=1 인 활성 인산화제만 (티온 P=S 자동 제외).
# warhead: 이탈기 제거 + [*]더미(=Ser-Oγ 부착점).
# ※ GNINA 공유 도킹 거동은 버전 의존적. 6CQT/6CQX VX 부가물 기하 재현으로 반드시 검증.
set -e
cd "$(dirname "$0")"
REC=receptor_H.pdb
CX=-5.337; CY=-41.735; CZ=29.543; SZ=18.0
OUT=out_covalent; mkdir -p "$OUT"
MAN=ligands/ligand_manifest.csv

tail -n +2 "$MAN" | while IFS=, read sid code cls role noncov cov leaving hascov; do
  [ "$hascov" = "1" ] || continue
  lig="ligands/$cov"
  [ -f "$lig" ] || { echo "  [건너뜀] $sid warhead 없음"; continue; }
  echo ">>> [cov] $sid (leaving=$leaving)"
  gnina -r "$REC" -l "$lig" \
    --covalent_rec_atom A:203:OG \
    --covalent_lig_atom_pattern "[#15]" \
    --covalent_optimize_lig --covalent_bond_order 1 \
    --center_x $CX --center_y $CY --center_z $CZ \
    --size_x $SZ --size_y $SZ --size_z $SZ \
    --exhaustiveness 16 --num_modes 9 \
    --cnn_scoring rescore \
    -o "$OUT/${sid}_cov.sdf.gz" --log "$OUT/${sid}_cov.log"
done
echo "공유 도킹 완료 -> $OUT/"
echo "※ [*]더미를 Ser-Oγ로 처리 안 하는 버전이면 warhead에서 [*] 제거 후 재시도."
