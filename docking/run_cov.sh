#!/bin/bash
cd "$(dirname "$0")"
REC=receptor_H.pdb
RECATOM="D:149:OG"
CX=-5.337; CY=-41.735; CZ=29.543; SZ=18.0
OUT=out_covalent; mkdir -p "$OUT"
sed 's/\r$//' ligands/ligand_manifest.csv | tail -n +2 | \
while IFS=, read -r sid code cls role noncov cov leaving hascov; do
  hascov="$(echo "$hascov" | tr -d ' \r\n')"
  [ "$hascov" = "1" ] || continue
  cov="$(echo "$cov" | tr -d ' \r\n')"
  lig="ligands/$cov"
  [ -f "$lig" ] || { echo "  [skip] $sid (no $lig)"; continue; }
  echo ">>> $sid"
  gnina -r "$REC" -l "$lig" \
    --covalent_rec_atom "$RECATOM" --covalent_lig_atom_pattern "[#15]" \
    --covalent_optimize_lig --covalent_bond_order 1 \
    --center_x $CX --center_y $CY --center_z $CZ \
    --size_x $SZ --size_y $SZ --size_z $SZ \
    --exhaustiveness 16 --num_modes 9 --cnn_scoring none \
    -o "$OUT/${sid}_cov.sdf.gz" --log "$OUT/${sid}_cov.log" \
    >/dev/null 2>&1 && echo "    ok $sid" || echo "    FAIL $sid ($?)"
done
echo "=== done. logs: $(ls "$OUT"/*_cov.log 2>/dev/null | wc -l) ==="
