#!/bin/bash
# covalent5 · 도킹 robustness/sensitivity (리뷰어 대응, WSL/GNINA)
# 단일 셋업 결론이 box 위치·random seed에 의존하지 않음을 보인다.
# 변형: 기준 box(±2 A 이동) + 대체 seed. 각 변형마다 전 리간드 비공유 도킹.
# 실행:  cd $COV5/docking && export ORCA... (불필요), gnina PATH 확인 후
#        bash run_sensitivity.sh 2>&1 | tee sens_run.log
# 결과:  out_sens_<label>/  → python parse_sensitivity.py 로 변형별 Spearman 산출
set -e
cd "$(dirname "$0")"
REC=receptor.pdbqt
CX=-5.337; CY=-41.735; CZ=29.543; SZ=22.5
MAN=ligands/ligand_manifest.csv

# label  dx  dy  dz  seed
VARIANTS=(
  "base    0   0   0   42"
  "boxxp2  2   0   0   42"
  "boxxm2 -2   0   0   42"
  "boxyp2  0   2   0   42"
  "boxym2  0  -2   0   42"
  "seed7   0   0   0    7"
  "seed123 0   0   0  123"
)

dock_variant () {
  local label="$1" dx="$2" dy="$3" dz="$4" seed="$5"
  local cx cy cz out
  cx=$(python3 -c "print($CX+$dx)"); cy=$(python3 -c "print($CY+$dy)"); cz=$(python3 -c "print($CZ+$dz)")
  out="out_sens_${label}"; mkdir -p "$out"
  echo "=== variant $label  center=($cx,$cy,$cz) seed=$seed ==="
  tail -n +2 "$MAN" | while IFS=, read sid code cls role noncov cov leaving hascov; do
    lig="ligands/$noncov"; [ -f "$lig" ] || continue
    gnina -r "$REC" -l "$lig" \
      --center_x $cx --center_y $cy --center_z $cz \
      --size_x $SZ --size_y $SZ --size_z $SZ \
      --exhaustiveness 16 --num_modes 9 --cnn_scoring rescore \
      --seed $seed \
      -o "$out/${sid}_nc.sdf.gz" --log "$out/${sid}_nc.log" >/dev/null 2>&1 \
      && echo "  ok $sid" || echo "  FAIL $sid"
  done
}

for v in "${VARIANTS[@]}"; do
  set -- $v
  dock_variant "$1" "$2" "$3" "$4" "$5"
done
echo "=== 완료. python parse_sensitivity.py 로 변형별 Spearman rho 산출 ==="
