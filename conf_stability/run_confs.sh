#!/bin/bash
set -e
ORCA="${ORCA_BIN:-orca}"
ROOT="/mnt/c/Users/kkaan/Downloads/claude_research/covalent5/conf_stability"
cd "$ROOT/VX_00_P3S/conf0"
"$ORCA" opt.inp > opt.out 2>&1
"$ORCA" sp.inp  > sp.out  2>&1
echo "DONE VX_00_P3S/conf0"
cd "$ROOT/VX_00_P3S/conf1"
"$ORCA" opt.inp > opt.out 2>&1
"$ORCA" sp.inp  > sp.out  2>&1
echo "DONE VX_00_P3S/conf1"
cd "$ROOT/VX_00_P3S/conf2"
"$ORCA" opt.inp > opt.out 2>&1
"$ORCA" sp.inp  > sp.out  2>&1
echo "DONE VX_00_P3S/conf2"
cd "$ROOT/VR_00_P8S/conf0"
"$ORCA" opt.inp > opt.out 2>&1
"$ORCA" sp.inp  > sp.out  2>&1
echo "DONE VR_00_P8S/conf0"
cd "$ROOT/VR_00_P8S/conf1"
"$ORCA" opt.inp > opt.out 2>&1
"$ORCA" sp.inp  > sp.out  2>&1
echo "DONE VR_00_P8S/conf1"
cd "$ROOT/VR_00_P8S/conf2"
"$ORCA" opt.inp > opt.out 2>&1
"$ORCA" sp.inp  > sp.out  2>&1
echo "DONE VR_00_P8S/conf2"
cd "$ROOT/DETAB_00_P3R/conf0"
"$ORCA" opt.inp > opt.out 2>&1
"$ORCA" sp.inp  > sp.out  2>&1
echo "DONE DETAB_00_P3R/conf0"
cd "$ROOT/DETAB_00_P3R/conf1"
"$ORCA" opt.inp > opt.out 2>&1
"$ORCA" sp.inp  > sp.out  2>&1
echo "DONE DETAB_00_P3R/conf1"
cd "$ROOT/DETAB_00_P3R/conf2"
"$ORCA" opt.inp > opt.out 2>&1
"$ORCA" sp.inp  > sp.out  2>&1
echo "DONE DETAB_00_P3R/conf2"
cd "$ROOT/FEN_00_P3R/conf0"
"$ORCA" opt.inp > opt.out 2>&1
"$ORCA" sp.inp  > sp.out  2>&1
echo "DONE FEN_00_P3R/conf0"
cd "$ROOT/FEN_00_P3R/conf1"
"$ORCA" opt.inp > opt.out 2>&1
"$ORCA" sp.inp  > sp.out  2>&1
echo "DONE FEN_00_P3R/conf1"
cd "$ROOT/FEN_00_P3R/conf2"
"$ORCA" opt.inp > opt.out 2>&1
"$ORCA" sp.inp  > sp.out  2>&1
echo "DONE FEN_00_P3R/conf2"
echo "=== 전체 완료. python parse_conf_stability.py ==="
