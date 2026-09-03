#!/bin/bash
# VX_00_P3S — opt+freq 후 single-point. ORCA는 병렬 시 풀 경로 필요.
# 실행 전: export ORCA_BIN=/home/k9/orca/orca
set -e
ORCA="${ORCA_BIN:-orca}"
"$ORCA" opt.inp > opt.out 2>&1
"$ORCA" sp.inp  > sp.out  2>&1
echo "DONE VX_00_P3S"
