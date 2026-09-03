#!/bin/bash
# covalent5 · 무인 도킹 일괄 실행 (DFT와 병행 — GPU 사용)
# 사용:
#   export COV5="/mnt/c/Users/kkaan/Downloads/claude_research/covalent5"
#   cd "$COV5/docking"
#   nohup bash run_all_docking.sh > docking_all.log 2>&1 &
# 진행: tail -f "$COV5/docking/docking_all.log"
set -e
cd "$(dirname "$0")"
# gnina(cuda12.8 빌드)가 conda의 cuDNN/cuBLAS를 찾도록 로더 경로 지정
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH}"
echo "===== [1/4] gnina 확인 $(date) ====="
gnina --version | head -2

echo "===== [2/4] 수용체 준비 $(date) ====="
if [ ! -f receptor.pdbqt ] || [ ! -f receptor_H.pdb ]; then
  bash prepare_receptor.sh
else
  echo "수용체 이미 있음, 건너뜀"
fi

echo "===== [3/4] 비공유 도킹 30종 $(date) ====="
bash run_gnina_noncovalent.sh

echo "===== [4/4] 공유 도킹 26종 $(date) ====="
bash run_gnina_covalent.sh

echo "===== 점수 파싱 $(date) ====="
python parse_docking.py

echo "===== 도킹 전부 완료 $(date) ====="
echo "결과: docking/docking_results.csv, out_noncovalent/, out_covalent/"
