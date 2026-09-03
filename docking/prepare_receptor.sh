#!/bin/bash
# covalent5 · Step 4-③ 수용체 준비
# 4EY7 chain A (도네페질/글리칸/물 제거 완료) -> 양성자화 -> 도킹용 포맷
# 요구 도구: openbabel(obabel). 더 정밀하게는 reduce 또는 pdb2pqr(pH 7.4) 권장.
set -e
cd "$(dirname "$0")"

IN=receptor_4EY7_chainA.pdb

# 1) 양성자화 (pH 7.4) — His/Glu 프로토네이션 상태 주의(촉매 삼중쌍 His447)
obabel "$IN" -O receptor_H.pdb -p 7.4

# 2) GNINA용 pdbqt (CNN 스코어링은 pdb도 가능하나 pdbqt 권장)
obabel receptor_H.pdb -O receptor.pdbqt -xr

echo "완료: receptor_H.pdb, receptor.pdbqt"
echo "도킹 박스 중심 = Ser203 Oγ: -5.337 -41.735 29.543 (size 22.5 Å)"
echo "※ 권장: 양성자화 후 His447 NE2/ND1, Glu334 상태를 수동 점검(촉매 삼중쌍)"
