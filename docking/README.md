# covalent5 · Step 4-③ GNINA 도킹 패키지

hAChE에 대한 (A) 비공유 도킹 — RQ1의 '잘못된 프록시' 축 — 과 (B) 공유 도킹(P→Ser203 Oγ)
설정. **계산은 Ku의 RTX 3090 머신(GNINA, GPU CNN 스코어링)에서 실행.**

## 수용체
- **4EY7** = 재조합 human AChE + 도네페질. 체인 A 단백질만 추출(`receptor_4EY7_chainA.pdb`,
  도네페질·글리칸·물·타체인 제거). 촉매 삼중쌍 Ser203/His447/Glu334 완비.
- 도킹 박스 중심 = **Ser203 Oγ (-5.337, -41.735, 29.543)**, gorge 깊이 고려 size 18–22.5 Å.

## 실행 순서
```bash
bash prepare_receptor.sh        # 양성자화 -> receptor.pdbqt / receptor_H.pdb
bash run_gnina_noncovalent.sh   # 비공유: intact 리간드 30종
bash run_gnina_covalent.sh      # 공유: 활성 인산화제 26종 (티온 제외)
python parse_docking.py         # 점수 -> docking_results.csv
```

## 리간드
- `ligands/noncovalent/<sid>.sdf` : intact 작용제·OP (30종)
- `ligands/covalent/<sid>.sdf` : 이탈기 제거 warhead + [*]더미(=Ser-Oγ 부착점). 26종.
  - 이탈기: G/Novichok = F, Tabun = CN, V = thiocholine(S), oxon = 아릴/비닐옥시.
  - **티온(P=S: parathion/chlorpyrifos/malathion)은 공유 warhead 미생성** — 직접 인산화 안 함(생체활성화 필요). 비공유만.

## 검증 (필수)
공유 도킹 결과는 **PDB 6CQT(VX⁻)/6CQX(VX⁺)** 의 실험적 부가물 기하와 대조:
- VX P3R/P3S warhead 공유 도킹 포즈 ↔ 6CQT/6CQX의 Ser203-P 결합 거리·각도 재현되는지 확인.
- 재현되면 Novichok 적용 신뢰. 안 되면 GNINA 버전별 covalent 플래그/[*] 처리 조정.

## RQ1 분석
`docking_results.csv`의 비공유 affinity 순위 ↔ QM 클러스터 ΔG‡(Step 4-④) 순위를
Spearman ρ로 비교. ρ 낮음 = 비공유 도킹이 비가역 억제제를 오정렬(H1).

## 주의
- GNINA covalent 거동은 버전 의존적(`--covalent_rec_atom`, `--covalent_lig_atom_pattern`).
  본 스크립트는 표준 플래그 기준. [*]더미 처리 방식은 설치본에 맞게 조정.
- 양성자화 후 His447 토토머 상태 수동 점검 권장.
