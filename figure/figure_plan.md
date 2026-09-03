# covalent5 · 피규어 구성안

타깃 Chem Sci(2단). 각 피규어는 1장씩 생성→피드백→다음 (파이프라인 규칙).
데이터 무결성: 측정/계산 확보된 것만 실제 플롯, 미확보는 레이아웃만(플레이스홀더).

| Fig | 제목 | 내용 | 필요 데이터 | 현재 가능? |
|-----|------|------|-------------|-----------|
| 1 | Concept/workflow | 닫힌 루프(DFT→도킹→QM ΔG‡→GPR→Ellman→갱신) 개념도 | — (스키마) | ✅ 지금 |
| 2 | 분자·입체화학 | A-230/232/234/242 구조 + P R/S × C=N E/Z 입체이성질체 맵 + 앵커/실험 OP | 구조(확보) | ✅ 지금 |
| 3 | 비공유 도킹 결과 | affinity/CNN by 화합물·입체이성질체, 좁은 범위·입체 near-degeneracy 강조 | docking_results.csv ✅ | ✅ 지금 |
| 4 | RQ1: 비공유 vs 공유 | 비공유 affinity 순위 ↔ QM 클러스터 ΔG‡ 순위 (Spearman ρ) | QM ΔG‡ (대기) | ⏳ QM 후 |
| 5 | RQ2: 입체선택성 | enantiomer별 ΔG‡ 차이 + VX RP/SP 실험 앵커(6CQT/6CQX) | QM ΔG‡ + 실험 | ⏳ |
| 6 | RQ3: 능동학습 | GPR hold-out 성능·보정 CI, baseline(RF/XGB) 대비, AL 곡선 | AL 결과 (대기) | ⏳ |
| 7 | 도킹 포즈 패널 | hAChE gorge 내 best pose(촉매 삼중쌍 표시), 대표 작용제 3~4종 | 포즈(.sdf.gz) ✅ | ▶ PyMOL 스크립트 |
| 8 | MD | 부가물/복합체 안정성(RMSD/RMSF), gorge 접촉 | MD (대기) | ⏳ |

**지금 만들 수 있는 것:** Fig 1(개념도), Fig 2(분자맵), Fig 3(도킹 데이터), Fig 7(포즈 이미지=PyMOL).
**대기:** Fig 4/5/6/8 (QM·AL·MD·실험 결과 확보 후).

## 포즈 이미지(Fig 7)
`render_poses.py`로 PyMOL에서 receptor + best pose + 촉매 삼중쌍 정렬 렌더 → `figure/poses/*.png`.
정렬: 모든 포즈를 Ser203 활성부위 기준 동일 시점으로 → 비교 가능한 패널.
