# covalent5 · 공유 부가물 MD 계획 (GROMACS + AmberTools)

목표: hAChE Ser203에 공유결합한 인산화 부가물의 입체이성질체별 안정성·gorge 상호작용 MD → H2 지지.

## 파일럿: VX 부가물 (검증 가능)
- 출발 구조: `vx_adduct_6CQX_chainA.pdb` — 6CQX 체인 A 단백질 + VX-Ser203 공유 부가물(P1·O1·O2·C1·C2·C3).
  실험적 P–Oγ 결합(1.32 Å) 포함 → 실험 구조라 셋업 검증에 이상적.

## 파라미터화 경로 (AmberTools, GAFF2/AM1-BCC + ff19SB)
접합부(Ser203 Oγ–P1)가 섬세 → 단계별 검증하며 진행.

### 단계 (설치 후 인터랙티브로 하나씩)
1. **구조 준비** — 단백질/부가물 분리, 부가물 프래그먼트 캡핑(Ser Oγ 측 처리)
2. **부가물 파라미터** — `antechamber`(GAFF2, AM1-BCC) + `parmchk2` → frcmod
3. **tleap 빌드** — ff19SB 로드, 부가물 라이브러리 로드, `bond` 로 P1–OG 공유 연결,
   TIP3P 용매화 + 중화 → prmtop/inpcrd
4. **GROMACS 변환** — `acpype`(또는 ParmEd) prmtop → `.top`/`.gro`
5. **실행** — EM → NVT(100 ps) → NPT(100 ps) → production(예: 100 ns), 3090 GPU
6. **분석** — RMSD/RMSF, P–Oγ 결합 안정성, gorge 잔기(Trp86, His447) 접촉, 입체이성질체 비교

## 확장
VX 파일럿 검증 후 → Novichok A-230/232/234/242 입체이성질체에 phosphyl 이식
(6CQX 부가물 기하 템플릿 사용) → 동일 파이프라인.

## 설치 (선행)
```bash
conda install -c conda-forge ambertools acpype -y
which antechamber tleap parmchk2 acpype   # 확인
```

## 주의
- 접합부 전하/원자타입 정합이 핵심 — 단계 2~3에서 검증.
- 인산화 Ser는 비표준 잔기 → tleap 라이브러리 수동 정의 필요할 수 있음.
- 본 MD는 보조 축(주축은 QM 클러스터 ΔG‡). DFT/도킹 우선.
