# covalent5 · Step 4-④ QM 클러스터 ΔG‡ 프로토콜 (설계)

목표: AChE Ser203 인산화 반응의 활성화 자유에너지 ΔG‡(±aging)를 입체이성질체별로 계산.
이것이 RQ1(비공유 도킹 순위 ≠ 공유 ΔG‡ 순위)·RQ2(입체선택성)의 **핵심 라벨**.

## 1. 메커니즘 (모델링 대상)
일반염기 촉매 인산화:
- His447이 Ser203 Oγ–H의 양성자를 받아 알콕사이드 활성화
- Ser203 Oγ⁻가 작용제 P를 친핵 공격 → **오각형 이중피라미드 TS**
- 이탈기 이탈 (G/Novichok=F⁻, Tabun=CN⁻, V=티올레이트)
- (선택) aging: P–O–C 탈알킬화 → 음이온성 포스포네이트
반응좌표 = **P···Oγ 거리** (3.0 → 1.6 Å). 옥시음이온 홀이 발달하는 음전하 안정화.

## 2. 클러스터 정의
**최소 모델 (1차):** 촉매 삼중쌍 + 기질
- Ser203 (사이드체인, Cβ 절단·H캡), His447 (이미다졸), Glu334 (카복실레이트), 작용제
**확장 모델 (정밀):** + 옥시음이온 홀 백본 아미드 (Gly121, Gly122, Ala204 N–H)
- 절단: 사이드체인은 Cα–Cβ 절단 후 Cβ에 H캡(Cα 방향, 1.09 Å)
- 제약: 캡 H(및 Cβ)를 **좌표 고정**(freeze)하여 단백질 골격 강성 모사

## 3. 출발 구조
- **VX 파일럿:** 6CQX(VX 부가물)에서 클러스터 추출 → 생성물(PC). 역방향 스캔으로 TS 탐색
  또는 비공유 도킹 포즈에서 반응물(Michaelis, RC) 구성. **6CQX 실험 구조로 검증** 가능.
- Novichok: VX 프로토콜 검증 후 동일 절단·동일 시점으로 이식.
- ※ 결정구조엔 H 없음 → 기질·잔기에 H 추가(반응물=중성 His, Ser–OH; 활성형은 His⁺/Ser–O⁻ 경로도 비교).

## 4. 계산 방법 (ORCA)
- RC/PC 최적화: `! wB97X-D3 def2-SVP def2/J RIJCOSX OPT CPCM(water) TightSCF defgrid2` + 제약
- P···Oγ relaxed scan: `%geom Scan B <iP> <iOG> = 3.0, 1.6, 15 end end`
- TS 최적화: `! OptTS NumFreq` + 스캔 최고점 guess + 활성원자 Hessian
- 진동수: TS는 **허수진동 1개**(P–Oγ 형성 모드) 확인 필수
- 단일점: `! wB97X-D3 def2-TZVP ... SP` (ε=4 단백질 또는 water 비교)
- ΔG‡ = G(TS) − G(RC), aging은 별도 스캔(P–O–C)

## 5. 검증 (필수)
- VX RP/SP의 ΔG‡ 차이가 실험적 입체특이 동역학(6CQT/6CQX + 문헌 kᵢ)과 **부호·크기 정합**하는지.
- 정합 시 Novichok 외삽 신뢰. 이게 RQ2의 계산-실험 교차검증.

## 6. 산출물
`qm_cluster/<sid>/{rc.inp, scan.inp, ts.inp, freq.inp, sp.inp}` + ΔG‡ 요약표
→ 비공유 도킹 순위와 Spearman ρ 비교(RQ1, Fig 4).

## 7. 단계 (DFT 완료 후 인터랙티브 — TS 탐색은 babysitting 필요)
1. VX 파일럿 클러스터 추출·H추가·시각 점검
2. RC 최적화 → P–Oγ 스캔 → TS opt → freq(허수 1개)
3. ΔG‡ 산출, 6CQX/kᵢ 검증
4. 검증 통과 → Novichok 4종 × 입체이성질체 확장
5. (선택) 옥시음이온 홀 확장 모델로 정밀화

## 주의
- TS 탐색은 자동화가 불안정 → 분자별 확인하며 진행(DFT처럼).
- 최소 모델로 경향성 먼저, 확장 모델로 정밀화.
