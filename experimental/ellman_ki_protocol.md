# covalent5 · Ellman anti-AChE kᵢ 측정 프로토콜 초안 (Tier 1–3)

> **범위·안전:** 본 문서는 **효소 저해 동역학 측정(분석·방어 연구)** 설계만 다룬다. 모든 작용제
> 취급·희석·폐기는 기관 **OPCW 인가 SOP** 내에서만 수행한다. **작용제 합성·제조 경로는 포함하지
> 않는다.** 측정 목적은 위협평가용 상대 효력(kᵢ) 산출.

## 0. 이 데이터가 논문에서 하는 일
- **Tier 1** 고전작용제+비규제 OP 자체측정 kᵢ → 문헌 이질성 제거, 균질 벤치마크(현 `data/anchor_ki.csv` 대체).
- **Tier 2** Novichok A-시리즈 racemic kᵢ → 물리모델 예측을 실제 타깃에서 검증("데이터 0" 해소).
- **Tier 3** enantiomer resolved kᵢ → H2(입체선택성) 직접검증 + 닫힌루프 AL(예측→측정→갱신).
→ 산출 kᵢ가 본문 numerical claim의 **재현 가능한 측정 source** ([[feedback_reproducibility_before_writing]] 준수).

## 1. 원리
Ellman 법: AChE가 acetylthiocholine(ATChI)를 가수분해 → thiocholine → DTNB와 반응해
5-thio-2-nitrobenzoate(TNB, 노란색) 생성. **412 nm 흡광 증가율 = 효소 활성.**
비가역 OP 억제제는 시간에 따라 활성을 소실시키므로, **잔존 활성의 시간·농도 의존성**에서
이분자 저해 속도상수 kᵢ를 구한다.

## 2. 재료·장비 (기관 보유 확인)
- 재조합 **human AChE** (인간 관련성; recombinant, e.g. erythrocyte/expressed)
- **DTNB** (Ellman's reagent), **acetylthiocholine iodide (ATChI)**
- 완충액: 0.1 M 인산완충액 pH 7.4 (또는 8.0), 37 °C
- 마이크로플레이트 리더(또는 분광광도계), 412 nm, 37 °C 항온
- (Tier 3) 키랄 분리: 키랄 HPLC 또는 GC 컬럼
- QC: **paraoxon-ethyl**(양성대조; Worek2004 kᵢ≈2.2×10⁶ M⁻¹min⁻¹ 기지값과 대조 → 방법 검증)

## 3. 사전 파라미터 (측정 전 확정)
1. **효소 활성 baseline**: AChE 농도별로 412 nm 증가율 선형구간 확인 → 작업 효소량 결정
   (초기속도가 5–10분 선형 유지되는 수준).
2. **ATChI Km**: ATChI 농도구간(예 0.05–2 mM)에서 Michaelis-Menten → Km 산출.
   저해 측정 시 기질보호(substrate protection) 보정을 위해 사용. DTNB ≈ 0.3 mM 고정.
3. **자발 가수분해 대조**: 효소 없이 ATChI+DTNB 흡광 증가(배경) → 모든 값에서 차감.

## 4. Tier 1–2: 이분자 kᵢ (비가역) — 표준 2가지 중 택1

### (A) 잔존활성법 (Aldridge / Hart & O'Brien)
1. 효소를 **여러 억제제 농도 [I]**(예 4–6점)와 37 °C 예비배양(별도 튜브, 기질 없음).
2. 정해진 **시간점 t**(예 0, 1, 2, 4, 6, 10 min)마다 분취 → ATChI/DTNB 함유 웰에 희석 첨가 →
   즉시 412 nm 초기속도(=잔존활성) 측정.
3. 각 [I]에서 ln(잔존활성 %) vs t 선형회귀 → 기울기 = −k_obs (유사1차).
4. **k_obs vs [I] 선형회귀 → 기울기 = kᵢ (M⁻¹min⁻¹).** 절편=자발/배경.
   ※ 예비배양 중 기질 없음 → substrate protection 무. [I]≫[E]로 유사1차 성립.

### (B) 연속법 (Worek continuous)
1. 큐벳/웰에 완충액+DTNB+ATChI+효소로 반응 시작 후 **억제제 첨가**, 412 nm를 연속 기록.
2. 진행성 저해곡선(속도 감소)을 1차 실효식으로 피팅 → k_obs → [I] 변화로 kᵢ.
   ※ 기질 존재 → substrate protection을 Km/[S]로 보정: kᵢ(true)=kᵢ(app)·(1+[S]/Km).

> **권장:** Tier 1 검증은 (A)로(직관적·문헌 비교 용이), 처리량 필요 시 (B) 병행.

### 농도·시간 설계(예시, 사전 Km/활성으로 조정)
- [I]: kᵢ 규모에 맞춰 반감기 t½=ln2/(kᵢ·[I])가 **1–10 min** 범위 되게 선정.
  (강력작용제 GD/VX/BSAR는 pM–nM, 약한 OP는 µM–mM 필요할 수 있음)
- 시간점 ≥5, 잔존활성 ~80%→~20% 구간 포착.
- **반복 n≥3**(독립 실험), 웰 triplicate.

## 5. Tier 3: 입체특이 kᵢ
1. racemic 작용제를 **키랄 HPLC/GC로 enantiomer 분리**(분석 수준; 합성 아님).
2. 각 enantiomer 분획으로 §4(A) kᵢ 측정 → Sp/Rp resolved kᵢ.
3. VX Sp/Rp 먼저(6CQT/6CQX 구조·문헌 대조 가능) → Novichok 확대.
4. **닫힌루프 AL:** GPR 예측 불확실도 큰 화합물 우선 측정 → 모델 갱신 → 반복.

## 6. 대조·QC (데이터 무결성)
- **양성대조 paraoxon-ethyl**: 측정 kᵢ가 Worek2004(2.2×10⁶) 오차범위 재현 → 방법 유효성 확인. **선행 필수.**
- 자발 가수분해 대조, DTNB 배경, 효소-only(비저해) 대조.
- **Aging 체크(선택, Novichok 중요):** 저해 후 oxime(2-PAM/HI-6) 재활성화 시도 →
  aged vs non-aged 구분(재활성화 안 되면 aging 진행). 위협평가 함의(의료대응).
- 단위 통일: **kᵢ = M⁻¹ min⁻¹**, pKᵢ=log₁₀(kᵢ). racemic vs resolved 명시.

## 7. 산출·데이터 포맷 (모델 입력)
`experimental/measured_ki.csv`:
`stereo_id, agent, config(Sp/Rp/racemic/achiral), ki, pKi, ki_SE, n_rep, method(A/B), source(self-measured), date`
→ 기존 `data/anchor_ki.csv`와 동일 스키마로 대체·병합. 모델 재학습·검증에 직접 투입.

## 8. 측정 우선순위(현실 목표)
1. QC: paraoxon-ethyl (방법검증) — **먼저**
2. Tier 1: GB, GD, GA, VX, VR + 비규제 OP(paraoxon/methyl-paraoxon/DFP 등)
3. Tier 2: **A-234 우선**(2025 CBI hAChE 문헌 대조 + Skripal 서사) → A-230/232/242
4. Tier 3: VX(Sp/Rp) → A-234 resolved

## 관련
[[project_covalent5]] [[JHM_framing_strategy]] [[feedback_reproducibility_before_writing]] [[feedback_active_learning_methodology]]
