#!/usr/bin/env python3
"""
covalent5 · Figure 4 (R2 재생성)
  (a) QM 클러스터 P···Oγ relaxed scan  — 실측 18점, "barrierless" 제거, 불연속 명시
  (b) 효소 P-입체선택성 (문헌)
  (c) VX 이성질체 실측 log k_i
  (d) [신규] leave-both-VX-out 외부검증 — 1차 캡션이 약속했으나 그림에 없던 패널
모든 수치 출처: qm_cluster/scan_alkoxide.out, data/anchor_ki.csv,
                data/anchor_ki_worek.csv, data/ml_robustness_out.md
"""
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams.update({
    'font.size': 8, 'axes.titlesize': 8.5, 'axes.labelsize': 8,
    'xtick.labelsize': 7, 'ytick.labelsize': 7, 'legend.fontsize': 6.8,
    'font.family': 'sans-serif', 'font.sans-serif': ['DejaVu Sans'],
    'axes.spines.top': False, 'axes.spines.right': False, 'axes.linewidth': 0.8,
    'axes.axisbelow': True, 'axes.grid': True, 'grid.linewidth': 0.5,
    'grid.alpha': 0.35, 'legend.frameon': False,
    'savefig.dpi': 300, 'savefig.bbox': 'tight', 'savefig.pad_inches': 0.06,
})
import matplotlib.pyplot as plt
import numpy as np

# 검증된 팔레트 (dataviz validate_palette.js 전 항목 PASS)
C_QM     = "#6F5BD5"   # QM 스캔
C_G      = "#4C7FE0"   # G 계열
C_V      = "#F4611A"   # V 계열 (VX)
C_PRED   = "#0F9D8F"   # 모델 예측
INK      = "#1f2328"
MUTED    = "#6B7280"

H2K = 627.5094740631

# ── (a) scan_alkoxide.out "Actual Energy" 표면 (Hartree) ─────────────────
scan = [(3.30, -2027.32227728), (3.20, -2027.32250943), (3.10, -2027.32255759),
        (3.00, -2027.32216380), (2.90, -2027.33065853), (2.80, -2027.33034085),
        (2.70, -2027.32961862), (2.60, -2027.32905612), (2.50, -2027.32883610),
        (2.40, -2027.32922961), (2.30, -2027.33043796), (2.20, -2027.33260031),
        (2.10, -2027.33574130), (2.00, -2027.33985789), (1.90, -2027.36017971),
        (1.80, -2027.37347600), (1.70, -2027.38311105), (1.60, -2027.38548685)]
r = np.array([p[0] for p in scan])
E = (np.array([p[1] for p in scan]) - scan[0][1]) * H2K

fig = plt.figure(figsize=(7.0, 5.6))
gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.30)

# ══ (a) ═══════════════════════════════════════════════════════════════════
ax = fig.add_subplot(gs[0, 0])
ax.axvspan(2.90, 3.00, color="#9CA3AF", alpha=0.16, lw=0, zorder=0)
ax.axhline(0, color=MUTED, lw=0.6, ls=(0, (4, 3)), zorder=1)
ax.plot(r, E, '-', color=C_QM, lw=1.6, zorder=3)
ax.plot(r, E, 'o', color=C_QM, ms=3.6, mec='white', mew=0.6, zorder=4)

ax.annotate("discontinuous 5.3 kcal mol$^{-1}$ step\n(geometry reorganisation)",
            xy=(2.95, -2.4), xytext=(2.60, 11.5), fontsize=6.2, color=INK, ha='center',
            arrowprops=dict(arrowstyle='-|>', color=INK, lw=0.7, shrinkA=0, shrinkB=2))
ax.annotate("+1.1 kcal mol$^{-1}$ rise\n(2.9 $\\to$ 2.5 Å)",
            xy=(2.50, -4.6), xytext=(2.18, -23.0), fontsize=6.2, color=INK, ha='center',
            arrowprops=dict(arrowstyle='-|>', color=INK, lw=0.7, shrinkA=0, shrinkB=2))
ax.annotate("$-$39.7 kcal mol$^{-1}$ adduct",
            xy=(1.60, -40.4), xytext=(1.72, -46.0), fontsize=6.2, color=C_QM, ha='center',
            arrowprops=dict(arrowstyle='-|>', color=C_QM, lw=0.7, shrinkA=0, shrinkB=2))

ax.set_xlim(3.42, 1.48)
ax.set_ylim(-49, 17)
ax.set_xlabel("P$\\cdots$O$_\\gamma$ (Å)")
ax.set_ylabel("Relative $E$ (kcal mol$^{-1}$)")
ax.set_title("(a) QM cluster approach coordinate", loc='left', fontweight='bold')
ax.text(0.02, 0.03, "70-atom cluster, charge $-$1\n"
        "VX + Ser203/His447/Glu334 side chains\n"
        r"$\omega$B97X-D3/def2-SVP, CPCM(water)" "\nrelaxed scan; no TS located",
        transform=ax.transAxes, ha='left', va='bottom', fontsize=5.6, color=MUTED,
        linespacing=1.35)

# ══ (b) ═══════════════════════════════════════════════════════════════════
ax = fig.add_subplot(gs[0, 1])
names = ["Soman (GD)", "Sarin (GB)", "VX"]
ratio = [40000.0, 4200.0, 1.4e8 / 1.2e6]          # VX = 116.7 (anchor_ki.csv)
lbl   = ["> 40,000$\\times$", "4,200$\\times$", "117$\\times$"]
cols  = [C_G, C_G, C_V]
y = np.arange(len(names))[::-1]
vals = np.log10(ratio)
ax.barh(y, vals, height=0.56, color=cols, zorder=3)
# Soman은 하한값 → 열린 끝 표시
ax.annotate("", xy=(vals[0] + 0.42, y[0]), xytext=(vals[0], y[0]),
            arrowprops=dict(arrowstyle='-|>', color=C_G, lw=1.4))
for yi, v, t in zip(y, vals, lbl):
    ax.text(v + (0.52 if t.startswith('>') else 0.12), yi, t, va='center',
            fontsize=6.8, color=INK, fontweight='bold')
ax.set_yticks(y); ax.set_yticklabels(names)
ax.set_xlim(0, 6.0)
ax.set_ylim(-1.30, 2.50)          # 하단 여백: 출처 주석이 막대·값라벨과 겹치지 않게
ax.set_xlabel("$\\log_{10}$ ( $k_i$ more-toxic / less-toxic P-isomer )")
ax.set_title("(b) Enzyme P-stereoselectivity (published $k_i$)", loc='left', fontweight='bold')
ax.grid(axis='y', visible=False)
ax.text(0.985, 0.025, "Rate constants replotted from\nWorek et al. 2004; Bester et al. 2018",
        transform=ax.transAxes, ha='right', va='bottom', fontsize=6.0, color=MUTED,
        linespacing=1.35)

# ══ (c) ═══════════════════════════════════════════════════════════════════
ax = fig.add_subplot(gs[1, 0])
meas = [8.146, 6.079]
xs = np.arange(2)
ax.bar(xs, meas, width=0.5, color=[C_V, "#F9B48D"], zorder=3)
for x, v, a in zip(xs, meas, ["1.4$\\times$10$^{8}$", "1.2$\\times$10$^{6}$"]):
    ax.text(x, v + 0.16, a, ha='center', fontsize=6.8, color=INK, fontweight='bold')
ax.set_xticks(xs); ax.set_xticklabels(["$S_\\mathrm{P}$ (toxic)", "$R_\\mathrm{P}$"])
ax.set_ylim(0, 9.6)
ax.set_ylabel("measured $\\log k_i$ (hAChE)")
ax.set_title("(c) VX enantiomers, measured", loc='left', fontweight='bold')
ax.grid(axis='x', visible=False)

# ══ (d) 신규 ══════════════════════════════════════════════════════════════
ax = fig.add_subplot(gs[1, 1])
pred = [8.50, 6.48]
sd   = [0.81, 2.33]
ax.bar(xs, pred, width=0.5, color=C_PRED, zorder=3, label="predicted (VX pair withheld)")
ax.errorbar(xs, pred, yerr=sd, fmt='none', ecolor=INK, elinewidth=0.9,
            capsize=3.2, capthick=0.9, zorder=4)
ax.plot(xs, meas, 'D', ms=5.2, mfc='white', mec=INK, mew=1.1, zorder=5,
        label="measured")
for x, p, e in zip(xs, pred, sd):
    ax.text(x, p + e + 0.45, f"{p:.2f}", va='bottom', ha='center', fontsize=6.8,
            color=INK, fontweight='bold')
ax.set_xticks(xs); ax.set_xticklabels(["$S_\\mathrm{P}$ (toxic)", "$R_\\mathrm{P}$"])
ax.set_xlim(-0.62, 1.62)
ax.set_ylim(0, 12.6)
ax.set_ylabel("$\\log k_i$ (hAChE)")
ax.set_title("(d) Leave-both-VX-out, out-of-sample", loc='left', fontweight='bold')
ax.grid(axis='x', visible=False)
ax.legend(loc='upper center', bbox_to_anchor=(0.52, 1.03), ncol=2, handlelength=1.2,
          columnspacing=1.1, borderaxespad=0.0)

fig.savefig("fig4_covalent_stereo_R2.png")
print("saved fig4_covalent_stereo_R2.png")
