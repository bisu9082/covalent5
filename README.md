# Benchmarking affinity-based docking for irreversible acetylcholinesterase inhibitors

Reproducibility package for the manuscript:

> **Benchmarking affinity-based docking for irreversible acetylcholinesterase
> inhibitors: a physics-informed reactivity framework for stereoselective potency
> prioritization**
> Jeongyun Kim, Jin Yoo, Doo-Hee Lee, Ku Kang.
> *ACS Omega* (under revision), Manuscript ID ao-2026-08355u.

This repository contains **all experimental anchors, quantum-chemical descriptors,
docking inputs/outputs, and analysis scripts** needed to regenerate every reported
statistic. The manuscript text and figure files are intentionally **not** included.

---

## What the study does

We benchmark whether conventional, affinity-based molecular docking can rank the
potency of **irreversible (covalent) organophosphorus inhibitors of acetylcholinesterase
(AChE)** across a stereochemically resolved set of 14 inhibitors (nerve agents and OP
pesticides), and we test whether **physics-informed reactivity descriptors** derived from
quantum chemistry recover a signal that docking scores miss. The headline finding is a
*negative benchmark*: docking scores do not reliably rank covalent potency, whereas a
mechanism-motivated descriptor gives a modest but reproducible correlation.

## Headline results (and where they come from)

| Quantity | Value | Regenerate with |
|---|---|---|
| Non-covalent docking vs pKi (Vina / GNINA-CNN) | Spearman rho ~ +0.12 / -0.40 | `docking/parse_docking.py` -> `docking/docking_results.csv` |
| Covalent docking vs pKi | rho ~ +0.23 | `docking/run_gnina_covalent.sh`, `data/covalent_docking.csv` |
| Physics-informed descriptor vs pKi | rho = -0.62 | `al_robustness.py`, `results/ml_robustness_out.md` |
| Docking box/seed sensitivity | Vina rho in [-0.03,+0.17]; CNN in [-0.33,-0.10] | `docking/run_sensitivity.sh` -> `docking/parse_sensitivity.py` |
| Active-site waters control | Vina +0.18 / CNN -0.13 | `docking/run_waters.sh` -> `docking/parse_waters.py` |
| Covalent-compatible (near-attack) pose scoring | ranking unchanged; 8/14 within 4.5 A of Ser203 Og | `docking/scardock_reanalysis.csv` |
| Conformer sensitivity of descriptors | qP range <= 0.007 e; gap <= 0.06 eV (VX 0.23 eV); dipole up to ~4.9 D | `conf_stability/parse_conf_stability.py` |
| Statistical robustness (permutation, bootstrap, leave-both-VX-out, AD) | see file | `al_robustness.py`, `data/ml_robustness.md` |

## Repository layout

```
data/                 Experimental anchors, descriptor tables, prepared ligand
                      structures (SDF/XYZ), and Ki kinetics provenance.
  anchor_ki.csv         14 inhibitors: Ki, pKi, stereo-configuration, source.
  descriptors.csv       DFT electronic descriptors (HOMO/LUMO/gap/dipole/qP).
  molecule_manifest.csv Canonical SMILES + stereo ids for every species.
  structures/           Prepared 3D geometries (44 SDF + 44 XYZ).
experimental/         Ellman assay Ki protocol.
code (repo root *.py) build_molecules.py, build_features.py,
                      generate_dft_inputs.py, extract_descriptors.py,
                      al_framework.py, al_rank.py, al_robustness.py,
                      render_molecules.py, plot_fig2.py, plot_fig4_full.py
dft/                  Per-molecule ORCA inputs (opt.inp/sp.inp), text outputs
                      (opt.out/sp.out), optimized geometries, job_manifest.csv.
                      Binary scratch (.gbw/.densities/...) is excluded by design.
conf_stability/       Conformer-sensitivity DFT jobs (inputs + text outputs)
                      and summary; gen_conf_dft.py, parse_conf_stability.py.
docking/              Receptor prep, ligand inputs, GNINA/Vina configs, pose
                      archives, logs, and all parse_* scripts + result CSVs.
qm_cluster/           Active-site cluster models (VX + Ser203/His447/Glu334 side
                      chains), the relaxed P...Ogamma scan, and ORCA templates.
figure/               Manuscript and SI figures with their generating scripts.

Analysis outputs are written into data/ (al_loo_pred.csv, al_predictions.csv,
ml_robustness_out.md). Earlier versions of this README named a results/
directory; it never existed and the reference has been removed.

verify_repo_layout.py checks this tree against the inventory above and exits
non-zero on any mismatch. Its output is deposited as repo_inventory.txt.
```

## Reproducing the analysis

1. **Environment.** Python >= 3.10 with the packages in `requirements.txt`
   (`rdkit` from conda-forge). External binaries: **ORCA 6.x**, **GNINA 1.x**,
   **Open Babel 3.x**.

2. **Ligand preparation.** `python build_molecules.py` enumerates stereoisomers
   (RDKit), embeds a single low-energy conformer (ETKDGv3, fixed seed; MMFF94),
   and writes `data/structures/`. Neutral (as-drawn) protonation is used.

3. **Quantum descriptors.** `python generate_dft_inputs.py` writes ORCA
   inputs (wB97X-D3/def2-SVP//def2-TZVP, CPCM(water)); run them, then
   `python extract_descriptors.py` -> `data/descriptors.csv`.

4. **Docking.** `docking/run_gnina_noncovalent.sh`, `run_gnina_covalent.sh`,
   `run_sensitivity.sh`, `run_waters.sh`; parse with the matching `parse_*.py`.

5. **Model + robustness.** `python al_robustness.py` reproduces the descriptor
   correlation, LOO cross-validation, label-permutation test, bootstrap
   intervals, leave-both-VX-enantiomers-out control, and applicability domain.

6. **Conformer sensitivity.** `cd conf_stability && python gen_conf_dft.py`,
   run `run_confs.sh`, then `python parse_conf_stability.py`.

## Data provenance & licensing

Experimental inhibition constants are compiled from the primary literature; the
per-inhibitor source is recorded in `data/anchor_ki.csv` and
`data/lit_anchor_kinetics.md`. Receptor structure: PDB **4EY7** (human AChE).

- **Code**: MIT License (see `LICENSE`).
- **Data**: CC-BY-4.0.

Full raw ORCA outputs and complete docking pose sets (binary/large) are archived
separately on Zenodo; see the Data and Code Availability statement in the paper.

## Citation

If you use this package, please cite the paper above. A `CITATION.cff` will be
added on acceptance with the final DOI.
