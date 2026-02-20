# 🌿 Network Pharmacology of Mahkota Dewa (*Phaleria macrocarpa*) for Diabetic Nephropathy

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

## 📖 About

This repository contains the **data, scripts, and computational results** for the network pharmacology study of Mahkota Dewa (*Phaleria macrocarpa*) against Diabetic Nephropathy (DN). This repository serves as the **data availability** companion to the associated manuscript.

### Study Overview

An integrated network pharmacology approach combined with molecular docking and molecular dynamics (MD) simulation was used to investigate the therapeutic mechanisms of Mahkota Dewa against Diabetic Nephropathy.

### Key Findings

- **13 intersection targets** identified between Mahkota Dewa compounds and DN-related genes
- **PPARG** emerged as the top hub gene (degree=12)
- **Phalerin** (avg. -8.56 kcal/mol) and **Mangiferin** (avg. -8.49 kcal/mol) are the top lead compounds
- Both outperformed approved drugs (Pioglitazone, Losartan, Atorvastatin) at multiple targets
- **100 ns MD simulations** validated binding stability:
  - Mangiferin-RELA: RMSD 1.48 ± 0.15 Å (exceptionally stable)
  - Phalerin-AGTR1: RMSD 2.90 ± 0.58 Å (stable, typical GPCR dynamics)

> **📄 Full research summary**: [`data/mahkota_dewa_dn/RESEARCH_SUMMARY.md`](data/mahkota_dewa_dn/RESEARCH_SUMMARY.md)

---

## 📁 Repository Structure

```
herbal-network-pharmacology/
├── data/mahkota_dewa_dn/
│   ├── RESEARCH_SUMMARY.md              # Full results & interpretation
│   ├── raw/                             # Source compound & disease gene data
│   │   ├── compounds.csv
│   │   ├── disease_genes.csv
│   │   └── source_references.csv
│   ├── processed/                       # Intersection targets, predictions
│   │   ├── intersection_targets.csv
│   │   ├── intersection_genes.txt
│   │   └── predicted_targets.csv
│   ├── results/                         # Network analysis & enrichment
│   │   ├── hub_genes.csv
│   │   ├── network_statistics.json
│   │   ├── network_centralities.csv
│   │   ├── kegg_enrichment.csv
│   │   ├── go_bp_enrichment.csv
│   │   ├── admet_predictions.csv
│   │   └── docking/                     # Docking receptor/ligand files & figures
│   └── reanalysis_v2/                   # Molecular docking & MD simulation
│       ├── ligands/                     # Prepared ligand files (SDF + PDBQT)
│       ├── receptors/                   # Prepared receptor files (PDB + PDBQT)
│       ├── controls/                    # Control drug files
│       ├── kaggle_docking_output/       # AutoDock Vina docking results
│       ├── results_moldock/             # Docking results (backup)
│       ├── kaggle_md/                   # MD simulation scripts & output
│       │   ├── Mangiferin_RELA/         # MD raw data
│       │   ├── analysis_output/         # MD analysis timeseries
│       │   ├── combined_plots/          # Publication-ready plots
│       │   └── *.py / *.ipynb           # Simulation & analysis notebooks
│       ├── md-mmgbsa-mangiferin-rela.log
│       └── md-phalerin-agtr1-mmgbsa.log
├── notebooks/                           # Kaggle GPU notebooks
│   └── kaggle/                          # MD simulation & docking notebooks
├── src/                                 # Analysis scripts (Python)
├── config/                              # Project configuration files
└── references.md                        # Literature references
```

---

## 🔬 Methods Summary

| Step | Method | Tool/Database |
|------|--------|---------------|
| Compound Collection | Literature review | KNApSAcK, PubChem |
| Target Prediction | Structure-based | SwissTargetPrediction (prob ≥ 0.1) |
| Disease Genes | Database mining | OpenTargets, DisGeNET |
| PPI Network | STRING interactions | STRING (score ≥ 0.7) |
| Pathway Enrichment | Over-representation | KEGG, GO via Enrichr |
| ADMET Screening | Lipinski + Veber rules | RDKit |
| Molecular Docking | AutoDock Vina | exhaustiveness=32, num_modes=9 |
| MD Simulation | 100 ns production | OpenMM 8.2 (ff14SB + GAFF2, GPU) |
| Binding Free Energy | MM-GBSA | OpenMM GB model |

---

## 📊 Key Data Files

| File | Description |
|------|-------------|
| `raw/compounds.csv` | 26 bioactive compounds with SMILES |
| `processed/intersection_targets.csv` | 13 drug-disease common targets |
| `results/hub_genes.csv` | Network centrality metrics (degree, betweenness) |
| `results/kegg_enrichment.csv` | KEGG pathway enrichment results |
| `reanalysis_v2/kaggle_docking_output/results/docking_results.csv` | Docking scores for all compound-target pairs |
| `reanalysis_v2/kaggle_md/Mangiferin_RELA/md_summary.json` | MD simulation statistics |
| `reanalysis_v2/md-mmgbsa-mangiferin-rela.log` | MM-GBSA binding free energy (Mangiferin-RELA) |
| `reanalysis_v2/md-phalerin-agtr1-mmgbsa.log` | MM-GBSA binding free energy (Phalerin-AGTR1) |

---

## 🖥️ Reproducibility

### Local Analysis (Network Pharmacology)

The network pharmacology pipeline uses Python scripts in `src/`:

```bash
python -m venv venv && venv\Scripts\activate
pip install -r requirements.txt
python -m src.cli -c config/projects/mahkota_dewa_dn.yaml run-all
```

### GPU Computation (Docking & MD)

Molecular docking and MD simulation were performed on **Kaggle GPU (T4 x2)**. Notebooks are available in `notebooks/kaggle/`.

---

## 📝 License

MIT License — see [LICENSE](LICENSE) for details.

## 📖 Citation

If you use data from this repository, please cite:

```bibtex
@misc{nandatama2025mahkotadewa,
  author = {Nandatama, Engki},
  title = {Network Pharmacology Study of Mahkota Dewa (Phaleria macrocarpa) for Diabetic Nephropathy},
  year = {2025},
  publisher = {GitHub},
  url = {https://github.com/engkinandatama/herbal-network-pharmacology}
}
```

---

*Developed by [Engki Nandatama](https://github.com/engkinandatama)*
