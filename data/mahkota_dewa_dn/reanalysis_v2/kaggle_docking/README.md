# Molecular Docking Dataset — Mahkota Dewa (Phaleria macrocarpa)

## Overview

Validated dataset for AutoDock Vina molecular docking.
Part of network pharmacology reanalysis.

## Contents

```
dataset/
├── ligands/           # 15 compound PDBQT files
│   ├── 264THM.pdbqt
│   ├── Apigenin.pdbqt
│   ├── ... (15 total)
├── controls/          # 4 positive control drugs
│   ├── Atorvastatin.pdbqt  → HMGCR
│   ├── Losartan.pdbqt      → AGTR1
│   ├── Pioglitazone.pdbqt  → PPARG
│   └── Sildenafil.pdbqt    → PDE5A
├── receptors/         # 5 protein targets (Meeko-prepared, with polar H)
│   ├── 1HWK.pdbqt  (HMGCR)
│   ├── 1TBF.pdbqt  (RELA)
│   ├── 3QXY.pdbqt  (AGTR1)
│   ├── 4YAY.pdbqt  (PPARG)
│   └── 6MS7.pdbqt  (PDE5A)
├── binding_sites.json  # Docking box coordinates
├── docking_config.json # Full configuration
└── dataset-metadata.json
```

## QC Status

- **29/29 validation checks PASSED**
- Ligands: CID verified vs PubChem, atom counts confirmed
- Receptors: Prepared with Meeko (gold standard), polar H present
- Binding sites: Centers verified within receptor bounding boxes

## Docking Jobs: 79 total

- 15 compounds × 5 targets = 75
- 4 controls × 1 target each = 4

## How to Use

1. Upload this `dataset/` folder as a Kaggle dataset
2. Create a new Kaggle notebook
3. Copy contents of `moldock_vina.py` into the notebook
4. Run all cells — results saved to `/kaggle/working/results/`
