# Reanalysis v2: Molecular Docking with Proper Ligand QC

**Date**: 2026-02-11  
**Reason**: Original docking used ligands without explicit hydrogens (264THM had 0 H-atoms).  
**Scope**: Re-dock all 13 drug-like compounds + 4 controls against 5 target proteins.

## Protocol

### Ligand Preparation

1. Download 3D SDF from PubChem (`/rest/pug/compound/cid/{CID}/SDF?record_type=3d`)
2. Convert to PDBQT: `obabel -isdf X.sdf -opdbqt -O X.pdbqt --gen3d -h -p 7.4`
3. QC: Verify H-atom count matches molecular formula

### Docking

- Tool: AutoDock Vina
- Exhaustiveness: 8
- Num modes: 9
- Binding sites: from `receptors/binding_sites.json`

### Compounds (13 drug-like)

| # | Name | PubChem CID |
|---|------|-------------|
| 1 | Gallic acid | 370 |
| 2 | Quercetin | 5280343 |
| 3 | Kaempferol | 5280863 |
| 4 | Phalerin | 189609 |
| 5 | 2,6,4'-trihydroxy-4-methoxybenzophenone | 10467773 |
| 6 | Swertianin | 5316833 |
| 7 | Apigenin | 5280443 |
| 8 | Luteolin | 5280445 |
| 9 | Caffeic acid | 689043 |
| 10 | Pinoresinol | 73399 |
| 11 | Lariciresinol | 101712 |
| 12 | Matairesinol | 119205 |
| 13 | Naringeninchalcone | 159654 |

### Controls (4)

| # | Name | PubChem CID | Target |
|---|------|-------------|--------|
| 1 | Pioglitazone | 4829 | PPARG (6MS7) |
| 2 | Atorvastatin | 60823 | HMGCR (1HWK) |
| 3 | Losartan | 3961 | AGTR1 (4YAY) |
| 4 | Sildenafil | 5212 | PDE5A (1TBF) |

### Target Proteins (5)

| Gene | PDB ID | Description |
|------|--------|-------------|
| PPARG | 6MS7 | PPAR-gamma |
| HMGCR | 1HWK | HMG-CoA reductase |
| AGTR1 | 4YAY | Angiotensin II receptor type 1 |
| PDE5A | 1TBF | Phosphodiesterase 5A |
| RELA | 3QXY | NF-κB p65 |
