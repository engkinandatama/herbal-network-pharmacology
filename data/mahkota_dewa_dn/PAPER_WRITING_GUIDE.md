# 📝 Paper Writing Guide: Network Pharmacology Study of Mahkota Dewa for Diabetic Nephropathy

## 1. Judul Paper (Title Options)

### Opsi A (Comprehensive)
>
> **"Integrated Network Pharmacology, Molecular Docking, and Molecular Dynamics Simulation Reveals Multi-Target Therapeutic Mechanism of *Phaleria macrocarpa* Against Diabetic Nephropathy"**

### Opsi B (Focus on Lead Compound)
>
> **"264-Trihydroxy-4-methoxybenzophenone from *Phaleria macrocarpa* as a Potential PPARG Modulator for Diabetic Nephropathy: A Network Pharmacology and Molecular Dynamics Study"**

### Opsi C (Indonesian Focus)
>
> **"An *In Silico* Investigation of Mahkota Dewa (*Phaleria macrocarpa*) Bioactive Compounds Against Diabetic Nephropathy Through Network Pharmacology Approach"**

**Rekomendasi**: Opsi A (paling komprehensif, mencakup semua metode)

---

## 2. Rumusan Masalah (Research Questions)

### Background Problem
>
> Diabetic Nephropathy (DN) adalah komplikasi mikrovaskular diabetes yang menjadi penyebab utama End-Stage Renal Disease (ESRD). Terapi konvensional seperti ACE inhibitor dan ARB tidak sepenuhnya efektif menghentikan progresi penyakit. Mahkota Dewa (*Phaleria macrocarpa*) secara empiris digunakan untuk manajemen diabetes di Indonesia, namun mekanisme molekulernya belum diketahui.

### Research Questions

1. Apa saja senyawa bioaktif dalam Mahkota Dewa yang berpotensi terhadap target DN?
2. Target protein apa yang menjadi hub genes dalam interaksi senyawa-target?
3. Jalur signaling apa yang terlibat dalam mekanisme aksi Mahkota Dewa terhadap DN?
4. Apakah prediksi docking didukung oleh stabilitas binding pada simulasi MD?

### Objectives

1. Mengidentifikasi senyawa bioaktif dan target molekuler Mahkota Dewa
2. Menganalisis jaringan interaksi protein-protein (PPI network)
3. Melakukan enrichment analysis untuk pathway signaling
4. Memvalidasi interaksi top compounds melalui molecular docking dan MD simulation

---

## 3. Flow Penulisan Paper

```
┌─────────────────────────────────────────────────────────────┐
│                        ABSTRACT                              │
│  Background → Methods → Key Results → Conclusion            │
└─────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────┐
│                     1. INTRODUCTION                          │
│  1.1 DN Background (epidemiology, current therapy gaps)     │
│  1.2 Mahkota Dewa (traditional use, phytochemistry)         │
│  1.3 Network Pharmacology approach                          │
│  1.4 Study objectives                                       │
└─────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────┐
│                   2. MATERIALS AND METHODS                   │
│  2.1 Compound collection                                    │
│  2.2 Target prediction (SwissTargetPrediction)              │
│  2.3 Disease gene collection (OpenTargets, DisGeNET)        │
│  2.4 PPI network construction (STRING)                      │
│  2.5 Pathway enrichment (KEGG, GO)                          │
│  2.6 ADMET screening                                        │
│  2.7 Molecular docking (AutoDock Vina)                      │
│  2.8 MD simulation (GROMACS, 50 ns)                         │
└─────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────┐
│                       3. RESULTS                             │
│  3.1 Compound-target network                                │
│  3.2 Intersection targets (Venn diagram)                    │
│  3.3 PPI network and hub genes                              │
│  3.4 KEGG/GO enrichment                                     │
│  3.5 Docking results (affinity heatmap)                     │
│  3.6 MD simulation (RMSD, RMSF, Rg)                         │
└─────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────┐
│                      4. DISCUSSION                           │
│  4.1 Multi-target mechanism                                 │
│  4.2 Key pathways in DN                                     │
│  4.3 Lead compound profiles                                 │
│  4.4 Comparison with approved drugs                         │
│  4.5 Limitations                                            │
└─────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────┐
│                      5. CONCLUSION                           │
│  Summary of findings + therapeutic implications             │
└─────────────────────────────────────────────────────────────┘
```

---

## 4. Hasil yang Akan Ditampilkan

### 4.1 Figures (Gambar)

| No | Figure | Deskripsi | Status |
|----|--------|-----------|--------|
| 1 | Venn Diagram | Overlap drug targets ∩ disease genes | 🔴 To create |
| 2 | PPI Network | Hub genes highlighted | 🔴 To create |
| 3 | Compound-Target Network | Bipartite network | 🔴 To create |
| 4 | KEGG Enrichment | Top 10 pathways bar chart | 🔴 To create |
| 5 | Docking Heatmap | Affinity matrix compounds vs targets | 🔴 To create |
| 6 | 3D Docking Poses | Top compound-target interactions | ✅ Done (3 poses) |
| 7 | MD Analysis - 264THM-PPARG | RMSD, RMSF, Rg plots | ✅ Done |
| 8 | MD Analysis - Luteolin-PDE5A | RMSD, RMSF, Rg plots | ✅ Done |

### 4.2 Tables (Tabel)

| No | Table | Konten | Status |
|----|-------|--------|--------|
| 1 | Compound List | 26 compounds with SMILES, class | ✅ Data ready |
| 2 | Intersection Targets | 14 targets with annotations | ✅ Data ready |
| 3 | Hub Genes | Top 10 by degree/betweenness | ✅ Data ready |
| 4 | KEGG Pathways | Top 10 with p-values | ✅ Data ready |
| 5 | Docking Results | Compounds vs 5 targets | ✅ Data ready |
| 6 | MD Summary | RMSD, Rg, RMSF comparison | ✅ Data ready |

---

## 5. Cara Menganalisis dan Menulis Diskusi

### 5.1 Template Diskusi per Subseksi

#### A. Multi-Target Mechanism

```markdown
Temuan utama: Mahkota Dewa bekerja melalui mekanisme multi-target.

Struktur paragraf:
1. "The network analysis revealed that [X compounds] from 
   P. macrocarpa interact with [Y targets] related to DN..."
2. "Among these, PPARG emerged as the central hub gene with 
   the highest degree centrality (12), suggesting..."
3. "This multi-target characteristic aligns with the 
   polypharmacology concept in traditional medicine..."
```

#### B. Key Pathways

```markdown
Temuan: AGE-RAGE dan Calcium signaling paling signifikan.

Struktur paragraf:
1. "KEGG enrichment analysis identified AGE-RAGE signaling 
   as the most relevant pathway (p = 3.34e-05)..."
2. "In DN pathogenesis, AGE-RAGE activation leads to..."
3. "Our findings suggest that Mahkota Dewa compounds may 
   attenuate this pathway through SERPINE1, AGTR1, and RELA..."
```

#### C. Lead Compound Analysis

```markdown
Temuan: 264THM terbaik secara docking DAN MD.

Struktur paragraf:
1. "264-trihydroxy-4-methoxybenzophenone demonstrated the 
   highest average binding affinity (-8.46 kcal/mol)..."
2. "Notably, this compound outperformed pioglitazone 
   (-8.48 kcal/mol) at PPARG..."
3. "MD simulation confirmed stable binding throughout 50 ns 
   with RMSD of 0.224 ± 0.056 nm..."
```

#### D. Comparison with Approved Drugs

```markdown
Key selling point: Compounds lebih baik dari obat approved!

Template:
| Target | Control Drug | Control Score | Best MD Compound | Score |
|--------|--------------|---------------|------------------|-------|
| PPARG  | Pioglitazone | -8.48         | 264THM           | -9.53 ✓ |
| PDE5A  | Sildenafil   | -8.64         | 264THM           | -9.14 ✓ |
```

#### E. Limitations

```markdown
Wajib ada! Template:
1. "This study has several limitations. First, the predictions 
   are computational and require experimental validation..."
2. "Second, bioavailability and pharmacokinetics of the 
   lead compounds have not been assessed..."
3. "Third, the 50 ns MD simulation, while demonstrating 
   stability, may not capture rare conformational events..."
```

---

## 6. Novelty Points (Keunggulan Riset)

### Apa yang membuat riset ini berbeda?

1. **First comprehensive network pharmacology study** of Mahkota Dewa for DN
2. **MD simulation validation** (not just docking) - jarang dilakukan
3. **Multi-target analysis** - 14 validated targets
4. **Comparison with approved drugs** - 264THM outperforms pioglitazone
5. **Indonesian traditional medicine** focus - underexplored

### Kalimat Novelty untuk Abstract/Introduction
>
> "To our knowledge, this is the first study integrating network pharmacology with molecular dynamics simulation to elucidate the multi-target mechanism of *Phaleria macrocarpa* against diabetic nephropathy."

---

## 7. Checklist Sebelum Submit

### Data & Figures

- [ ] Venn diagram (intersection targets)
- [ ] PPI network visualization
- [ ] Compound-target network
- [ ] KEGG enrichment bar chart
- [ ] Docking heatmap
- [ ] 3D docking poses
- [ ] MD plots (264THM-PPARG, Luteolin-PDE5A)

### Writing

- [ ] Abstract (250-300 words)
- [ ] Introduction (3-4 paragraphs)
- [ ] Methods (detailed, reproducible)
- [ ] Results (figures + narrative)
- [ ] Discussion (interpret, compare, limitations)
- [ ] Conclusion (2-3 sentences)
- [ ] References (40-60 citations)

### Supplementary

- [ ] Full compound list with SMILES
- [ ] All docking results table
- [ ] KEGG/GO enrichment full tables

---

## 8. Target Jurnal

| Jurnal | Impact Factor | Kesesuaian | Review Time |
|--------|---------------|------------|-------------|
| **Frontiers in Pharmacology** | 5.6 | ⭐⭐⭐⭐⭐ | 2-3 months |
| **Journal of Ethnopharmacology** | 5.4 | ⭐⭐⭐⭐⭐ | 2-4 months |
| **Molecules (MDPI)** | 4.6 | ⭐⭐⭐⭐ | 2-3 weeks |
| **Phytomedicine** | 7.0 | ⭐⭐⭐⭐ | 3-4 months |
| **Computers in Biology and Medicine** | 7.7 | ⭐⭐⭐⭐ | 2-3 months |

**Recommendation**: Frontiers in Pharmacology (open access, good IF, accepts network pharmacology)

---

*Guide Created: January 10, 2026*
