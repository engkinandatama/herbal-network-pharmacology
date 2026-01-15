# Standard Operating Procedure (SOP): Molecular Dynamics Preparation & QA
>
> [!IMPORTANT]
> **"Garbage In, Garbage Out"**. Validasi di tahap awal adalah fase paling kritikal. Kesalahan di sini (seperti atom H yang hilang) akan membuat simulasi berhari-hari menjadi sia-sia.

## 1. Post-Mortem Evaluasi: Kasus 264THM

**Masalah:** Ligand 264THM kehilangan semua atom Hidrogen (29 atoms vs 31 atoms expected).
**Penyebab:** Konversi dari format 2D (SDF/PubChem) ke 3D (Mol2) tanpa flag eksplisit untuk penambahan hidrogen, atau kegagalan software (OpenBabel) mendeteksi valensi kosong.
**Dampak:**

- Massa molekul salah.
- Geometri kacau (atom C yang harusnya tetrahedral jadi planar/undefined).
- Interaksi elektrostatik & VdW salah total.
- **H-bond analysis gagal** (karena tidak ada donor H).

---

## 2. Pre-Simulation Checklist (The "Golden Rules")

### Phase A: Ligand Preparation (CRITICAL)

1. **Source Validation:**
   - [ ] Download struktur 3D (SDF 3D) jika tersedia, jangan andalkan 2D.
   - [ ] Jika dari SMILES/2D, **WAJIB** generate 3D coordinates + Energy Minimization awal (e.g., pakai Avogadro/Chimera/OpenBabel `--gen3d`).

2. **Protonation State (pH 7.4):**
   - [ ] Cek pKa ligand. Pastikan protonasi benar pada pH fisiologis (7.4).
   - [ ] Tool: `obabel -p 7.4` atau server seperti SwissParam/Avogadro.

3. **Hydrogen Check (THE "264THM" RULE):**
   - [ ] **HITUNG JUMLAH ATOM.** Bandingkan formula kimia (misal C14H12O5) dengan file input.
   - [ ] Pastikan C-H, O-H, N-H lengkap.
   - [ ] Command check: `grep "H" ligand.mol2 | wc -l` (harus > 0 untuk organik).

4. **Topology Generation (ACPYPE/LigParGen):**
   - [ ] Pastikan `Net Charge` ligand masuk akal (bulat, biasanya 0, +1, -1).
   - [ ] Cek output log ACPYPE: apakah ada warning "Atom type not found" atau "Guessing parameter"?

### Phase B: Protein Preparation

1. **Clean PDB:**
   - [ ] Hapus atom non-standar yang tidak perlu (air kristal, ion buffer).
   - [ ] Pastikan tidak ada "Missing Atoms" atau "Missing Residues" di active site.
2. **Protonation (H++ / GROMACS `pdb2gmx`):**
   - [ ] Cek residu Histidine (HIE/HID/HIP) - status protonasi penting untuk binding.

### Phase C: System Building (GROMACS)

1. **Complex Assembly:**
   - [ ] Visualisasi (VMD/Chimera) setelah menggabungkan Protein + Ligand (`complex.gro`).
   - [ ] **CLASH CHECK:** Pastikan ligand tidak "menabrak" atom protein (jarak < 1.0 Å).

2. **Solvation & Neutralization:**
   - [ ] Box size cukup? (Min distance 1.0 nm ke dinding box).
   - [ ] Ion netralisasi ditambahkan? (Total charge sistem harus 0).

3. **Energy Minimization (EM):**
   - [ ] **Fmax < 1000 kJ/mol/nm?** (Syarat mutlak sebelum lanjut).
   - [ ] Jika gagal konvergen, Cek "Broken Molecule" atau steric clash parah.

### Phase D: Equilibration (NVT/NPT)

1. **Temperature & Pressure:**
   - [ ] NVT: Temperature stabil di 300K? (Fluktuasi wajar, tapi rata-rata harus tepat).
   - [ ] NPT: Density stabil di ~1000 kg/m^3 (untuk air)? Pressure stabil di 1 bar?

---

## 3. Automation Script Idea

Simpan script ini sebagai `check_ligand.sh` untuk dijalankan sebelum MD:

```bash
#!/bin/bash
file=$1
# Cek jumlah atom
natoms=$(grep -c "ATOM" $file)
nhydrogens=$(grep -i "H" $file | grep "ATOM" | wc -l)

echo "File: $file"
echo "Total Atoms Found: $natoms"
echo "Total Hydrogens Found: $nhydrogens"

if [ "$nhydrogens" -eq "0" ]; then
    echo "⚠️  CRITICAL WARNING: NO HYDROGENS DETECTED! DO NOT PROCEED!"
    exit 1
else
    echo "✅ Hydrogens present. Please verify count with chemical formula."
fi
```
