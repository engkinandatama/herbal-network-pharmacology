# MD Simulation Bug Report: Missing Hydrogens in 264THM Ligand

## Issue Summary

- **Date Discovered:** 2026-01-13
- **Affected Simulation:** 264THM-PPARG (50ns MD simulation)
- **Status:** ⚠️ NEEDS RE-RUN

## Problem Description

The 264THM ligand topology is missing **all hydrogen atoms**. This causes:

1. **H-bond analysis underestimation** - Detected only 0.69 H-bonds (should be higher)
2. **Incorrect molecular geometry** - Carbon valences not properly filled
3. **Potential force field errors** - Missing non-bonded interactions

### Evidence

**264THM in md.gro file:**

- Atoms: C1-C24 (24 carbons) + O1-O5 (5 oxygens) = **29 atoms**
- **NO HYDROGEN ATOMS**

**Expected (from PubChem CID 10467773):**

- Molecular formula: **C₁₄H₁₂O₅**
- Should have: 14 C + 12 H + 5 O = **31 atoms**

**Luteolin (for comparison):**

- Has explicit H atoms (H7, H8, H9, H10) ✅
- H-bond analysis: 3.65 H-bonds (valid)

## Root Cause

Likely during ligand preparation:

1. SDF file downloaded from PubChem was 2D format (no explicit H)
2. OpenBabel conversion did not add hydrogens (`-h` flag missing)
3. ACPYPE topology generation did not detect missing H

## Impact on Results

| Analysis | 264THM-PPARG | Luteolin-PDE5A |
|----------|--------------|----------------|
| H-bond | ❌ Underestimated (0.69) | ✅ Valid (3.65) |
| RMSD | ⚠️ May be affected | ✅ Valid |
| Rg | ⚠️ May be affected | ✅ Valid |
| RMSF | ⚠️ May be affected | ✅ Valid |

## Fix Required

### Step 1: Regenerate ligand with explicit hydrogens

```bash
# Download 3D SDF from PubChem
wget "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/10467773/SDF?record_type=3d" -O 264THM_3d.sdf

# Convert with explicit hydrogens
obabel 264THM_3d.sdf -O 264THM.mol2 -h --gen3d

# Verify hydrogen count
grep -c "^.*H.*" 264THM.mol2
```

### Step 2: Regenerate ACPYPE topology

```bash
acpype -i 264THM.mol2 -c bcc -n 0
```

### Step 3: Re-run MD simulation

- Use corrected topology files
- Estimated time: ~8 hours on Kaggle GPU

## Files to Update

- [ ] `notebooks/kaggle/md_simulation_264THM_PPARG_gpu.ipynb` - Add hydrogen verification step
- [ ] Ligand preparation workflow - Add explicit check for H atoms
- [ ] `data/mahkota_dewa_dn/RESEARCH_SUMMARY.md` - Note limitation

## Priority

**HIGH** - This affects the validity of the 264THM-PPARG MD simulation results.

The Luteolin-PDE5A simulation is valid and can be used in the paper.

---

*Note: Created automatically during H-bond analysis debugging session*
