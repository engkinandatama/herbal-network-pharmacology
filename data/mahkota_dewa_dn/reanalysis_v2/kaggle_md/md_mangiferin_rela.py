# %% [markdown]
# # MD Simulation: Mangiferin -> RELA (100 ns)
#
# **System:** Mangiferin (ligand from *Phaleria macrocarpa*) bound to RELA (NF-kB p65 subunit)
#
# **PDB:** 1TBF | **Docking Affinity:** -10.22 kcal/mol (better than control Losartan -8.99)
#
# **Protocol:** AMBER ff14SB (protein) + GAFF2/AM1-BCC (ligand) | TIP3P | 310K | 100 ns
#
# **Analysis:** RMSD, RMSF, Rg, H-bonds, MM-GBSA binding free energy

# %%
# ============================================================
# CELL 0: Install Dependencies (Kaggle GPU)
# ============================================================
# Kaggle has NO conda/mamba and uses system Python 3.12 with
# NumPy 2.0.2. We install Miniforge3 to get mamba, then install
# scientific packages via conda + pip.
#
# CRITICAL: Conda packages (rdkit, scipy) are compiled with
# newer GCC needing GLIBCXX_3.4.31 / CXXABI_1.3.15 which
# Kaggle's system libstdc++ lacks. These cannot be replaced
# mid-process. FIX: Install rdkit+scipy via pip instead
# (compiled against system libs).
# ============================================================
import subprocess, sys, os, glob, time, shutil

t_start = time.time()

MINIFORGE_DIR = "/tmp/miniforge"
MAMBA = f"{MINIFORGE_DIR}/bin/mamba"

PY_VER = f"{sys.version_info.major}.{sys.version_info.minor}"
print(f"Kernel Python: {PY_VER} ({sys.executable})")

# ?????????????????????????????????????????????????????????
# STEP 1: Install Miniforge3
# ?????????????????????????????????????????????????????????
print("=" * 60)
print("STEP 1: Installing Miniforge3...")
print("=" * 60)

if not os.path.exists(MAMBA):
    url = "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
    print("  Downloading...")
    r = subprocess.run(f"wget -q {url} -O /tmp/miniforge.sh",
                       shell=True, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"Cannot download Miniforge3: {r.stderr}")

    print(f"  Installing to {MINIFORGE_DIR}...")
    r = subprocess.run(f"bash /tmp/miniforge.sh -b -p {MINIFORGE_DIR}",
                       shell=True, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"Install failed: {r.stderr[-500:]}")

    os.remove("/tmp/miniforge.sh")
    print(f"  [OK] Miniforge3 installed ({time.time()-t_start:.0f}s)")
else:
    print(f"  [OK] Already present")

# ?????????????????????????????????????????????????????????
# STEP 2: Install conda packages + cleanup
# ?????????????????????????????????????????????????????????
print("\n" + "=" * 60)
print("STEP 2: Mamba install...")
print("=" * 60)

# rdkit and scipy are NOT installed here -- they come from pip
# to avoid GLIBCXX/CXXABI conflicts with system libstdc++.
packages = f"python={PY_VER} numpy openmm cudatoolkit openmmforcefields openff-toolkit ambertools pdbfixer mdtraj parmed"
print(f"  -> {packages}")
t_install = time.time()
r = subprocess.run(
    f"{MAMBA} install -y -c conda-forge {packages}",
    shell=True, capture_output=True, text=True, timeout=600
)
if r.stdout:
    for line in r.stdout.strip().split('\n')[-10:]:
        print(f"    {line}")
if r.returncode != 0:
    print(f"\n  [WARN] mamba failed (exit {r.returncode})")
    if r.stderr:
        for line in r.stderr.strip().split('\n')[-10:]:
            print(f"    STDERR: {line}")
else:
    print(f"\n  [OK] Packages installed ({time.time()-t_install:.0f}s)")

# Show conda numpy version
r_np = subprocess.run(
    f"{MINIFORGE_DIR}/bin/python -c 'import numpy; print(numpy.__version__)'",
    shell=True, capture_output=True, text=True
)
print(f"  Conda numpy: {r_np.stdout.strip() if r_np.returncode == 0 else '?'}")

# Remove conda rdkit+scipy (auto-installed as openff-toolkit deps).
# These need GLIBCXX_3.4.31 / CXXABI_1.3.15 -> use pip versions.
print("  Cleaning conda rdkit+scipy (will use pip versions)...")
subprocess.run(f"{MAMBA} remove -y --force rdkit rdkit-core scipy 2>/dev/null",
               shell=True, capture_output=True, text=True)
removed = 0
for pattern in [f"{MINIFORGE_DIR}/lib/libRDKit*",
                f"{MINIFORGE_DIR}/lib/python*/site-packages/rdkit*",
                f"{MINIFORGE_DIR}/lib/python*/site-packages/scipy*"]:
    for f in glob.glob(pattern):
        if os.path.isdir(f):
            shutil.rmtree(f, ignore_errors=True)
        else:
            os.remove(f)
        removed += 1
print(f"  [OK] Removed conda rdkit+scipy ({removed} files)")

# ?????????????????????????????????????????????????????????
# STEP 2.5: Install pip packages (BEFORE removing system numpy)
# ?????????????????????????????????????????????????????????
print("\n" + "=" * 60)
print("STEP 2.5: Pip install (rdkit, scipy, meeko, gemmi, vina)...")
print("=" * 60)
for pkg in ["rdkit", "scipy", "meeko", "gemmi", "vina"]:
    r = subprocess.run([sys.executable, "-m", "pip", "install", "-q", pkg],
                       capture_output=True, text=True, timeout=300)
    if r.returncode == 0:
        print(f"  [OK] {pkg}")
    else:
        print(f"  [FAIL] {pkg}: {r.stderr[-300:]}")

# ?????????????????????????????????????????????????????????
# STEP 3: Remove system numpy & patch paths
# ?????????????????????????????????????????????????????????
print("\n" + "=" * 60)
print("STEP 3: Remove system numpy & patch paths...")
print("=" * 60)

# 3a: Remove system numpy (conflicts with conda numpy)
print("  Removing system numpy...")
r = subprocess.run([sys.executable, "-m", "pip", "uninstall", "-y", "numpy"],
                   capture_output=True, text=True)
if r.returncode == 0:
    print("  [OK] System numpy removed")
else:
    for np_dir in glob.glob("/usr/local/lib/python3.12/dist-packages/numpy*"):
        shutil.rmtree(np_dir, ignore_errors=True)
    print("  [OK] System numpy force-removed")

# 3b: Flush cached numpy modules
mods_to_flush = [k for k in sys.modules if 'numpy' in k.lower()]
for mod in mods_to_flush:
    del sys.modules[mod]
print(f"  [OK] Flushed {len(mods_to_flush)} cached modules")

# 3c: Add miniforge to sys.path (FIRST position)
site_dirs = glob.glob(f"{MINIFORGE_DIR}/lib/python{PY_VER}/site-packages")
site_dirs += glob.glob(f"{MINIFORGE_DIR}/lib/python*/site-packages")
for sp in sorted(set(site_dirs)):
    if os.path.isdir(sp) and sp not in sys.path:
        sys.path.insert(0, sp)
        print(f"  [OK] sys.path += {sp}")

# 3d: Environment variables
os.environ["PATH"] = f"{MINIFORGE_DIR}/bin:" + os.environ.get("PATH", "")
os.environ["LD_LIBRARY_PATH"] = f"{MINIFORGE_DIR}/lib:" + os.environ.get("LD_LIBRARY_PATH", "")
os.environ["AMBERHOME"] = MINIFORGE_DIR
print(f"  [OK] PATH, LD_LIBRARY_PATH, AMBERHOME set")

# ?????????????????????????????????????????????????????????
# STEP 4: Verify critical imports
# ?????????????????????????????????????????????????????????
print("\n" + "=" * 60)
print("STEP 4: Verifying critical imports...")
print("=" * 60)

# Verify numpy
try:
    import numpy as np
    print(f"  [OK] numpy {np.__version__}")
except Exception as e:
    print(f"  [FAIL] numpy FAILED: {e}")

# Verify OpenMM + CUDA
try:
    import openmm
    print(f"  [OK] OpenMM {openmm.__version__}")
    platforms = [openmm.Platform.getPlatform(i).getName()
                 for i in range(openmm.Platform.getNumPlatforms())]
    print(f"    Platforms: {platforms}")
    if "CUDA" in platforms:
        print(f"  [OK] CUDA platform ready")
    else:
        print(f"  [WARN] CUDA not found -- will use CPU (much slower!)")
except ImportError as e:
    print(f"  [FAIL] OpenMM import FAILED: {e}")
    raise

# All other packages
all_ok = True
for mod_name in ["pdbfixer", "mdtraj", "rdkit", "parmed",
                 "openff.toolkit", "openmmforcefields", "scipy",
                 "meeko", "vina"]:
    try:
        __import__(mod_name)
        print(f"  [OK] {mod_name}")
    except ImportError as e:
        print(f"  [FAIL] {mod_name}: {e}")
        all_ok = False

elapsed = time.time() - t_start
print(f"\n{'=' * 60}")
print(f"Setup complete in {elapsed:.0f}s")
if all_ok:
    print("[OK] All dependencies verified!")
else:
    print("[WARN] Some packages missing -- check logs above")
print(f"{'=' * 60}")

# %%
# ============================================================
# CELL 1: Configuration
# ============================================================
import json
from pathlib import Path

# --- System Definition ---
SYSTEM_NAME = "Mangiferin_RELA"
RECEPTOR_PDB_ID = "1TBF"
LIGAND_NAME = "Mangiferin"

# --- Paths (Kaggle) ---
# Check for dedicated dataset first
DATASET_DIR = Path("/kaggle/input/md-mangiferin-rela-input")
if not DATASET_DIR.exists():
    DATASET_DIR = Path("/kaggle/input/mahkota-dewa-docking")
OUTPUT_DIR = Path("/kaggle/working/md_results") / SYSTEM_NAME
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Input files from dataset
# PDB: Full protein structure WITH hydrogens (same file that was converted to PDBQT for docking)
# SDF: Ligand with full bond-order info (needed for AM1-BCC charge assignment)
# PDBQT: For re-docking only (AutoDock Vina format)
# NOTE: Support both flat layout (all files in root) and subdirectory layout (ligands/, receptors/)
def _find_file(dataset_dir, filename, subdirs=None):
    """Find file in dataset, trying root first, then subdirectories."""
    # Try root first (flat layout)
    p = dataset_dir / filename
    if p.exists():
        return p
    # Try subdirectories
    if subdirs:
        for sub in subdirs:
            p = dataset_dir / sub / filename
            if p.exists():
                return p
    # Default to root path (will fail with clear error)
    return dataset_dir / filename

RECEPTOR_PDB = _find_file(DATASET_DIR, f"{RECEPTOR_PDB_ID}_withH.pdb", ["receptors"])
RECEPTOR_PDBQT = _find_file(DATASET_DIR, f"{RECEPTOR_PDB_ID}.pdbqt", ["receptors"])
LIGAND_SDF = _find_file(DATASET_DIR, f"{LIGAND_NAME}.sdf", ["ligands"])
LIGAND_PDBQT = _find_file(DATASET_DIR, f"{LIGAND_NAME}.pdbqt", ["ligands"])
# --- Binding Site (hardcoded from docking_config.json for 1TBF/RELA) ---
# Hardcoded to avoid FileNotFoundError if dataset version is still processing.
BINDING_SITE = {
    "center": (28.73, 30.17, 65.13),
    "size": (20.71, 20.0, 20.0),
}

# Diagnostic: list dataset root contents
print("Dataset contents:")
if DATASET_DIR.exists():
    for item in sorted(DATASET_DIR.iterdir()):
        tag = "DIR" if item.is_dir() else f"{item.stat().st_size:,}B"
        print(f"  [{tag}] {item.name}")
else:
    print(f"  [WARN] Dataset not found at {DATASET_DIR}")

# --- MD Parameters ---
MD_PARAMS = {
    "temperature_K": 310.15,        # 37?C physiological
    "pressure_atm": 1.0,
    "ionic_strength_M": 0.15,       # physiological NaCl
    "timestep_fs": 2.0,             # with SHAKE/HBonds constraints
    "nonbonded_cutoff_nm": 1.0,
    "padding_nm": 1.2,              # solvation box padding
    # Equilibration
    "minimization_steps": 5000,
    "nvt_steps": 50000,             # 100 ps NVT heating
    "npt_steps": 500000,            # 1 ns NPT equilibration
    # Production
    "production_ns": 100,
    "save_interval_ps": 10,         # save frame every 10 ps
    "log_interval_ps": 10,          # log energy every 10 ps
    "checkpoint_interval_ns": 10,   # checkpoint every 10 ns
}

# Derived values
_dt = MD_PARAMS["timestep_fs"]
MD_PARAMS["production_steps"] = int(MD_PARAMS["production_ns"] * 1e6 / _dt)
MD_PARAMS["save_interval_steps"] = int(MD_PARAMS["save_interval_ps"] * 1000 / _dt)
MD_PARAMS["log_interval_steps"] = int(MD_PARAMS["log_interval_ps"] * 1000 / _dt)
MD_PARAMS["checkpoint_interval_steps"] = int(MD_PARAMS["checkpoint_interval_ns"] * 1e6 / _dt)

print("=" * 60)
print(f"System: {SYSTEM_NAME}")
print(f"Receptor PDB: {RECEPTOR_PDB}")
print(f"Ligand SDF: {LIGAND_SDF}")
print(f"Binding site center: {BINDING_SITE['center']}")
print(f"Production: {MD_PARAMS['production_ns']} ns = {MD_PARAMS['production_steps']:,} steps")
print(f"Save every: {MD_PARAMS['save_interval_ps']} ps ({MD_PARAMS['save_interval_steps']} steps)")
print(f"Checkpoint every: {MD_PARAMS['checkpoint_interval_ns']} ns")
print(f"Expected frames: {MD_PARAMS['production_steps'] // MD_PARAMS['save_interval_steps']:,}")
print(f"Output: {OUTPUT_DIR}")
print("=" * 60)

# Verify ALL input files exist before proceeding
print("\nVerifying input files:")
all_ok = True
for f in [RECEPTOR_PDB, RECEPTOR_PDBQT, LIGAND_SDF, LIGAND_PDBQT]:
    if f.exists():
        print(f"  [OK] {f.name} ({f.stat().st_size:,} bytes)")
    else:
        print(f"  [FAIL] MISSING: {f}")
        all_ok = False

assert all_ok, "Some input files are missing! Check dataset upload."

# %%
# ============================================================
# CELL 2: Prepare Protein (PDBFixer)
# ============================================================
# The dataset contains 1TBF_withH.pdb, preprocessed by `reduce`
# for docking. For MD we need to clean it because:
#   1. `reduce` added H to BOTH alternate conformations, then
#      stripped alt-loc indicators -> duplicate heavy atoms
#   2. The H placement was for docking, not MD (wrong pH model)
#
# Strategy (gold standard for MD prep):
#   (a) Strip ALL hydrogens (PDBFixer will re-add at pH 7.4)
#   (b) Remove duplicate heavy atoms (keep first occurrence = conformer A)
#   (c) PDBFixer: fix missing residues/atoms -> add H at pH 7.4
#
# Binding site geometry is preserved because heavy atom coords
# are from the SAME crystal structure (1TBF), and the ligand
# pose comes from re-docking (Cell 4) using the clean PDBQT.

from pdbfixer import PDBFixer
from openmm.app import PDBFile
import openmm.unit as unit

print("Preparing receptor from dataset PDB...")

# Step 1: Strip H + dedup heavy atoms
clean_pdb = OUTPUT_DIR / "receptor_clean.pdb"
n_h_removed = 0
n_dupes = 0
seen_atoms = set()
with open(RECEPTOR_PDB) as fin, open(clean_pdb, "w") as fout:
    for line in fin:
        if line.startswith(("ATOM", "HETATM")):
            # (a) Skip hydrogens -- element is in columns 77-78
            element = line[76:78].strip() if len(line) > 77 else ''
            if element == 'H':
                n_h_removed += 1
                continue

            # (b) Skip alt-loc conformers other than A/blank
            alt_loc = line[16]
            if alt_loc not in (' ', '', 'A'):
                n_dupes += 1
                continue

            # (c) Dedup by atom identity (name + residue + chain + seq)
            atom_key = (line[12:16], line[17:20], line[21], line[22:27])
            if atom_key in seen_atoms:
                n_dupes += 1
                continue
            seen_atoms.add(atom_key)

            # Clear alt-loc indicator
            line = line[:16] + ' ' + line[17:]
        fout.write(line)

print(f"  Stripped {n_h_removed} hydrogens (will re-add at pH 7.4)")
print(f"  Removed {n_dupes} duplicate/alt-conformation atoms")
print(f"  Heavy atoms retained: {len(seen_atoms)}")

# Step 2: PDBFixer
fixer = PDBFixer(str(clean_pdb))

# Remove heterogens (crystal waters, co-crystallized ligands)
fixer.removeHeterogens(keepWater=False)

n_chains = len(list(fixer.topology.chains()))
n_res = len(list(fixer.topology.residues()))
n_atoms_before = len(list(fixer.topology.atoms()))
print(f"  Chains: {n_chains}, Residues: {n_res}, Atoms: {n_atoms_before}")

# Find and fix missing residues
fixer.findMissingResidues()
missing_res = dict(fixer.missingResidues)

if missing_res:
    # Remove terminal missing residues (disordered, not needed for MD)
    keys_to_remove = []
    for key in missing_res:
        chain_idx, res_idx = key
        chain = list(fixer.topology.chains())[chain_idx]
        n_res_chain = len(list(chain.residues()))
        if res_idx == 0 or res_idx >= n_res_chain:
            keys_to_remove.append(key)
    for key in keys_to_remove:
        del fixer.missingResidues[key]
    n_internal = len(fixer.missingResidues)
    print(f"  Missing residues: {len(missing_res)} total, {n_internal} internal (adding these)")

# Find and add missing atoms
fixer.findMissingAtoms()
n_missing_atoms = sum(len(v) for v in fixer.missingAtoms.values())
n_missing_terminals = sum(len(v) for v in fixer.missingTerminals.values())
print(f"  Missing atoms: {n_missing_atoms} in {len(fixer.missingAtoms)} residues")
print(f"  Missing terminals: {n_missing_terminals}")

fixer.addMissingAtoms()
print("  [OK] Missing atoms added")

# Add ALL hydrogens at pH 7.4 (physiological)
fixer.addMissingHydrogens(7.4)

n_atoms_final = len(list(fixer.topology.atoms()))
n_res_final = len(list(fixer.topology.residues()))
print(f"  Final: {n_atoms_final} atoms, {n_res_final} residues")

# Save fixed protein
protein_pdb = OUTPUT_DIR / "protein_fixed.pdb"
with open(protein_pdb, "w") as f:
    PDBFile.writeFile(fixer.topology, fixer.positions, f)
print(f"\n[OK] Protein saved: {protein_pdb}")

# %%
# ============================================================
# CELL 3: Prepare Ligand (from SDF)
# ============================================================
# WHY SDF (not PDBQT)?
#   SDF contains full bond-order info (single/double/aromatic bonds).
#   This is REQUIRED for correct AM1-BCC charge calculation and
#   GAFF2 atom typing. PDBQT loses this info.
#   The SDF coordinates come from PubChem 3D conformer (same source
#   used to generate the PDBQT for docking).

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

print(f"Loading ligand from {LIGAND_SDF}...")
supplier = Chem.SDMolSupplier(str(LIGAND_SDF), removeHs=False)
mol = next(supplier)
assert mol is not None, "Failed to load ligand SDF!"

# Check if hydrogens are present
n_total = mol.GetNumAtoms()
mol_noH = Chem.RemoveAllHs(mol)
n_heavy = mol_noH.GetNumAtoms()
has_H = (n_total > n_heavy)

if not has_H:
    print("  Adding hydrogens (not present in SDF)...")
    mol = Chem.AddHs(mol, addCoords=True)
    # Re-generate 3D coords with H
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
    n_total = mol.GetNumAtoms()

formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
mw = Descriptors.ExactMolWt(mol)

print(f"  Name: {LIGAND_NAME}")
print(f"  Formula: {formula}")
print(f"  MW: {mw:.2f} Da")
print(f"  Heavy atoms: {n_heavy}")
print(f"  Total atoms (with H): {n_total}")
print(f"  Has 3D coords: {mol.GetNumConformers() > 0}")

# Save prepared ligand
ligand_sdf_prep = OUTPUT_DIR / f"{LIGAND_NAME}_prepared.sdf"
writer = Chem.SDWriter(str(ligand_sdf_prep))
writer.write(mol)
writer.close()
print(f"\n[OK] Ligand saved: {ligand_sdf_prep}")

# %%
# ============================================================
# CELL 4: Re-dock Ligand to Get Binding Pose
# ============================================================
# We re-dock because our Kaggle docking output only has affinity
# scores in CSV, not the actual 3D pose coordinates needed for MD.
# Using the SAME PDBQT files and SAME binding site as original docking
# ensures the pose is consistent.

from vina import Vina

print(f"Re-docking {LIGAND_NAME} into {RECEPTOR_PDB_ID}...")
print(f"  Receptor PDBQT: {RECEPTOR_PDBQT}")
print(f"  Ligand PDBQT: {LIGAND_PDBQT}")
print(f"  Center: {BINDING_SITE['center']}")
print(f"  Size: {BINDING_SITE['size']}")

v = Vina(sf_name="vina")
v.set_receptor(str(RECEPTOR_PDBQT))
v.set_ligand_from_file(str(LIGAND_PDBQT))

cx, cy, cz = BINDING_SITE["center"]
sx, sy, sz = BINDING_SITE["size"]
v.compute_vina_maps(center=[cx, cy, cz], box_size=[sx, sy, sz])

v.dock(exhaustiveness=32, n_poses=10)

# Get best pose
energies = v.energies()
best_affinity = energies[0][0]
print(f"\n  Best affinity: {best_affinity:.2f} kcal/mol")
print(f"  All poses: {[f'{e[0]:.2f}' for e in energies]}")

# Save best pose
docked_pdbqt = OUTPUT_DIR / f"{LIGAND_NAME}_docked.pdbqt"
v.write_poses(str(docked_pdbqt), n_poses=1, overwrite=True)

# Convert docked PDBQT -> RDKit Mol for SDF/PDB output
# Try multiple methods since Meeko can fail silently (returns None)

docked_mol = None

# --- Method 1: Meeko ---
try:
    from meeko import PDBQTMolecule, RDKitMolCreate
    pdbqt_mol = PDBQTMolecule.from_file(str(docked_pdbqt), skip_typing=True)
    rdkit_mols = RDKitMolCreate.from_pdbqt_mol(pdbqt_mol)
    # Filter out None entries
    valid_mols = [m for m in rdkit_mols if m is not None]
    if valid_mols:
        docked_mol = valid_mols[0]
        print("  Converted via Meeko [OK]")
    else:
        print("  Meeko returned None, trying fallback...")
except Exception as e:
    print(f"  Meeko failed: {e}, trying fallback...")

# --- Method 2: obabel ---
if docked_mol is None:
    try:
        import subprocess as sp
        docked_sdf_tmp = OUTPUT_DIR / f"{LIGAND_NAME}_docked_tmp.sdf"
        r = sp.run(
            ["obabel", str(docked_pdbqt), "-O", str(docked_sdf_tmp), "-h"],
            capture_output=True, text=True, timeout=60
        )
        if r.returncode == 0 and docked_sdf_tmp.exists():
            supplier = Chem.SDMolSupplier(str(docked_sdf_tmp), removeHs=False)
            for m in supplier:
                if m is not None:
                    docked_mol = m
                    print("  Converted via obabel [OK]")
                    break
        if docked_mol is None:
            print("  obabel conversion failed, trying manual method...")
    except Exception as e:
        print(f"  obabel failed: {e}, trying manual method...")

# --- Method 3: Map docked coords onto original mol ---
if docked_mol is None:
    print("  Using manual PDBQT coordinate extraction...")
    import numpy as np

    # Parse PDBQT heavy atom coordinates
    docked_coords = []
    with open(docked_pdbqt) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                element = line[77:79].strip()
                if element != 'H' and element != 'HD':
                    docked_coords.append((x, y, z))

    # Get prepared mol from Cell 3 (has correct bond orders)
    docked_mol = Chem.RWMol(mol)  # Copy from Cell 3's prepared mol
    docked_mol = Chem.RemoveAllHs(docked_mol)  # Work with heavy atoms

    n_heavy = docked_mol.GetNumAtoms()
    if len(docked_coords) == n_heavy:
        conf = Chem.Conformer(n_heavy)
        for i, (x, y, z) in enumerate(docked_coords):
            conf.SetAtomPosition(i, (x, y, z))
        docked_mol.RemoveAllConformers()
        docked_mol.AddConformer(conf, assignId=True)

        # Re-add hydrogens with 3D coords
        docked_mol = Chem.AddHs(docked_mol, addCoords=True)
        AllChem.MMFFOptimizeMolecule(docked_mol, maxIters=200,
                                     nonBondedThresh=100.0)
        print(f"  Mapped {n_heavy} heavy atom coords [OK]")
    else:
        raise RuntimeError(
            f"Atom count mismatch: PDBQT has {len(docked_coords)} "
            f"heavy atoms, mol has {n_heavy}"
        )

assert docked_mol is not None, "All PDBQT conversion methods failed!"

# Save as SDF (for OpenFF parameterization)
docked_sdf = OUTPUT_DIR / f"{LIGAND_NAME}_docked.sdf"
writer = Chem.SDWriter(str(docked_sdf))
writer.write(docked_mol)
writer.close()

# Save as PDB (for merging with protein in OpenMM)
docked_pdb = OUTPUT_DIR / f"{LIGAND_NAME}_docked.pdb"
Chem.MolToPDBFile(docked_mol, str(docked_pdb))

print(f"\n[OK] Docked ligand saved:")
print(f"  SDF: {docked_sdf}")
print(f"  PDB: {docked_pdb}")

# %%
# ============================================================
# CELL 5: Build Solvated Complex
# ============================================================
from openmm.app import (
    ForceField, Modeller, PDBFile, PME, HBonds, NoCutoff,
    Simulation, StateDataReporter, DCDReporter
)
from openmm import (
    LangevinMiddleIntegrator, MonteCarloBarostat, Platform
)
from openmmforcefields.generators import GAFFTemplateGenerator
from openff.toolkit import Molecule as OFFMolecule
import openmm.unit as unit

print("Building protein-ligand complex...")

# --- Load protein ---
protein_pdb_obj = PDBFile(str(protein_pdb))
n_prot_atoms = protein_pdb_obj.topology.getNumAtoms()
print(f"  Protein atoms: {n_prot_atoms}")

# --- Load docked ligand for OpenFF parameterization ---
off_mol = OFFMolecule.from_file(str(docked_sdf))
n_lig_atoms = off_mol.n_atoms
print(f"  Ligand atoms: {n_lig_atoms}")

# --- Setup ForceField with GAFF2 for ligand ---
ff = ForceField("amber/ff14SB.xml", "amber/tip3p_standard.xml")

gaff_gen = GAFFTemplateGenerator(
    molecules=[off_mol],
    forcefield="gaff-2.11"
)
ff.registerTemplateGenerator(gaff_gen.generator)
print("  Registered GAFF2 template generator (AM1-BCC charges)")

# --- Merge protein + ligand ---
modeller = Modeller(protein_pdb_obj.topology, protein_pdb_obj.positions)

lig_pdb = PDBFile(str(docked_pdb))
modeller.add(lig_pdb.topology, lig_pdb.positions)
n_complex = modeller.topology.getNumAtoms()
print(f"  Complex atoms (protein+ligand): {n_complex}")

# --- Solvate with TIP3P + NaCl ---
modeller.addSolvent(
    ff,
    model="tip3p",
    padding=MD_PARAMS["padding_nm"] * unit.nanometers,
    ionicStrength=MD_PARAMS["ionic_strength_M"] * unit.molar,
    positiveIon="Na+",
    negativeIon="Cl-",
)

n_atoms_total = modeller.topology.getNumAtoms()
n_residues_total = modeller.topology.getNumResidues()

# Count components (OpenMM uses HOH, NA, CL as residue names)
n_water = sum(1 for r in modeller.topology.residues() if r.name == "HOH")
n_na = sum(1 for r in modeller.topology.residues() if r.name in ["NA", "Na+"])
n_cl = sum(1 for r in modeller.topology.residues() if r.name in ["CL", "Cl-"])

print(f"\n  Solvated system:")
print(f"    Total atoms:     {n_atoms_total:,}")
print(f"    Total residues:  {n_residues_total:,}")
print(f"    Protein atoms:   {n_prot_atoms:,}")
print(f"    Ligand atoms:    {n_lig_atoms}")
print(f"    Water molecules: {n_water:,}")
print(f"    Na+ ions:        {n_na}")
print(f"    Cl- ions:        {n_cl}")

# Save solvated complex PDB (used as topology reference for trajectory loading)
solvated_pdb = OUTPUT_DIR / "complex_solvated.pdb"
with open(solvated_pdb, "w") as f:
    PDBFile.writeFile(modeller.topology, modeller.positions, f)
print(f"\n[OK] Solvated complex saved: {solvated_pdb}")

# %%
# ============================================================
# CELL 6: Energy Minimization + Equilibration
# ============================================================
import time
import numpy as np

print("=" * 60)
print("ENERGY MINIMIZATION + EQUILIBRATION")
print("=" * 60)

# --- Create System ---
system = ff.createSystem(
    modeller.topology,
    nonbondedMethod=PME,
    nonbondedCutoff=MD_PARAMS["nonbonded_cutoff_nm"] * unit.nanometers,
    constraints=HBonds,
    rigidWater=True,
    hydrogenMass=1.5 * unit.amu,  # HMR for stability with 2fs timestep
)

# --- Select Platform ---
try:
    platform = Platform.getPlatformByName("CUDA")
    properties = {"DeviceIndex": "0", "Precision": "mixed"}
    print(f"  Using CUDA platform (mixed precision)")
except Exception:
    platform = Platform.getPlatformByName("CPU")
    properties = {}
    print(f"  [WARN] Using CPU platform (CUDA not available)")

# --- Integrator ---
integrator = LangevinMiddleIntegrator(
    MD_PARAMS["temperature_K"] * unit.kelvin,
    1.0 / unit.picoseconds,
    MD_PARAMS["timestep_fs"] * unit.femtoseconds,
)

# --- Create Simulation ---
simulation = Simulation(modeller.topology, system, integrator, platform, properties)
simulation.context.setPositions(modeller.positions)

# --- [1/3] Energy Minimization ---
print(f"\n[1/3] Energy Minimization ({MD_PARAMS['minimization_steps']} steps)...")
state = simulation.context.getState(getEnergy=True)
e_before = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
print(f"  Energy before: {e_before:.1f} kJ/mol")

t0 = time.time()
simulation.minimizeEnergy(maxIterations=MD_PARAMS["minimization_steps"])

state = simulation.context.getState(getEnergy=True, getPositions=True)
e_after = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
print(f"  Energy after:  {e_after:.1f} kJ/mol")
print(f"  DeltaE = {e_after - e_before:.1f} kJ/mol ({time.time()-t0:.1f}s)")

# Save minimized
min_pdb = OUTPUT_DIR / "complex_minimized.pdb"
with open(min_pdb, "w") as f:
    PDBFile.writeFile(simulation.topology, state.getPositions(), f)
print(f"  Saved: {min_pdb}")

# --- [2/3] NVT Heating (100 ps) ---
print(f"\n[2/3] NVT Heating to {MD_PARAMS['temperature_K']}K ({MD_PARAMS['nvt_steps']} steps = 100 ps)...")
simulation.context.setVelocitiesToTemperature(50 * unit.kelvin)

equil_log = OUTPUT_DIR / "equilibration.log"
simulation.reporters.append(
    StateDataReporter(
        str(equil_log), 5000,
        step=True, time=True, potentialEnergy=True,
        temperature=True, speed=True,
    )
)

t0 = time.time()
simulation.step(MD_PARAMS["nvt_steps"])
print(f"  NVT done ({time.time()-t0:.1f}s)")

# --- [3/3] NPT Equilibration (1 ns) ---
print(f"\n[3/3] NPT Equilibration ({MD_PARAMS['npt_steps']} steps = 1 ns)...")
system.addForce(MonteCarloBarostat(
    MD_PARAMS["pressure_atm"] * unit.atmospheres,
    MD_PARAMS["temperature_K"] * unit.kelvin,
    25,
))
simulation.context.reinitialize(preserveState=True)

t0 = time.time()
simulation.step(MD_PARAMS["npt_steps"])
print(f"  NPT done ({time.time()-t0:.1f}s)")

# Save equilibrated structure (THIS is the topology reference for analysis)
state = simulation.context.getState(getPositions=True, getEnergy=True)
equil_pdb = OUTPUT_DIR / "complex_equilibrated.pdb"
with open(equil_pdb, "w") as f:
    PDBFile.writeFile(simulation.topology, state.getPositions(), f)

# Save equilibration checkpoint (for resume if production crashes)
equil_chk = OUTPUT_DIR / "equilibrated.chk"
simulation.saveCheckpoint(str(equil_chk))

simulation.reporters.clear()

print(f"\n[OK] Equilibration complete!")
print(f"  Equilibrated PDB: {equil_pdb}")
print(f"  Checkpoint: {equil_chk}")

# %%
# ============================================================
# CELL 7: Production MD (100 ns) -- MONOLITH RUN
# ============================================================
# DESIGN DECISIONS:
#   1. Single continuous loop (NOT chunked restarts) -> clean RMSD
#   2. CheckpointReporter every 10ns -> crash recovery
#   3. Progress reporting via manual print (not a separate reporter)
#   4. If kernel crashes, re-run this cell: it loads checkpoint
#      and APPENDS to existing trajectory (no data loss)
#   5. TIME SAFETY: Auto-stop before Kaggle 12h limit to ensure
#      analysis cells run and output is saved

from openmm.app import CheckpointReporter

# ============================================================
# WALL-CLOCK TIME SAFETY
# Kaggle limit = 12h. We reserve 1h for analysis/export.
# ============================================================
MAX_WALL_SECONDS = 12.0 * 3600  # 12 hours total notebook limit
ANALYSIS_BUFFER_SECONDS = 0.25 * 3600  # 15 min for saving final state
MAX_MD_WALL = MAX_WALL_SECONDS - ANALYSIS_BUFFER_SECONDS  # 11.75h for MD
wall_start = t_start  # t_start defined in Cell 0 (line 28)

elapsed_so_far = time.time() - wall_start
print(f"[TIME]  Time elapsed so far: {elapsed_so_far/3600:.2f}h")
print(f"  MD time budget: {MAX_MD_WALL/3600:.1f}h")
print(f"  Remaining for MD: {(MAX_MD_WALL - elapsed_so_far)/3600:.1f}h")

print("=" * 60)
print(f"PRODUCTION MD: {MD_PARAMS['production_ns']} ns")
print(f"  Steps: {MD_PARAMS['production_steps']:,}")
print(f"  Timestep: {MD_PARAMS['timestep_fs']} fs")
print(f"  Temperature: {MD_PARAMS['temperature_K']} K")
print(f"  Save interval: {MD_PARAMS['save_interval_ps']} ps")
print("=" * 60)

# --- Output files ---
traj_dcd = OUTPUT_DIR / "production.dcd"
prod_log = OUTPUT_DIR / "production.log"
checkpoint_file = OUTPUT_DIR / "checkpoint.chk"
final_checkpoint = OUTPUT_DIR / "final.chk"

# --- Check for existing checkpoint (crash recovery) ---
# If the kernel crashed mid-production, we can resume from the last
# checkpoint. The DCD and log files will be appended to.
append_mode = False
if checkpoint_file.exists():
    try:
        simulation.loadCheckpoint(str(checkpoint_file))
        start_step = simulation.currentStep
        if start_step > 0:
            remaining = MD_PARAMS["production_steps"] - start_step
            ns_done = start_step * MD_PARAMS["timestep_fs"] / 1e6
            print(f"\n[WARN]  RESUMING from checkpoint!")
            print(f"  Step: {start_step:,} ({ns_done:.1f} ns)")
            print(f"  Remaining: {remaining:,} steps")
            append_mode = True
        else:
            print("\n  Checkpoint exists but at step 0, starting fresh...")
    except Exception as e:
        print(f"\n  Checkpoint load failed ({e}), starting fresh...")
        # Reset from equilibrated state
        if equil_chk.exists():
            simulation.loadCheckpoint(str(equil_chk))

if not append_mode:
    print("\n  Starting fresh production run...")
    start_step = simulation.currentStep

# --- Setup Reporters ---
simulation.reporters.clear()

# DCD trajectory
simulation.reporters.append(
    DCDReporter(str(traj_dcd), MD_PARAMS["save_interval_steps"], append=append_mode)
)

# Energy/temperature log
simulation.reporters.append(
    StateDataReporter(
        str(prod_log), MD_PARAMS["log_interval_steps"],
        step=True, time=True, potentialEnergy=True,
        kineticEnergy=True, totalEnergy=True,
        temperature=True, volume=True,
        density=True, speed=True,
        append=append_mode,
    )
)

# Checkpoint (every 10 ns for crash recovery)
simulation.reporters.append(
    CheckpointReporter(str(checkpoint_file), MD_PARAMS["checkpoint_interval_steps"])
)

# --- RUN PRODUCTION ---
print(f"\n[START] Starting at {time.strftime('%Y-%m-%d %H:%M:%S')}...")

run_start_time = time.time()
run_start_step = simulation.currentStep
total_target = MD_PARAMS["production_steps"]

# Report progress every 100 ps (50,000 steps at 2fs)
PROGRESS_CHUNK = 50000
time_stopped_early = False

while simulation.currentStep < total_target:
    # ===== TIME SAFETY CHECK =====
    wall_elapsed = time.time() - wall_start
    if wall_elapsed > MAX_MD_WALL:
        ns_now = simulation.currentStep * MD_PARAMS["timestep_fs"] / 1e6
        print(f"\n  [TIME]  TIME LIMIT REACHED after {wall_elapsed/3600:.1f}h wall time")
        print(f"  Stopping production at {ns_now:.1f} ns to ensure analysis runs")
        print(f"  Use continuation notebook to complete the remaining simulation")
        time_stopped_early = True
        break
    # =============================

    steps_to_run = min(PROGRESS_CHUNK, total_target - simulation.currentStep)
    simulation.step(steps_to_run)

    # Progress report
    current = simulation.currentStep
    elapsed = time.time() - run_start_time
    done_in_run = current - run_start_step
    if done_in_run > 0:
        rate = done_in_run / elapsed
        remaining = total_target - current
        eta_sec = remaining / rate if rate > 0 else 0
        pct = current / total_target * 100
        ns_done = current * MD_PARAMS["timestep_fs"] / 1e6
        wall_total = time.time() - wall_start
        wall_remaining = MAX_MD_WALL - wall_total
        print(f"  [{pct:5.1f}%] {ns_done:.1f}/{MD_PARAMS['production_ns']} ns | "
              f"{rate:.0f} steps/s | ETA: {eta_sec/3600:.1f}h | "
              f"Wall: {wall_total/3600:.1f}h (limit in {wall_remaining/3600:.1f}h)", flush=True)

# Save final state (ALWAYS, even if stopped early)
simulation.saveCheckpoint(str(final_checkpoint))

state = simulation.context.getState(getPositions=True, getVelocities=True, getEnergy=True)
final_pdb = OUTPUT_DIR / "production_final.pdb"
with open(final_pdb, "w") as f:
    PDBFile.writeFile(simulation.topology, state.getPositions(), f)

total_elapsed = time.time() - run_start_time
ACTUAL_NS = simulation.currentStep * MD_PARAMS["timestep_fs"] / 1e6

print(f"\n{'=' * 60}")
if time_stopped_early:
    print(f"[TIME]  Production MD stopped early (time safety)")
    print(f"  Achieved: {ACTUAL_NS:.1f} / {MD_PARAMS['production_ns']} ns")
    print(f"  [PIN] Save output as dataset, then run continuation notebook")
else:
    print(f"[OK] Production MD complete!")
    print(f"  Total time: {MD_PARAMS['production_ns']} ns")
print(f"  Wall time:  {total_elapsed/3600:.2f} hours")
print(f"  Trajectory: {traj_dcd}")
print(f"  Final PDB:  {final_pdb}")
print(f"  Checkpoint: {final_checkpoint}")
print(f"  Finished at {time.strftime('%Y-%m-%d %H:%M:%S')}")
print(f"{'=' * 60}")

