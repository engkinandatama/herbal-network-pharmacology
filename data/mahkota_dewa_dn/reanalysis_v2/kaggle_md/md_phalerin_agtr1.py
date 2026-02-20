# %% [markdown]
# # MD Simulation: Phalerin → AGTR1 (100 ns)
#
# **System:** Phalerin (ligand from *Phaleria macrocarpa*) bound to AGTR1 (Angiotensin II Type 1 Receptor)
#
# **PDB:** 3QXY | **Docking Affinity:** -9.96 kcal/mol (better than control Losartan -8.99)
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

# ─────────────────────────────────────────────────────────
# STEP 1: Install Miniforge3
# ─────────────────────────────────────────────────────────
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
    print(f"  ✅ Miniforge3 installed ({time.time()-t_start:.0f}s)")
else:
    print(f"  ✅ Already present")

# ─────────────────────────────────────────────────────────
# STEP 2: Install conda packages + cleanup
# ─────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("STEP 2: Mamba install...")
print("=" * 60)

# rdkit and scipy are NOT installed here — they come from pip
# to avoid GLIBCXX/CXXABI conflicts with system libstdc++.
packages = f"python={PY_VER} numpy openmm cudatoolkit openmmforcefields openff-toolkit ambertools pdbfixer mdtraj parmed"
print(f"  → {packages}")
t_install = time.time()
r = subprocess.run(
    f"{MAMBA} install -y -c conda-forge {packages}",
    shell=True, capture_output=True, text=True, timeout=600
)
if r.stdout:
    for line in r.stdout.strip().split('\n')[-10:]:
        print(f"    {line}")
if r.returncode != 0:
    print(f"\n  ⚠️ mamba failed (exit {r.returncode})")
    if r.stderr:
        for line in r.stderr.strip().split('\n')[-10:]:
            print(f"    STDERR: {line}")
else:
    print(f"\n  ✅ Packages installed ({time.time()-t_install:.0f}s)")

# Show conda numpy version
r_np = subprocess.run(
    f"{MINIFORGE_DIR}/bin/python -c 'import numpy; print(numpy.__version__)'",
    shell=True, capture_output=True, text=True
)
print(f"  Conda numpy: {r_np.stdout.strip() if r_np.returncode == 0 else '?'}")

# Remove conda rdkit+scipy (auto-installed as openff-toolkit deps).
# These need GLIBCXX_3.4.31 / CXXABI_1.3.15 → use pip versions.
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
print(f"  ✅ Removed conda rdkit+scipy ({removed} files)")

# ─────────────────────────────────────────────────────────
# STEP 2.5: Install pip packages (BEFORE removing system numpy)
# ─────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("STEP 2.5: Pip install (rdkit, scipy, meeko, gemmi, vina)...")
print("=" * 60)
for pkg in ["rdkit", "scipy", "meeko", "gemmi", "vina"]:
    r = subprocess.run([sys.executable, "-m", "pip", "install", "-q", pkg],
                       capture_output=True, text=True, timeout=300)
    if r.returncode == 0:
        print(f"  ✅ {pkg}")
    else:
        print(f"  ❌ {pkg}: {r.stderr[-300:]}")

# ─────────────────────────────────────────────────────────
# STEP 3: Remove system numpy & patch paths
# ─────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("STEP 3: Remove system numpy & patch paths...")
print("=" * 60)

# 3a: Remove system numpy (conflicts with conda numpy)
print("  Removing system numpy...")
r = subprocess.run([sys.executable, "-m", "pip", "uninstall", "-y", "numpy"],
                   capture_output=True, text=True)
if r.returncode == 0:
    print("  ✅ System numpy removed")
else:
    for np_dir in glob.glob("/usr/local/lib/python3.12/dist-packages/numpy*"):
        shutil.rmtree(np_dir, ignore_errors=True)
    print("  ✅ System numpy force-removed")

# 3b: Flush cached numpy modules
mods_to_flush = [k for k in sys.modules if 'numpy' in k.lower()]
for mod in mods_to_flush:
    del sys.modules[mod]
print(f"  ✅ Flushed {len(mods_to_flush)} cached modules")

# 3c: Add miniforge to sys.path (FIRST position)
site_dirs = glob.glob(f"{MINIFORGE_DIR}/lib/python{PY_VER}/site-packages")
site_dirs += glob.glob(f"{MINIFORGE_DIR}/lib/python*/site-packages")
for sp in sorted(set(site_dirs)):
    if os.path.isdir(sp) and sp not in sys.path:
        sys.path.insert(0, sp)
        print(f"  ✅ sys.path += {sp}")

# 3d: Environment variables
os.environ["PATH"] = f"{MINIFORGE_DIR}/bin:" + os.environ.get("PATH", "")
os.environ["LD_LIBRARY_PATH"] = f"{MINIFORGE_DIR}/lib:" + os.environ.get("LD_LIBRARY_PATH", "")
os.environ["AMBERHOME"] = MINIFORGE_DIR
print(f"  ✅ PATH, LD_LIBRARY_PATH, AMBERHOME set")

# ─────────────────────────────────────────────────────────
# STEP 4: Verify critical imports
# ─────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("STEP 4: Verifying critical imports...")
print("=" * 60)

# Verify numpy
try:
    import numpy as np
    print(f"  ✅ numpy {np.__version__}")
except Exception as e:
    print(f"  ❌ numpy FAILED: {e}")

# Verify OpenMM + CUDA
try:
    import openmm
    print(f"  ✅ OpenMM {openmm.__version__}")
    platforms = [openmm.Platform.getPlatform(i).getName()
                 for i in range(openmm.Platform.getNumPlatforms())]
    print(f"    Platforms: {platforms}")
    if "CUDA" in platforms:
        print(f"  ✅ CUDA platform ready")
    else:
        print(f"  ⚠️ CUDA not found — will use CPU (much slower!)")
except ImportError as e:
    print(f"  ❌ OpenMM import FAILED: {e}")
    raise

# All other packages
all_ok = True
for mod_name in ["pdbfixer", "mdtraj", "rdkit", "parmed",
                 "openff.toolkit", "openmmforcefields", "scipy",
                 "meeko", "vina"]:
    try:
        __import__(mod_name)
        print(f"  ✅ {mod_name}")
    except ImportError as e:
        print(f"  ❌ {mod_name}: {e}")
        all_ok = False

elapsed = time.time() - t_start
print(f"\n{'=' * 60}")
print(f"Setup complete in {elapsed:.0f}s")
if all_ok:
    print("✅ All dependencies verified!")
else:
    print("⚠️ Some packages missing — check logs above")
print(f"{'=' * 60}")

# %%
# ============================================================
# CELL 1: Configuration
# ============================================================
import json
from pathlib import Path

# --- System Definition ---
SYSTEM_NAME = "Phalerin_AGTR1"
RECEPTOR_PDB_ID = "3QXY"
LIGAND_NAME = "Phalerin"

# --- Paths (Kaggle) ---
DATASET_DIR = Path("/kaggle/input/mahkota-dewa-docking")
OUTPUT_DIR = Path("/kaggle/working/md_results") / SYSTEM_NAME
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Input files from dataset
# PDB: Full protein structure WITH hydrogens (same file that was converted to PDBQT for docking)
# SDF: Ligand with full bond-order info (needed for AM1-BCC charge assignment)
# PDBQT: For re-docking only (AutoDock Vina format)
RECEPTOR_PDB = DATASET_DIR / "receptors" / f"{RECEPTOR_PDB_ID}_withH.pdb"
RECEPTOR_PDBQT = DATASET_DIR / "receptors" / f"{RECEPTOR_PDB_ID}.pdbqt"
LIGAND_SDF = DATASET_DIR / "ligands" / f"{LIGAND_NAME}.sdf"
LIGAND_PDBQT = DATASET_DIR / "ligands" / f"{LIGAND_NAME}.pdbqt"
BINDING_SITES_JSON = DATASET_DIR / "binding_sites.json"

# --- Binding Site (from docking config) ---
with open(BINDING_SITES_JSON) as f:
    all_sites = json.load(f)
    bs = all_sites[RECEPTOR_PDB_ID]

BINDING_SITE = {
    "center": (bs["center_x"], bs["center_y"], bs["center_z"]),
    "size": (bs["size_x"], bs["size_y"], bs["size_z"]),
}

# --- MD Parameters ---
MD_PARAMS = {
    "temperature_K": 310.15,        # 37°C physiological
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
for f in [RECEPTOR_PDB, RECEPTOR_PDBQT, LIGAND_SDF, LIGAND_PDBQT, BINDING_SITES_JSON]:
    if f.exists():
        print(f"  ✅ {f.name} ({f.stat().st_size:,} bytes)")
    else:
        print(f"  ❌ MISSING: {f}")
        all_ok = False

assert all_ok, "Some input files are missing! Check dataset upload."

# %%
# ============================================================
# CELL 2: Prepare Protein (PDBFixer)
# ============================================================
# The dataset contains 3QXY_withH.pdb, preprocessed by `reduce`
# for docking. For MD we need to clean it because:
#   1. `reduce` added H to BOTH alternate conformations, then
#      stripped alt-loc indicators → duplicate heavy atoms
#   2. The H placement was for docking, not MD (wrong pH model)
#
# Strategy (gold standard for MD prep):
#   (a) Strip ALL hydrogens (PDBFixer will re-add at pH 7.4)
#   (b) Remove duplicate heavy atoms (keep first occurrence = conformer A)
#   (c) PDBFixer: fix missing residues/atoms → add H at pH 7.4
#
# Binding site geometry is preserved because heavy atom coords
# are from the SAME crystal structure (3QXY), and the ligand
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
            # (a) Skip hydrogens — element is in columns 77-78
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
print("  ✅ Missing atoms added")

# Add ALL hydrogens at pH 7.4 (physiological)
fixer.addMissingHydrogens(7.4)

n_atoms_final = len(list(fixer.topology.atoms()))
n_res_final = len(list(fixer.topology.residues()))
print(f"  Final: {n_atoms_final} atoms, {n_res_final} residues")

# Save fixed protein
protein_pdb = OUTPUT_DIR / "protein_fixed.pdb"
with open(protein_pdb, "w") as f:
    PDBFile.writeFile(fixer.topology, fixer.positions, f)
print(f"\n✅ Protein saved: {protein_pdb}")

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
print(f"\n✅ Ligand saved: {ligand_sdf_prep}")

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

# Convert docked PDBQT → RDKit Mol for SDF/PDB output
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
        print("  Converted via Meeko ✅")
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
                    print("  Converted via obabel ✅")
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
        print(f"  Mapped {n_heavy} heavy atom coords ✅")
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

print(f"\n✅ Docked ligand saved:")
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
print(f"\n✅ Solvated complex saved: {solvated_pdb}")

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
    print(f"  ⚠️ Using CPU platform (CUDA not available)")

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
print(f"  ΔE = {e_after - e_before:.1f} kJ/mol ({time.time()-t0:.1f}s)")

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

print(f"\n✅ Equilibration complete!")
print(f"  Equilibrated PDB: {equil_pdb}")
print(f"  Checkpoint: {equil_chk}")

# %%
# ============================================================
# CELL 7: Production MD (100 ns) — MONOLITH RUN
# ============================================================
# DESIGN DECISIONS:
#   1. Single continuous loop (NOT chunked restarts) → clean RMSD
#   2. CheckpointReporter every 10ns → crash recovery
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
MAX_WALL_SECONDS = 11.0 * 3600  # 11 hours total notebook limit
ANALYSIS_BUFFER_SECONDS = 1.0 * 3600  # 1 hour for analysis
MAX_MD_WALL = MAX_WALL_SECONDS - ANALYSIS_BUFFER_SECONDS  # 10h for MD
wall_start = t_start  # t_start defined in Cell 0 (line 28)

elapsed_so_far = time.time() - wall_start
print(f"⏱️  Time elapsed so far: {elapsed_so_far/3600:.2f}h")
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
            print(f"\n⚠️  RESUMING from checkpoint!")
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
print(f"\n🚀 Starting at {time.strftime('%Y-%m-%d %H:%M:%S')}...")

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
        print(f"\n  ⏱️  TIME LIMIT REACHED after {wall_elapsed/3600:.1f}h wall time")
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
    print(f"⏱️  Production MD stopped early (time safety)")
    print(f"  Achieved: {ACTUAL_NS:.1f} / {MD_PARAMS['production_ns']} ns")
    print(f"  📌 Save output as dataset, then run continuation notebook")
else:
    print(f"✅ Production MD complete!")
    print(f"  Total time: {MD_PARAMS['production_ns']} ns")
print(f"  Wall time:  {total_elapsed/3600:.2f} hours")
print(f"  Trajectory: {traj_dcd}")
print(f"  Final PDB:  {final_pdb}")
print(f"  Checkpoint: {final_checkpoint}")
print(f"  Finished at {time.strftime('%Y-%m-%d %H:%M:%S')}")
print(f"{'=' * 60}")

# %%
# ============================================================
# CELL 8: Trajectory Analysis (10 metrics)
# ============================================================
# Metrics: RMSD, RMSF, Rg, H-bonds, Contact Frequency,
#          Key Distances, SASA, Equilibration Check,
#          Representative Frames
import mdtraj as md
import numpy as np
import warnings
warnings.filterwarnings("ignore")

print("Loading trajectory...")

# Load trajectory using equilibrated PDB as topology reference
traj = md.load(str(traj_dcd), top=str(equil_pdb))
print(f"  Frames: {traj.n_frames}")
print(f"  Atoms: {traj.n_atoms}")
print(f"  Time: {traj.time[0]:.1f} - {traj.time[-1]:.1f} ps")

# --- [0/10] PBC Reimaging (fix ligand jumping artifacts) ---
print("\n[0/10] Reimaging trajectory (fixing PBC artifacts)...")
# This centers the protein and ensures the ligand is placed
# adjacent to it, fixing visual "jumping" across periodic boundaries.
# Does NOT change the physics — only coordinate representation.
try:
    protein_atoms_set = set(traj.topology.select("protein").tolist())
    anchor = [protein_atoms_set]
    traj = traj.image_molecules(inplace=False, anchor_molecules=anchor)
    print("  ✅ Reimaging complete — PBC artifacts fixed")
except Exception as e:
    print(f"  ⚠️ Reimaging failed ({e}), continuing without it")
    print("  (Visualization may show ligand jumping, but analysis is unaffected)")

# --- Atom Selections ---
protein_atoms = traj.topology.select("protein")
ca_atoms = traj.topology.select("protein and name CA")
backbone_atoms = traj.topology.select("protein and (name CA or name C or name N or name O)")

# Robust ligand detection: everything that's NOT protein, water, or ions
water_atoms = set(traj.topology.select("water"))
# Handle both naming conventions (NA/CL and Na+/Cl-)
try:
    ion_atoms = set(traj.topology.select("resname NA CL"))
except:
    ion_atoms = set()

all_known = set(protein_atoms) | water_atoms | ion_atoms
ligand_atoms = np.array(sorted([i for i in range(traj.n_atoms) if i not in all_known]))

print(f"  Protein atoms: {len(protein_atoms)}")
print(f"  Backbone atoms: {len(backbone_atoms)}")
print(f"  CA atoms: {len(ca_atoms)}")
print(f"  Ligand atoms: {len(ligand_atoms)}")

assert len(ligand_atoms) > 0, "No ligand atoms detected! Check topology."
assert len(ca_atoms) > 0, "No CA atoms detected!"

# --- Time axis ---
time_ns = traj.time / 1000.0  # ps → ns

# --- [1/10] RMSD (protein backbone) ---
print("\n[1/10] Calculating backbone RMSD...")
traj_bb = traj.atom_slice(backbone_atoms)
rmsd_protein = md.rmsd(traj_bb, traj_bb, frame=0) * 10  # nm → Å
print(f"  Mean: {np.mean(rmsd_protein):.2f} ± {np.std(rmsd_protein):.2f} Å")

# --- [2/10] RMSD (ligand, after protein alignment) ---
print("[2/10] Calculating ligand RMSD...")
traj_aligned = traj.superpose(traj, frame=0, atom_indices=ca_atoms)
lig_traj = traj_aligned.atom_slice(ligand_atoms)
rmsd_ligand = md.rmsd(lig_traj, lig_traj, frame=0) * 10  # nm → Å
print(f"  Mean: {np.mean(rmsd_ligand):.2f} ± {np.std(rmsd_ligand):.2f} Å")

# --- [3/10] RMSF (Cα, full trajectory) ---
print("[3/10] Calculating RMSF (Cα)...")
ca_traj = traj_aligned.atom_slice(ca_atoms)
rmsf = md.rmsf(ca_traj, ca_traj, frame=0) * 10  # nm → Å
residue_ids = [traj.topology.atom(a).residue.resSeq for a in ca_atoms]
print(f"  Max RMSF: {np.max(rmsf):.2f} Å (residue {residue_ids[np.argmax(rmsf)]})")

# --- [4/10] Radius of Gyration ---
print("[4/10] Calculating Radius of Gyration...")
rg = md.compute_rg(traj.atom_slice(protein_atoms)) * 10  # nm → Å
print(f"  Mean: {np.mean(rg):.2f} ± {np.std(rg):.2f} Å")

# --- [5/10] Hydrogen Bonds (OPTIMIZED: sample every Nth frame) ---
print("\n[5/10] Calculating protein-ligand H-bonds (sampled)...")

lig_set = set(ligand_atoms.tolist())
prot_set = set(protein_atoms.tolist())

# Sample every 10th frame to avoid O(n²) slowness
# For 10K frames → 1K sampled → still statistically valid
sample_step = max(1, traj.n_frames // 1000)
sampled_indices = list(range(0, traj.n_frames, sample_step))

hbond_counts_sampled = []
for idx in sampled_indices:
    frame = traj[idx]
    try:
        hbs = md.baker_hubbard(frame, freq=0.0)
        count = 0
        for hb in hbs:
            d, h, a = hb
            if (d in lig_set and a in prot_set) or (d in prot_set and a in lig_set):
                count += 1
        hbond_counts_sampled.append(count)
    except Exception:
        hbond_counts_sampled.append(0)

# Interpolate to full frame count for plotting
hbond_counts_sampled = np.array(hbond_counts_sampled)
hbond_counts = np.interp(
    np.arange(traj.n_frames),
    sampled_indices,
    hbond_counts_sampled
)

print(f"  Sampled {len(sampled_indices)} frames (every {sample_step})")
print(f"  Mean H-bonds: {np.mean(hbond_counts):.1f} ± {np.std(hbond_counts):.1f}")

# --- [6/10] Protein-Ligand Contact Frequency ---
print("\n[6/10] Calculating protein-ligand contact frequency...")
# Find residues within 4.5 Å of ligand in each sampled frame
CONTACT_CUTOFF_NM = 0.45  # 4.5 Å in nm
residue_contact_counts = {}  # resSeq → count

# Get ligand atom indices as set for fast lookup
for idx in sampled_indices:
    frame = traj[idx]
    for ca_idx in ca_atoms:
        res = traj.topology.atom(ca_idx).residue
        # Compute min distance between this residue's atoms and any ligand atom
        res_atoms = [a.index for a in res.atoms]
        pairs = [(r, l) for r in res_atoms for l in ligand_atoms]
        if len(pairs) == 0:
            continue
        distances = md.compute_distances(frame, pairs)[0]
        if np.min(distances) < CONTACT_CUTOFF_NM:
            key = (res.resSeq, res.name)
            residue_contact_counts[key] = residue_contact_counts.get(key, 0) + 1

# Convert to frequency (fraction of sampled frames)
n_sampled = len(sampled_indices)
contact_freq = {k: v / n_sampled for k, v in residue_contact_counts.items()}
contact_freq_sorted = sorted(contact_freq.items(), key=lambda x: -x[1])

print(f"  Residues contacting ligand (>10% of time):")
top_contacts = [(k, v) for k, v in contact_freq_sorted if v > 0.1]
for (resSeq, resName), freq in top_contacts[:15]:
    print(f"    {resName}{resSeq}: {freq*100:.1f}%")

# --- [7/10] Key Distance Monitoring ---
print("\n[7/10] Monitoring key protein-ligand distances...")
# Track center-of-mass distance between ligand and top-3 contacting residues
lig_com_distances = {}  # resSeq → array of distances over time

# Get top 3 contacting residues
top3_residues = [(resSeq, resName) for (resSeq, resName), _ in contact_freq_sorted[:3]]

for resSeq, resName in top3_residues:
    # Find all atoms of this residue
    res_atom_indices = []
    for a in traj.topology.atoms:
        if a.residue.resSeq == resSeq and a.residue.name == resName:
            res_atom_indices.append(a.index)
    if len(res_atom_indices) == 0:
        continue

    # Compute min distance between residue and ligand for each frame
    # Sample for speed
    dists = []
    for idx in range(0, traj.n_frames, sample_step):
        frame = traj[idx]
        pairs = [(r, l) for r in res_atom_indices for l in ligand_atoms]
        d = md.compute_distances(frame, pairs)[0]
        dists.append(np.min(d) * 10)  # nm → Å

    # Interpolate to full frame count
    dists_full = np.interp(
        np.arange(traj.n_frames),
        list(range(0, traj.n_frames, sample_step)),
        dists
    )
    lig_com_distances[f"{resName}{resSeq}"] = dists_full
    print(f"  {resName}{resSeq}: mean={np.mean(dists_full):.2f} ± {np.std(dists_full):.2f} Å")

# --- [8/10] SASA (Solvent Accessible Surface Area) ---
print("\n[8/10] Calculating SASA...")
# Compute SASA for protein only (sampled for speed)
sasa_sample_step = max(1, traj.n_frames // 500)
sasa_indices = list(range(0, traj.n_frames, sasa_sample_step))

prot_traj = traj.atom_slice(protein_atoms)
sasa_sampled = []
for idx in sasa_indices:
    frame_sasa = md.shrake_rupley(prot_traj[idx], mode='residue')
    sasa_sampled.append(np.sum(frame_sasa))  # total SASA in nm²

sasa_sampled = np.array(sasa_sampled) * 100  # nm² → Å²
sasa_full = np.interp(np.arange(traj.n_frames), sasa_indices, sasa_sampled)
print(f"  Mean total SASA: {np.mean(sasa_full):.1f} ± {np.std(sasa_full):.1f} Å²")
print(f"  (sampled {len(sasa_indices)} frames)")

# --- [9/10] Equilibration Check ---
print("\n[9/10] Checking equilibration...")
# Split trajectory into 5 blocks and check if RMSD mean is stable
n_blocks = 5
block_size = len(rmsd_protein) // n_blocks
block_means = []
for i in range(n_blocks):
    start = i * block_size
    end = start + block_size
    block_means.append(np.mean(rmsd_protein[start:end]))

print(f"  RMSD block means (each ~{MD_PARAMS['production_ns']/n_blocks:.0f} ns):")
for i, bm in enumerate(block_means):
    print(f"    Block {i+1}: {bm:.2f} Å")

# Check: if first block deviates >50% from last 3 blocks, equilibration is <20ns
last3_mean = np.mean(block_means[2:])
first_dev = abs(block_means[0] - last3_mean) / last3_mean * 100
if first_dev > 50:
    equil_ns = MD_PARAMS['production_ns'] / n_blocks  # skip first block
    print(f"  ⚠️ First block deviates {first_dev:.0f}% — recommend skipping first {equil_ns:.0f} ns")
else:
    equil_ns = 0
    print(f"  ✅ System appears well-equilibrated (deviation: {first_dev:.0f}%)")

# Drift check on last 3 blocks
from scipy import stats as scipy_stats
last3_time = np.arange(len(block_means[2:]))
slope, _, _, p_val, _ = scipy_stats.linregress(last3_time, block_means[2:])
if p_val < 0.05:
    print(f"  ⚠️ Possible drift in last 60ns (slope={slope:.4f} Å/block, p={p_val:.3f})")
else:
    print(f"  ✅ No significant drift (p={p_val:.3f})")

# --- [10/10] Representative Frame Extraction ---
print("\n[10/10] Extracting representative frames...")
# Save: first frame, middle frame, last frame, and frame closest to mean RMSD
# These are useful for making figures in PyMOL/ChimeraX

rep_frames = {
    "initial": 0,
    "middle": traj.n_frames // 2,
    "final": traj.n_frames - 1,
    "equilibrated_mean": int(np.argmin(np.abs(rmsd_protein - np.mean(rmsd_protein[traj.n_frames//5:])))),
}

for label, fidx in rep_frames.items():
    # Save solute only (no water/ions) for visualization
    solute_atoms_arr = np.concatenate([protein_atoms, ligand_atoms])
    frame = traj[fidx].atom_slice(solute_atoms_arr)
    out_path = OUTPUT_DIR / f"frame_{label}.pdb"
    frame.save(str(out_path))
    rmsd_val = rmsd_protein[fidx]
    print(f"  {label}: frame {fidx} (t={time_ns[fidx]:.1f} ns, RMSD={rmsd_val:.2f} Å) → {out_path.name}")

# --- SAVE intermediate results NOW (before risky MM-GBSA) ---
print("\n💾 Saving intermediate analysis results...")
import pandas as pd

ts_df = pd.DataFrame({
    "Time_ns": time_ns,
    "RMSD_Protein_A": rmsd_protein,
    "RMSD_Ligand_A": rmsd_ligand,
    "Rg_A": rg,
    "Hbonds": hbond_counts,
    "SASA_A2": sasa_full,
})
ts_csv = OUTPUT_DIR / "timeseries.csv"
ts_df.to_csv(ts_csv, index=False)

rmsf_df = pd.DataFrame({"Residue": residue_ids, "RMSF_A": rmsf})
rmsf_csv = OUTPUT_DIR / "rmsf_per_residue.csv"
rmsf_df.to_csv(rmsf_csv, index=False)

# Save contact frequency
contact_df = pd.DataFrame([
    {"Residue": f"{rn}{rs}", "ResSeq": rs, "ResName": rn, "Contact_Frequency": freq}
    for (rs, rn), freq in contact_freq_sorted if freq > 0.05
])
contact_csv = OUTPUT_DIR / "contact_frequency.csv"
contact_df.to_csv(contact_csv, index=False)

# Save key distances
if lig_com_distances:
    dist_df = pd.DataFrame({"Time_ns": time_ns, **lig_com_distances})
    dist_csv = OUTPUT_DIR / "key_distances.csv"
    dist_df.to_csv(dist_csv, index=False)

print(f"  ✅ timeseries.csv ({len(ts_df)} rows)")
print(f"  ✅ rmsf_per_residue.csv ({len(rmsf_df)} rows)")
print(f"  ✅ contact_frequency.csv ({len(contact_df)} rows)")
if lig_com_distances:
    print(f"  ✅ key_distances.csv")

print(f"\n{'=' * 60}")
print("ANALYSIS SUMMARY:")
print(f"  RMSD protein:  {np.mean(rmsd_protein):.2f} ± {np.std(rmsd_protein):.2f} Å")
print(f"  RMSD ligand:   {np.mean(rmsd_ligand):.2f} ± {np.std(rmsd_ligand):.2f} Å")
print(f"  Rg:            {np.mean(rg):.2f} ± {np.std(rg):.2f} Å")
print(f"  H-bonds:       {np.mean(hbond_counts):.1f} ± {np.std(hbond_counts):.1f}")
print(f"  SASA:          {np.mean(sasa_full):.1f} ± {np.std(sasa_full):.1f} Å²")
print(f"  Top contacts:  {', '.join([f'{rn}{rs}' for (rs, rn), _ in contact_freq_sorted[:5]])}")
print(f"  Equil. check:  {'OK' if first_dev <= 50 else f'skip first {equil_ns:.0f}ns'}")
print(f"  Rep. frames:   {len(rep_frames)} PDBs saved")
print(f"{'=' * 60}")

# %%
# ============================================================
# CELL 9: MM-GBSA Binding Free Energy
# ============================================================
# ΔG_bind = G_complex - G_receptor - G_ligand
# Using OpenMM implicit solvent (GB-OBC2) for the solvation energy.
#
# OPTIMIZATION: Create ForceField and System ONCE, reuse for all frames.
# Previous version created them 300x (3 per frame × 100 frames) = very slow.

print("=" * 60)
print("MM-GBSA BINDING FREE ENERGY ANALYSIS")
print("=" * 60)

N_FRAMES_MMGBSA = 100
# Skip first 20% as equilibration period
start_frame = int(traj.n_frames * 0.2)
frame_indices = np.linspace(start_frame, traj.n_frames - 1, N_FRAMES_MMGBSA, dtype=int)

print(f"  Analyzing {N_FRAMES_MMGBSA} frames (from {start_frame} to {traj.n_frames-1})")

# --- Pre-create stripped trajectory components (no water/ions) ---
solute_atoms = np.concatenate([protein_atoms, ligand_atoms])
prot_local_idx = np.arange(len(protein_atoms))  # indices within solute
lig_local_idx = np.arange(len(protein_atoms), len(protein_atoms) + len(ligand_atoms))

# --- Create ForceField ONCE with implicit solvent ---
print("  Setting up implicit solvent force field (one-time)...")
ff_gb = ForceField("amber/ff14SB.xml", "implicit/obc2.xml")
gaff_gen_gb = GAFFTemplateGenerator(molecules=[off_mol], forcefield="gaff-2.11")
ff_gb.registerTemplateGenerator(gaff_gen_gb.generator)

# --- Pre-build topologies and systems for each component ---
# Extract one frame to build topology templates
ref_frame = traj[frame_indices[0]]
ref_solute = ref_frame.atom_slice(solute_atoms)
ref_receptor = ref_frame.atom_slice(protein_atoms)
ref_ligand = ref_frame.atom_slice(ligand_atoms)

# Save reference PDBs for topology
tmp_dir = OUTPUT_DIR / "_tmp_mmgbsa"
tmp_dir.mkdir(exist_ok=True)

tmp_complex_pdb = tmp_dir / "complex.pdb"
tmp_receptor_pdb = tmp_dir / "receptor.pdb"
tmp_ligand_pdb = tmp_dir / "ligand.pdb"

ref_solute.save(str(tmp_complex_pdb))
ref_receptor.save(str(tmp_receptor_pdb))
ref_ligand.save(str(tmp_ligand_pdb))

# Build OpenMM systems ONCE
print("  Building OpenMM systems for complex/receptor/ligand...")

def build_gb_system(pdb_path):
    """Build an implicit-solvent system from PDB."""
    pdb = PDBFile(str(pdb_path))
    sys_gb = ff_gb.createSystem(pdb.topology, nonbondedMethod=NoCutoff, constraints=HBonds)
    integ = LangevinMiddleIntegrator(300*unit.kelvin, 1/unit.picoseconds, 2*unit.femtoseconds)
    sim = Simulation(pdb.topology, sys_gb, integ)
    return sim

try:
    sim_complex = build_gb_system(tmp_complex_pdb)
    sim_receptor = build_gb_system(tmp_receptor_pdb)
    sim_ligand = build_gb_system(tmp_ligand_pdb)
    print("  ✅ All three systems built successfully")
    mmgbsa_ready = True
except Exception as e:
    print(f"  ❌ System build failed: {e}")
    print("  Skipping MM-GBSA analysis.")
    mmgbsa_ready = False

# --- Calculate energies per frame ---
dG_values = []
errors = 0

if mmgbsa_ready:
    print(f"\n  Computing {N_FRAMES_MMGBSA} frames...")
    t0 = time.time()

    for i, fidx in enumerate(frame_indices):
        try:
            frame = traj[fidx]
            frame_solute = frame.atom_slice(solute_atoms)

            # Get positions for each component (convert nm → OpenMM units)
            pos_complex = frame_solute.xyz[0]  # nm, shape (N_solute, 3)
            pos_receptor = pos_complex[prot_local_idx]
            pos_ligand = pos_complex[lig_local_idx]

            # Calculate energy: complex
            sim_complex.context.setPositions(pos_complex)
            E_complex = sim_complex.context.getState(getEnergy=True).getPotentialEnergy()
            E_complex = E_complex.value_in_unit(unit.kilocalories_per_mole)

            # Calculate energy: receptor
            sim_receptor.context.setPositions(pos_receptor)
            E_receptor = sim_receptor.context.getState(getEnergy=True).getPotentialEnergy()
            E_receptor = E_receptor.value_in_unit(unit.kilocalories_per_mole)

            # Calculate energy: ligand
            sim_ligand.context.setPositions(pos_ligand)
            E_ligand = sim_ligand.context.getState(getEnergy=True).getPotentialEnergy()
            E_ligand = E_ligand.value_in_unit(unit.kilocalories_per_mole)

            dG = E_complex - E_receptor - E_ligand
            dG_values.append(dG)

        except Exception as e:
            errors += 1
            if errors <= 3:
                print(f"    Error at frame {fidx}: {e}")

        if (i + 1) % 25 == 0:
            elapsed = time.time() - t0
            eta = elapsed / (i + 1) * (N_FRAMES_MMGBSA - i - 1)
            print(f"  [{i+1}/{N_FRAMES_MMGBSA}] {len(dG_values)} OK, {errors} errors | ETA: {eta:.0f}s")

dG_values = np.array(dG_values)

# Cleanup temp files
import shutil
if tmp_dir.exists():
    shutil.rmtree(tmp_dir)

if len(dG_values) > 0:
    dG_mean = np.mean(dG_values)
    dG_std = np.std(dG_values)
    dG_sem = dG_std / np.sqrt(len(dG_values))

    print(f"\n{'=' * 60}")
    print(f"MM-GBSA Results ({len(dG_values)}/{N_FRAMES_MMGBSA} frames):")
    print(f"  ΔG_bind = {dG_mean:.2f} ± {dG_sem:.2f} kcal/mol")
    print(f"  Std dev = {dG_std:.2f} kcal/mol")
    print(f"  Range: [{np.min(dG_values):.2f}, {np.max(dG_values):.2f}] kcal/mol")
    print(f"{'=' * 60}")
else:
    print("\n⚠️ MM-GBSA failed for all frames. Skipping.")
    dG_mean, dG_std, dG_sem = 0.0, 0.0, 0.0

# %%
# ============================================================
# CELL 10: Publication Figures
# ============================================================
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 10,
    "axes.labelsize": 12,
    "axes.titlesize": 13,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
    "figure.dpi": 300,
    "axes.linewidth": 1.2,
})

# --- Figure 1: 4-Panel MD Analysis ---
fig = plt.figure(figsize=(14, 10))
gs = GridSpec(2, 2, hspace=0.35, wspace=0.3)

window = max(1, len(rmsd_protein) // 100)

# Panel A: Protein RMSD
ax1 = fig.add_subplot(gs[0, 0])
ax1.plot(time_ns, rmsd_protein, color="#2196F3", linewidth=0.3, alpha=0.4)
rmsd_smooth = np.convolve(rmsd_protein, np.ones(window)/window, mode="valid")
time_smooth = time_ns[:len(rmsd_smooth)]
ax1.plot(time_smooth, rmsd_smooth, color="#1565C0", linewidth=1.5, label="Running avg")
ax1.set_xlabel("Time (ns)")
ax1.set_ylabel("RMSD (Å)")
ax1.set_title(f"(A) Protein Backbone RMSD\nmean={np.mean(rmsd_protein):.2f}±{np.std(rmsd_protein):.2f} Å")
ax1.legend()
ax1.set_xlim(0, MD_PARAMS["production_ns"])

# Panel B: Ligand RMSD
ax2 = fig.add_subplot(gs[0, 1])
ax2.plot(time_ns, rmsd_ligand, color="#FF9800", linewidth=0.3, alpha=0.4)
rmsd_lig_smooth = np.convolve(rmsd_ligand, np.ones(window)/window, mode="valid")
ax2.plot(time_smooth, rmsd_lig_smooth, color="#E65100", linewidth=1.5, label="Running avg")
ax2.set_xlabel("Time (ns)")
ax2.set_ylabel("RMSD (Å)")
ax2.set_title(f"(B) Ligand RMSD\nmean={np.mean(rmsd_ligand):.2f}±{np.std(rmsd_ligand):.2f} Å")
ax2.legend()
ax2.set_xlim(0, MD_PARAMS["production_ns"])

# Panel C: RMSF
ax3 = fig.add_subplot(gs[1, 0])
ax3.plot(residue_ids, rmsf, color="#4CAF50", linewidth=1.0)
ax3.fill_between(residue_ids, 0, rmsf, color="#4CAF50", alpha=0.15)
ax3.set_xlabel("Residue Number")
ax3.set_ylabel("RMSF (Å)")
ax3.set_title("(C) Cα RMSF per Residue")
threshold = np.mean(rmsf) + 2 * np.std(rmsf)
ax3.axhline(threshold, color="red", linestyle="--", alpha=0.5,
            label=f"Mean+2σ ({threshold:.1f} Å)")
ax3.legend()

# Panel D: Radius of Gyration
ax4 = fig.add_subplot(gs[1, 1])
ax4.plot(time_ns, rg, color="#9C27B0", linewidth=0.3, alpha=0.4)
rg_smooth = np.convolve(rg, np.ones(window)/window, mode="valid")
ax4.plot(time_smooth, rg_smooth, color="#6A1B9A", linewidth=1.5, label="Running avg")
ax4.set_xlabel("Time (ns)")
ax4.set_ylabel("Rg (Å)")
ax4.set_title(f"(D) Radius of Gyration\nmean={np.mean(rg):.2f}±{np.std(rg):.2f} Å")
ax4.legend()
ax4.set_xlim(0, MD_PARAMS["production_ns"])

fig.suptitle(f"MD Simulation: {LIGAND_NAME} – {RECEPTOR_PDB_ID} (AGTR1) | {MD_PARAMS['production_ns']} ns",
             fontsize=14, fontweight="bold", y=0.98)

fig_path_1 = OUTPUT_DIR / "md_analysis_4panel.png"
plt.savefig(fig_path_1, dpi=300, bbox_inches="tight", facecolor="white")
plt.close()
print(f"✅ Figure 1 saved: {fig_path_1}")

# --- Figure 2: H-bond + MM-GBSA ---
fig2, (ax5, ax6) = plt.subplots(1, 2, figsize=(14, 5))

# Panel E: H-bonds
ax5.plot(time_ns, hbond_counts, color="#F44336", linewidth=0.3, alpha=0.4)
hb_smooth = np.convolve(hbond_counts, np.ones(window)/window, mode="valid")
ax5.plot(time_smooth, hb_smooth, color="#B71C1C", linewidth=1.5, label="Running avg")
ax5.set_xlabel("Time (ns)")
ax5.set_ylabel("Number of H-bonds")
ax5.set_title(f"(E) Protein–Ligand H-bonds\nmean={np.mean(hbond_counts):.1f}±{np.std(hbond_counts):.1f}")
ax5.legend()
ax5.set_xlim(0, MD_PARAMS["production_ns"])

# Panel F: MM-GBSA
if len(dG_values) > 0:
    ax6.hist(dG_values, bins=25, color="#00BCD4", edgecolor="black", alpha=0.7)
    ax6.axvline(dG_mean, color="red", linewidth=2, linestyle="--",
                label=f"Mean: {dG_mean:.2f} ± {dG_sem:.2f}")
    ax6.set_xlabel("ΔG_bind (kcal/mol)")
    ax6.set_ylabel("Count")
    ax6.set_title(f"(F) MM-GBSA Binding Free Energy\nΔG = {dG_mean:.2f} ± {dG_sem:.2f} kcal/mol")
    ax6.legend()
else:
    ax6.text(0.5, 0.5, "MM-GBSA\nnot available", ha="center", va="center",
             transform=ax6.transAxes, fontsize=14)
    ax6.set_title("(F) MM-GBSA Binding Free Energy")

fig2.suptitle(f"H-bond & Binding Energy: {LIGAND_NAME} – {RECEPTOR_PDB_ID}",
              fontsize=14, fontweight="bold")
plt.tight_layout()

fig_path_2 = OUTPUT_DIR / "md_hbond_mmgbsa.png"
plt.savefig(fig_path_2, dpi=300, bbox_inches="tight", facecolor="white")
plt.close()
print(f"✅ Figure 2 saved: {fig_path_2}")

# --- Figure 3: SASA, Contact Frequency, Key Distances ---
fig3 = plt.figure(figsize=(16, 5))
gs3 = GridSpec(1, 3, wspace=0.35)

# Panel G: SASA over time
ax7 = fig3.add_subplot(gs3[0, 0])
ax7.plot(time_ns, sasa_full, color="#607D8B", linewidth=0.3, alpha=0.4)
sasa_smooth = np.convolve(sasa_full, np.ones(window)/window, mode="valid")
ax7.plot(time_smooth, sasa_smooth, color="#37474F", linewidth=1.5, label="Running avg")
ax7.set_xlabel("Time (ns)")
ax7.set_ylabel("SASA (Å²)")
ax7.set_title(f"(G) Protein SASA\nmean={np.mean(sasa_full):.0f}±{np.std(sasa_full):.0f} Å²")
ax7.legend()
ax7.set_xlim(0, MD_PARAMS["production_ns"])

# Panel H: Contact frequency (bar chart)
ax8 = fig3.add_subplot(gs3[0, 1])
top_n = min(12, len(top_contacts))
if top_n > 0:
    labels_c = [f"{rn}{rs}" for (rs, rn), _ in top_contacts[:top_n]]
    freqs_c = [freq * 100 for _, freq in top_contacts[:top_n]]
    colors_c = plt.cm.YlOrRd(np.linspace(0.3, 0.9, top_n))
    bars = ax8.barh(range(top_n), freqs_c, color=colors_c, edgecolor="black", linewidth=0.5)
    ax8.set_yticks(range(top_n))
    ax8.set_yticklabels(labels_c)
    ax8.set_xlabel("Contact Frequency (%)")
    ax8.set_title("(H) Protein-Ligand Contacts")
    ax8.invert_yaxis()
else:
    ax8.text(0.5, 0.5, "No contacts\ndetected", ha="center", va="center",
             transform=ax8.transAxes, fontsize=14)

# Panel I: Key distances over time
ax9 = fig3.add_subplot(gs3[0, 2])
dist_colors = ["#E91E63", "#3F51B5", "#009688"]
for i, (label, dists) in enumerate(lig_com_distances.items()):
    c = dist_colors[i % len(dist_colors)]
    ax9.plot(time_ns, dists, color=c, linewidth=0.3, alpha=0.3)
    d_smooth = np.convolve(dists, np.ones(window)/window, mode="valid")
    ax9.plot(time_smooth, d_smooth, color=c, linewidth=1.5, label=label)
ax9.axhline(4.5, color="gray", linestyle=":", alpha=0.5, label="Contact cutoff")
ax9.set_xlabel("Time (ns)")
ax9.set_ylabel("Min Distance (Å)")
ax9.set_title("(I) Key Residue–Ligand Distances")
ax9.legend(fontsize=7, loc="upper right")
ax9.set_xlim(0, MD_PARAMS["production_ns"])

fig3.suptitle(f"Extended Analysis: {LIGAND_NAME} – {RECEPTOR_PDB_ID}",
              fontsize=14, fontweight="bold")
plt.tight_layout()

fig_path_3 = OUTPUT_DIR / "md_extended_analysis.png"
plt.savefig(fig_path_3, dpi=300, bbox_inches="tight", facecolor="white")
plt.close()
print(f"✅ Figure 3 saved: {fig_path_3}")

# %%
# ============================================================
# CELL 11: Export Summary Statistics
# ============================================================
import json as json_mod

summary = {
    "System": SYSTEM_NAME,
    "Receptor": RECEPTOR_PDB_ID,
    "Ligand": LIGAND_NAME,
    "Simulation_Time_ns": MD_PARAMS["production_ns"],
    "N_Frames": int(traj.n_frames),
    "N_Atoms_Total": int(n_atoms_total),
    "Docking_Affinity_kcal_mol": float(best_affinity),
    "RMSD_Protein_Mean_A": float(np.mean(rmsd_protein)),
    "RMSD_Protein_Std_A": float(np.std(rmsd_protein)),
    "RMSD_Ligand_Mean_A": float(np.mean(rmsd_ligand)),
    "RMSD_Ligand_Std_A": float(np.std(rmsd_ligand)),
    "Rg_Mean_A": float(np.mean(rg)),
    "Rg_Std_A": float(np.std(rg)),
    "Hbond_Mean": float(np.mean(hbond_counts)),
    "Hbond_Std": float(np.std(hbond_counts)),
    "SASA_Mean_A2": float(np.mean(sasa_full)),
    "SASA_Std_A2": float(np.std(sasa_full)),
    "Top_Contacts": ", ".join([f"{rn}{rs}" for (rs, rn), _ in contact_freq_sorted[:5]]),
    "Equilibration_OK": first_dev <= 50,
    "MMGBSA_dG_Mean_kcal_mol": float(dG_mean) if len(dG_values) > 0 else None,
    "MMGBSA_dG_SEM_kcal_mol": float(dG_sem) if len(dG_values) > 0 else None,
    "MMGBSA_dG_Std_kcal_mol": float(dG_std) if len(dG_values) > 0 else None,
    "MMGBSA_N_Frames": int(len(dG_values)),
}

# Save CSV
csv_path = OUTPUT_DIR / "md_summary.csv"
pd.DataFrame([summary]).to_csv(csv_path, index=False)

# Save JSON
json_path = OUTPUT_DIR / "md_summary.json"
with open(json_path, "w") as f:
    json_mod.dump(summary, f, indent=2, default=str)

print(f"✅ Summary CSV: {csv_path}")
print(f"✅ Summary JSON: {json_path}")

print(f"\n{'=' * 60}")
print(f"SUMMARY: {SYSTEM_NAME}")
print(f"{'=' * 60}")
for k, v in summary.items():
    if isinstance(v, float):
        print(f"  {k}: {v:.4f}")
    else:
        print(f"  {k}: {v}")

# %%
# ============================================================
# CELL 12: Package Results for Download
# ============================================================
import zipfile

print("Packaging results for download...")

zip_path = Path("/kaggle/working") / f"MD_{SYSTEM_NAME}_results.zip"

with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
    for f in sorted(OUTPUT_DIR.glob("*")):
        # Skip large binary files (DCD ~GB, checkpoints)
        if f.suffix in [".dcd", ".chk"]:
            size_mb = f.stat().st_size / 1e6
            print(f"  Skip (large): {f.name} ({size_mb:.1f} MB)")
            continue
        if f.name.startswith("_tmp"):
            continue
        if f.is_dir():
            continue
        zf.write(f, f"{SYSTEM_NAME}/{f.name}")
        print(f"  Added: {f.name}")

print(f"\n✅ Download: {zip_path}")
print(f"   Size: {zip_path.stat().st_size / 1e6:.1f} MB")

# %%
# ============================================================
# CELL 13: Final Report
# ============================================================
print("=" * 60)
print(f"📋 FINAL REPORT: {SYSTEM_NAME}")
print("=" * 60)
print(f"""
System: {LIGAND_NAME} bound to {RECEPTOR_PDB_ID} (AGTR1)
Simulation: {MD_PARAMS['production_ns']} ns @ {MD_PARAMS['temperature_K']}K

Docking:
  Re-dock affinity: {best_affinity:.2f} kcal/mol

Stability:
  Protein RMSD: {np.mean(rmsd_protein):.2f} ± {np.std(rmsd_protein):.2f} Å
  Ligand RMSD:  {np.mean(rmsd_ligand):.2f} ± {np.std(rmsd_ligand):.2f} Å
  Rg:           {np.mean(rg):.2f} ± {np.std(rg):.2f} Å
  SASA:         {np.mean(sasa_full):.0f} ± {np.std(sasa_full):.0f} Å²

Interactions:
  H-bonds:     {np.mean(hbond_counts):.1f} ± {np.std(hbond_counts):.1f}
  Top contacts: {', '.join([f'{rn}{rs}' for (rs, rn), _ in contact_freq_sorted[:5]])}

Binding Energy:
  MM-GBSA ΔG:  {dG_mean:.2f} ± {dG_sem:.2f} kcal/mol ({len(dG_values)} frames)

Files:
  Figures:  md_analysis_4panel.png, md_hbond_mmgbsa.png, md_extended_analysis.png
  Data:     md_summary.csv/json, timeseries.csv, rmsf_per_residue.csv
            contact_frequency.csv, key_distances.csv
  Frames:   frame_initial.pdb, frame_middle.pdb, frame_final.pdb, frame_equilibrated_mean.pdb
  Download: MD_{SYSTEM_NAME}_results.zip
""")
print("🎉 MD Simulation Complete!")
