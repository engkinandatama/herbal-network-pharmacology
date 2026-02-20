# %% [markdown]
# # MM-GBSA Binding Free Energy — Mangiferin–RELA (100 ns) — GPU
# Computes MM-GBSA ΔG_bind from 100 ns trajectory.
#
# **GPU REQUIRED**: Uses OpenMM for energy calculations on CUDA.
#
# Run this in **parallel** with the CPU analysis notebook.
#
# **Input**: Dataset `md-mangiferin-rela-run2-output`

# %%
# ============================================================
# CELL 0: Install Dependencies (needs OpenMM + CUDA)
# ============================================================
import subprocess, sys, os, glob, time, shutil

t_start = time.time()
PY_VER = f"{sys.version_info.major}.{sys.version_info.minor}"
print(f"Kernel Python: {PY_VER} ({sys.executable})")
print("=" * 60)
print("STEP 1: Installing Miniforge3...")
print("=" * 60)

MINIFORGE_DIR = "/tmp/miniforge"
MAMBA = f"{MINIFORGE_DIR}/bin/mamba"

if not os.path.exists(MAMBA):
    url = "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
    print("  Downloading...")
    subprocess.run(f"wget -q {url} -O /tmp/miniforge.sh", shell=True, check=True)
    r = subprocess.run(f"bash /tmp/miniforge.sh -b -p {MINIFORGE_DIR}",
                       shell=True, capture_output=True, text=True)
    os.remove("/tmp/miniforge.sh")
    print(f"  ✅ Miniforge3 installed ({time.time()-t_start:.0f}s)")

print("\n" + "=" * 60)
print("STEP 2: Mamba install...")
print("=" * 60)
# rdkit NOT installed via mamba — conda rdkit needs GLIBCXX_3.4.31 which Kaggle lacks
packages = f"python={PY_VER} numpy openmm cudatoolkit openmmforcefields openff-toolkit ambertools pdbfixer parmed"
print(f"  → {packages}")
r = subprocess.run(
    f"{MAMBA} install -y -c conda-forge {packages}",
    shell=True, capture_output=True, text=True, timeout=600
)
if r.returncode != 0:
    print(f"\n  ⚠️ mamba failed (exit {r.returncode})")
    if r.stderr:
        for line in r.stderr.strip().split('\n')[-10:]:
            print(f"    STDERR: {line}")
else:
    print(f"  ✅ Packages installed ({time.time()-t_start:.0f}s)")

# Remove conda rdkit+scipy (auto-installed as openff-toolkit deps).
# Conda rdkit compiled with newer GCC → needs GLIBCXX_3.4.31 → crashes on Kaggle.
# Solution: use pip versions (compiled against system libs).
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

subprocess.run(f"{MAMBA} clean -afy 2>/dev/null", shell=True,
               capture_output=True, text=True)

# Pip install rdkit+scipy (BEFORE module flush — compiled against system libs)
print("  Installing rdkit+scipy via pip...")
for pkg in ["rdkit", "scipy"]:
    r_pip = subprocess.run([sys.executable, "-m", "pip", "install", "-q", pkg],
                           capture_output=True, text=True, timeout=300)
    if r_pip.returncode == 0:
        print(f"  ✅ pip: {pkg}")
    else:
        print(f"  ❌ pip {pkg}: {r_pip.stderr[-200:]}")

print("\n" + "=" * 60)
print("STEP 3: Remove system numpy & patch paths...")
print("=" * 60)

subprocess.run([sys.executable, "-m", "pip", "uninstall", "-y", "numpy"],
               capture_output=True, text=True)
for np_dir in glob.glob("/usr/local/lib/python3.12/dist-packages/numpy*"):
    shutil.rmtree(np_dir, ignore_errors=True)
    print(f"  ✅ System numpy removed")

mods_to_flush = [k for k in sys.modules if any(x in k.lower() for x in
                 ['numpy', 'openmm', 'mdtraj', 'openff', 'parmed', 'rdkit'])]
for mod in mods_to_flush:
    del sys.modules[mod]
print(f"  ✅ Flushed {len(mods_to_flush)} cached modules")

for sp in glob.glob(f"{MINIFORGE_DIR}/lib/python*/site-packages"):
    if sp not in sys.path:
        sys.path.insert(0, sp)
        print(f"  ✅ sys.path += {sp}")

os.environ["PATH"] = f"{MINIFORGE_DIR}/bin:" + os.environ.get("PATH", "")
os.environ["LD_LIBRARY_PATH"] = f"{MINIFORGE_DIR}/lib:" + os.environ.get("LD_LIBRARY_PATH", "")
os.environ["AMBERHOME"] = MINIFORGE_DIR
print("  ✅ Environment set")

print("\n" + "=" * 60)
print("STEP 4: Verifying...")
print("=" * 60)
try:
    import numpy as np; print(f"  ✅ numpy {np.__version__}")
except Exception as e:
    print(f"  ❌ numpy FAILED: {e}")
try:
    import openmm; print(f"  ✅ OpenMM {openmm.__version__}")
    plats = [openmm.Platform.getPlatform(i).getName() for i in range(openmm.Platform.getNumPlatforms())]
    print(f"  ✅ CUDA: {'ready' if 'CUDA' in plats else 'NOT FOUND — ' + str(plats)}")
except Exception as e:
    print(f"  ❌ OpenMM FAILED: {e}")
try:
    import openff.toolkit; print(f"  ✅ openff.toolkit")
except Exception as e:
    print(f"  ❌ openff FAILED: {e}")
try:
    import openmmforcefields; print(f"  ✅ openmmforcefields")
except Exception as e:
    print(f"  ❌ openmmforcefields FAILED: {e}")
try:
    import mdtraj as md; print(f"  ✅ mdtraj {md.__version__}")
except Exception as e:
    print(f"  ❌ mdtraj FAILED: {e}")

print(f"\n  Setup: {time.time()-t_start:.0f}s")
print("=" * 60)

# %%
# ============================================================
# CELL 1: Load Trajectory & Prepare
# ============================================================
import numpy as np
import mdtraj as md
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

SYSTEM_NAME = "Mangiferin_RELA"
LIGAND_NAME = "Mangiferin"
RECEPTOR_PDB_ID = "1NFI"
SAVE_INTERVAL_PS = 10

# --- Find dataset (robust search) ---
print("Searching for dataset...")
INPUT_ROOT = Path("/kaggle/input")
RESULTS_DIR = None

# Strategy 1: Check known dataset slugs
for name in ["md-mangiferin-rela-run2-output", "md-mangiferin-rela-output"]:
    for pattern in [
        INPUT_ROOT / name / "md_results" / SYSTEM_NAME,
        INPUT_ROOT / name / SYSTEM_NAME,
        INPUT_ROOT / name / "Mangiferin_RELA",
    ]:
        if pattern.exists() and list(pattern.glob("production_run*.dcd")):
            RESULTS_DIR = pattern
            break
    if RESULTS_DIR:
        break

# Strategy 2: Recursive search for DCD files
if RESULTS_DIR is None:
    print("  Known slugs not found, searching recursively...")
    for dcd in INPUT_ROOT.rglob("production_run1.dcd"):
        RESULTS_DIR = dcd.parent
        break

# Strategy 3: Search for Mangiferin_RELA folder
if RESULTS_DIR is None:
    for d in INPUT_ROOT.rglob("Mangiferin_RELA"):
        if d.is_dir():
            RESULTS_DIR = d
            break

if RESULTS_DIR is None:
    print("  ❌ Dataset not found! Listing /kaggle/input contents:")
    for p in sorted(INPUT_ROOT.rglob("*")):
        if p.is_file() and p.stat().st_size > 1e6:
            print(f"    {p.relative_to(INPUT_ROOT)} ({p.stat().st_size/1e6:.1f} MB)")
        elif p.is_dir():
            print(f"    {p.relative_to(INPUT_ROOT)}/ (dir)")
    raise FileNotFoundError("Cannot find DCD files in /kaggle/input!")

OUTPUT_DIR = Path("/kaggle/working/mmgbsa") / SYSTEM_NAME
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

print("=" * 60)
print(f"MM-GBSA: {SYSTEM_NAME}")
print(f"Dataset: {RESULTS_DIR}")
print("=" * 60)

# --- Topology ---
topo_file = RESULTS_DIR / "complex_solvated.pdb"
if not topo_file.exists():
    topo_file = RESULTS_DIR / "complex_equilibrated.pdb"

# --- Load all DCDs with stride ---
STRIDE = 2
EFFECTIVE_DT_PS = SAVE_INTERVAL_PS * STRIDE  # 20 ps

dcd_files = sorted(RESULTS_DIR.glob("production_run*.dcd"))
print(f"\nFound {len(dcd_files)} DCD files (stride={STRIDE}):")

import gc
trajs = []
for dcd in dcd_files:
    t = md.load(str(dcd), top=str(topo_file), stride=STRIDE)
    sz_gb = dcd.stat().st_size / 1e9
    print(f"  {dcd.name}: {t.n_frames} frames ({sz_gb:.2f} GB on disk)")
    trajs.append(t)
    gc.collect()

# --- Override time axis ---
if len(trajs) > 0:
    dt_ps = EFFECTIVE_DT_PS
    trajs[0].time = np.arange(trajs[0].n_frames) * dt_ps
    print(f"  {dcd_files[0].name} time → 0-{trajs[0].time[-1]/1000:.1f} ns")
    for i in range(1, len(trajs)):
        t_start_ps = trajs[i-1].time[-1] + dt_ps
        trajs[i].time = np.arange(trajs[i].n_frames) * dt_ps + t_start_ps
        print(f"  {dcd_files[i].name} time → {trajs[i].time[0]/1000:.1f}-{trajs[i].time[-1]/1000:.1f} ns")

if len(trajs) > 1:
    traj = md.join(trajs)
else:
    traj = trajs[0]

del trajs
gc.collect()

# --- Trim to 100 ns ---
TARGET_NS = 100.0
mask = traj.time <= TARGET_NS * 1000
if not np.all(mask):
    n_before = traj.n_frames
    traj = traj[mask]
    gc.collect()
    print(f"  ✂️  Trimmed {n_before} → {traj.n_frames} frames")

total_ns = (traj.time[-1] - traj.time[0]) / 1000
mem_gb = traj.n_frames * traj.n_atoms * 3 * 4 / 1e9
print(f"\n📊 Combined: {traj.n_frames:,} frames, {total_ns:.1f} ns (~{mem_gb:.1f} GB)")

# --- Atom selections ---
protein_atoms = traj.topology.select("protein")
ca_atoms = traj.topology.select("protein and name CA")
water_atoms = set(traj.topology.select("water"))
try:
    ion_atoms = set(traj.topology.select("resname NA CL"))
except:
    ion_atoms = set()

all_known = set(protein_atoms) | water_atoms | ion_atoms
ligand_atoms = np.array(sorted([i for i in range(traj.n_atoms) if i not in all_known]))

print(f"  Protein: {len(protein_atoms)} atoms")
print(f"  Ligand: {len(ligand_atoms)} atoms")
assert len(ligand_atoms) > 0, "No ligand atoms!"

# --- Find ligand SDF ---
sdf_candidates = [
    Path(f"/kaggle/input/mahkota-dewa-docking/ligands/{LIGAND_NAME}.sdf"),
    Path(f"/kaggle/input/mahkota-dewa-docking/{LIGAND_NAME}.sdf"),
    RESULTS_DIR / f"{LIGAND_NAME}.sdf",
]
LIGAND_SDF = None
for p in sdf_candidates:
    if p.exists():
        LIGAND_SDF = p
        break

if LIGAND_SDF is None:
    # Search recursively
    for p in Path("/kaggle/input").rglob(f"{LIGAND_NAME}.sdf"):
        LIGAND_SDF = p
        break

assert LIGAND_SDF is not None, f"Cannot find {LIGAND_NAME}.sdf!"
print(f"  Ligand SDF: {LIGAND_SDF}")

# %%
# ============================================================
# CELL 2: MM-GBSA Binding Free Energy (GPU)
# ============================================================
from openmm.app import ForceField, PDBFile, NoCutoff, HBonds, Simulation
from openmm import LangevinMiddleIntegrator, Platform
import openmm.unit as unit
from openff.toolkit.topology import Molecule as OFFMolecule
from openmmforcefields.generators import GAFFTemplateGenerator

print("=" * 60)
print("MM-GBSA BINDING FREE ENERGY")
print("=" * 60)

# --- Load ligand molecule (FIX for UNL residue bug) ---
print("\nLoading ligand molecule for GAFF template...")
off_mol = OFFMolecule.from_file(str(LIGAND_SDF), allow_undefined_stereo=True)
print(f"  ✅ Loaded {LIGAND_NAME}: {off_mol.n_atoms} atoms")

# --- Setup implicit solvent FF with GAFF ---
print("Setting up ForceField + GAFF2...")
ff_gb = ForceField("amber/ff14SB.xml", "implicit/obc2.xml")
gaff_gen = GAFFTemplateGenerator(molecules=[off_mol], forcefield="gaff-2.11")
ff_gb.registerTemplateGenerator(gaff_gen.generator)
print("  ✅ GAFF2 template registered (fixes UNL residue bug)")

# --- Interior dielectric screening (ε=4.0) ---
# Mangiferin is highly polar (11 O, 7 OH groups) causing GB solvation
# overestimation with ε=1. Using ε=4 reduces Coulomb by 4x and adjusts
# GB polar solvation, giving more realistic binding free energies.
import math
SOLUTE_DIELECTRIC = 4.0
print(f"  ⚡ Interior dielectric: ε = {SOLUTE_DIELECTRIC} (polar ligand screening)")

# --- Select frames ---
N_FRAMES = 100
start_frame = int(traj.n_frames * 0.2)  # Skip first 20% (equilibration)
frame_indices = np.linspace(start_frame, traj.n_frames - 1, N_FRAMES, dtype=int)
print(f"\n  Analyzing {N_FRAMES} frames (from frame {start_frame} to {traj.n_frames-1})")
print(f"  Time range: {traj.time[start_frame]/1000:.1f} - {traj.time[-1]/1000:.1f} ns")

# --- Prepare atom indices ---
solute_atoms = np.concatenate([protein_atoms, ligand_atoms])
prot_local_idx = np.arange(len(protein_atoms))
lig_local_idx = np.arange(len(protein_atoms), len(protein_atoms) + len(ligand_atoms))

# --- Write reference PDBs ---
ref_frame = traj[frame_indices[0]]
ref_solute = ref_frame.atom_slice(solute_atoms)
ref_receptor = ref_frame.atom_slice(protein_atoms)

tmp_dir = OUTPUT_DIR / "_tmp_mmgbsa"
tmp_dir.mkdir(exist_ok=True)

ref_receptor.save(str(tmp_dir / "receptor.pdb"))
print("  ✅ Reference receptor PDB written")

# --- Build ligand topology from OpenFF molecule (correct bonds from SDF) ---
# CRITICAL: Do NOT use MDTraj to save ligand or complex PDB!
# MDTraj guesses bonds from distances, which for complex molecules like
# Mangiferin produces WRONG bond topology → GAFF assigns different charges
# → ΔG offset of -416,000 kcal/mol instead of -27 kcal/mol.
def _build_ligand_openmm_topology(off_mol):
    """Create OpenMM Topology for ligand from OpenFF molecule.
    This preserves the exact bond structure from the SDF file."""
    from openmm.app import Topology, Element
    topo = Topology()
    chain = topo.addChain()
    res = topo.addResidue("UNL", chain)
    omm_atoms = []
    for atom in off_mol.atoms:
        elem = Element.getByAtomicNumber(atom.atomic_number)
        omm_atoms.append(topo.addAtom(atom.name if atom.name else f"{atom.symbol}{atom.molecule_atom_index}",
                                       elem, res))
    for bond in off_mol.bonds:
        topo.addBond(omm_atoms[bond.atom1_index], omm_atoms[bond.atom2_index])
    return topo

ligand_omm_topology = _build_ligand_openmm_topology(off_mol)
print(f"  ✅ Ligand OpenMM topology built from SDF ({off_mol.n_atoms} atoms, {off_mol.n_bonds} bonds)")

# --- Build COMPLEX topology by combining receptor PDB + OpenFF ligand ---
# This ensures ligand bonds in complex match standalone ligand exactly.
def _build_combined_topology(receptor_pdb_path, off_mol):
    """Build combined OpenMM topology: receptor (from PDB) + ligand (from SDF).
    
    This is the KEY FIX: instead of saving complex.pdb via MDTraj (which
    guesses ligand bonds from distances — WRONG for Mangiferin), we combine:
    - Receptor topology from PDB (correct protein bonds)
    - Ligand topology from OpenFF molecule (correct bonds from SDF)
    
    This guarantees GAFF assigns identical charges to ligand atoms in both
    the complex and standalone ligand systems."""
    from openmm.app import PDBFile, Topology, Element
    
    rec_pdb = PDBFile(str(receptor_pdb_path))
    topo = Topology()
    
    # Copy receptor topology (chains, residues, atoms, bonds)
    atom_map = {}
    for chain in rec_pdb.topology.chains():
        new_chain = topo.addChain(chain.id)
        for residue in chain.residues():
            new_res = topo.addResidue(residue.name, new_chain, residue.id)
            for atom in residue.atoms():
                new_atom = topo.addAtom(atom.name, atom.element, new_res)
                atom_map[atom] = new_atom
    for bond in rec_pdb.topology.bonds():
        topo.addBond(atom_map[bond[0]], atom_map[bond[1]])
    
    # Add ligand from OpenFF molecule (correct bonds from SDF)
    lig_chain = topo.addChain()
    lig_res = topo.addResidue("UNL", lig_chain)
    lig_atoms = []
    for atom in off_mol.atoms:
        elem = Element.getByAtomicNumber(atom.atomic_number)
        name = atom.name if atom.name else f"{atom.symbol}{atom.molecule_atom_index}"
        lig_atoms.append(topo.addAtom(name, elem, lig_res))
    for bond in off_mol.bonds:
        topo.addBond(lig_atoms[bond.atom1_index], lig_atoms[bond.atom2_index])
    
    n_rec = sum(1 for _ in rec_pdb.topology.atoms())
    n_lig = off_mol.n_atoms
    print(f"  ✅ Combined topology: {n_rec} receptor + {n_lig} ligand = {n_rec + n_lig} atoms")
    return topo

combined_complex_topology = _build_combined_topology(tmp_dir / "receptor.pdb", off_mol)

# --- Build GB systems with force groups for energy decomposition ---
def _build_gb_system(topology, name):
    """Create System, assign force groups, and build Simulation."""
    sys_gb = ff_gb.createSystem(topology, nonbondedMethod=NoCutoff, constraints=HBonds)
    
    # Apply interior dielectric screening (ε=4.0)
    # 1. NonbondedForce: scale charges by 1/√ε → Coulomb divided by ε, VdW unchanged
    # 2. CustomGBForce: set soluteDielectric → GB solvation uses (1/ε - 1/78.5)
    from openmm import NonbondedForce, CustomGBForce as CGBForce
    charge_scale = 1.0 / math.sqrt(SOLUTE_DIELECTRIC)
    
    for force in sys_gb.getForces():
        if isinstance(force, NonbondedForce):
            for idx in range(force.getNumParticles()):
                q, sig, eps = force.getParticleParameters(idx)
                force.setParticleParameters(idx, q * charge_scale, sig, eps)
            for idx in range(force.getNumExceptions()):
                p1, p2, cprod, sig, eps = force.getExceptionParameters(idx)
                force.setExceptionParameters(idx, p1, p2, cprod * charge_scale * charge_scale, sig, eps)
        elif isinstance(force, CGBForce):
            for gp_idx in range(force.getNumGlobalParameters()):
                if force.getGlobalParameterName(gp_idx) == "soluteDielectric":
                    force.setGlobalParameterDefaultValue(gp_idx, SOLUTE_DIELECTRIC)
    
    # Assign force groups for energy decomposition
    force_map = {}
    for i, force in enumerate(sys_gb.getForces()):
        fname = force.__class__.__name__
        force.setForceGroup(i)
        force_map[i] = fname
    integ = LangevinMiddleIntegrator(300*unit.kelvin, 1/unit.picoseconds, 2*unit.femtoseconds)
    try:
        platform = Platform.getPlatformByName("CUDA")
        props = {"Precision": "mixed"}
        sim = Simulation(topology, sys_gb, integ, platform, props)
        print(f"  → {name}: CUDA ✅")
    except:
        sim = Simulation(topology, sys_gb, integ)
        print(f"  → {name}: CPU (CUDA not available)")
    return sim, force_map

def get_energy_components(sim, positions, force_map):
    """Get total + decomposed energies using force groups."""
    sim.context.setPositions(positions)
    E_total = sim.context.getState(getEnergy=True).getPotentialEnergy()
    E_total = E_total.value_in_unit(unit.kilocalories_per_mole)
    
    # Decompose by force group
    E_gas = 0.0  # gas-phase (bonded + nonbonded)
    E_solv = 0.0  # solvation (GB)
    E_vdw_ele = 0.0  # NonbondedForce specifically
    
    for grp, fname in force_map.items():
        E_grp = sim.context.getState(getEnergy=True, groups={grp}).getPotentialEnergy()
        E_grp = E_grp.value_in_unit(unit.kilocalories_per_mole)
        if "GB" in fname or "CustomGB" in fname or "GBSA" in fname:
            E_solv += E_grp
        else:
            E_gas += E_grp
            if "Nonbonded" in fname:
                E_vdw_ele = E_grp
    
    return E_total, E_gas, E_solv, E_vdw_ele

print("\nBuilding 3 GB systems (complex, receptor, ligand)...")
try:
    # Complex: built from combined topology (receptor PDB + OpenFF ligand) — KEY FIX
    sim_complex, fmap_complex = _build_gb_system(combined_complex_topology, "complex (combined)")
    # Receptor: built from receptor PDB
    rec_pdb_obj = PDBFile(str(tmp_dir / "receptor.pdb"))
    sim_receptor, fmap_receptor = _build_gb_system(rec_pdb_obj.topology, "receptor")
    # Ligand: built from OpenFF topology (correct bonds from SDF)
    sim_ligand, fmap_ligand = _build_gb_system(ligand_omm_topology, "ligand (OpenFF)")
    print("  ✅ All three systems built (complex uses combined topology)")
    print(f"  Force groups (complex): {fmap_complex}")
    mmgbsa_ready = True
except Exception as e:
    print(f"  ❌ System build failed: {e}")
    import traceback
    traceback.print_exc()
    mmgbsa_ready = False

# --- Compute energies with decomposition ---
results = []  # list of dicts per frame
errors = 0

if mmgbsa_ready:
    print(f"\n🚀 Computing {N_FRAMES} frames with energy decomposition...")
    t0 = time.time()

    for i, fidx in enumerate(frame_indices):
        try:
            frame = traj[fidx]
            frame_solute = frame.atom_slice(solute_atoms)
            pos_complex = frame_solute.xyz[0]
            pos_receptor = pos_complex[prot_local_idx]
            pos_ligand = pos_complex[lig_local_idx]

            E_tot_c, E_gas_c, E_solv_c, E_nb_c = get_energy_components(sim_complex, pos_complex, fmap_complex)
            E_tot_r, E_gas_r, E_solv_r, E_nb_r = get_energy_components(sim_receptor, pos_receptor, fmap_receptor)
            E_tot_l, E_gas_l, E_solv_l, E_nb_l = get_energy_components(sim_ligand, pos_ligand, fmap_ligand)

            dG_total = E_tot_c - E_tot_r - E_tot_l
            dE_gas = E_gas_c - E_gas_r - E_gas_l  # gas-phase (VdW + Ele + bonded)
            dG_solv = E_solv_c - E_solv_r - E_solv_l  # solvation (GB polar)
            dE_vdw_ele = E_nb_c - E_nb_r - E_nb_l  # NonbondedForce (VdW + Ele)

            results.append({
                "frame_idx": int(fidx),
                "time_ns": float(traj.time[fidx] / 1000),
                "dG_total": dG_total,
                "dE_gas": dE_gas,
                "dG_solv": dG_solv,
                "dE_vdw_ele": dE_vdw_ele,
            })
        except Exception as e:
            errors += 1
            if errors <= 5:
                print(f"    ⚠️ Error at frame {fidx}: {e}")

        if (i + 1) % 10 == 0:
            elapsed = time.time() - t0
            eta = elapsed / (i + 1) * (N_FRAMES - i - 1)
            n_ok = len(results)
            if n_ok > 0:
                last_dG = results[-1]["dG_total"]
                print(f"  [{i+1}/{N_FRAMES}] {n_ok} OK, {errors} err | last ΔG={last_dG:.1f} | ETA: {eta:.0f}s")
            else:
                print(f"  [{i+1}/{N_FRAMES}] {n_ok} OK, {errors} err | ETA: {eta:.0f}s")

    total_time = time.time() - t0
    print(f"\n  ⏱️  MM-GBSA completed in {total_time:.1f}s ({total_time/60:.1f} min)")

# --- Cleanup ---
import shutil
if tmp_dir.exists():
    shutil.rmtree(tmp_dir)

# Free GPU/simulation memory before plotting
try: del sim_complex, sim_receptor, sim_ligand
except: pass
try: del traj
except: pass
gc.collect()
print("  ✅ GPU memory freed")

# --- Results ---
import pandas as pd
results_df = pd.DataFrame(results)

if len(results_df) > 0:
    dG_total = results_df["dG_total"].values
    dE_gas = results_df["dE_gas"].values
    dG_solv = results_df["dG_solv"].values
    dE_vdw_ele = results_df["dE_vdw_ele"].values

    # Statistics
    stats = {}
    for col in ["dG_total", "dE_gas", "dG_solv", "dE_vdw_ele"]:
        vals = results_df[col].values
        stats[col] = {
            "mean": float(np.mean(vals)),
            "std": float(np.std(vals)),
            "sem": float(np.std(vals) / np.sqrt(len(vals))),
        }

    print(f"\n{'=' * 60}")
    print(f"📊 MM-GBSA RESULTS ({len(results_df)}/{N_FRAMES} frames)")
    print(f"{'=' * 60}")
    print(f"  ΔG_total   = {stats['dG_total']['mean']:8.2f} ± {stats['dG_total']['sem']:.2f} kcal/mol")
    print(f"  ΔE_gas     = {stats['dE_gas']['mean']:8.2f} ± {stats['dE_gas']['sem']:.2f} kcal/mol  (VdW + Ele + bonded)")
    print(f"  ΔE_vdw+ele = {stats['dE_vdw_ele']['mean']:8.2f} ± {stats['dE_vdw_ele']['sem']:.2f} kcal/mol  (NonbondedForce only)")
    print(f"  ΔG_solv    = {stats['dG_solv']['mean']:8.2f} ± {stats['dG_solv']['sem']:.2f} kcal/mol  (GB solvation)")
    print(f"  Errors:    {errors}")
    print(f"{'=' * 60}")


    # --- Save CSV/JSON immediately (BEFORE cell boundary to avoid OOM kernel death) ---
    import json as json_mod
    
    results_df.to_csv(OUTPUT_DIR / "mmgbsa_per_frame.csv", index=False)
    
    summary = {
        "System": SYSTEM_NAME,
        "Receptor": RECEPTOR_PDB_ID,
        "Ligand": LIGAND_NAME,
        "Simulation_Time_ns": float(total_ns),
        "N_Frames_Analyzed": int(len(results_df)),
        "N_Frames_Errors": int(errors),
        "Equilibration_Skip_Percent": 20,
        "Stride": STRIDE,
        "dG_total_mean": stats["dG_total"]["mean"],
        "dG_total_std": stats["dG_total"]["std"],
        "dG_total_sem": stats["dG_total"]["sem"],
        "dG_total_min": float(np.min(dG_total)),
        "dG_total_max": float(np.max(dG_total)),
        "dE_gas_mean": stats["dE_gas"]["mean"],
        "dE_gas_std": stats["dE_gas"]["std"],
        "dE_gas_sem": stats["dE_gas"]["sem"],
        "dE_vdw_ele_mean": stats["dE_vdw_ele"]["mean"],
        "dE_vdw_ele_std": stats["dE_vdw_ele"]["std"],
        "dE_vdw_ele_sem": stats["dE_vdw_ele"]["sem"],
        "dG_solv_mean": stats["dG_solv"]["mean"],
        "dG_solv_std": stats["dG_solv"]["std"],
        "dG_solv_sem": stats["dG_solv"]["sem"],
    }
    
    with open(OUTPUT_DIR / "mmgbsa_summary.json", "w") as f:
        json_mod.dump(summary, f, indent=2, default=str)
    pd.DataFrame([summary]).to_csv(OUTPUT_DIR / "mmgbsa_summary.csv", index=False)
    print("✅ CSV and JSON saved")
    
    # --- Final Report ---
    print(f"\n{'=' * 60}")
    print(f"📋 MM-GBSA FINAL REPORT: {SYSTEM_NAME}")
    print(f"{'=' * 60}")
    print(f"""
System: {LIGAND_NAME} bound to {RECEPTOR_PDB_ID} (RELA/NF-κB p65)
Simulation: {total_ns:.0f} ns | Stride: {STRIDE} | Frames analyzed: {len(results_df)}

MM-GBSA Energy Decomposition (kcal/mol):
  ΔE_gas     = {stats['dE_gas']['mean']:8.2f} ± {stats['dE_gas']['sem']:.2f}  (gas-phase: VdW + Ele + bonded)
  ΔE_vdw+ele = {stats['dE_vdw_ele']['mean']:8.2f} ± {stats['dE_vdw_ele']['sem']:.2f}  (NonbondedForce only)
  ΔG_solv    = {stats['dG_solv']['mean']:8.2f} ± {stats['dG_solv']['sem']:.2f}  (GB polar solvation)
  ─────────────────────────────────────
  ΔG_total   = {stats['dG_total']['mean']:8.2f} ± {stats['dG_total']['sem']:.2f}

Interpretation:
  {'✅ FAVORABLE binding (ΔG < -10 kcal/mol)' if stats['dG_total']['mean'] < -10 else '⚠️ MODERATE binding affinity' if stats['dG_total']['mean'] < 0 else '❌ UNFAVORABLE binding (may indicate GB solvation overestimation for highly polar ligand)'}
  Gas-phase {'dominates' if abs(stats['dE_gas']['mean']) > abs(stats['dG_solv']['mean']) else 'is offset by'} solvation penalty

NOTE: Plots will be generated locally from downloaded CSV files.

Output files:
  - mmgbsa_per_frame.csv (with decomposition)
  - mmgbsa_summary.json
  - mmgbsa_summary.csv
""")
else:
    print("❌ No MM-GBSA results to save.")

print(f"\n⏱️  Total time: {(time.time()-t_start)/60:.1f} min")
print("🎉 MM-GBSA Analysis Complete!")

