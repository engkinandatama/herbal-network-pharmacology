# %% [markdown]
# # MM-GBSA Binding Free Energy — Phalerin–AGTR1 (100 ns) — GPU
# Computes MM-GBSA ΔG_bind from 100 ns trajectory.
#
# **GPU REQUIRED**: Uses OpenMM for energy calculations on CUDA.
#
# Run this in **parallel** with the CPU analysis notebook.
#
# **Input**: Dataset `md-phalerin-agtr1-final-output`

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

SYSTEM_NAME = "Phalerin_AGTR1"
LIGAND_NAME = "Phalerin"
RECEPTOR_PDB_ID = "6OS1"
SAVE_INTERVAL_PS = 10

# --- Find dataset (robust search) ---
print("Searching for dataset...")
INPUT_ROOT = Path("/kaggle/input")
RESULTS_DIR = None

# Strategy 1: Check known dataset slugs
for name in ["md-phalerin-agtr1-final-output", "md-phalerin-agtr1-run3"]:
    for pattern in [
        INPUT_ROOT / name / "md_results" / SYSTEM_NAME,
        INPUT_ROOT / name / SYSTEM_NAME,
        INPUT_ROOT / name / "Phalerin_AGTR1",
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

# Strategy 3: Search for Phalerin_AGTR1 folder
if RESULTS_DIR is None:
    for d in INPUT_ROOT.rglob("Phalerin_AGTR1"):
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
    RESULTS_DIR / f"{LIGAND_NAME}.sdf",
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

ref_solute.save(str(tmp_dir / "complex.pdb"))
ref_receptor.save(str(tmp_dir / "receptor.pdb"))
# NOTE: Do NOT save ligand.pdb via MDTraj — it guesses bonds from distances
# which produce bond topology different from the GAFF template (from SDF).
# Instead, we build the ligand system directly from OpenFF molecule topology.
print("  ✅ Reference PDBs written (complex + receptor)")

# --- Build ligand topology from OpenFF molecule (correct bonds from SDF) ---
def _build_ligand_openmm_topology(off_mol):
    """Create OpenMM Topology for ligand from OpenFF molecule.
    This preserves the exact bond structure from the SDF file,
    avoiding MDTraj's distance-based bond guessing that causes UNL errors."""
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

# --- Build GB systems with force groups for energy decomposition ---
def build_gb_system_from_pdb(pdb_path):
    """Build GB system from a PDB file (for complex and receptor)."""
    pdb_obj = PDBFile(str(pdb_path))
    return _build_gb_system(pdb_obj.topology, pdb_path.name)

def build_gb_system_from_topology(topology, name="ligand"):
    """Build GB system from an OpenMM topology (for ligand with correct bonds)."""
    return _build_gb_system(topology, name)

def _build_gb_system(topology, name):
    """Shared: create System, assign force groups, and build Simulation."""
    sys_gb = ff_gb.createSystem(topology, nonbondedMethod=NoCutoff, constraints=HBonds)
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
    sim_complex, fmap_complex = build_gb_system_from_pdb(tmp_dir / "complex.pdb")
    sim_receptor, fmap_receptor = build_gb_system_from_pdb(tmp_dir / "receptor.pdb")
    # Build ligand system from OpenFF topology (NOT from MDTraj-saved PDB)
    sim_ligand, fmap_ligand = build_gb_system_from_topology(ligand_omm_topology, "ligand (OpenFF)")
    print("  ✅ All three systems built")
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
else:
    print("\n⚠️ MM-GBSA failed for all frames.")
    stats = {}

# %%
# ============================================================
# CELL 3: Save Results & Plot
# ============================================================
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import json as json_mod

if len(results_df) > 0:
    # --- Save CSV (per-frame with decomposition) ---
    results_df.to_csv(OUTPUT_DIR / "mmgbsa_per_frame.csv", index=False)

    # --- Save summary JSON ---
    summary = {
        "System": SYSTEM_NAME,
        "Receptor": RECEPTOR_PDB_ID,
        "Ligand": LIGAND_NAME,
        "Simulation_Time_ns": float(total_ns),
        "N_Frames_Analyzed": int(len(results_df)),
        "N_Frames_Errors": int(errors),
        "Equilibration_Skip_Percent": 20,
        "Stride": STRIDE,
        # Total binding free energy
        "dG_total_mean": stats["dG_total"]["mean"],
        "dG_total_std": stats["dG_total"]["std"],
        "dG_total_sem": stats["dG_total"]["sem"],
        "dG_total_min": float(np.min(dG_total)),
        "dG_total_max": float(np.max(dG_total)),
        # Energy decomposition
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

    # --- Plot: 2x2 (Histogram, Time series, Components bar, Component time) ---
    plt.rcParams.update({
        "font.family": "sans-serif", "font.size": 10,
        "axes.labelsize": 12, "axes.titlesize": 13,
        "figure.dpi": 300, "axes.linewidth": 1.2,
    })

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # (A) Histogram of total ΔG
    ax = axes[0, 0]
    dG_m = stats["dG_total"]["mean"]
    dG_s = stats["dG_total"]["sem"]
    dG_sd = stats["dG_total"]["std"]
    ax.hist(dG_total, bins=25, color="#00BCD4", edgecolor="black", alpha=0.7)
    ax.axvline(dG_m, color="red", linewidth=2, linestyle="--",
              label=f"Mean: {dG_m:.2f} ± {dG_s:.2f}")
    ax.axvline(dG_m - dG_sd, color="orange", linewidth=1, linestyle=":")
    ax.axvline(dG_m + dG_sd, color="orange", linewidth=1, linestyle=":",
              label=f"±1σ: {dG_sd:.2f}")
    ax.set_xlabel("ΔG_bind (kcal/mol)")
    ax.set_ylabel("Count")
    ax.set_title(f"(A) MM-GBSA Distribution\nΔG = {dG_m:.2f} ± {dG_s:.2f} kcal/mol")
    ax.legend(fontsize=8)

    # (B) Time series
    ax = axes[0, 1]
    t_ns = results_df["time_ns"].values
    ax.scatter(t_ns, dG_total, s=10, alpha=0.5, color="#3F51B5", label="Per-frame")
    if len(dG_total) >= 10:
        w = max(1, len(dG_total) // 10)
        dg_smooth = np.convolve(dG_total, np.ones(w)/w, mode="valid")
        ax.plot(t_ns[:len(dg_smooth)], dg_smooth, color="#E91E63", linewidth=2,
                label=f"Running avg (w={w})")
    ax.axhline(dG_m, color="red", linestyle="--", alpha=0.5)
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("ΔG_bind (kcal/mol)")
    ax.set_title("(B) MM-GBSA over Time")
    ax.legend(fontsize=8)

    # (C) Energy decomposition bar chart
    ax = axes[1, 0]
    components = ["ΔE_gas\n(VdW+Ele+Bonded)", "ΔG_solv\n(GB Polar)", "ΔG_total"]
    means = [stats["dE_gas"]["mean"], stats["dG_solv"]["mean"], stats["dG_total"]["mean"]]
    sems = [stats["dE_gas"]["sem"], stats["dG_solv"]["sem"], stats["dG_total"]["sem"]]
    colors = ["#2196F3" if m < 0 else "#FF5722" for m in means]
    bars = ax.bar(components, means, yerr=sems, capsize=5, color=colors,
                  edgecolor="black", linewidth=0.5, alpha=0.85)
    ax.axhline(0, color="black", linewidth=0.5)
    for bar, m, s in zip(bars, means, sems):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height(),
                f"{m:.1f}±{s:.1f}", ha="center", va="bottom" if m >= 0 else "top",
                fontsize=9, fontweight="bold")
    ax.set_ylabel("Energy (kcal/mol)")
    ax.set_title("(C) Energy Decomposition")

    # (D) Component time series
    ax = axes[1, 1]
    ax.scatter(t_ns, dE_gas, s=8, alpha=0.3, color="#2196F3", label="ΔE_gas")
    ax.scatter(t_ns, dG_solv, s=8, alpha=0.3, color="#FF9800", label="ΔG_solv")
    if len(dG_total) >= 10:
        gas_smooth = np.convolve(dE_gas, np.ones(w)/w, mode="valid")
        solv_smooth = np.convolve(dG_solv, np.ones(w)/w, mode="valid")
        ax.plot(t_ns[:len(gas_smooth)], gas_smooth, color="#1565C0", linewidth=2)
        ax.plot(t_ns[:len(solv_smooth)], solv_smooth, color="#E65100", linewidth=2)
    ax.axhline(0, color="black", linewidth=0.5, linestyle=":")
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Energy (kcal/mol)")
    ax.set_title("(D) Energy Components over Time")
    ax.legend(fontsize=8)

    fig.suptitle(f"MM-GBSA Binding Free Energy: {LIGAND_NAME} – {RECEPTOR_PDB_ID} (AGTR1)",
                 fontsize=14, fontweight="bold")
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "mmgbsa_plot.png", dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()
    print("✅ Plot saved (2x2 with decomposition)")

    # --- Final Report ---
    print(f"\n{'=' * 60}")
    print(f"📋 MM-GBSA FINAL REPORT: {SYSTEM_NAME}")
    print(f"{'=' * 60}")
    print(f"""
System: {LIGAND_NAME} bound to {RECEPTOR_PDB_ID} (AGTR1)
Simulation: {total_ns:.0f} ns | Stride: {STRIDE} | Frames analyzed: {len(results_df)}

MM-GBSA Energy Decomposition (kcal/mol):
  ΔE_gas     = {stats['dE_gas']['mean']:8.2f} ± {stats['dE_gas']['sem']:.2f}  (gas-phase: VdW + Ele + bonded)
  ΔE_vdw+ele = {stats['dE_vdw_ele']['mean']:8.2f} ± {stats['dE_vdw_ele']['sem']:.2f}  (NonbondedForce only)
  ΔG_solv    = {stats['dG_solv']['mean']:8.2f} ± {stats['dG_solv']['sem']:.2f}  (GB polar solvation)
  ─────────────────────────────────────
  ΔG_total   = {stats['dG_total']['mean']:8.2f} ± {stats['dG_total']['sem']:.2f}

Interpretation:
  {'✅ FAVORABLE binding (ΔG < -10 kcal/mol)' if stats['dG_total']['mean'] < -10 else '⚠️ MODERATE binding affinity' if stats['dG_total']['mean'] < 0 else '❌ UNFAVORABLE binding'}
  Gas-phase {'dominates' if abs(stats['dE_gas']['mean']) > abs(stats['dG_solv']['mean']) else 'is offset by'} solvation penalty

Output files:
  - mmgbsa_per_frame.csv (with decomposition)
  - mmgbsa_summary.json
  - mmgbsa_summary.csv
  - mmgbsa_plot.png (2x2 panel)
""")
else:
    print("❌ No MM-GBSA results to save.")

print(f"\n⏱️  Total time: {(time.time()-t_start)/60:.1f} min")
print("🎉 MM-GBSA Analysis Complete!")
