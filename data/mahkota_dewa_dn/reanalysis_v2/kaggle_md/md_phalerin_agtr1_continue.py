# %% [markdown]
# # MD Simulation — Phalerin–AGTR1 (CONTINUATION)
# Resume production MD from Run 1 checkpoint.
#
# **Prerequisites:**
#   1. Run 1 output saved as Kaggle dataset (e.g. "md-phalerin-agtr1-run1")
#   2. Original docking dataset still attached
#
# **What this notebook does:**
#   - Install same dependencies
#   - Rebuild system from saved topology
#   - Load checkpoint → resume production
#   - Combine trajectories (run1 + run2)
#   - Full analysis on combined 100ns trajectory

# %%
# ============================================================
# CELL 0: Install Dependencies (same as Run 1)
# ============================================================
import subprocess, sys, os, glob, time, shutil

t_start = time.time()  # Global timer for time-safety

PY_VER = f"{sys.version_info.major}.{sys.version_info.minor}"
print(f"Kernel Python: {PY_VER} ({sys.executable})")

# ─────────────────────────────────────────────────────────
# STEP 1: Install Miniforge3
# ─────────────────────────────────────────────────────────
print("=" * 60)
print("STEP 1: Installing Miniforge3...")
print("=" * 60)

MINIFORGE_DIR = "/tmp/miniforge"
MAMBA = f"{MINIFORGE_DIR}/bin/mamba"

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
# CELL 1: Configuration (CONTINUATION MODE)
# ============================================================
from pathlib import Path
import json

# --- Paths ---
# Run 1 output (saved as dataset)
# IMPORTANT: Change this to your actual dataset name!
RUN1_DATASET = Path("/kaggle/input/md-phalerin-agtr1-run1")

# Check multiple possible dataset names
possible_names = [
    "md-phalerin-agtr1-run1",
    "md-phalerin-run1",
    "md-run1-output",
]
for name in possible_names:
    candidate = Path(f"/kaggle/input/{name}")
    if candidate.exists():
        RUN1_DATASET = candidate
        break

# The Run 1 output structure
SYSTEM_NAME = "Phalerin_AGTR1"
RUN1_RESULTS = RUN1_DATASET / "md_results" / SYSTEM_NAME

# If dataset saved output directly (no md_results subfolder)
if not RUN1_RESULTS.exists():
    RUN1_RESULTS = RUN1_DATASET / SYSTEM_NAME
if not RUN1_RESULTS.exists():
    # Try flat structure
    RUN1_RESULTS = RUN1_DATASET
    print(f"⚠️ Using flat dataset structure: {RUN1_RESULTS}")

# Original docking dataset
DOCKING_DATASET = Path("/kaggle/input/mahkota-dewa-docking")
LIGAND_SDF = DOCKING_DATASET / "ligands" / "Phalerin.sdf"

# Output for this run
OUTPUT_DIR = Path("/kaggle/working/md_results") / SYSTEM_NAME
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# MD parameters (MUST match Run 1)
RECEPTOR_PDB_ID = "3QXY"
LIGAND_NAME = "Phalerin"

MD_PARAMS = {
    "temperature_K": 310.15,
    "pressure_atm": 1.0,
    "ionic_strength_M": 0.15,
    "timestep_fs": 2.0,
    "nonbonded_cutoff_nm": 1.0,
    "padding_nm": 1.2,
    "minimization_steps": 5000,
    "nvt_steps": 50000,
    "npt_steps": 500000,
    "production_ns": 100,
    "save_interval_ps": 10,
    "log_interval_ps": 10,
    "checkpoint_interval_ns": 10,
}

_dt = MD_PARAMS["timestep_fs"]
MD_PARAMS["production_steps"] = int(MD_PARAMS["production_ns"] * 1e6 / _dt)
MD_PARAMS["save_interval_steps"] = int(MD_PARAMS["save_interval_ps"] * 1000 / _dt)
MD_PARAMS["log_interval_steps"] = int(MD_PARAMS["log_interval_ps"] * 1000 / _dt)
MD_PARAMS["checkpoint_interval_steps"] = int(MD_PARAMS["checkpoint_interval_ns"] * 1e6 / _dt)

# --- Verify Run 1 files ---
print("=" * 60)
print(f"CONTINUATION MODE: {SYSTEM_NAME}")
print(f"Run 1 results: {RUN1_RESULTS}")
print("=" * 60)

required_files = {
    "checkpoint": None,  # Will search for it
    "topology": None,    # Will search for it
    "trajectory": None,  # DCD from run 1
    "docked_sdf": None,  # For GAFF2 parameterization
}

# Find checkpoint (prefer checkpoint.chk, fallback to final.chk)
for name in ["checkpoint.chk", "final.chk", "equilibrated.chk"]:
    p = RUN1_RESULTS / name
    if p.exists():
        required_files["checkpoint"] = p
        break

# Find topology PDB
for name in ["complex_solvated.pdb", "complex_equilibrated.pdb", "complex_minimized.pdb"]:
    p = RUN1_RESULTS / name
    if p.exists():
        required_files["topology"] = p
        break

# Find trajectory
p = RUN1_RESULTS / "production.dcd"
if p.exists():
    required_files["trajectory"] = p

# Find docked SDF
p = RUN1_RESULTS / f"{LIGAND_NAME}_docked.sdf"
if p.exists():
    required_files["docked_sdf"] = p

print("\nRun 1 files found:")
all_ok = True
for key, path in required_files.items():
    if path and path.exists():
        size_mb = path.stat().st_size / 1e6
        print(f"  ✅ {key}: {path.name} ({size_mb:.1f} MB)")
    else:
        print(f"  ❌ {key}: NOT FOUND")
        if key in ["checkpoint", "topology"]:
            all_ok = False

assert all_ok, "Critical Run 1 files missing! Check dataset upload."

# Also check docked SDF
if required_files["docked_sdf"] is None:
    print("  ⚠️ Docked SDF not found in Run 1, will use original SDF")
    required_files["docked_sdf"] = LIGAND_SDF

docked_sdf = required_files["docked_sdf"]
print(f"\n  Using ligand SDF: {docked_sdf}")

# %%
# ============================================================
# CELL 2: Rebuild System & Resume Production MD
# ============================================================
# SAFETY FEATURES:
#   1. Anti-loop assertion: verifies checkpoint step > 0
#   2. Wall-clock time limit: auto-stops before Kaggle 12h
#   3. Explicit step verification prints
# ============================================================
from openmm.app import (
    ForceField, Modeller, PDBFile, PME, HBonds, NoCutoff,
    Simulation, StateDataReporter, DCDReporter, CheckpointReporter
)
from openmm import (
    LangevinMiddleIntegrator, MonteCarloBarostat, Platform
)
from openmmforcefields.generators import GAFFTemplateGenerator
from openff.toolkit import Molecule as OFFMolecule
import openmm.unit as unit
from rdkit import Chem

# ============================================================
# WALL-CLOCK TIME SAFETY
# Kaggle limit = 12h (43200s). We reserve 1.5h for analysis.
# If we approach this limit, we STOP MD and proceed to analysis.
# This ensures output is ALWAYS saved (no more lost runs).
# ============================================================
t_notebook_start = time.time()
MAX_WALL_SECONDS = 11.0 * 3600  # 11 hours max for entire notebook
ANALYSIS_BUFFER_SECONDS = 1.5 * 3600  # 1.5 hours reserved for analysis
MAX_MD_SECONDS = MAX_WALL_SECONDS - ANALYSIS_BUFFER_SECONDS  # ~9.5h for MD

print(f"⏱️  TIME SAFETY: Notebook started at {time.strftime('%Y-%m-%d %H:%M:%S')}")
print(f"  Max wall time: {MAX_WALL_SECONDS/3600:.1f}h")
print(f"  Reserved for analysis: {ANALYSIS_BUFFER_SECONDS/3600:.1f}h")
print(f"  Max time for MD: {MAX_MD_SECONDS/3600:.1f}h")

print("\nRebuilding simulation system from Run 1 topology...")

# Load topology
topology_pdb = required_files["topology"]
pdb = PDBFile(str(topology_pdb))
print(f"  Topology: {pdb.topology.getNumAtoms()} atoms, "
      f"{pdb.topology.getNumResidues()} residues")

# Load ligand for GAFF2 parameterization
off_mol = OFFMolecule.from_file(str(docked_sdf))
print(f"  Ligand: {off_mol.n_atoms} atoms")

# Setup ForceField (must match Run 1 exactly)
ff = ForceField("amber/ff14SB.xml", "amber/tip3p_standard.xml")
gaff_gen = GAFFTemplateGenerator(molecules=[off_mol], forcefield="gaff-2.11")
ff.registerTemplateGenerator(gaff_gen.generator)
print("  ✅ ForceField + GAFF2 ready")

# Create System (must match Run 1)
system = ff.createSystem(
    pdb.topology,
    nonbondedMethod=PME,
    nonbondedCutoff=MD_PARAMS["nonbonded_cutoff_nm"] * unit.nanometers,
    constraints=HBonds,
    rigidWater=True,
    hydrogenMass=1.5 * unit.amu,
)

# Add barostat (was added during NPT in Run 1)
system.addForce(MonteCarloBarostat(
    MD_PARAMS["pressure_atm"] * unit.atmospheres,
    MD_PARAMS["temperature_K"] * unit.kelvin,
    25,
))

# Select Platform
try:
    platform = Platform.getPlatformByName("CUDA")
    properties = {"DeviceIndex": "0", "Precision": "mixed"}
    print("  ✅ Using CUDA platform")
except Exception:
    platform = Platform.getPlatformByName("CPU")
    properties = {}
    print("  ⚠️ Using CPU platform")

# Integrator
integrator = LangevinMiddleIntegrator(
    MD_PARAMS["temperature_K"] * unit.kelvin,
    1.0 / unit.picoseconds,
    MD_PARAMS["timestep_fs"] * unit.femtoseconds,
)

# Create Simulation
simulation = Simulation(pdb.topology, system, integrator, platform, properties)
simulation.context.setPositions(pdb.positions)

# ============================================================
# CHECKPOINT LOADING + ANTI-LOOP VERIFICATION
# ============================================================
checkpoint_file = required_files["checkpoint"]

# Step BEFORE loading checkpoint (should be 0 or from setPositions)
step_before = simulation.currentStep
print(f"\n{'=' * 60}")
print(f"CHECKPOINT VERIFICATION")
print(f"{'=' * 60}")
print(f"  Step BEFORE loading checkpoint: {step_before:,}")

# Load checkpoint
print(f"  Loading checkpoint: {checkpoint_file.name}...")
simulation.loadCheckpoint(str(checkpoint_file))

# Step AFTER loading checkpoint (should be ~22M for 44ns)
step_after = simulation.currentStep
ns_done = step_after * MD_PARAMS["timestep_fs"] / 1e6
remaining_steps = MD_PARAMS["production_steps"] - step_after
remaining_ns = remaining_steps * MD_PARAMS["timestep_fs"] / 1e6

print(f"  Step AFTER loading checkpoint: {step_after:,}")
print(f"  = {ns_done:.2f} ns completed in Run 1")
print(f"  Remaining: {remaining_steps:,} steps ({remaining_ns:.1f} ns)")

# ============================================================
# ANTI-LOOP ASSERTIONS
# These will ABORT the notebook if checkpoint didn't load properly
# ============================================================
MIN_EXPECTED_NS = 10.0  # We expect at least 10 ns from Run 1
min_expected_steps = int(MIN_EXPECTED_NS * 1e6 / MD_PARAMS["timestep_fs"])

assert step_after > 0, (
    f"❌ ABORT: Checkpoint loaded but currentStep is 0! "
    f"This would restart the simulation from the beginning. "
    f"Checkpoint file may be corrupt: {checkpoint_file}"
)

assert step_after > min_expected_steps, (
    f"❌ ABORT: Checkpoint step ({step_after:,} = {ns_done:.1f} ns) "
    f"is less than minimum expected ({min_expected_steps:,} = {MIN_EXPECTED_NS} ns). "
    f"This suggests the checkpoint is from an early/wrong run."
)

assert step_after < MD_PARAMS["production_steps"], (
    f"❌ ABORT: Checkpoint step ({step_after:,}) >= target ({MD_PARAMS['production_steps']:,}). "
    f"Simulation is already complete! No need to continue."
)

assert step_after != step_before, (
    f"❌ ABORT: Step counter didn't change after loading checkpoint! "
    f"Before: {step_before}, After: {step_after}. Checkpoint may not have loaded."
)

print(f"\n  ✅ ALL CHECKS PASSED — checkpoint is valid")
print(f"  ✅ NOT a loop — step counter jumped from {step_before:,} → {step_after:,}")
print(f"  ✅ Will resume from {ns_done:.1f} ns (not from 0)")

# --- Copy Run 1 DCD to output for later concatenation ---
import shutil

run1_dcd = required_files["trajectory"]
if run1_dcd and run1_dcd.exists():
    run1_dcd_copy = OUTPUT_DIR / "production_run1.dcd"
    if not run1_dcd_copy.exists():
        shutil.copy2(run1_dcd, run1_dcd_copy)
        print(f"  Copied Run 1 DCD: {run1_dcd_copy.name} ({run1_dcd.stat().st_size/1e6:.1f} MB)")

# Copy topology files too
for name in ["complex_solvated.pdb", "complex_equilibrated.pdb", "protein_fixed.pdb"]:
    src = RUN1_RESULTS / name
    dst = OUTPUT_DIR / name
    if src.exists() and not dst.exists():
        shutil.copy2(src, dst)

# --- Setup Reporters for continuation ---
run2_dcd = OUTPUT_DIR / "production_run2.dcd"
run2_log = OUTPUT_DIR / "production_run2.log"
run2_checkpoint = OUTPUT_DIR / "checkpoint.chk"

simulation.reporters.clear()

simulation.reporters.append(
    DCDReporter(str(run2_dcd), MD_PARAMS["save_interval_steps"], append=False)
)

simulation.reporters.append(
    StateDataReporter(
        str(run2_log), MD_PARAMS["log_interval_steps"],
        step=True, time=True, potentialEnergy=True,
        kineticEnergy=True, totalEnergy=True,
        temperature=True, volume=True,
        density=True, speed=True,
    )
)

simulation.reporters.append(
    CheckpointReporter(str(run2_checkpoint), MD_PARAMS["checkpoint_interval_steps"])
)

# ============================================================
# RUN PRODUCTION (continuation with TIME SAFETY)
# ============================================================
print(f"\n{'=' * 60}")
print(f"PRODUCTION MD CONTINUATION: {remaining_ns:.1f} ns remaining")
print(f"  From: {ns_done:.1f} ns → To: {MD_PARAMS['production_ns']} ns")
print(f"  Steps: {remaining_steps:,}")
print(f"  ⏱️  Time limit: {MAX_MD_SECONDS/3600:.1f}h for MD portion")
print(f"{'=' * 60}")
print(f"\n🚀 Resuming at {time.strftime('%Y-%m-%d %H:%M:%S')}...")

run_start_time = time.time()
run_start_step = step_after
total_target = MD_PARAMS["production_steps"]

PROGRESS_CHUNK = 50000  # Report every 100 ps
time_stopped_early = False

while simulation.currentStep < total_target:
    # ===== TIME SAFETY CHECK =====
    wall_elapsed = time.time() - t_notebook_start
    if wall_elapsed > MAX_MD_SECONDS:
        ns_now = simulation.currentStep * MD_PARAMS["timestep_fs"] / 1e6
        print(f"\n  ⏱️  TIME LIMIT REACHED after {wall_elapsed/3600:.1f}h wall time")
        print(f"  Stopping production at {ns_now:.1f} ns to ensure analysis runs")
        print(f"  (Another continuation run can pick up from here)")
        time_stopped_early = True
        break
    # =============================

    steps_to_run = min(PROGRESS_CHUNK, total_target - simulation.currentStep)
    simulation.step(steps_to_run)

    current = simulation.currentStep
    elapsed = time.time() - run_start_time
    done_in_run = current - run_start_step
    if done_in_run > 0:
        rate = done_in_run / elapsed
        remaining = total_target - current
        eta_sec = remaining / rate if rate > 0 else 0
        pct = current / total_target * 100
        ns_now = current * MD_PARAMS["timestep_fs"] / 1e6
        wall_total = time.time() - t_notebook_start
        wall_remaining = MAX_MD_SECONDS - wall_total
        print(f"  [{pct:5.1f}%] {ns_now:.1f}/{MD_PARAMS['production_ns']} ns | "
              f"{rate:.0f} steps/s | ETA: {eta_sec/3600:.1f}h | "
              f"Wall: {wall_total/3600:.1f}h (limit in {wall_remaining/3600:.1f}h)", flush=True)

# Save final state (ALWAYS, even if stopped early)
final_checkpoint = OUTPUT_DIR / "final.chk"
simulation.saveCheckpoint(str(final_checkpoint))

state = simulation.context.getState(getPositions=True, getVelocities=True, getEnergy=True)
final_pdb = OUTPUT_DIR / "production_final.pdb"
with open(final_pdb, "w") as f:
    PDBFile.writeFile(simulation.topology, state.getPositions(), f)

total_elapsed = time.time() - run_start_time
ns_completed = simulation.currentStep * MD_PARAMS["timestep_fs"] / 1e6
# Store actual ns for use in analysis cells
ACTUAL_NS = ns_completed

print(f"\n{'=' * 60}")
if time_stopped_early:
    print(f"⏱️  Production MD stopped early (time safety)")
    print(f"  Simulated in this run: {(simulation.currentStep - run_start_step) * MD_PARAMS['timestep_fs'] / 1e6:.1f} ns")
    print(f"  Total so far: {ns_completed:.1f} / {MD_PARAMS['production_ns']} ns")
    print(f"  📌 To reach 100 ns, run another continuation with this output as dataset")
else:
    print(f"✅ Production MD complete!")
    print(f"  Total simulated: {ns_completed:.1f} ns")
print(f"  Wall time (this run): {total_elapsed/3600:.2f} hours")
print(f"  Run 2 DCD: {run2_dcd}")
print(f"  Final checkpoint: {final_checkpoint}")
print(f"  Finished at {time.strftime('%Y-%m-%d %H:%M:%S')}")
print(f"{'=' * 60}")

# %%
# ============================================================
# CELL 3: Combine Trajectories
# ============================================================
import mdtraj as md
import numpy as np

print("Combining Run 1 + Run 2 trajectories...")

# Use the equilibrated PDB as topology reference
equil_pdb = OUTPUT_DIR / "complex_equilibrated.pdb"
if not equil_pdb.exists():
    equil_pdb = OUTPUT_DIR / "complex_solvated.pdb"

# Load both DCDs
run1_dcd_path = OUTPUT_DIR / "production_run1.dcd"
run2_dcd_path = OUTPUT_DIR / "production_run2.dcd"

trajs = []
if run1_dcd_path.exists():
    t1 = md.load(str(run1_dcd_path), top=str(equil_pdb))
    print(f"  Run 1: {t1.n_frames} frames ({t1.time[0]:.1f} - {t1.time[-1]:.1f} ps)")
    trajs.append(t1)
if run2_dcd_path.exists():
    t2 = md.load(str(run2_dcd_path), top=str(equil_pdb))
    print(f"  Run 2: {t2.n_frames} frames ({t2.time[0]:.1f} - {t2.time[-1]:.1f} ps)")
    trajs.append(t2)

if len(trajs) == 2:
    # Fix time axis for Run 2 (mdtraj may reset time to 0)
    # Run 2 starts where Run 1 ended
    if t2.time[0] < t1.time[-1]:
        dt_ps = MD_PARAMS["save_interval_ps"]
        t2_start = t1.time[-1] + dt_ps
        t2.time = np.arange(t2.n_frames) * dt_ps + t2_start
        print(f"  Fixed Run 2 time: {t2.time[0]:.1f} - {t2.time[-1]:.1f} ps")
    traj = md.join([t1, t2])
elif len(trajs) == 1:
    traj = trajs[0]
else:
    raise RuntimeError("No trajectory files found!")

print(f"\n  Combined trajectory:")
print(f"    Frames: {traj.n_frames}")
print(f"    Time: {traj.time[0]:.1f} - {traj.time[-1]:.1f} ps")
print(f"    Duration: {(traj.time[-1] - traj.time[0]) / 1000:.1f} ns")
print(f"    Atoms: {traj.n_atoms}")

# %%
# ============================================================
# CELL 4: Full Trajectory Analysis (10 metrics)
# ============================================================
import warnings
warnings.filterwarnings("ignore")

print("Analyzing combined trajectory...")

# --- [0/10] PBC Reimaging ---
print("\n[0/10] Reimaging trajectory (fixing PBC artifacts)...")
try:
    protein_atoms_set = set(traj.topology.select("protein").tolist())
    anchor = [protein_atoms_set]
    traj = traj.image_molecules(inplace=False, anchor_molecules=anchor)
    print("  ✅ Reimaging complete")
except Exception as e:
    print(f"  ⚠️ Reimaging failed ({e}), continuing without it")

# --- Atom Selections ---
protein_atoms = traj.topology.select("protein")
ca_atoms = traj.topology.select("protein and name CA")
backbone_atoms = traj.topology.select("protein and (name CA or name C or name N or name O)")

water_atoms = set(traj.topology.select("water"))
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

assert len(ligand_atoms) > 0, "No ligand atoms detected!"
assert len(ca_atoms) > 0, "No CA atoms detected!"

# --- Time axis ---
time_ns = traj.time / 1000.0

# --- [1/10] RMSD (protein backbone) ---
print("\n[1/10] Calculating backbone RMSD...")
traj_bb = traj.atom_slice(backbone_atoms)
rmsd_protein = md.rmsd(traj_bb, traj_bb, frame=0) * 10
print(f"  Mean: {np.mean(rmsd_protein):.2f} ± {np.std(rmsd_protein):.2f} Å")

# --- [2/10] RMSD (ligand) ---
print("[2/10] Calculating ligand RMSD...")
traj_aligned = traj.superpose(traj, frame=0, atom_indices=ca_atoms)
lig_traj = traj_aligned.atom_slice(ligand_atoms)
rmsd_ligand = md.rmsd(lig_traj, lig_traj, frame=0) * 10
print(f"  Mean: {np.mean(rmsd_ligand):.2f} ± {np.std(rmsd_ligand):.2f} Å")

# --- [3/10] RMSF ---
print("[3/10] Calculating RMSF (Cα)...")
ca_traj = traj_aligned.atom_slice(ca_atoms)
rmsf = md.rmsf(ca_traj, ca_traj, frame=0) * 10
residue_ids = [traj.topology.atom(a).residue.resSeq for a in ca_atoms]
print(f"  Max RMSF: {np.max(rmsf):.2f} Å (residue {residue_ids[np.argmax(rmsf)]})")

# --- [4/10] Rg ---
print("[4/10] Calculating Radius of Gyration...")
rg = md.compute_rg(traj.atom_slice(protein_atoms)) * 10
print(f"  Mean: {np.mean(rg):.2f} ± {np.std(rg):.2f} Å")

# --- [5/10] H-bonds ---
print("\n[5/10] Calculating protein-ligand H-bonds (sampled)...")
lig_set = set(ligand_atoms.tolist())
prot_set = set(protein_atoms.tolist())
sample_step = max(1, traj.n_frames // 1000)
sampled_indices = list(range(0, traj.n_frames, sample_step))

hbond_counts_sampled = []
for idx in sampled_indices:
    frame = traj[idx]
    try:
        hbs = md.baker_hubbard(frame, freq=0.0)
        count = sum(1 for hb in hbs
                    if (hb[0] in lig_set and hb[2] in prot_set) or
                       (hb[0] in prot_set and hb[2] in lig_set))
        hbond_counts_sampled.append(count)
    except Exception:
        hbond_counts_sampled.append(0)

hbond_counts_sampled = np.array(hbond_counts_sampled)
hbond_counts = np.interp(np.arange(traj.n_frames), sampled_indices, hbond_counts_sampled)
print(f"  Sampled {len(sampled_indices)} frames")
print(f"  Mean H-bonds: {np.mean(hbond_counts):.1f} ± {np.std(hbond_counts):.1f}")

# --- [6/10] Contact Frequency ---
print("\n[6/10] Calculating protein-ligand contact frequency...")
CONTACT_CUTOFF_NM = 0.45
residue_contact_counts = {}

for idx in sampled_indices:
    frame = traj[idx]
    for ca_idx in ca_atoms:
        res = traj.topology.atom(ca_idx).residue
        res_atoms = [a.index for a in res.atoms]
        pairs = [(r, l) for r in res_atoms for l in ligand_atoms]
        if len(pairs) == 0:
            continue
        distances = md.compute_distances(frame, pairs)[0]
        if np.min(distances) < CONTACT_CUTOFF_NM:
            key = (res.resSeq, res.name)
            residue_contact_counts[key] = residue_contact_counts.get(key, 0) + 1

n_sampled = len(sampled_indices)
contact_freq = {k: v / n_sampled for k, v in residue_contact_counts.items()}
contact_freq_sorted = sorted(contact_freq.items(), key=lambda x: -x[1])

print(f"  Residues contacting ligand (>10% of time):")
top_contacts = [(k, v) for k, v in contact_freq_sorted if v > 0.1]
for (resSeq, resName), freq in top_contacts[:15]:
    print(f"    {resName}{resSeq}: {freq*100:.1f}%")

# --- [7/10] Key Distances ---
print("\n[7/10] Monitoring key protein-ligand distances...")
lig_com_distances = {}
top3_residues = [(resSeq, resName) for (resSeq, resName), _ in contact_freq_sorted[:3]]

for resSeq, resName in top3_residues:
    res_atom_indices = [a.index for a in traj.topology.atoms
                        if a.residue.resSeq == resSeq and a.residue.name == resName]
    if len(res_atom_indices) == 0:
        continue
    dists = []
    for idx in range(0, traj.n_frames, sample_step):
        frame = traj[idx]
        pairs = [(r, l) for r in res_atom_indices for l in ligand_atoms]
        d = md.compute_distances(frame, pairs)[0]
        dists.append(np.min(d) * 10)

    dists_full = np.interp(np.arange(traj.n_frames),
                           list(range(0, traj.n_frames, sample_step)), dists)
    lig_com_distances[f"{resName}{resSeq}"] = dists_full
    print(f"  {resName}{resSeq}: mean={np.mean(dists_full):.2f} ± {np.std(dists_full):.2f} Å")

# --- [8/10] SASA ---
print("\n[8/10] Calculating SASA...")
sasa_sample_step = max(1, traj.n_frames // 500)
sasa_indices = list(range(0, traj.n_frames, sasa_sample_step))

prot_traj = traj.atom_slice(protein_atoms)
sasa_sampled = []
for idx in sasa_indices:
    frame_sasa = md.shrake_rupley(prot_traj[idx], mode='residue')
    sasa_sampled.append(np.sum(frame_sasa))

sasa_sampled = np.array(sasa_sampled) * 100
sasa_full = np.interp(np.arange(traj.n_frames), sasa_indices, sasa_sampled)
print(f"  Mean total SASA: {np.mean(sasa_full):.1f} ± {np.std(sasa_full):.1f} Å²")

# --- [9/10] Equilibration Check ---
print("\n[9/10] Checking equilibration...")
n_blocks = 5
block_size = len(rmsd_protein) // n_blocks
block_means = [np.mean(rmsd_protein[i*block_size:(i+1)*block_size]) for i in range(n_blocks)]

print(f"  RMSD block means:")
for i, bm in enumerate(block_means):
    ns_block = MD_PARAMS['production_ns'] / n_blocks
    print(f"    Block {i+1} ({i*ns_block:.0f}-{(i+1)*ns_block:.0f} ns): {bm:.2f} Å")

last3_mean = np.mean(block_means[2:])
first_dev = abs(block_means[0] - last3_mean) / last3_mean * 100
if first_dev > 50:
    equil_ns = MD_PARAMS['production_ns'] / n_blocks
    print(f"  ⚠️ First block deviates {first_dev:.0f}% — skip first {equil_ns:.0f} ns")
else:
    equil_ns = 0
    print(f"  ✅ Well-equilibrated (deviation: {first_dev:.0f}%)")

from scipy import stats as scipy_stats
last3_time = np.arange(len(block_means[2:]))
slope, _, _, p_val, _ = scipy_stats.linregress(last3_time, block_means[2:])
if p_val < 0.05:
    print(f"  ⚠️ Possible drift (slope={slope:.4f}, p={p_val:.3f})")
else:
    print(f"  ✅ No significant drift (p={p_val:.3f})")

# --- [10/10] Representative Frames ---
print("\n[10/10] Extracting representative frames...")
rep_frames = {
    "initial": 0,
    "middle": traj.n_frames // 2,
    "final": traj.n_frames - 1,
    "equilibrated_mean": int(np.argmin(np.abs(rmsd_protein - np.mean(rmsd_protein[traj.n_frames//5:])))),
}

for label, fidx in rep_frames.items():
    solute_atoms_arr = np.concatenate([protein_atoms, ligand_atoms])
    frame = traj[fidx].atom_slice(solute_atoms_arr)
    out_path = OUTPUT_DIR / f"frame_{label}.pdb"
    frame.save(str(out_path))
    print(f"  {label}: frame {fidx} (t={time_ns[fidx]:.1f} ns, RMSD={rmsd_protein[fidx]:.2f} Å)")

# --- Save intermediate results ---
print("\n💾 Saving analysis results...")
import pandas as pd

ts_df = pd.DataFrame({
    "Time_ns": time_ns, "RMSD_Protein_A": rmsd_protein,
    "RMSD_Ligand_A": rmsd_ligand, "Rg_A": rg,
    "Hbonds": hbond_counts, "SASA_A2": sasa_full,
})
ts_df.to_csv(OUTPUT_DIR / "timeseries.csv", index=False)

rmsf_df = pd.DataFrame({"Residue": residue_ids, "RMSF_A": rmsf})
rmsf_df.to_csv(OUTPUT_DIR / "rmsf_per_residue.csv", index=False)

contact_df = pd.DataFrame([
    {"Residue": f"{rn}{rs}", "ResSeq": rs, "ResName": rn, "Contact_Frequency": freq}
    for (rs, rn), freq in contact_freq_sorted if freq > 0.05
])
contact_df.to_csv(OUTPUT_DIR / "contact_frequency.csv", index=False)

if lig_com_distances:
    dist_df = pd.DataFrame({"Time_ns": time_ns, **lig_com_distances})
    dist_df.to_csv(OUTPUT_DIR / "key_distances.csv", index=False)

print(f"  ✅ All CSVs saved")

# %%
# ============================================================
# CELL 5: MM-GBSA Binding Free Energy
# ============================================================
from openmm.app import ForceField, PDBFile, NoCutoff, HBonds, Simulation
from openmm import LangevinMiddleIntegrator
import openmm.unit as unit

print("=" * 60)
print("MM-GBSA BINDING FREE ENERGY ANALYSIS")
print("=" * 60)

N_FRAMES_MMGBSA = 100
start_frame = int(traj.n_frames * 0.2)
frame_indices = np.linspace(start_frame, traj.n_frames - 1, N_FRAMES_MMGBSA, dtype=int)

print(f"  Analyzing {N_FRAMES_MMGBSA} frames (from {start_frame} to {traj.n_frames-1})")

solute_atoms = np.concatenate([protein_atoms, ligand_atoms])
prot_local_idx = np.arange(len(protein_atoms))
lig_local_idx = np.arange(len(protein_atoms), len(protein_atoms) + len(ligand_atoms))

print("  Setting up implicit solvent force field...")
ff_gb = ForceField("amber/ff14SB.xml", "implicit/obc2.xml")
gaff_gen_gb = GAFFTemplateGenerator(molecules=[off_mol], forcefield="gaff-2.11")
ff_gb.registerTemplateGenerator(gaff_gen_gb.generator)

ref_frame = traj[frame_indices[0]]
ref_solute = ref_frame.atom_slice(solute_atoms)
ref_receptor = ref_frame.atom_slice(protein_atoms)
ref_ligand = ref_frame.atom_slice(ligand_atoms)

tmp_dir = OUTPUT_DIR / "_tmp_mmgbsa"
tmp_dir.mkdir(exist_ok=True)

ref_solute.save(str(tmp_dir / "complex.pdb"))
ref_receptor.save(str(tmp_dir / "receptor.pdb"))
ref_ligand.save(str(tmp_dir / "ligand.pdb"))

def build_gb_system(pdb_path):
    pdb_obj = PDBFile(str(pdb_path))
    sys_gb = ff_gb.createSystem(pdb_obj.topology, nonbondedMethod=NoCutoff, constraints=HBonds)
    integ = LangevinMiddleIntegrator(300*unit.kelvin, 1/unit.picoseconds, 2*unit.femtoseconds)
    return Simulation(pdb_obj.topology, sys_gb, integ)

try:
    sim_complex = build_gb_system(tmp_dir / "complex.pdb")
    sim_receptor = build_gb_system(tmp_dir / "receptor.pdb")
    sim_ligand = build_gb_system(tmp_dir / "ligand.pdb")
    print("  ✅ All three systems built")
    mmgbsa_ready = True
except Exception as e:
    print(f"  ❌ System build failed: {e}")
    mmgbsa_ready = False

dG_values = []
errors = 0

if mmgbsa_ready:
    print(f"\n  Computing {N_FRAMES_MMGBSA} frames...")
    t0 = time.time()
    for i, fidx in enumerate(frame_indices):
        try:
            frame = traj[fidx]
            frame_solute = frame.atom_slice(solute_atoms)
            pos_complex = frame_solute.xyz[0]
            pos_receptor = pos_complex[prot_local_idx]
            pos_ligand = pos_complex[lig_local_idx]

            sim_complex.context.setPositions(pos_complex)
            E_complex = sim_complex.context.getState(getEnergy=True).getPotentialEnergy()
            E_complex = E_complex.value_in_unit(unit.kilocalories_per_mole)

            sim_receptor.context.setPositions(pos_receptor)
            E_receptor = sim_receptor.context.getState(getEnergy=True).getPotentialEnergy()
            E_receptor = E_receptor.value_in_unit(unit.kilocalories_per_mole)

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
    print(f"{'=' * 60}")
else:
    print("\n⚠️ MM-GBSA failed for all frames.")
    dG_mean, dG_std, dG_sem = 0.0, 0.0, 0.0

# %%
# ============================================================
# CELL 6: Publication Figures
# ============================================================
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 10,
    "axes.labelsize": 12, "axes.titlesize": 13,
    "xtick.labelsize": 9, "ytick.labelsize": 9,
    "legend.fontsize": 9, "figure.dpi": 300, "axes.linewidth": 1.2,
})

window = max(1, len(rmsd_protein) // 100)

# --- Figure 1: 4-Panel ---
fig = plt.figure(figsize=(14, 10))
gs = GridSpec(2, 2, hspace=0.35, wspace=0.3)

ax1 = fig.add_subplot(gs[0, 0])
ax1.plot(time_ns, rmsd_protein, color="#2196F3", linewidth=0.3, alpha=0.4)
rmsd_smooth = np.convolve(rmsd_protein, np.ones(window)/window, mode="valid")
time_smooth = time_ns[:len(rmsd_smooth)]
ax1.plot(time_smooth, rmsd_smooth, color="#1565C0", linewidth=1.5, label="Running avg")
ax1.set_xlabel("Time (ns)"); ax1.set_ylabel("RMSD (Å)")
ax1.set_title(f"(A) Protein Backbone RMSD\nmean={np.mean(rmsd_protein):.2f}±{np.std(rmsd_protein):.2f} Å")
ax1.legend(); ax1.set_xlim(0, time_ns[-1])

ax2 = fig.add_subplot(gs[0, 1])
ax2.plot(time_ns, rmsd_ligand, color="#FF9800", linewidth=0.3, alpha=0.4)
rmsd_lig_smooth = np.convolve(rmsd_ligand, np.ones(window)/window, mode="valid")
ax2.plot(time_smooth, rmsd_lig_smooth, color="#E65100", linewidth=1.5, label="Running avg")
ax2.set_xlabel("Time (ns)"); ax2.set_ylabel("RMSD (Å)")
ax2.set_title(f"(B) Ligand RMSD\nmean={np.mean(rmsd_ligand):.2f}±{np.std(rmsd_ligand):.2f} Å")
ax2.legend(); ax2.set_xlim(0, time_ns[-1])

ax3 = fig.add_subplot(gs[1, 0])
ax3.plot(residue_ids, rmsf, color="#4CAF50", linewidth=1.0)
ax3.fill_between(residue_ids, 0, rmsf, color="#4CAF50", alpha=0.15)
ax3.set_xlabel("Residue Number"); ax3.set_ylabel("RMSF (Å)")
ax3.set_title("(C) Cα RMSF per Residue")
threshold = np.mean(rmsf) + 2 * np.std(rmsf)
ax3.axhline(threshold, color="red", linestyle="--", alpha=0.5, label=f"Mean+2σ ({threshold:.1f} Å)")
ax3.legend()

ax4 = fig.add_subplot(gs[1, 1])
ax4.plot(time_ns, rg, color="#9C27B0", linewidth=0.3, alpha=0.4)
rg_smooth = np.convolve(rg, np.ones(window)/window, mode="valid")
ax4.plot(time_smooth, rg_smooth, color="#6A1B9A", linewidth=1.5, label="Running avg")
ax4.set_xlabel("Time (ns)"); ax4.set_ylabel("Rg (Å)")
ax4.set_title(f"(D) Radius of Gyration\nmean={np.mean(rg):.2f}±{np.std(rg):.2f} Å")
ax4.legend(); ax4.set_xlim(0, time_ns[-1])

actual_ns = time_ns[-1]
fig.suptitle(f"MD Simulation: {LIGAND_NAME} – {RECEPTOR_PDB_ID} (AGTR1) | {actual_ns:.0f} ns",
             fontsize=14, fontweight="bold", y=0.98)

fig_path_1 = OUTPUT_DIR / "md_analysis_4panel.png"
plt.savefig(fig_path_1, dpi=300, bbox_inches="tight", facecolor="white")
plt.close()
print(f"✅ Figure 1 saved: {fig_path_1}")

# --- Figure 2: H-bond + MM-GBSA ---
fig2, (ax5, ax6) = plt.subplots(1, 2, figsize=(14, 5))

ax5.plot(time_ns, hbond_counts, color="#F44336", linewidth=0.3, alpha=0.4)
hb_smooth = np.convolve(hbond_counts, np.ones(window)/window, mode="valid")
ax5.plot(time_smooth, hb_smooth, color="#B71C1C", linewidth=1.5, label="Running avg")
ax5.set_xlabel("Time (ns)"); ax5.set_ylabel("Number of H-bonds")
ax5.set_title(f"(E) Protein–Ligand H-bonds\nmean={np.mean(hbond_counts):.1f}±{np.std(hbond_counts):.1f}")
ax5.legend(); ax5.set_xlim(0, time_ns[-1])

if len(dG_values) > 0:
    ax6.hist(dG_values, bins=25, color="#00BCD4", edgecolor="black", alpha=0.7)
    ax6.axvline(dG_mean, color="red", linewidth=2, linestyle="--",
                label=f"Mean: {dG_mean:.2f} ± {dG_sem:.2f}")
    ax6.set_xlabel("ΔG_bind (kcal/mol)"); ax6.set_ylabel("Count")
    ax6.set_title(f"(F) MM-GBSA ΔG = {dG_mean:.2f} ± {dG_sem:.2f} kcal/mol")
    ax6.legend()
else:
    ax6.text(0.5, 0.5, "MM-GBSA\nnot available", ha="center", va="center",
             transform=ax6.transAxes, fontsize=14)

fig2.suptitle(f"H-bond & Binding Energy: {LIGAND_NAME} – {RECEPTOR_PDB_ID}",
              fontsize=14, fontweight="bold")
plt.tight_layout()
plt.savefig(OUTPUT_DIR / "md_hbond_mmgbsa.png", dpi=300, bbox_inches="tight", facecolor="white")
plt.close()
print(f"✅ Figure 2 saved")

# --- Figure 3: SASA + Contacts + Distances ---
fig3 = plt.figure(figsize=(16, 5))
gs3 = GridSpec(1, 3, wspace=0.35)

ax7 = fig3.add_subplot(gs3[0, 0])
ax7.plot(time_ns, sasa_full, color="#607D8B", linewidth=0.3, alpha=0.4)
sasa_smooth = np.convolve(sasa_full, np.ones(window)/window, mode="valid")
ax7.plot(time_smooth, sasa_smooth, color="#37474F", linewidth=1.5, label="Running avg")
ax7.set_xlabel("Time (ns)"); ax7.set_ylabel("SASA (Å²)")
ax7.set_title(f"(G) Protein SASA\nmean={np.mean(sasa_full):.0f}±{np.std(sasa_full):.0f} Å²")
ax7.legend(); ax7.set_xlim(0, time_ns[-1])

ax8 = fig3.add_subplot(gs3[0, 1])
top_n = min(12, len(top_contacts))
if top_n > 0:
    labels_c = [f"{rn}{rs}" for (rs, rn), _ in top_contacts[:top_n]]
    freqs_c = [freq * 100 for _, freq in top_contacts[:top_n]]
    colors_c = plt.cm.YlOrRd(np.linspace(0.3, 0.9, top_n))
    ax8.barh(range(top_n), freqs_c, color=colors_c, edgecolor="black", linewidth=0.5)
    ax8.set_yticks(range(top_n)); ax8.set_yticklabels(labels_c)
    ax8.set_xlabel("Contact Frequency (%)"); ax8.set_title("(H) Protein-Ligand Contacts")
    ax8.invert_yaxis()

ax9 = fig3.add_subplot(gs3[0, 2])
dist_colors = ["#E91E63", "#3F51B5", "#009688"]
for i, (label, dists) in enumerate(lig_com_distances.items()):
    c = dist_colors[i % len(dist_colors)]
    ax9.plot(time_ns, dists, color=c, linewidth=0.3, alpha=0.3)
    d_smooth = np.convolve(dists, np.ones(window)/window, mode="valid")
    ax9.plot(time_smooth, d_smooth, color=c, linewidth=1.5, label=label)
ax9.axhline(4.5, color="gray", linestyle=":", alpha=0.5, label="Contact cutoff")
ax9.set_xlabel("Time (ns)"); ax9.set_ylabel("Min Distance (Å)")
ax9.set_title("(I) Key Residue–Ligand Distances")
ax9.legend(fontsize=7); ax9.set_xlim(0, time_ns[-1])

fig3.suptitle(f"Extended Analysis: {LIGAND_NAME} – {RECEPTOR_PDB_ID}", fontsize=14, fontweight="bold")
plt.tight_layout()
plt.savefig(OUTPUT_DIR / "md_extended_analysis.png", dpi=300, bbox_inches="tight", facecolor="white")
plt.close()
print(f"✅ Figure 3 saved")

# %%
# ============================================================
# CELL 7: Summary + Package
# ============================================================
import json as json_mod

# Get best affinity from Run 1 log or set manually
best_affinity = -9.96  # from Run 1

summary = {
    "System": SYSTEM_NAME,
    "Receptor": RECEPTOR_PDB_ID,
    "Ligand": LIGAND_NAME,
    "Simulation_Time_ns": float(actual_ns),
    "N_Frames": int(traj.n_frames),
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
    "MMGBSA_N_Frames": int(len(dG_values)),
}

pd.DataFrame([summary]).to_csv(OUTPUT_DIR / "md_summary.csv", index=False)
with open(OUTPUT_DIR / "md_summary.json", "w") as f:
    json_mod.dump(summary, f, indent=2, default=str)

print(f"✅ Summary saved")

# --- Package for download ---
import zipfile

zip_path = Path("/kaggle/working") / f"MD_{SYSTEM_NAME}_results.zip"
with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
    for f in sorted(OUTPUT_DIR.glob("*")):
        if f.suffix in [".dcd", ".chk"]:
            print(f"  Skip: {f.name} ({f.stat().st_size/1e6:.1f} MB)")
            continue
        if f.name.startswith("_tmp") or f.is_dir():
            continue
        zf.write(f, f"{SYSTEM_NAME}/{f.name}")
        print(f"  Added: {f.name}")

print(f"\n✅ Download: {zip_path} ({zip_path.stat().st_size/1e6:.1f} MB)")

# --- Final Report ---
print(f"\n{'=' * 60}")
print(f"📋 FINAL REPORT: {SYSTEM_NAME}")
print(f"{'=' * 60}")
print(f"""
System: {LIGAND_NAME} bound to {RECEPTOR_PDB_ID} (AGTR1)
Simulation: {actual_ns:.0f} ns @ {MD_PARAMS['temperature_K']}K

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
""")
print("🎉 MD Simulation Complete!")
