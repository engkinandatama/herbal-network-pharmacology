# %% [markdown]
# # MD Simulation -- Mangiferin-RELA (RUN 3: MD-ONLY)
# Resume from Run 2 checkpoint (75.6 ns) -> finish at 100 ns.
#
# **MD-ONLY MODE**: No analysis in this notebook.
# Analysis will be done separately after 100 ns is reached.
#
# **Prerequisites:**
#   1. Run 2 output saved as dataset (e.g. "md-mangiferin-rela-run2")
#   2. Original docking dataset still attached

# %%
# ============================================================
# CELL 0: Install Dependencies
# ============================================================
import subprocess, sys, os, glob, time, shutil

t_start = time.time()

PY_VER = f"{sys.version_info.major}.{sys.version_info.minor}"
print(f"Kernel Python: {PY_VER} ({sys.executable})")

# ?????????????????????????????????????????????????????????
# STEP 1: Install Miniforge3
# ?????????????????????????????????????????????????????????
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
    print(f"  [OK] Miniforge3 installed ({time.time()-t_start:.0f}s)")
else:
    print(f"  [OK] Already present")

# ?????????????????????????????????????????????????????????
# STEP 2: Mamba install (ONLY what's needed for MD)
# ?????????????????????????????????????????????????????????
print("\n" + "=" * 60)
print("STEP 2: Mamba install...")
print("=" * 60)

packages = f"python={PY_VER} numpy openmm cudatoolkit openmmforcefields openff-toolkit ambertools pdbfixer parmed"
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

# Remove conda rdkit (not needed for MD-only, but pulled as dep)
print("  Cleaning conda rdkit+scipy...")
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
print(f"  [OK] Cleaned ({removed} files)")

# Pip install rdkit only (needed for openff ligand loading)
print("  Installing rdkit via pip...")
subprocess.run([sys.executable, "-m", "pip", "install", "-q", "rdkit"],
               capture_output=True, text=True, timeout=300)

# ?????????????????????????????????????????????????????????
# STEP 3: Patch paths
# ?????????????????????????????????????????????????????????
print("\n" + "=" * 60)
print("STEP 3: Remove system numpy & patch paths...")
print("=" * 60)

subprocess.run([sys.executable, "-m", "pip", "uninstall", "-y", "numpy"],
               capture_output=True, text=True)
for np_dir in glob.glob("/usr/local/lib/python3.12/dist-packages/numpy*"):
    shutil.rmtree(np_dir, ignore_errors=True)
print("  [OK] System numpy removed")

mods_to_flush = [k for k in sys.modules if 'numpy' in k.lower()]
for mod in mods_to_flush:
    del sys.modules[mod]
print(f"  [OK] Flushed {len(mods_to_flush)} cached modules")

site_dirs = glob.glob(f"{MINIFORGE_DIR}/lib/python{PY_VER}/site-packages")
site_dirs += glob.glob(f"{MINIFORGE_DIR}/lib/python*/site-packages")
for sp in sorted(set(site_dirs)):
    if os.path.isdir(sp) and sp not in sys.path:
        sys.path.insert(0, sp)
        print(f"  [OK] sys.path += {sp}")

os.environ["PATH"] = f"{MINIFORGE_DIR}/bin:" + os.environ.get("PATH", "")
os.environ["LD_LIBRARY_PATH"] = f"{MINIFORGE_DIR}/lib:" + os.environ.get("LD_LIBRARY_PATH", "")
os.environ["AMBERHOME"] = MINIFORGE_DIR
print("  [OK] Environment set")

# ?????????????????????????????????????????????????????????
# STEP 4: Quick verify
# ?????????????????????????????????????????????????????????
print("\n" + "=" * 60)
print("STEP 4: Verifying...")
print("=" * 60)

import numpy as np; print(f"  [OK] numpy {np.__version__}")
import openmm; print(f"  [OK] OpenMM {openmm.__version__}")
from openmm import Platform
cuda_ok = False
for i in range(Platform.getNumPlatforms()):
    if Platform.getPlatform(i).getName() == "CUDA":
        cuda_ok = True
        break
print(f"  [OK] CUDA: {'ready' if cuda_ok else '[WARN] NOT FOUND'}")

from openff.toolkit import Molecule as OFFMolecule; print("  [OK] openff.toolkit")
import openmmforcefields; print("  [OK] openmmforcefields")

elapsed = time.time() - t_start
print(f"\n  Setup: {elapsed:.0f}s")
print("=" * 60)

# %%
# ============================================================
# CELL 1: Configuration + Find Files
# ============================================================
from pathlib import Path

SYSTEM_NAME = "Mangiferin_RELA"

# --- Locate Run 2 output dataset ---
RUN2_DATASET = Path("/kaggle/input/md-mangiferin-rela-run2")

for name in ["md-mangiferin-rela-run2", "md-phalerin-run2",
             "md-mangiferin-rela-continue"]:
    candidate = Path(f"/kaggle/input/{name}")
    if candidate.exists():
        RUN2_DATASET = candidate
        break

if not RUN2_DATASET.exists():
    print("[WARN] Run 2 dataset not found at expected paths.")
    # Fallback to checking any md-results folder
    md_results = list(Path("/kaggle/input").glob("**/md_results"))
    if md_results:
        RUN2_DATASET = md_results[0].parent
        print(f"  Found potential dataset: {RUN2_DATASET}")

# Find the results folder
RUN2_RESULTS = RUN2_DATASET / "md_results" / SYSTEM_NAME
if not RUN2_RESULTS.exists():
    RUN2_RESULTS = RUN2_DATASET / SYSTEM_NAME
if not RUN2_RESULTS.exists():
    RUN2_RESULTS = RUN2_DATASET
    print(f"[WARN] Using flat structure: {RUN2_RESULTS}")

# Original docking dataset (or dedicated input)
DOCKING_DATASET = Path("/kaggle/input/mahkota-dewa-docking")
if Path("/kaggle/input/md-mangiferin-rela-input").exists():
    DOCKING_DATASET = Path("/kaggle/input/md-mangiferin-rela-input")

# Update Run 1 dataset check
RUN1_DATASET = None
for name in ["md-mangiferin-rela-run1-output", "md-mangiferin-rela-run1", "md-run1-output"]:
    candidate = Path(f"/kaggle/input/{name}")
    if candidate.exists():
        RUN1_DATASET = candidate
        break

# Output for this run
OUTPUT_DIR = Path("/kaggle/working/md_results") / SYSTEM_NAME
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# --- MD Parameters (MUST match previous runs) ---
MD_PARAMS = {
    "production_ns": 100,
    "timestep_fs": 2.0,
    "temperature_K": 300,
    "pressure_atm": 1.0,
    "nonbonded_cutoff_nm": 1.0,
    "save_interval_ps": 10,
    "log_interval_ps": 100,
    "checkpoint_interval_ns": 10,
}
MD_PARAMS["production_steps"] = int(MD_PARAMS["production_ns"] * 1e6 / MD_PARAMS["timestep_fs"])
MD_PARAMS["save_interval_steps"] = int(MD_PARAMS["save_interval_ps"] * 1000 / MD_PARAMS["timestep_fs"])
MD_PARAMS["log_interval_steps"] = int(MD_PARAMS["log_interval_ps"] * 1000 / MD_PARAMS["timestep_fs"])
MD_PARAMS["checkpoint_interval_steps"] = int(MD_PARAMS["checkpoint_interval_ns"] * 1e6 / MD_PARAMS["timestep_fs"])

# --- Find required files ---
print("=" * 60)
print(f"RUN 3 (MD-ONLY): {SYSTEM_NAME}")
print(f"Run 2 results: {RUN2_RESULTS}")
print("=" * 60)

# Prefer final.chk (75.6 ns) over checkpoint.chk (70 ns)
final_chk = RUN2_RESULTS / "final.chk"
regular_chk = RUN2_RESULTS / "checkpoint.chk"

if final_chk.exists():
    checkpoint_file = final_chk
    print(f"\n  [OK] Using final.chk ({final_chk.stat().st_size/1e6:.1f} MB) -- 75.6 ns")
elif regular_chk.exists():
    checkpoint_file = regular_chk
    print(f"\n  [WARN] final.chk not found, using checkpoint.chk ({regular_chk.stat().st_size/1e6:.1f} MB) -- 70 ns")
else:
    raise FileNotFoundError("No checkpoint found in Run 2 output!")

# Topology
topology_file = RUN2_RESULTS / "complex_solvated.pdb"
if not topology_file.exists():
    raise FileNotFoundError(f"Topology not found: {topology_file}")
print(f"  [OK] Topology: {topology_file.name} ({topology_file.stat().st_size/1e6:.1f} MB)")

# Ligand SDF (check multiple locations)
docked_sdf = None
for sdf_candidate in [
    RUN2_RESULTS / "Mangiferin_docked.sdf",
    RUN2_RESULTS / f"{SYSTEM_NAME.split('_')[0]}_docked.sdf",
    DOCKING_DATASET / "ligands" / "Mangiferin.sdf" if DOCKING_DATASET.exists() else None,
]:
    if sdf_candidate and sdf_candidate.exists():
        docked_sdf = sdf_candidate
        break

# If not found in Run 2 output, check Run 1
if docked_sdf is None and RUN1_DATASET:
    run1_results = RUN1_DATASET / "md_results" / SYSTEM_NAME
    for sdf_candidate in [
        run1_results / "Mangiferin_docked.sdf",
        run1_results / f"{SYSTEM_NAME.split('_')[0]}_docked.sdf",
    ]:
        if sdf_candidate.exists():
            docked_sdf = sdf_candidate
            break

if docked_sdf is None:
    raise FileNotFoundError("Ligand SDF not found anywhere!")
print(f"  [OK] Ligand SDF: {docked_sdf}")

# List existing DCD files in Run 2 output for later reference
print(f"\n  DCD files in Run 2 output:")
for dcd in sorted(RUN2_RESULTS.glob("*.dcd")):
    print(f"    {dcd.name}: {dcd.stat().st_size/1e9:.2f} GB")

# %%
# ============================================================
# CELL 2: Rebuild System -> Load Checkpoint -> Run MD to 100 ns
# ============================================================
# MD-ONLY: No analysis. Just reach 100 ns and save.
# Time limit: use nearly full 12h for MD.
# ============================================================
from openmm.app import (
    ForceField, PDBFile, PME, HBonds,
    Simulation, StateDataReporter, DCDReporter, CheckpointReporter
)
from openmm import LangevinMiddleIntegrator, MonteCarloBarostat, Platform
from openmmforcefields.generators import GAFFTemplateGenerator
from openff.toolkit import Molecule as OFFMolecule
import openmm.unit as unit

# ============================================================
# TIME SAFETY -- Use nearly all 12h for MD
# Only 15 minutes reserved for saving final state
# ============================================================
MAX_WALL_SECONDS = 11.5 * 3600   # 11.5 hours total
SAVE_BUFFER = 0.25 * 3600        # 15 min for saving
MAX_MD_SECONDS = MAX_WALL_SECONDS - SAVE_BUFFER

print(f"[TIME]  TIME SAFETY: Started at {time.strftime('%Y-%m-%d %H:%M:%S')}")
print(f"  Max MD time: {MAX_MD_SECONDS/3600:.1f}h")
print(f"  Save buffer: {SAVE_BUFFER/60:.0f} min")

# --- Rebuild system ---
print("\nRebuilding simulation system...")

pdb = PDBFile(str(topology_file))
print(f"  Topology: {pdb.topology.getNumAtoms()} atoms")

off_mol = OFFMolecule.from_file(str(docked_sdf))
print(f"  Ligand: {off_mol.n_atoms} atoms")

ff = ForceField("amber/ff14SB.xml", "amber/tip3p_standard.xml")
gaff_gen = GAFFTemplateGenerator(molecules=[off_mol], forcefield="gaff-2.11")
ff.registerTemplateGenerator(gaff_gen.generator)
print("  [OK] ForceField + GAFF2 ready")

system = ff.createSystem(
    pdb.topology,
    nonbondedMethod=PME,
    nonbondedCutoff=MD_PARAMS["nonbonded_cutoff_nm"] * unit.nanometers,
    constraints=HBonds,
    rigidWater=True,
    hydrogenMass=1.5 * unit.amu,
)
system.addForce(MonteCarloBarostat(
    MD_PARAMS["pressure_atm"] * unit.atmospheres,
    MD_PARAMS["temperature_K"] * unit.kelvin, 25,
))

try:
    platform = Platform.getPlatformByName("CUDA")
    properties = {"DeviceIndex": "0", "Precision": "mixed"}
    print("  [OK] CUDA platform")
except Exception:
    platform = Platform.getPlatformByName("CPU")
    properties = {}
    print("  [WARN] CPU platform (slow!)")

integrator = LangevinMiddleIntegrator(
    MD_PARAMS["temperature_K"] * unit.kelvin,
    1.0 / unit.picoseconds,
    MD_PARAMS["timestep_fs"] * unit.femtoseconds,
)

simulation = Simulation(pdb.topology, system, integrator, platform, properties)
simulation.context.setPositions(pdb.positions)

# --- Load checkpoint ---
step_before = simulation.currentStep
print(f"\n{'=' * 60}")
print("CHECKPOINT VERIFICATION")
print(f"{'=' * 60}")
print(f"  Step BEFORE: {step_before:,}")
print(f"  Loading: {checkpoint_file.name}...")

simulation.loadCheckpoint(str(checkpoint_file))

step_after = simulation.currentStep
ns_done = step_after * MD_PARAMS["timestep_fs"] / 1e6
remaining_steps = MD_PARAMS["production_steps"] - step_after
remaining_ns = remaining_steps * MD_PARAMS["timestep_fs"] / 1e6

print(f"  Step AFTER: {step_after:,}")
print(f"  = {ns_done:.2f} ns completed")
print(f"  Remaining: {remaining_steps:,} steps ({remaining_ns:.1f} ns)")

# Anti-loop assertions
assert step_after > 0, f"[FAIL] ABORT: Step is 0! Checkpoint corrupt."
assert step_after > 5_000_000, f"[FAIL] ABORT: Step too low ({step_after:,}). Wrong checkpoint?"
assert step_after < MD_PARAMS["production_steps"], f"[FAIL] Already at 100 ns!"
assert step_after != step_before, f"[FAIL] Checkpoint didn't load!"

print(f"\n  [OK] Checkpoint valid ({ns_done:.1f} ns)")
print(f"  [OK] Need {remaining_ns:.1f} ns more to reach 100 ns")

# Estimate time needed
estimated_hours = remaining_ns / 3.75  # ~3.75 ns/h at 521 steps/s
print(f"  [TIME]  Estimated time: {estimated_hours:.1f}h (budget: {MAX_MD_SECONDS/3600:.1f}h)")
if estimated_hours > MAX_MD_SECONDS / 3600:
    print(f"  [WARN] May not finish! Will save checkpoint for next run.")

# --- Copy topology + existing DCDs to output ---
import shutil

for name in ["complex_solvated.pdb", "complex_equilibrated.pdb", "protein_fixed.pdb"]:
    src = RUN2_RESULTS / name
    dst = OUTPUT_DIR / name
    if src.exists() and not dst.exists():
        shutil.copy2(src, dst)

# Copy ligand SDF to output (needed for future analysis)
ligand_dst = OUTPUT_DIR / docked_sdf.name
if not ligand_dst.exists():
    shutil.copy2(docked_sdf, ligand_dst)

# Copy existing DCD files (for final trajectory combination)
for dcd_name in ["production_run1.dcd", "production_run2.dcd"]:
    src = RUN2_RESULTS / dcd_name
    dst = OUTPUT_DIR / dcd_name
    if src.exists() and not dst.exists():
        sz = src.stat().st_size / 1e9
        print(f"  Copying {dcd_name} ({sz:.2f} GB)...")
        shutil.copy2(src, dst)

# --- Setup reporters ---
run3_dcd = OUTPUT_DIR / "production_run3.dcd"
run3_log = OUTPUT_DIR / "production_run3.log"
run3_checkpoint = OUTPUT_DIR / "checkpoint.chk"

simulation.reporters.clear()

simulation.reporters.append(
    DCDReporter(str(run3_dcd), MD_PARAMS["save_interval_steps"], append=False)
)
simulation.reporters.append(
    StateDataReporter(
        str(run3_log), MD_PARAMS["log_interval_steps"],
        step=True, time=True, potentialEnergy=True,
        kineticEnergy=True, totalEnergy=True,
        temperature=True, volume=True,
        density=True, speed=True,
    )
)
simulation.reporters.append(
    CheckpointReporter(str(run3_checkpoint), MD_PARAMS["checkpoint_interval_steps"])
)

# ============================================================
# RUN MD (MD-ONLY -- no analysis after this)
# ============================================================
print(f"\n{'=' * 60}")
print(f"PRODUCTION MD -- RUN 3 (MD-ONLY)")
print(f"  From: {ns_done:.1f} ns -> To: {MD_PARAMS['production_ns']} ns")
print(f"  Remaining: {remaining_ns:.1f} ns ({remaining_steps:,} steps)")
print(f"  [TIME]  Time budget: {MAX_MD_SECONDS/3600:.1f}h")
print(f"{'=' * 60}")
print(f"\n[START] Resuming at {time.strftime('%Y-%m-%d %H:%M:%S')}...")

run_start_time = time.time()
run_start_step = step_after
total_target = MD_PARAMS["production_steps"]

PROGRESS_CHUNK = 50000  # every 100 ps
time_stopped_early = False

while simulation.currentStep < total_target:
    # Time safety check
    wall_elapsed = time.time() - t_start
    if wall_elapsed > MAX_MD_SECONDS:
        ns_now = simulation.currentStep * MD_PARAMS["timestep_fs"] / 1e6
        print(f"\n  [TIME]  TIME LIMIT after {wall_elapsed/3600:.1f}h")
        print(f"  Stopped at {ns_now:.1f} ns")
        time_stopped_early = True
        break

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
        wall_total = time.time() - t_start
        wall_left = MAX_MD_SECONDS - wall_total
        print(f"  [{pct:5.1f}%] {ns_now:.1f}/{MD_PARAMS['production_ns']} ns | "
              f"{rate:.0f} steps/s | ETA: {eta_sec/3600:.1f}h | "
              f"Wall: {wall_total/3600:.1f}h (limit in {wall_left/3600:.1f}h)", flush=True)

# ============================================================
# SAVE FINAL STATE (ALWAYS)
# ============================================================
final_chk = OUTPUT_DIR / "final.chk"
simulation.saveCheckpoint(str(final_chk))

state = simulation.context.getState(getPositions=True, getVelocities=True, getEnergy=True)
final_pdb = OUTPUT_DIR / "production_final.pdb"
with open(final_pdb, "w") as f:
    PDBFile.writeFile(simulation.topology, state.getPositions(), f)

total_elapsed = time.time() - run_start_time
ns_final = simulation.currentStep * MD_PARAMS["timestep_fs"] / 1e6

print(f"\n{'=' * 60}")
if time_stopped_early:
    ns_this_run = (simulation.currentStep - run_start_step) * MD_PARAMS["timestep_fs"] / 1e6
    print(f"[TIME]  MD stopped early (time safety)")
    print(f"  This run: +{ns_this_run:.1f} ns")
    print(f"  Total: {ns_final:.1f} / {MD_PARAMS['production_ns']} ns")
    print(f"  [PIN] Save output as dataset, run again to finish")
else:
    print(f"? 100 NS REACHED!")
    print(f"  Total simulated: {ns_final:.1f} ns")
    print(f"  [PIN] Now run analysis notebook on this output")
print(f"  Wall time: {total_elapsed/3600:.2f}h")
print(f"  Run 3 DCD: {run3_dcd}")
print(f"  Final checkpoint: {final_chk}")
print(f"  Finished: {time.strftime('%Y-%m-%d %H:%M:%S')}")
print(f"{'=' * 60}")

# List all output files
print(f"\n? Output files:")
total_size = 0
for f in sorted(OUTPUT_DIR.iterdir()):
    sz = f.stat().st_size
    total_size += sz
    unit_str = "GB" if sz > 1e9 else "MB" if sz > 1e6 else "kB"
    val = sz / 1e9 if sz > 1e9 else sz / 1e6 if sz > 1e6 else sz / 1e3
    print(f"  {f.name}: {val:.2f} {unit_str}")
print(f"  TOTAL: {total_size/1e9:.2f} GB")
