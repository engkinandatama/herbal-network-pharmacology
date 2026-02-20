# %% [markdown]
# # Supplementary: Core RMSD Analysis — Phalerin–AGTR1
# Computes RMSD for secondary structure core (Helix+Sheet) and binding site residues.
# 
# **Purpose**: Show that protein core remains stable despite flexible loop regions.
#
# **Input**: Dataset `md-phalerin-agtr1-final-output`
#
# **Runtime**: ~5-10 minutes (trajectory loading only, no heavy analysis)

# %%
# ============================================================
# CELL 0: Setup (minimal — only MDTraj needed)
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
    print(f"  Done ({time.time()-t_start:.0f}s)")

print("\n" + "=" * 60)
print("STEP 2: Mamba install (minimal)...")
print("=" * 60)
# Only need mdtraj + numpy for this lightweight analysis
packages = f"python={PY_VER} numpy mdtraj"
print(f"  -> {packages}")
r = subprocess.run(
    f"{MAMBA} install -y -c conda-forge {packages}",
    shell=True, capture_output=True, text=True, timeout=600
)
if r.returncode != 0:
    print(f"  mamba failed (exit {r.returncode})")
    for line in (r.stderr or "").strip().split('\n')[-5:]:
        print(f"    {line}")
else:
    print(f"  Installed ({time.time()-t_start:.0f}s)")

subprocess.run(f"{MAMBA} clean -afy 2>/dev/null", shell=True,
               capture_output=True, text=True)

print("\n" + "=" * 60)
print("STEP 3: Patch paths...")
print("=" * 60)

subprocess.run([sys.executable, "-m", "pip", "uninstall", "-y", "numpy"],
               capture_output=True, text=True)
for np_dir in glob.glob("/usr/local/lib/python3.12/dist-packages/numpy*"):
    shutil.rmtree(np_dir, ignore_errors=True)

mods_to_flush = [k for k in sys.modules if any(x in k.lower() for x in
                 ['numpy', 'mdtraj'])]
for mod in mods_to_flush:
    del sys.modules[mod]
print(f"  Flushed {len(mods_to_flush)} cached modules")

for sp in glob.glob(f"{MINIFORGE_DIR}/lib/python*/site-packages"):
    if sp not in sys.path:
        sys.path.insert(0, sp)
        print(f"  sys.path += {sp}")

os.environ["PATH"] = f"{MINIFORGE_DIR}/bin:" + os.environ.get("PATH", "")
os.environ["LD_LIBRARY_PATH"] = f"{MINIFORGE_DIR}/lib:" + os.environ.get("LD_LIBRARY_PATH", "")
print("  Environment set")

import numpy as np
import mdtraj as md
print(f"  numpy {np.__version__}, mdtraj {md.__version__}")
print(f"\n  Setup: {time.time()-t_start:.0f}s")

# %%
# ============================================================
# CELL 1: Load Trajectory
# ============================================================
from pathlib import Path
import warnings, gc
warnings.filterwarnings("ignore")

SYSTEM_NAME = "Phalerin_AGTR1"
LIGAND_NAME = "Phalerin"
SAVE_INTERVAL_PS = 10

# --- Find dataset ---
print("Searching for dataset...")
INPUT_ROOT = Path("/kaggle/input")
RESULTS_DIR = None

for name in ["md-phalerin-agtr1-final-output", "md-phalerin-agtr1-run3"]:
    for pattern in [
        INPUT_ROOT / name / "md_results" / SYSTEM_NAME,
        INPUT_ROOT / name / SYSTEM_NAME,
    ]:
        if pattern.exists() and list(pattern.glob("production_run*.dcd")):
            RESULTS_DIR = pattern
            break
    if RESULTS_DIR:
        break

if RESULTS_DIR is None:
    for dcd in INPUT_ROOT.rglob("production_run1.dcd"):
        RESULTS_DIR = dcd.parent
        break

assert RESULTS_DIR is not None, "Dataset not found!"

OUTPUT_DIR = Path("/kaggle/working/analysis") / SYSTEM_NAME
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

print(f"Dataset: {RESULTS_DIR}")
print("=" * 60)

# --- Load trajectory ---
topo_file = RESULTS_DIR / "complex_solvated.pdb"
if not topo_file.exists():
    topo_file = RESULTS_DIR / "complex_equilibrated.pdb"

STRIDE = 2
EFFECTIVE_DT_PS = SAVE_INTERVAL_PS * STRIDE

dcd_files = sorted(RESULTS_DIR.glob("production_run*.dcd"))
print(f"Found {len(dcd_files)} DCD files (stride={STRIDE}):")

trajs = []
for dcd in dcd_files:
    t = md.load(str(dcd), top=str(topo_file), stride=STRIDE)
    print(f"  {dcd.name}: {t.n_frames} frames")
    trajs.append(t)
    gc.collect()

# Override time axis
if len(trajs) > 0:
    dt_ps = EFFECTIVE_DT_PS
    trajs[0].time = np.arange(trajs[0].n_frames) * dt_ps
    for i in range(1, len(trajs)):
        t_start_ps = trajs[i-1].time[-1] + dt_ps
        trajs[i].time = np.arange(trajs[i].n_frames) * dt_ps + t_start_ps

traj = md.join(trajs) if len(trajs) > 1 else trajs[0]
del trajs; gc.collect()

# Trim to 100 ns
TARGET_NS = 100.0
mask = traj.time <= TARGET_NS * 1000
if not np.all(mask):
    traj = traj[mask]
    gc.collect()

time_ns = traj.time / 1000.0
print(f"\nLoaded: {traj.n_frames:,} frames, {time_ns[-1]:.1f} ns")

# --- Atom selections ---
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

print(f"  Protein: {len(protein_atoms)} atoms")
print(f"  Backbone: {len(backbone_atoms)} atoms")
print(f"  Ligand: {len(ligand_atoms)} atoms")

# %%
# ============================================================
# CELL 2: Core RMSD Analysis
# ============================================================
print("=" * 60)
print("CORE RMSD ANALYSIS")
print("=" * 60)

# --- 1. Standard backbone RMSD (for comparison) ---
print("\n1. Computing standard backbone RMSD...")
prot_traj = traj.atom_slice(protein_atoms)
rmsd_backbone = md.rmsd(prot_traj, prot_traj, frame=0, atom_indices=None) * 10  # nm -> A
# Use backbone atoms relative to protein-only trajectory
bb_local = prot_traj.topology.select("backbone")
rmsd_backbone_proper = md.rmsd(prot_traj, prot_traj, frame=0, atom_indices=bb_local) * 10
print(f"  Backbone RMSD: {np.mean(rmsd_backbone_proper):.2f} +/- {np.std(rmsd_backbone_proper):.2f} A")

# --- 2. DSSP: Identify secondary structure ---
print("\n2. Computing DSSP secondary structure assignment...")
# Compute DSSP on first frame (reference) to define core residues
dssp_ref = md.compute_dssp(prot_traj[0], simplified=True)  # 'H', 'E', 'C'
dssp_ref = dssp_ref[0]  # first frame

n_helix = np.sum(dssp_ref == 'H')
n_sheet = np.sum(dssp_ref == 'E')
n_coil = np.sum(dssp_ref == 'C')
n_total = len(dssp_ref)
print(f"  Total residues: {n_total}")
print(f"  Helix (H): {n_helix} ({100*n_helix/n_total:.1f}%)")
print(f"  Sheet (E): {n_sheet} ({100*n_sheet/n_total:.1f}%)")
print(f"  Coil  (C): {n_coil} ({100*n_coil/n_total:.1f}%)")

# Get residue indices for core (H+E)
core_residue_indices = np.where((dssp_ref == 'H') | (dssp_ref == 'E'))[0]
print(f"  Core residues (H+E): {len(core_residue_indices)}")

# Map to atom indices (backbone atoms of core residues)
core_backbone_atoms = []
for res_idx in core_residue_indices:
    res = prot_traj.topology.residue(res_idx)
    for atom in res.atoms:
        if atom.name in ['CA', 'C', 'N', 'O']:
            core_backbone_atoms.append(atom.index)
core_backbone_atoms = np.array(core_backbone_atoms)
print(f"  Core backbone atoms: {len(core_backbone_atoms)}")

# --- 3. RMSD on core residues only ---
print("\n3. Computing Core RMSD (Helix + Sheet only)...")
rmsd_core = md.rmsd(prot_traj, prot_traj, frame=0, atom_indices=core_backbone_atoms) * 10
print(f"  Core RMSD: {np.mean(rmsd_core):.2f} +/- {np.std(rmsd_core):.2f} A")

# --- 4. Binding site RMSD ---
print("\n4. Computing Binding Site RMSD...")
# Top contact residues from analysis (resSeq numbers)
binding_site_resseq = [30, 31, 105, 106, 107, 108, 110, 111, 112, 113, 177, 181]
print(f"  Binding site residues (resSeq): {binding_site_resseq}")

# Map resSeq to residue indices + backbone atoms
binding_backbone_atoms = []
binding_residues_found = []
for res in prot_traj.topology.residues:
    if res.resSeq in binding_site_resseq:
        binding_residues_found.append(f"{res.name}{res.resSeq}")
        for atom in res.atoms:
            if atom.name in ['CA', 'C', 'N', 'O']:
                binding_backbone_atoms.append(atom.index)

binding_backbone_atoms = np.array(binding_backbone_atoms)
print(f"  Found: {', '.join(binding_residues_found)} ({len(binding_backbone_atoms)} backbone atoms)")

rmsd_binding = md.rmsd(prot_traj, prot_traj, frame=0, atom_indices=binding_backbone_atoms) * 10
print(f"  Binding Site RMSD: {np.mean(rmsd_binding):.2f} +/- {np.std(rmsd_binding):.2f} A")

# --- 5. Also compute DSSP over time for secondary structure stability ---
print("\n5. Computing DSSP over time (every 10th frame for speed)...")
dssp_stride = 10
dssp_frames = range(0, prot_traj.n_frames, dssp_stride)
dssp_all = md.compute_dssp(prot_traj[list(dssp_frames)], simplified=True)
# Count H, E, C per frame
ss_time_ns = time_ns[list(dssp_frames)]
n_helix_t = np.sum(dssp_all == 'H', axis=1)
n_sheet_t = np.sum(dssp_all == 'E', axis=1)
n_coil_t = np.sum(dssp_all == 'C', axis=1)
print(f"  DSSP computed for {len(dssp_frames)} frames")
print(f"  Helix range: {np.min(n_helix_t)}-{np.max(n_helix_t)} residues")
print(f"  Sheet range: {np.min(n_sheet_t)}-{np.max(n_sheet_t)} residues")

# --- 6. Summary ---
print(f"\n{'='*60}")
print(f"CORE RMSD SUMMARY")
print(f"{'='*60}")
print(f"  Full Backbone RMSD:    {np.mean(rmsd_backbone_proper):.2f} +/- {np.std(rmsd_backbone_proper):.2f} A")
print(f"  Core (H+E) RMSD:      {np.mean(rmsd_core):.2f} +/- {np.std(rmsd_core):.2f} A")
print(f"  Binding Site RMSD:     {np.mean(rmsd_binding):.2f} +/- {np.std(rmsd_binding):.2f} A")
print(f"  Reduction (core):      {(1 - np.mean(rmsd_core)/np.mean(rmsd_backbone_proper))*100:.1f}% lower")
print(f"  Reduction (binding):   {(1 - np.mean(rmsd_binding)/np.mean(rmsd_backbone_proper))*100:.1f}% lower")

# --- Free trajectory memory (~4.7 GB) before plotting ---
# All RMSD data is already in numpy arrays, trajectory no longer needed
import gc
try: del traj
except: pass
try: del prot_traj
except: pass
gc.collect()
print("\nTrajectory freed from memory for plotting")

# %%
# ============================================================
# CELL 3: Save CSV & Generate Plots
# ============================================================
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --- Save core RMSD timeseries ---
df_rmsd = pd.DataFrame({
    "Time_ns": time_ns,
    "RMSD_Backbone_A": rmsd_backbone_proper,
    "RMSD_Core_HE_A": rmsd_core,
    "RMSD_BindingSite_A": rmsd_binding,
})
df_rmsd.to_csv(OUTPUT_DIR / "core_rmsd_timeseries.csv", index=False)
print(f"Saved: core_rmsd_timeseries.csv ({len(df_rmsd)} frames)")

# Save DSSP over time
df_dssp = pd.DataFrame({
    "Time_ns": ss_time_ns,
    "N_Helix": n_helix_t,
    "N_Sheet": n_sheet_t,
    "N_Coil": n_coil_t,
})
df_dssp.to_csv(OUTPUT_DIR / "dssp_over_time.csv", index=False)
print(f"Saved: dssp_over_time.csv ({len(df_dssp)} frames)")

# Save binding site residue list
df_binding = pd.DataFrame({
    "Residue": binding_residues_found,
    "ResSeq": binding_site_resseq[:len(binding_residues_found)],
})
df_binding.to_csv(OUTPUT_DIR / "binding_site_residues.csv", index=False)

# Save core residue info with DSSP
core_res_info = []
for i, res_idx in enumerate(range(n_total)):
    res = prot_traj.topology.residue(res_idx)
    core_res_info.append({
        "Residue": f"{res.name}{res.resSeq}",
        "ResSeq": res.resSeq,
        "ResName": res.name,
        "DSSP": dssp_ref[res_idx],
        "Is_Core": dssp_ref[res_idx] in ['H', 'E'],
        "Is_BindingSite": res.resSeq in binding_site_resseq,
    })
df_core_info = pd.DataFrame(core_res_info)
df_core_info.to_csv(OUTPUT_DIR / "residue_classification.csv", index=False)
print(f"Saved: residue_classification.csv ({len(df_core_info)} residues)")

# --- Plot ---
plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 10,
    "axes.labelsize": 12, "axes.titlesize": 13,
    "xtick.labelsize": 9, "ytick.labelsize": 9,
    "legend.fontsize": 9, "figure.dpi": 300, "axes.linewidth": 1.2,
})

window = max(1, len(time_ns) // 100)

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# (A) RMSD Comparison: Backbone vs Core vs Binding Site
ax = axes[0, 0]
for data, label, color, alpha_raw in [
    (rmsd_backbone_proper, "Full Backbone", "#2196F3", 0.15),
    (rmsd_core, "Core (Helix+Sheet)", "#4CAF50", 0.2),
    (rmsd_binding, "Binding Site", "#FF9800", 0.25),
]:
    ax.plot(time_ns, data, color=color, linewidth=0.3, alpha=alpha_raw)
    smooth = np.convolve(data, np.ones(window)/window, mode="valid")
    ax.plot(time_ns[:len(smooth)], smooth, color=color, linewidth=2.0, label=label)

ax.set_xlabel("Time (ns)"); ax.set_ylabel("RMSD (A)")
ax.set_title("(A) RMSD Comparison")
ax.legend(fontsize=8)
ax.set_xlim(0, time_ns[-1])

# (B) Bar chart: mean RMSD comparison
ax = axes[0, 1]
categories = ["Full\nBackbone", "Core\n(Helix+Sheet)", "Binding\nSite"]
means = [np.mean(rmsd_backbone_proper), np.mean(rmsd_core), np.mean(rmsd_binding)]
stds = [np.std(rmsd_backbone_proper), np.std(rmsd_core), np.std(rmsd_binding)]
colors = ["#2196F3", "#4CAF50", "#FF9800"]

bars = ax.bar(categories, means, yerr=stds, capsize=5, color=colors,
              edgecolor="black", linewidth=0.5, alpha=0.85)
for bar, m, s in zip(bars, means, stds):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + s + 0.05,
            f"{m:.2f}+/-{s:.2f}", ha="center", va="bottom", fontsize=9, fontweight="bold")
ax.set_ylabel("RMSD (A)")
ax.set_title("(B) Mean RMSD Comparison")

# (C) Secondary structure over time
ax = axes[1, 0]
ax.fill_between(ss_time_ns, 0, n_helix_t, color="#E91E63", alpha=0.6, label="Helix")
ax.fill_between(ss_time_ns, n_helix_t, n_helix_t + n_sheet_t, color="#3F51B5", alpha=0.6, label="Sheet")
ax.fill_between(ss_time_ns, n_helix_t + n_sheet_t, n_helix_t + n_sheet_t + n_coil_t,
                color="#9E9E9E", alpha=0.4, label="Coil")
ax.set_xlabel("Time (ns)"); ax.set_ylabel("Number of Residues")
ax.set_title("(C) Secondary Structure over Time")
ax.legend(fontsize=8)
ax.set_xlim(0, time_ns[-1])

# (D) DSSP per-residue map (reference frame)
ax = axes[1, 1]
dssp_numeric = np.zeros(n_total)
dssp_numeric[dssp_ref == 'H'] = 2
dssp_numeric[dssp_ref == 'E'] = 1
dssp_numeric[dssp_ref == 'C'] = 0

resseq = [prot_traj.topology.residue(i).resSeq for i in range(n_total)]
from matplotlib.colors import ListedColormap
cmap = ListedColormap(["#9E9E9E", "#3F51B5", "#E91E63"])
ax.bar(resseq, 1, bottom=0, color=[cmap(int(v)) for v in dssp_numeric], width=1.0)
# Mark binding site residues
for rs in binding_site_resseq:
    if rs in resseq:
        ax.axvline(rs, color="orange", linewidth=0.8, alpha=0.7, linestyle="--")
ax.set_xlabel("Residue Number"); ax.set_ylabel("")
ax.set_yticks([])
ax.set_title("(D) DSSP Map (ref frame) | Orange = Binding Site")

# Legend for DSSP
from matplotlib.patches import Patch
ax.legend(handles=[
    Patch(facecolor="#E91E63", label="Helix"),
    Patch(facecolor="#3F51B5", label="Sheet"),
    Patch(facecolor="#9E9E9E", label="Coil"),
], fontsize=8, loc="upper right")

PRODUCTION_NS = time_ns[-1]
fig.suptitle(f"Core RMSD Analysis: {LIGAND_NAME} - AGTR1 | {PRODUCTION_NS:.0f} ns",
             fontsize=14, fontweight="bold", y=0.98)
plt.tight_layout()
plt.savefig(OUTPUT_DIR / "core_rmsd_analysis.png", dpi=300, bbox_inches="tight", facecolor="white")
plt.close()
print("Saved: core_rmsd_analysis.png")

# --- Final report ---
print(f"\n{'='*60}")
print(f"SUPPLEMENTARY ANALYSIS COMPLETE: {SYSTEM_NAME}")
print(f"{'='*60}")
print(f"""
Results:
  Full Backbone RMSD:  {np.mean(rmsd_backbone_proper):.2f} +/- {np.std(rmsd_backbone_proper):.2f} A
  Core (H+E) RMSD:    {np.mean(rmsd_core):.2f} +/- {np.std(rmsd_core):.2f} A
  Binding Site RMSD:   {np.mean(rmsd_binding):.2f} +/- {np.std(rmsd_binding):.2f} A

Secondary Structure (reference frame):
  Helix: {n_helix} residues ({100*n_helix/n_total:.1f}%)
  Sheet: {n_sheet} residues ({100*n_sheet/n_total:.1f}%)
  Coil:  {n_coil} residues ({100*n_coil/n_total:.1f}%)

Output files:
  - core_rmsd_timeseries.csv
  - dssp_over_time.csv
  - residue_classification.csv
  - binding_site_residues.csv
  - core_rmsd_analysis.png

Total time: {(time.time()-t_start)/60:.1f} min
""")
