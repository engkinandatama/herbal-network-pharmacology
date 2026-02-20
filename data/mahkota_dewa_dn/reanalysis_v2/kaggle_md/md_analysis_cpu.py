# %% [markdown]
# # MD Analysis — Phalerin–AGTR1 (100 ns) — CPU Only
# Full trajectory analysis from 3 production runs.
#
# **CPU-ONLY**: No GPU needed. MM-GBSA runs in separate GPU notebook.
#
# **Input**: Dataset `md-phalerin-agtr1-final-output`
#   - production_run1.dcd (0→40 ns)
#   - production_run2.dcd (40→75.6 ns)
#   - production_run3.dcd (75.6→100 ns)

# %%
# ============================================================
# CELL 0: Install Dependencies (CPU only — no CUDA needed)
# ============================================================
import subprocess, sys, os, glob, time, shutil

t_start = time.time()
PY_VER = f"{sys.version_info.major}.{sys.version_info.minor}"
print(f"Kernel Python: {PY_VER}")

# --- Miniforge ---
MINIFORGE_DIR = "/tmp/miniforge"
MAMBA = f"{MINIFORGE_DIR}/bin/mamba"

if not os.path.exists(MAMBA):
    url = "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
    print("Downloading Miniforge3...")
    subprocess.run(f"wget -q {url} -O /tmp/miniforge.sh", shell=True, check=True)
    subprocess.run(f"bash /tmp/miniforge.sh -b -p {MINIFORGE_DIR}",
                   shell=True, capture_output=True, text=True, check=True)
    os.remove("/tmp/miniforge.sh")
    print(f"  ✅ Miniforge3 installed ({time.time()-t_start:.0f}s)")

# --- Install ONLY analysis deps (no openmm/cuda) ---
print("Installing analysis packages...")
packages = f"python={PY_VER} numpy mdtraj parmed"
r = subprocess.run(f"{MAMBA} install -y -c conda-forge {packages}",
                   shell=True, capture_output=True, text=True, timeout=600)
if r.returncode != 0:
    print(f"⚠️ mamba failed: {r.stderr[-300:]}")
else:
    print(f"  ✅ Conda packages installed")

# Remove conda scipy (use pip version)
subprocess.run(f"{MAMBA} remove -y --force scipy 2>/dev/null",
               shell=True, capture_output=True, text=True)
for p in glob.glob(f"{MINIFORGE_DIR}/lib/python*/site-packages/scipy*"):
    shutil.rmtree(p, ignore_errors=True)

# Pip installs
for pkg in ["scipy", "matplotlib", "pandas"]:
    subprocess.run([sys.executable, "-m", "pip", "install", "-q", pkg],
                   capture_output=True, text=True, timeout=120)
    print(f"  ✅ pip: {pkg}")

# --- Patch paths ---
subprocess.run([sys.executable, "-m", "pip", "uninstall", "-y", "numpy"],
               capture_output=True, text=True)
for np_dir in glob.glob("/usr/local/lib/python3.12/dist-packages/numpy*"):
    shutil.rmtree(np_dir, ignore_errors=True)
mods_to_flush = [k for k in sys.modules if 'numpy' in k.lower()]
for mod in mods_to_flush:
    del sys.modules[mod]

for sp in glob.glob(f"{MINIFORGE_DIR}/lib/python*/site-packages"):
    if sp not in sys.path:
        sys.path.insert(0, sp)

os.environ["PATH"] = f"{MINIFORGE_DIR}/bin:" + os.environ.get("PATH", "")
os.environ["LD_LIBRARY_PATH"] = f"{MINIFORGE_DIR}/lib:" + os.environ.get("LD_LIBRARY_PATH", "")

# --- Verify ---
import numpy as np; print(f"  ✅ numpy {np.__version__}")
import mdtraj as md; print(f"  ✅ mdtraj {md.__version__}")
from scipy import stats; print(f"  ✅ scipy")
import matplotlib; print(f"  ✅ matplotlib")
import pandas as pd; print(f"  ✅ pandas")
print(f"\nSetup: {time.time()-t_start:.0f}s ✅")

# %%
# ============================================================
# CELL 1: Load & Combine Trajectories (100 ns)
# ============================================================
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

SYSTEM_NAME = "Phalerin_AGTR1"
LIGAND_NAME = "Phalerin"
RECEPTOR_PDB_ID = "6OS1"

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

# Strategy 2: Recursive search for DCD files under any input folder
if RESULTS_DIR is None:
    print("  Known slugs not found, searching recursively...")
    for dcd in INPUT_ROOT.rglob("production_run1.dcd"):
        RESULTS_DIR = dcd.parent
        break

# Strategy 3: Search for Phalerin_AGTR1 folder anywhere
if RESULTS_DIR is None:
    for d in INPUT_ROOT.rglob("Phalerin_AGTR1"):
        if d.is_dir():
            RESULTS_DIR = d
            break

if RESULTS_DIR is None:
    # Print debug info
    print("  ❌ Dataset not found! Listing /kaggle/input contents:")
    for p in sorted(INPUT_ROOT.rglob("*")):
        if p.is_file() and p.stat().st_size > 1e6:
            print(f"    {p.relative_to(INPUT_ROOT)} ({p.stat().st_size/1e6:.1f} MB)")
        elif p.is_dir():
            print(f"    {p.relative_to(INPUT_ROOT)}/ (dir)")
    raise FileNotFoundError("Cannot find DCD trajectory files in /kaggle/input!")

OUTPUT_DIR = Path("/kaggle/working/analysis") / SYSTEM_NAME
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# MD parameters (must match simulation)
SAVE_INTERVAL_PS = 10  # DCD save frequency (actual)
STRIDE = 2  # Load every Nth frame to save memory (Kaggle CPU = 13 GB RAM)
EFFECTIVE_DT_PS = SAVE_INTERVAL_PS * STRIDE  # 20 ps between loaded frames

print("=" * 60)
print(f"ANALYSIS: {SYSTEM_NAME} (100 ns)")
print(f"Dataset: {RESULTS_DIR}")
print(f"Memory optimization: stride={STRIDE} (every {STRIDE}th frame)")
print("=" * 60)

# --- Find topology ---
topo_file = RESULTS_DIR / "complex_solvated.pdb"
if not topo_file.exists():
    topo_file = RESULTS_DIR / "complex_equilibrated.pdb"
print(f"\nTopology: {topo_file.name}")

# --- Load all DCDs with stride ---
dcd_files = sorted(RESULTS_DIR.glob("production_run*.dcd"))
print(f"Found {len(dcd_files)} DCD files:")

import gc
trajs = []
for dcd in dcd_files:
    t = md.load(str(dcd), top=str(topo_file), stride=STRIDE)
    sz_gb = dcd.stat().st_size / 1e9
    print(f"  {dcd.name}: {t.n_frames} frames (stride={STRIDE}) ({sz_gb:.2f} GB on disk)")
    trajs.append(t)
    gc.collect()

# --- Override time axis (DCD header often stores wrong timestep) ---
if len(trajs) > 0:
    dt_ps = EFFECTIVE_DT_PS
    trajs[0].time = np.arange(trajs[0].n_frames) * dt_ps
    print(f"  {dcd_files[0].name} time → 0-{trajs[0].time[-1]/1000:.1f} ns")
    for i in range(1, len(trajs)):
        t_start_ps = trajs[i-1].time[-1] + dt_ps
        trajs[i].time = np.arange(trajs[i].n_frames) * dt_ps + t_start_ps
        print(f"  {dcd_files[i].name} time → {trajs[i].time[0]/1000:.1f}-{trajs[i].time[-1]/1000:.1f} ns")

# --- Join ---
if len(trajs) > 1:
    traj = md.join(trajs)
else:
    traj = trajs[0]

# Free memory from individual trajectories
del trajs
gc.collect()

# --- Trim to exactly 100 ns ---
TARGET_NS = 100.0
target_ps = TARGET_NS * 1000
mask_100ns = traj.time <= target_ps
if not np.all(mask_100ns):
    n_before = traj.n_frames
    traj = traj[mask_100ns]
    gc.collect()
    print(f"  ✂️  Trimmed {n_before} → {traj.n_frames} frames (capped at {TARGET_NS:.0f} ns)")

total_ns = (traj.time[-1] - traj.time[0]) / 1000
mem_est_gb = traj.n_frames * traj.n_atoms * 3 * 4 / 1e9
print(f"\n📊 Combined trajectory:")
print(f"  Frames: {traj.n_frames:,}")
print(f"  Time: {traj.time[0]/1000:.1f} - {traj.time[-1]/1000:.1f} ns")
print(f"  Duration: {total_ns:.1f} ns")
print(f"  Atoms: {traj.n_atoms:,}")
print(f"  Memory est: ~{mem_est_gb:.1f} GB")

# %%
# ============================================================
# CELL 2: Full Trajectory Analysis (10 metrics)
# ============================================================
print("=" * 60)
print("TRAJECTORY ANALYSIS")
print("=" * 60)

# --- [0/10] PBC Reimaging ---
print("\n[0/10] Reimaging trajectory...")
try:
    # Find molecules in topology — protein is usually the largest
    molecules = traj.topology.find_molecules()
    # Sort by size (largest first = protein)
    molecules_sorted = sorted(molecules, key=len, reverse=True)
    # Anchor = protein (largest molecule)
    anchor = [molecules_sorted[0]]
    traj = traj.image_molecules(inplace=False, anchor_molecules=anchor)
    print("  ✅ Reimaging complete")
except Exception as e:
    print(f"  ⚠️ Reimaging skipped ({e}) — analysis still valid")

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
print(f"  Cα: {len(ca_atoms)} atoms")
print(f"  Ligand: {len(ligand_atoms)} atoms")

assert len(ligand_atoms) > 0, "No ligand atoms!"
assert len(ca_atoms) > 0, "No CA atoms!"

time_ns = traj.time / 1000.0

# --- [1/10] Backbone RMSD ---
print("\n[1/10] Backbone RMSD...")
t1 = time.time()
traj_bb = traj.atom_slice(backbone_atoms)
rmsd_protein = md.rmsd(traj_bb, traj_bb, frame=0) * 10  # nm → Å
print(f"  Mean: {np.mean(rmsd_protein):.2f} ± {np.std(rmsd_protein):.2f} Å ({time.time()-t1:.1f}s)")

# --- [2/10] Ligand RMSD ---
print("[2/10] Ligand RMSD...")
t1 = time.time()
traj_aligned = traj.superpose(traj, frame=0, atom_indices=ca_atoms)
lig_traj = traj_aligned.atom_slice(ligand_atoms)
rmsd_ligand = md.rmsd(lig_traj, lig_traj, frame=0) * 10
print(f"  Mean: {np.mean(rmsd_ligand):.2f} ± {np.std(rmsd_ligand):.2f} Å ({time.time()-t1:.1f}s)")

# --- [3/10] RMSF ---
print("[3/10] Cα RMSF...")
t1 = time.time()
ca_traj = traj_aligned.atom_slice(ca_atoms)
rmsf = md.rmsf(ca_traj, ca_traj, frame=0) * 10
residue_ids = [traj.topology.atom(a).residue.resSeq for a in ca_atoms]
residue_names = [traj.topology.atom(a).residue.name for a in ca_atoms]
top_rmsf_idx = np.argmax(rmsf)
print(f"  Max: {np.max(rmsf):.2f} Å (residue {residue_names[top_rmsf_idx]}{residue_ids[top_rmsf_idx]}) ({time.time()-t1:.1f}s)")

# --- [4/10] Radius of Gyration ---
print("[4/10] Radius of Gyration...")
t1 = time.time()
rg = md.compute_rg(traj.atom_slice(protein_atoms)) * 10
print(f"  Mean: {np.mean(rg):.2f} ± {np.std(rg):.2f} Å ({time.time()-t1:.1f}s)")

# --- [5/10] H-bonds ---
print("\n[5/10] Protein-ligand H-bonds...")
t1 = time.time()
lig_set = set(ligand_atoms.tolist())
prot_set = set(protein_atoms.tolist())
sample_step = max(1, traj.n_frames // 1000)
sampled_indices = list(range(0, traj.n_frames, sample_step))

hbond_counts_sampled = []
for idx in sampled_indices:
    try:
        hbs = md.baker_hubbard(traj[idx], freq=0.0)
        count = sum(1 for hb in hbs
                    if (hb[0] in lig_set and hb[2] in prot_set) or
                       (hb[0] in prot_set and hb[2] in lig_set))
        hbond_counts_sampled.append(count)
    except:
        hbond_counts_sampled.append(0)

hbond_counts = np.interp(np.arange(traj.n_frames), sampled_indices,
                         np.array(hbond_counts_sampled))
print(f"  Sampled {len(sampled_indices)} frames")
print(f"  Mean: {np.mean(hbond_counts):.1f} ± {np.std(hbond_counts):.1f} ({time.time()-t1:.1f}s)")

# --- [6/10] Contact Frequency ---
print("\n[6/10] Contact frequency...")
t1 = time.time()
CONTACT_CUTOFF_NM = 0.45
residue_contact_counts = {}

for idx in sampled_indices:
    frame = traj[idx]
    for ca_idx in ca_atoms:
        res = traj.topology.atom(ca_idx).residue
        res_atoms = [a.index for a in res.atoms]
        pairs = [(r, l) for r in res_atoms for l in ligand_atoms]
        if not pairs:
            continue
        distances = md.compute_distances(frame, pairs)[0]
        if np.min(distances) < CONTACT_CUTOFF_NM:
            key = (res.resSeq, res.name)
            residue_contact_counts[key] = residue_contact_counts.get(key, 0) + 1

n_sampled = len(sampled_indices)
contact_freq = {k: v / n_sampled for k, v in residue_contact_counts.items()}
contact_freq_sorted = sorted(contact_freq.items(), key=lambda x: -x[1])
top_contacts = [(k, v) for k, v in contact_freq_sorted if v > 0.1]

print(f"  Top residues (>10% contact):")
for (resSeq, resName), freq in top_contacts[:15]:
    print(f"    {resName}{resSeq}: {freq*100:.1f}%")
print(f"  ({time.time()-t1:.1f}s)")

# --- [7/10] Key Distances ---
print("\n[7/10] Key residue-ligand distances...")
t1 = time.time()
lig_com_distances = {}
top3_residues = [(rs, rn) for (rs, rn), _ in contact_freq_sorted[:3]]

for resSeq, resName in top3_residues:
    res_atom_indices = [a.index for a in traj.topology.atoms
                        if a.residue.resSeq == resSeq and a.residue.name == resName]
    if not res_atom_indices:
        continue
    dists = []
    for idx in range(0, traj.n_frames, sample_step):
        pairs = [(r, l) for r in res_atom_indices for l in ligand_atoms]
        d = md.compute_distances(traj[idx], pairs)[0]
        dists.append(np.min(d) * 10)  # nm → Å

    dists_full = np.interp(np.arange(traj.n_frames),
                           list(range(0, traj.n_frames, sample_step)), dists)
    lig_com_distances[f"{resName}{resSeq}"] = dists_full
    print(f"  {resName}{resSeq}: {np.mean(dists_full):.2f} ± {np.std(dists_full):.2f} Å")
print(f"  ({time.time()-t1:.1f}s)")

# --- [8/10] SASA ---
print("\n[8/10] Protein SASA...")
t1 = time.time()
sasa_sample_step = max(1, traj.n_frames // 500)
sasa_indices = list(range(0, traj.n_frames, sasa_sample_step))

prot_traj = traj.atom_slice(protein_atoms)
sasa_sampled = []
for idx in sasa_indices:
    frame_sasa = md.shrake_rupley(prot_traj[idx], mode='residue')
    sasa_sampled.append(np.sum(frame_sasa))

sasa_sampled = np.array(sasa_sampled) * 100  # nm² → Å²
sasa_full = np.interp(np.arange(traj.n_frames), sasa_indices, sasa_sampled)
print(f"  Mean: {np.mean(sasa_full):.1f} ± {np.std(sasa_full):.1f} Å² ({time.time()-t1:.1f}s)")

# --- [9/10] Equilibration Check ---
print("\n[9/10] Equilibration check...")
PRODUCTION_NS = time_ns[-1]
n_blocks = 5
block_size = len(rmsd_protein) // n_blocks
block_means = [np.mean(rmsd_protein[i*block_size:(i+1)*block_size]) for i in range(n_blocks)]

print(f"  RMSD block means:")
for i, bm in enumerate(block_means):
    ns_block = PRODUCTION_NS / n_blocks
    print(f"    Block {i+1} ({i*ns_block:.0f}-{(i+1)*ns_block:.0f} ns): {bm:.2f} Å")

last3_mean = np.mean(block_means[2:])
first_dev = abs(block_means[0] - last3_mean) / last3_mean * 100
if first_dev > 50:
    equil_ns = PRODUCTION_NS / n_blocks
    print(f"  ⚠️ First block deviates {first_dev:.0f}% — consider skipping first {equil_ns:.0f} ns")
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

# --- Save CSVs ---
print("\n💾 Saving data...")
import pandas as pd

ts_df = pd.DataFrame({
    "Time_ns": time_ns, "RMSD_Protein_A": rmsd_protein,
    "RMSD_Ligand_A": rmsd_ligand, "Rg_A": rg,
    "Hbonds": hbond_counts, "SASA_A2": sasa_full,
})
ts_df.to_csv(OUTPUT_DIR / "timeseries.csv", index=False)

rmsf_df = pd.DataFrame({"Residue": residue_ids, "ResName": residue_names, "RMSF_A": rmsf})
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
print(f"\n⏱️  Analysis completed in {(time.time()-t_start)/60:.1f} min")

# %%
# ============================================================
# CELL 3: Publication Figures
# ============================================================

# Free trajectory memory (~4.7 GB) — all data is already in arrays
import gc
try: del traj
except: pass
try: del prot_traj
except: pass
gc.collect()
print("🧹 Trajectory freed from memory for plotting")

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

# ===============================
# Figure 1: 4-Panel (RMSD, RMSF, Rg)
# ===============================
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

fig.suptitle(f"MD Simulation: {LIGAND_NAME} – {RECEPTOR_PDB_ID} (AGTR1) | {PRODUCTION_NS:.0f} ns",
             fontsize=14, fontweight="bold", y=0.98)
plt.savefig(OUTPUT_DIR / "md_analysis_4panel.png", dpi=300, bbox_inches="tight", facecolor="white")
plt.close()
print("✅ Figure 1: 4-panel saved")

# ===============================
# Figure 2: H-bond + SASA
# ===============================
fig2, (ax5, ax6) = plt.subplots(1, 2, figsize=(14, 5))

ax5.plot(time_ns, hbond_counts, color="#F44336", linewidth=0.3, alpha=0.4)
hb_smooth = np.convolve(hbond_counts, np.ones(window)/window, mode="valid")
ax5.plot(time_smooth, hb_smooth, color="#B71C1C", linewidth=1.5, label="Running avg")
ax5.set_xlabel("Time (ns)"); ax5.set_ylabel("Number of H-bonds")
ax5.set_title(f"(E) Protein–Ligand H-bonds\nmean={np.mean(hbond_counts):.1f}±{np.std(hbond_counts):.1f}")
ax5.legend(); ax5.set_xlim(0, time_ns[-1])

ax6.plot(time_ns, sasa_full, color="#607D8B", linewidth=0.3, alpha=0.4)
sasa_smooth = np.convolve(sasa_full, np.ones(window)/window, mode="valid")
ax6.plot(time_smooth, sasa_smooth, color="#37474F", linewidth=1.5, label="Running avg")
ax6.set_xlabel("Time (ns)"); ax6.set_ylabel("SASA (Å²)")
ax6.set_title(f"(F) Protein SASA\nmean={np.mean(sasa_full):.0f}±{np.std(sasa_full):.0f} Å²")
ax6.legend(); ax6.set_xlim(0, time_ns[-1])

fig2.suptitle(f"H-bond & SASA: {LIGAND_NAME} – {RECEPTOR_PDB_ID}", fontsize=14, fontweight="bold")
plt.tight_layout()
plt.savefig(OUTPUT_DIR / "md_hbond_sasa.png", dpi=300, bbox_inches="tight", facecolor="white")
plt.close()
print("✅ Figure 2: H-bond + SASA saved")

# ===============================
# Figure 3: Contacts + Distances
# ===============================
fig3 = plt.figure(figsize=(14, 5))
gs3 = GridSpec(1, 2, wspace=0.35)

ax7 = fig3.add_subplot(gs3[0, 0])
top_n = min(15, len(top_contacts))
if top_n > 0:
    labels_c = [f"{rn}{rs}" for (rs, rn), _ in top_contacts[:top_n]]
    freqs_c = [freq * 100 for _, freq in top_contacts[:top_n]]
    colors_c = plt.cm.YlOrRd(np.linspace(0.3, 0.9, top_n))
    ax7.barh(range(top_n), freqs_c, color=colors_c, edgecolor="black", linewidth=0.5)
    ax7.set_yticks(range(top_n)); ax7.set_yticklabels(labels_c)
    ax7.set_xlabel("Contact Frequency (%)"); ax7.set_title("(G) Protein-Ligand Contacts")
    ax7.invert_yaxis()

ax8 = fig3.add_subplot(gs3[0, 1])
dist_colors = ["#E91E63", "#3F51B5", "#009688"]
for i, (label, dists) in enumerate(lig_com_distances.items()):
    c = dist_colors[i % len(dist_colors)]
    ax8.plot(time_ns, dists, color=c, linewidth=0.3, alpha=0.3)
    d_smooth = np.convolve(dists, np.ones(window)/window, mode="valid")
    ax8.plot(time_smooth, d_smooth, color=c, linewidth=1.5, label=label)
ax8.axhline(4.5, color="gray", linestyle=":", alpha=0.5, label="Contact cutoff")
ax8.set_xlabel("Time (ns)"); ax8.set_ylabel("Min Distance (Å)")
ax8.set_title("(H) Key Residue–Ligand Distances")
ax8.legend(fontsize=8); ax8.set_xlim(0, time_ns[-1])

fig3.suptitle(f"Contact Analysis: {LIGAND_NAME} – {RECEPTOR_PDB_ID}", fontsize=14, fontweight="bold")
plt.tight_layout()
plt.savefig(OUTPUT_DIR / "md_contacts_distances.png", dpi=300, bbox_inches="tight", facecolor="white")
plt.close()
print("✅ Figure 3: Contacts + Distances saved")

# %%
# ============================================================
# CELL 4: Summary & Package
# ============================================================
import json as json_mod

best_affinity = -9.96  # from docking

summary = {
    "System": SYSTEM_NAME,
    "Receptor": RECEPTOR_PDB_ID,
    "Ligand": LIGAND_NAME,
    "Simulation_Time_ns": float(PRODUCTION_NS),
    "N_Frames": int(traj.n_frames),
    "N_Runs": len(dcd_files),
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
    "MMGBSA_dG_Mean_kcal_mol": None,  # Will be computed in GPU notebook
    "MMGBSA_dG_SEM_kcal_mol": None,
}

pd.DataFrame([summary]).to_csv(OUTPUT_DIR / "md_summary.csv", index=False)
with open(OUTPUT_DIR / "md_summary.json", "w") as f:
    json_mod.dump(summary, f, indent=2, default=str)

print(f"✅ Summary saved")

# --- ZIP for download ---
import zipfile

zip_path = Path("/kaggle/working") / f"MD_{SYSTEM_NAME}_analysis.zip"
with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
    for f in sorted(OUTPUT_DIR.glob("*")):
        if f.is_dir():
            continue
        zf.write(f, f"{SYSTEM_NAME}/{f.name}")
        print(f"  📦 {f.name}")

print(f"\n✅ Download: {zip_path} ({zip_path.stat().st_size/1e6:.1f} MB)")

# --- Final Report ---
print(f"\n{'=' * 60}")
print(f"📋 FINAL REPORT: {SYSTEM_NAME}")
print(f"{'=' * 60}")
print(f"""
System: {LIGAND_NAME} bound to {RECEPTOR_PDB_ID} (AGTR1)
Simulation: {PRODUCTION_NS:.0f} ns @ 300K ({traj.n_frames:,} frames, {len(dcd_files)} runs)

Docking:
  Affinity: {best_affinity:.2f} kcal/mol

Stability:
  Protein RMSD: {np.mean(rmsd_protein):.2f} ± {np.std(rmsd_protein):.2f} Å
  Ligand RMSD:  {np.mean(rmsd_ligand):.2f} ± {np.std(rmsd_ligand):.2f} Å
  Rg:           {np.mean(rg):.2f} ± {np.std(rg):.2f} Å
  SASA:         {np.mean(sasa_full):.0f} ± {np.std(sasa_full):.0f} Å²
  Equilibrated: {'Yes' if first_dev <= 50 else 'No'}

Interactions:
  H-bonds:      {np.mean(hbond_counts):.1f} ± {np.std(hbond_counts):.1f}
  Top contacts:  {', '.join([f'{rn}{rs}' for (rs, rn), _ in contact_freq_sorted[:5]])}

MM-GBSA: Pending (separate GPU notebook)
""")
print("🎉 CPU Analysis Complete!")
