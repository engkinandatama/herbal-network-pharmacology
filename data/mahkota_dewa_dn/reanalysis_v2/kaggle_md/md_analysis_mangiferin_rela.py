# %% [markdown]
# # MD Analysis — Mangiferin–RELA (100 ns) — CPU Only
# Full trajectory analysis from 2 production runs + Core RMSD analysis.
#
# **CPU-ONLY**: No GPU needed. MM-GBSA runs in separate GPU notebook.
#
# **Input**: Dataset `md-mangiferin-rela-run2-output`
#   - production_run1.dcd (0→58.7 ns from Run 1)
#   - production_run2.dcd (58.7→100 ns from Run 2)
#
# **Analyses**:
#   - 10-metric standard analysis (RMSD, RMSF, Rg, H-bonds, etc.)
#   - Publication-quality figures (3 standard + 1 core RMSD)
#   - Core RMSD (Helix+Sheet), Binding Site RMSD, DSSP

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
    print(f"  [OK] Miniforge3 installed ({time.time()-t_start:.0f}s)")

# --- Install ONLY analysis deps (no openmm/cuda) ---
print("Installing analysis packages...")
packages = f"python={PY_VER} numpy mdtraj parmed"
r = subprocess.run(f"{MAMBA} install -y -c conda-forge {packages}",
                   shell=True, capture_output=True, text=True, timeout=600)
if r.returncode != 0:
    print(f"[WARN] mamba failed: {r.stderr[-300:]}")
else:
    print(f"  [OK] Conda packages installed")

# Remove conda scipy (use pip version)
subprocess.run(f"{MAMBA} remove -y --force scipy 2>/dev/null",
               shell=True, capture_output=True, text=True)
for p in glob.glob(f"{MINIFORGE_DIR}/lib/python*/site-packages/scipy*"):
    shutil.rmtree(p, ignore_errors=True)

# Pip installs
for pkg in ["scipy", "matplotlib", "pandas"]:
    subprocess.run([sys.executable, "-m", "pip", "install", "-q", pkg],
                   capture_output=True, text=True, timeout=120)
    print(f"  [OK] pip: {pkg}")

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
import numpy as np; print(f"  [OK] numpy {np.__version__}")
import mdtraj as md; print(f"  [OK] mdtraj {md.__version__}")
from scipy import stats; print(f"  [OK] scipy")
import matplotlib; print(f"  [OK] matplotlib")
import pandas as pd; print(f"  [OK] pandas")
print(f"\nSetup: {time.time()-t_start:.0f}s [OK]")

# %%
# ============================================================
# CELL 1: Load & Combine Trajectories (100 ns)
# ============================================================
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

SYSTEM_NAME = "Mangiferin_RELA"
LIGAND_NAME = "Mangiferin"
RECEPTOR_PDB_ID = "1NFI"

# --- Find dataset (robust search) ---
print("Searching for dataset...")
INPUT_ROOT = Path("/kaggle/input")
RESULTS_DIR = None

# Strategy 1: Check known dataset slugs
for name in ["md-mangiferin-rela-run2-output", "md-mangiferin-rela-final-output",
             "md-mangiferin-rela-run2-fixed"]:
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

# Strategy 2: Recursive search for DCD files under any input folder
if RESULTS_DIR is None:
    print("  Known slugs not found, searching recursively...")
    for dcd in INPUT_ROOT.rglob("production_run1.dcd"):
        RESULTS_DIR = dcd.parent
        break

# Strategy 3: Search for Mangiferin_RELA folder anywhere
if RESULTS_DIR is None:
    for d in INPUT_ROOT.rglob("Mangiferin_RELA"):
        if d.is_dir():
            RESULTS_DIR = d
            break

if RESULTS_DIR is None:
    # Print debug info
    print("  [FAIL] Dataset not found! Listing /kaggle/input contents:")
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
    print(f"  {dcd_files[0].name} time -> 0-{trajs[0].time[-1]/1000:.1f} ns")
    for i in range(1, len(trajs)):
        t_start_ps = trajs[i-1].time[-1] + dt_ps
        trajs[i].time = np.arange(trajs[i].n_frames) * dt_ps + t_start_ps
        print(f"  {dcd_files[i].name} time -> {trajs[i].time[0]/1000:.1f}-{trajs[i].time[-1]/1000:.1f} ns")

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
    print(f"  Trimmed {n_before} -> {traj.n_frames} frames (capped at {TARGET_NS:.0f} ns)")

total_ns = (traj.time[-1] - traj.time[0]) / 1000
mem_est_gb = traj.n_frames * traj.n_atoms * 3 * 4 / 1e9
print(f"\nCombined trajectory:")
print(f"  Frames: {traj.n_frames:,}")
print(f"  Time: {traj.time[0]/1000:.1f} - {traj.time[-1]/1000:.1f} ns")
print(f"  Duration: {total_ns:.1f} ns")
print(f"  Atoms: {traj.n_atoms:,}")
print(f"  Memory est: ~{mem_est_gb:.1f} GB")

# %%
# ============================================================
# CELL 2: Reimaging + Atom Selections
# ============================================================
print("=" * 60)
print("TRAJECTORY PREPROCESSING")
print("=" * 60)

# --- [0/10] PBC Reimaging ---
print("\n[0/10] Reimaging trajectory...")
try:
    # Find molecules in topology -- protein is usually the largest
    molecules = traj.topology.find_molecules()
    # Sort by size (largest first = protein)
    molecules_sorted = sorted(molecules, key=len, reverse=True)
    # Anchor = protein (largest molecule)
    anchor = [molecules_sorted[0]]
    traj = traj.image_molecules(inplace=False, anchor_molecules=anchor)
    print("  [OK] Reimaging complete")
except Exception as e:
    print(f"  [WARN] Reimaging skipped ({e}) -- analysis still valid")

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
print(f"  CA: {len(ca_atoms)} atoms")
print(f"  Ligand: {len(ligand_atoms)} atoms")

assert len(ligand_atoms) > 0, "No ligand atoms!"
assert len(ca_atoms) > 0, "No CA atoms!"

time_ns = traj.time / 1000.0

# %%
# ============================================================
# CELL 3: Full Trajectory Analysis (10 metrics)
# ============================================================
print("=" * 60)
print("TRAJECTORY ANALYSIS (10 METRICS)")
print("=" * 60)

# --- [1/10] Backbone RMSD ---
print("\n[1/10] Backbone RMSD...")
t1 = time.time()
traj_bb = traj.atom_slice(backbone_atoms)
rmsd_protein = md.rmsd(traj_bb, traj_bb, frame=0) * 10  # nm -> A
print(f"  Mean: {np.mean(rmsd_protein):.2f} +/- {np.std(rmsd_protein):.2f} A ({time.time()-t1:.1f}s)")

# --- [2/10] Ligand RMSD ---
print("[2/10] Ligand RMSD...")
t1 = time.time()
traj_aligned = traj.superpose(traj, frame=0, atom_indices=ca_atoms)
lig_traj = traj_aligned.atom_slice(ligand_atoms)
rmsd_ligand = md.rmsd(lig_traj, lig_traj, frame=0) * 10
print(f"  Mean: {np.mean(rmsd_ligand):.2f} +/- {np.std(rmsd_ligand):.2f} A ({time.time()-t1:.1f}s)")

# --- [3/10] RMSF ---
print("[3/10] CA RMSF...")
t1 = time.time()
ca_traj = traj_aligned.atom_slice(ca_atoms)
rmsf = md.rmsf(ca_traj, ca_traj, frame=0) * 10
residue_ids = [traj.topology.atom(a).residue.resSeq for a in ca_atoms]
residue_names = [traj.topology.atom(a).residue.name for a in ca_atoms]
top_rmsf_idx = np.argmax(rmsf)
print(f"  Max: {np.max(rmsf):.2f} A (residue {residue_names[top_rmsf_idx]}{residue_ids[top_rmsf_idx]}) ({time.time()-t1:.1f}s)")

# --- [4/10] Radius of Gyration ---
print("[4/10] Radius of Gyration...")
t1 = time.time()
rg = md.compute_rg(traj.atom_slice(protein_atoms)) * 10
print(f"  Mean: {np.mean(rg):.2f} +/- {np.std(rg):.2f} A ({time.time()-t1:.1f}s)")

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
print(f"  Mean: {np.mean(hbond_counts):.1f} +/- {np.std(hbond_counts):.1f} ({time.time()-t1:.1f}s)")

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
        dists.append(np.min(d) * 10)  # nm -> A

    dists_full = np.interp(np.arange(traj.n_frames),
                           list(range(0, traj.n_frames, sample_step)), dists)
    lig_com_distances[f"{resName}{resSeq}"] = dists_full
    print(f"  {resName}{resSeq}: {np.mean(dists_full):.2f} +/- {np.std(dists_full):.2f} A")
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

sasa_sampled = np.array(sasa_sampled) * 100  # nm^2 -> A^2
sasa_full = np.interp(np.arange(traj.n_frames), sasa_indices, sasa_sampled)
print(f"  Mean: {np.mean(sasa_full):.1f} +/- {np.std(sasa_full):.1f} A^2 ({time.time()-t1:.1f}s)")

# --- [9/10] Equilibration Check ---
print("\n[9/10] Equilibration check...")
PRODUCTION_NS = time_ns[-1]
n_blocks = 5
block_size = len(rmsd_protein) // n_blocks
block_means = [np.mean(rmsd_protein[i*block_size:(i+1)*block_size]) for i in range(n_blocks)]

print(f"  RMSD block means:")
for i, bm in enumerate(block_means):
    ns_block = PRODUCTION_NS / n_blocks
    print(f"    Block {i+1} ({i*ns_block:.0f}-{(i+1)*ns_block:.0f} ns): {bm:.2f} A")

last3_mean = np.mean(block_means[2:])
first_dev = abs(block_means[0] - last3_mean) / last3_mean * 100
if first_dev > 50:
    equil_ns = PRODUCTION_NS / n_blocks
    print(f"  [WARN] First block deviates {first_dev:.0f}% -- consider skipping first {equil_ns:.0f} ns")
else:
    equil_ns = 0
    print(f"  [OK] Well-equilibrated (deviation: {first_dev:.0f}%)")

from scipy import stats as scipy_stats
last3_time = np.arange(len(block_means[2:]))
slope, _, _, p_val, _ = scipy_stats.linregress(last3_time, block_means[2:])
if p_val < 0.05:
    print(f"  [WARN] Possible drift (slope={slope:.4f}, p={p_val:.3f})")
else:
    print(f"  [OK] No significant drift (p={p_val:.3f})")

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
    print(f"  {label}: frame {fidx} (t={time_ns[fidx]:.1f} ns, RMSD={rmsd_protein[fidx]:.2f} A)")

# --- Save CSVs ---
print("\nSaving analysis data...")

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

print(f"  [OK] All CSVs saved")
print(f"\nAnalysis completed in {(time.time()-t_start)/60:.1f} min")

# %%
# ============================================================
# CELL 4: Core RMSD Analysis (DSSP + Binding Site)
# ============================================================
print("\n" + "=" * 60)
print("CORE RMSD ANALYSIS")
print("=" * 60)

# --- 1. Standard backbone RMSD (protein-only traj) ---
print("\n1. Computing standard backbone RMSD on protein-only traj...")
# prot_traj already exists from SASA calculation
bb_local = prot_traj.topology.select("backbone")
rmsd_backbone_proper = md.rmsd(prot_traj, prot_traj, frame=0, atom_indices=bb_local) * 10
print(f"  Backbone RMSD: {np.mean(rmsd_backbone_proper):.2f} +/- {np.std(rmsd_backbone_proper):.2f} A")

# --- 2. DSSP: Identify secondary structure ---
print("\n2. Computing DSSP secondary structure assignment...")
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
# Top contact residues from Run 2 analysis (resSeq numbers)
# These are the residues with >70% contact frequency from the Run 2 log:
# LEU233, PHE254, LEU193, THR229, PHE288, TYR79, HIS80, ILE292, TYR196, THR191, ASP232
binding_site_resseq = [233, 254, 193, 229, 288, 79, 80, 292, 196, 191, 232]

# If contact frequency data is available, use top residues dynamically
if contact_freq_sorted:
    dynamic_binding = [rs for (rs, rn), freq in contact_freq_sorted if freq > 0.70]
    if len(dynamic_binding) >= 5:
        binding_site_resseq = dynamic_binding
        print(f"  Using dynamic binding site (>70% contact): {binding_site_resseq}")
    else:
        print(f"  Using predefined binding site residues: {binding_site_resseq}")
else:
    print(f"  Using predefined binding site residues: {binding_site_resseq}")

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

# --- 5. DSSP over time ---
print("\n5. Computing DSSP over time (every 10th frame)...")
dssp_stride = 10
dssp_frames = range(0, prot_traj.n_frames, dssp_stride)
dssp_all = md.compute_dssp(prot_traj[list(dssp_frames)], simplified=True)
ss_time_ns = time_ns[list(dssp_frames)]
n_helix_t = np.sum(dssp_all == 'H', axis=1)
n_sheet_t = np.sum(dssp_all == 'E', axis=1)
n_coil_t = np.sum(dssp_all == 'C', axis=1)
print(f"  DSSP computed for {len(list(dssp_frames))} frames")
print(f"  Helix range: {np.min(n_helix_t)}-{np.max(n_helix_t)} residues")
print(f"  Sheet range: {np.min(n_sheet_t)}-{np.max(n_sheet_t)} residues")

# --- 6. Core RMSD Summary ---
print(f"\n{'='*60}")
print(f"CORE RMSD SUMMARY")
print(f"{'='*60}")
print(f"  Full Backbone RMSD:    {np.mean(rmsd_backbone_proper):.2f} +/- {np.std(rmsd_backbone_proper):.2f} A")
print(f"  Core (H+E) RMSD:      {np.mean(rmsd_core):.2f} +/- {np.std(rmsd_core):.2f} A")
print(f"  Binding Site RMSD:     {np.mean(rmsd_binding):.2f} +/- {np.std(rmsd_binding):.2f} A")
print(f"  Reduction (core):      {(1 - np.mean(rmsd_core)/np.mean(rmsd_backbone_proper))*100:.1f}% lower")
print(f"  Reduction (binding):   {(1 - np.mean(rmsd_binding)/np.mean(rmsd_backbone_proper))*100:.1f}% lower")

# --- Save Core RMSD CSVs ---
df_rmsd_core = pd.DataFrame({
    "Time_ns": time_ns,
    "RMSD_Backbone_A": rmsd_backbone_proper,
    "RMSD_Core_HE_A": rmsd_core,
    "RMSD_BindingSite_A": rmsd_binding,
})
df_rmsd_core.to_csv(OUTPUT_DIR / "core_rmsd_timeseries.csv", index=False)

df_dssp = pd.DataFrame({
    "Time_ns": ss_time_ns,
    "N_Helix": n_helix_t,
    "N_Sheet": n_sheet_t,
    "N_Coil": n_coil_t,
})
df_dssp.to_csv(OUTPUT_DIR / "dssp_over_time.csv", index=False)

df_binding = pd.DataFrame({
    "Residue": binding_residues_found,
    "ResSeq": binding_site_resseq[:len(binding_residues_found)],
})
df_binding.to_csv(OUTPUT_DIR / "binding_site_residues.csv", index=False)

# Save core residue info with DSSP
core_res_info = []
for res_idx in range(n_total):
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

print(f"  [OK] Core RMSD CSVs saved")

# --- Free trajectory memory BEFORE plotting cell ---
# (~3.2 GB trajectory + sub-trajectories must be freed here
# to prevent OOM kernel crash at cell boundary)
print("\nFreeing trajectory memory...")
for _var in ['traj', 'prot_traj', 'traj_aligned', 'lig_traj', 'ca_traj', 'traj_bb',
             'dssp_all', 'dssp_matrix']:
    try: exec(f"del {_var}")
    except: pass
gc.collect()
mem_after = gc.get_count()
print(f"  [OK] Trajectory freed (gc objects: {mem_after})")

# %%
# ============================================================
# CELL 5: Plotting skipped (done locally from CSV data)
# ============================================================
# Plotting removed to prevent OOM kernel crash on Kaggle.
# Download the CSV outputs and run plotting locally.
print("Cell 5: Plotting skipped (will be done locally from CSV data)")
print("  Available CSVs: timeseries.csv, rmsf_per_residue.csv, contact_frequency.csv,")
print("  key_distances.csv, core_rmsd_timeseries.csv, dssp_over_time.csv, residue_classification.csv")

# %%
# ============================================================
# CELL 6: Summary & Package
# ============================================================
import json as json_mod

best_affinity = -10.27  # Mangiferin best docking affinity

summary = {
    "System": SYSTEM_NAME,
    "Receptor": RECEPTOR_PDB_ID,
    "Ligand": LIGAND_NAME,
    "Simulation_Time_ns": float(PRODUCTION_NS),
    "N_Frames_Analyzed": len(rmsd_protein),
    "Stride": STRIDE,
    "Docking_Affinity_kcal_mol": float(best_affinity),
    "RMSD_Protein_Mean_A": float(np.mean(rmsd_protein)),
    "RMSD_Protein_Std_A": float(np.std(rmsd_protein)),
    "RMSD_Ligand_Mean_A": float(np.mean(rmsd_ligand)),
    "RMSD_Ligand_Std_A": float(np.std(rmsd_ligand)),
    "RMSD_Core_Mean_A": float(np.mean(rmsd_core)),
    "RMSD_Core_Std_A": float(np.std(rmsd_core)),
    "RMSD_BindingSite_Mean_A": float(np.mean(rmsd_binding)),
    "RMSD_BindingSite_Std_A": float(np.std(rmsd_binding)),
    "Rg_Mean_A": float(np.mean(rg)),
    "Rg_Std_A": float(np.std(rg)),
    "Hbond_Mean": float(np.mean(hbond_counts)),
    "Hbond_Std": float(np.std(hbond_counts)),
    "SASA_Mean_A2": float(np.mean(sasa_full)),
    "SASA_Std_A2": float(np.std(sasa_full)),
    "Top_Contacts": ", ".join([f"{rn}{rs}" for (rs, rn), _ in contact_freq_sorted[:5]]),
    "Equilibration_OK": first_dev <= 50,
    "SS_Helix_Pct": float(100*n_helix/n_total),
    "SS_Sheet_Pct": float(100*n_sheet/n_total),
    "SS_Coil_Pct": float(100*n_coil/n_total),
    "MMGBSA_dG_Mean_kcal_mol": None,  # Will be computed in GPU notebook
    "MMGBSA_dG_SEM_kcal_mol": None,
}

pd.DataFrame([summary]).to_csv(OUTPUT_DIR / "md_summary.csv", index=False)
with open(OUTPUT_DIR / "md_summary.json", "w") as f:
    json_mod.dump(summary, f, indent=2, default=str)

print(f"[OK] Summary saved")

# --- ZIP for download ---
import zipfile

zip_path = Path("/kaggle/working") / f"MD_{SYSTEM_NAME}_analysis.zip"
with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
    for f in sorted(OUTPUT_DIR.glob("*")):
        if f.is_dir():
            continue
        zf.write(f, f"{SYSTEM_NAME}/{f.name}")
        print(f"  {f.name}")

print(f"\n[OK] Download: {zip_path} ({zip_path.stat().st_size/1e6:.1f} MB)")

# --- Final Report ---
print(f"\n{'=' * 60}")
print(f"FINAL REPORT: {SYSTEM_NAME}")
print(f"{'=' * 60}")
print(f"""
System: {LIGAND_NAME} bound to {RECEPTOR_PDB_ID} (RELA / NF-kB p65)
Simulation: {PRODUCTION_NS:.0f} ns @ 310K ({len(rmsd_protein):,} frames analyzed, stride={STRIDE})

Docking:
  Affinity: {best_affinity:.2f} kcal/mol

Stability:
  Protein RMSD:      {np.mean(rmsd_protein):.2f} +/- {np.std(rmsd_protein):.2f} A
  Core (H+E) RMSD:   {np.mean(rmsd_core):.2f} +/- {np.std(rmsd_core):.2f} A
  Binding Site RMSD:  {np.mean(rmsd_binding):.2f} +/- {np.std(rmsd_binding):.2f} A
  Ligand RMSD:        {np.mean(rmsd_ligand):.2f} +/- {np.std(rmsd_ligand):.2f} A
  Rg:                 {np.mean(rg):.2f} +/- {np.std(rg):.2f} A
  SASA:               {np.mean(sasa_full):.0f} +/- {np.std(sasa_full):.0f} A^2
  Equilibrated:       {'Yes' if first_dev <= 50 else 'No'}

Interactions:
  H-bonds:            {np.mean(hbond_counts):.1f} +/- {np.std(hbond_counts):.1f}
  Top contacts:       {', '.join([f'{rn}{rs}' for (rs, rn), _ in contact_freq_sorted[:5]])}

Secondary Structure:
  Helix: {n_helix} res ({100*n_helix/n_total:.1f}%)
  Sheet: {n_sheet} res ({100*n_sheet/n_total:.1f}%)
  Coil:  {n_coil} res ({100*n_coil/n_total:.1f}%)

Output Files:
  - timeseries.csv (RMSD, Rg, H-bonds, SASA over time)
  - rmsf_per_residue.csv
  - contact_frequency.csv
  - key_distances.csv
  - core_rmsd_timeseries.csv (Backbone/Core/BindingSite RMSD)
  - dssp_over_time.csv
  - residue_classification.csv
  - md_summary.csv / .json
  - frame_*.pdb (representative structures)
  - Plots: generate locally from CSV data

MM-GBSA: Pending (separate GPU notebook)

Total time: {(time.time()-t_start)/60:.1f} min
""")
print("CPU Analysis Complete!")
