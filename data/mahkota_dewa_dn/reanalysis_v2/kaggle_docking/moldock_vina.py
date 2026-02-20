"""
================================================================================
Molecular Docking: Mahkota Dewa (Phaleria macrocarpa) Network Pharmacology
================================================================================
Kaggle Notebook — AutoDock Vina Batch Docking

Dataset: 15 compounds x 5 targets + 4 controls = 79 docking jobs
Receptors prepared with Meeko (gold standard, polar H included)
Ligands prepared with Open Babel (3D, protonated at pH 7.4)

Targets:
  1HWK (HMGCR)  — HMG-CoA Reductase
  1TBF (RELA)   — NF-kB p65
  3QXY (AGTR1)  — Angiotensin II Receptor Type 1
  4YAY (PPARG)  — PPAR-gamma
  6MS7 (PDE5A)  — Phosphodiesterase 5A

Controls:
  Atorvastatin  → HMGCR (statin)
  Losartan      → AGTR1 (ARB)
  Pioglitazone  → PPARG (TZD)
  Sildenafil    → PDE5A (PDE5i)
================================================================================
"""

# %% [markdown]
# # 1. Setup & Installation

# %%
import subprocess, sys

# Install AutoDock Vina
subprocess.check_call([sys.executable, "-m", "pip", "install", "vina", "-q"])

# Verify
from vina import Vina
print(f"AutoDock Vina installed successfully")

# %%
import json, os, time, csv
from pathlib import Path
from vina import Vina

# ---- Paths ----
# When running on Kaggle, the dataset is mounted at /kaggle/input/
# Adjust DATASET_DIR based on your Kaggle dataset name
DATASET_DIR = Path("/kaggle/input/mahkota-dewa-docking")

# If running locally for testing, use relative path
if not DATASET_DIR.exists():
    DATASET_DIR = Path("dataset")

OUTPUT_DIR = Path("/kaggle/working/results")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

print(f"Dataset: {DATASET_DIR}")
print(f"Output: {OUTPUT_DIR}")

# ---- Load config ----
with open(DATASET_DIR / "docking_config.json") as f:
    config = json.load(f)

with open(DATASET_DIR / "binding_sites.json") as f:
    binding_sites = json.load(f)

print(f"\nReceptors: {len(config['receptors'])}")
print(f"Compounds: {len(config['compounds'])}")
print(f"Controls: {len(config['controls'])}")
print(f"Total jobs: {config['total_jobs']}")

# %% [markdown]
# # 2. Verify Dataset Integrity

# %%
print("=" * 60)
print("DATASET INTEGRITY CHECK")
print("=" * 60)

errors = []

# Check ligands
for compound in config["compounds"]:
    path = DATASET_DIR / "ligands" / f"{compound}.pdbqt"
    if not path.exists():
        errors.append(f"Missing ligand: {compound}")
    elif path.stat().st_size < 100:
        errors.append(f"Empty ligand: {compound}")
    else:
        print(f"  [OK] Ligand: {compound}")

# Check controls
for control, target in config["controls"].items():
    path = DATASET_DIR / "controls" / f"{control}.pdbqt"
    if not path.exists():
        errors.append(f"Missing control: {control}")
    else:
        print(f"  [OK] Control: {control} -> {target}")

# Check receptors
for pdb_id, info in config["receptors"].items():
    path = DATASET_DIR / "receptors" / f"{pdb_id}.pdbqt"
    if not path.exists():
        errors.append(f"Missing receptor: {pdb_id}")
    else:
        # Quick check for HD atoms
        content = path.read_text()
        hd_count = content.count(" HD\n") + content.count(" HD\r")
        print(f"  [OK] Receptor: {pdb_id} ({info['gene']}) — HD={hd_count}")

if errors:
    print(f"\n*** {len(errors)} ERRORS FOUND ***")
    for e in errors:
        print(f"  ! {e}")
    raise RuntimeError("Dataset integrity check failed!")
else:
    print(f"\nAll {19 + 5} files verified OK")

# %% [markdown]
# # 3. Run Molecular Docking

# %%
def run_docking(receptor_pdbqt, ligand_pdbqt, center, box_size, 
                exhaustiveness=32, n_poses=9, energy_range=3):
    """
    Run AutoDock Vina docking and return results.
    
    Returns:
        list of tuples: [(mode, affinity, rmsd_lb, rmsd_ub), ...]
    """
    v = Vina(sf_name='vina')
    v.set_receptor(str(receptor_pdbqt))
    v.set_ligand_from_file(str(ligand_pdbqt))
    v.compute_vina_maps(center=center, box_size=box_size)
    
    v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses, 
           min_rmsd=1.0, max_evals=0)
    
    energies = v.energies()
    return energies

# %%
# === RUN ALL DOCKING JOBS ===

results = []
job_count = 0
total_jobs = config["total_jobs"]
start_time = time.time()

exhaustiveness = config["vina_params"]["exhaustiveness"]
n_poses = config["vina_params"]["num_modes"]
energy_range = config["vina_params"]["energy_range"]

print("=" * 70)
print(f"STARTING BATCH DOCKING — {total_jobs} jobs")
print(f"Exhaustiveness: {exhaustiveness}, Poses: {n_poses}")
print("=" * 70)

# ---- 3a. Compound docking (15 compounds × 5 targets = 75 jobs) ----
for pdb_id, rec_info in config["receptors"].items():
    gene = rec_info["gene"]
    receptor_file = DATASET_DIR / "receptors" / f"{pdb_id}.pdbqt"
    site = binding_sites[pdb_id]
    
    center = [site["center_x"], site["center_y"], site["center_z"]]
    box_size = [site["size_x"], site["size_y"], site["size_z"]]
    
    print(f"\n--- {pdb_id} ({gene}) ---")
    
    for compound in config["compounds"]:
        ligand_file = DATASET_DIR / "ligands" / f"{compound}.pdbqt"
        job_count += 1
        
        try:
            energies = run_docking(
                receptor_file, ligand_file, center, box_size,
                exhaustiveness=exhaustiveness, n_poses=n_poses,
                energy_range=energy_range
            )
            
            best_affinity = energies[0][0]  # First pose, first energy value
            n_found = len(energies)
            
            results.append({
                "receptor": pdb_id,
                "gene": gene,
                "ligand": compound,
                "type": "compound",
                "best_affinity_kcal": round(best_affinity, 2),
                "poses_found": n_found,
                "status": "success"
            })
            
            elapsed = time.time() - start_time
            rate = job_count / elapsed if elapsed > 0 else 0
            eta = (total_jobs - job_count) / rate if rate > 0 else 0
            
            print(f"  [{job_count:>2d}/{total_jobs}] {compound:<25s} → {best_affinity:>6.1f} kcal/mol "
                  f"({n_found} poses) [ETA: {eta:.0f}s]")
            
        except Exception as e:
            results.append({
                "receptor": pdb_id,
                "gene": gene,
                "ligand": compound,
                "type": "compound",
                "best_affinity_kcal": None,
                "poses_found": 0,
                "status": f"error: {str(e)[:100]}"
            })
            print(f"  [{job_count:>2d}/{total_jobs}] {compound:<25s} → ERROR: {e}")

# ---- 3b. Control docking (4 controls × 1 target each = 4 jobs) ----
print(f"\n--- CONTROL DRUGS ---")

for control, target_pdb in config["controls"].items():
    gene = config["receptors"][target_pdb]["gene"]
    receptor_file = DATASET_DIR / "receptors" / f"{target_pdb}.pdbqt"
    ligand_file = DATASET_DIR / "controls" / f"{control}.pdbqt"
    site = binding_sites[target_pdb]
    
    center = [site["center_x"], site["center_y"], site["center_z"]]
    box_size = [site["size_x"], site["size_y"], site["size_z"]]
    
    job_count += 1
    
    try:
        energies = run_docking(
            receptor_file, ligand_file, center, box_size,
            exhaustiveness=exhaustiveness, n_poses=n_poses,
            energy_range=energy_range
        )
        
        best_affinity = energies[0][0]
        n_found = len(energies)
        
        results.append({
            "receptor": target_pdb,
            "gene": gene,
            "ligand": control,
            "type": "control",
            "best_affinity_kcal": round(best_affinity, 2),
            "poses_found": n_found,
            "status": "success"
        })
        
        print(f"  [{job_count:>2d}/{total_jobs}] {control:<25s} → {target_pdb} ({gene}) "
              f"→ {best_affinity:>6.1f} kcal/mol ({n_found} poses)")
        
    except Exception as e:
        results.append({
            "receptor": target_pdb,
            "gene": gene,
            "ligand": control,
            "type": "control",
            "best_affinity_kcal": None,
            "poses_found": 0,
            "status": f"error: {str(e)[:100]}"
        })
        print(f"  [{job_count:>2d}/{total_jobs}] {control:<25s} → ERROR: {e}")

total_time = time.time() - start_time
print(f"\n{'=' * 70}")
print(f"COMPLETED: {job_count}/{total_jobs} jobs in {total_time:.0f}s "
      f"({total_time/60:.1f} min)")
print(f"{'=' * 70}")

# %% [markdown]
# # 4. Save Results

# %%
# Save as CSV
csv_path = OUTPUT_DIR / "docking_results.csv"
with open(csv_path, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=[
        "receptor", "gene", "ligand", "type", 
        "best_affinity_kcal", "poses_found", "status"
    ])
    writer.writeheader()
    writer.writerows(results)
print(f"Saved: {csv_path}")

# Save as JSON
json_path = OUTPUT_DIR / "docking_results.json"
with open(json_path, "w") as f:
    json.dump({
        "metadata": {
            "total_jobs": total_jobs,
            "completed": sum(1 for r in results if r["status"] == "success"),
            "failed": sum(1 for r in results if r["status"] != "success"),
            "exhaustiveness": exhaustiveness,
            "runtime_seconds": round(total_time, 1)
        },
        "results": results
    }, f, indent=2)
print(f"Saved: {json_path}")

# %% [markdown]
# # 5. Quick Analysis

# %%
import pandas as pd

df = pd.DataFrame(results)
df_success = df[df["status"] == "success"].copy()

if len(df_success) == 0:
    print("No successful docking results to analyze!")
else:
    print("=" * 70)
    print("DOCKING RESULTS SUMMARY")
    print("=" * 70)
    
    # Best affinity per compound per target
    pivot = df_success[df_success["type"] == "compound"].pivot_table(
        values="best_affinity_kcal", 
        index="ligand", 
        columns="gene",
        aggfunc="first"
    )
    
    print("\nBinding Affinity (kcal/mol) — Compounds vs Targets:")
    print("-" * 70)
    print(pivot.to_string())
    
    # Control comparison
    print("\n\nControl Drug Affinities:")
    print("-" * 40)
    controls = df_success[df_success["type"] == "control"]
    for _, row in controls.iterrows():
        print(f"  {row['ligand']:<20s} → {row['gene']:<6s}: {row['best_affinity_kcal']:.1f} kcal/mol")
    
    # Top compounds per target
    print("\n\nTop 3 Compounds per Target (most negative = strongest binding):")
    print("-" * 70)
    for gene in config["receptors"].values():
        gene_name = gene["gene"]
        gene_data = df_success[(df_success["gene"] == gene_name) & 
                               (df_success["type"] == "compound")]
        if len(gene_data) > 0:
            top3 = gene_data.nsmallest(3, "best_affinity_kcal")
            control_data = controls[controls["gene"] == gene_name]
            control_val = control_data["best_affinity_kcal"].values[0] if len(control_data) > 0 else None
            
            print(f"\n  {gene_name} (control: {control_val} kcal/mol):")
            for _, row in top3.iterrows():
                better = " *** BETTER THAN CONTROL" if control_val and row["best_affinity_kcal"] <= control_val else ""
                print(f"    {row['ligand']:<25s}: {row['best_affinity_kcal']:.1f} kcal/mol{better}")
    
    # Save pivot table
    pivot.to_csv(OUTPUT_DIR / "affinity_matrix.csv")
    print(f"\nSaved: {OUTPUT_DIR / 'affinity_matrix.csv'}")

# %%
# Final summary for download
print("\n" + "=" * 70)
print("OUTPUT FILES (download from /kaggle/working/results/):")
print("=" * 70)
for f in sorted(OUTPUT_DIR.glob("*")):
    print(f"  {f.name} ({f.stat().st_size:,} bytes)")
