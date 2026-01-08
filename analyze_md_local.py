#!/usr/bin/env python3
"""
MD Analysis Script - Local (With PBC Fix)
Analyzes 264THM_PPARG 50ns trajectory using MDAnalysis
Includes: Unwrap PBC, Center protein, Align to reference

Usage: python analyze_md_local.py
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Check if MDAnalysis is installed
try:
    import MDAnalysis as mda
    from MDAnalysis.analysis import rms, align
    from MDAnalysis.analysis.rms import RMSD, RMSF
    from MDAnalysis import transformations
except ImportError:
    print("MDAnalysis not installed. Installing...")
    os.system("pip install MDAnalysis")
    import MDAnalysis as mda
    from MDAnalysis.analysis import rms, align
    from MDAnalysis.analysis.rms import RMSD, RMSF
    from MDAnalysis import transformations

# Configuration
BASE_DIR = Path(__file__).parent / "kaggle_output_verification" / "264THM_PPARG"
TPR_FILE = BASE_DIR / "md" / "md.tpr"
CHECKPOINTS = [10, 20, 30, 40, 50]
OUTPUT_DIR = BASE_DIR / "analysis_local_fixed"

def running_average(data, window=50):
    """Calculate running average with specified window size"""
    if len(data) < window:
        return data
    cumsum = np.cumsum(np.insert(data, 0, 0))
    return (cumsum[window:] - cumsum[:-window]) / window

def main():
    print("=" * 60)
    print("MD TRAJECTORY ANALYSIS - 264THM_PPARG (50 ns)")
    print("With PBC Fix: Unwrap + Center + Align")
    print("=" * 60)
    
    # Create output directory
    OUTPUT_DIR.mkdir(exist_ok=True)
    
    # Find all XTC files
    xtc_files = []
    for ns in CHECKPOINTS:
        xtc = BASE_DIR / "checkpoints" / f"checkpoint_{ns}ns" / "md.xtc"
        if xtc.exists():
            xtc_files.append(str(xtc))
            print(f"✓ Found: checkpoint_{ns}ns/md.xtc")
        else:
            print(f"✗ Missing: checkpoint_{ns}ns/md.xtc")
    
    if not xtc_files:
        print("ERROR: No XTC files found!")
        return
    
    print(f"\nLoading {len(xtc_files)} trajectory files...")
    
    # Load universe with all trajectories
    u = mda.Universe(str(TPR_FILE), xtc_files)
    
    print(f"Total frames: {len(u.trajectory)}")
    print(f"Time per frame: {u.trajectory.dt} ps")
    print(f"Total time: {len(u.trajectory) * u.trajectory.dt / 1000:.1f} ns")
    
    # Select atoms for analysis
    protein = u.select_atoms("protein")
    backbone = u.select_atoms("protein and backbone")
    ca_atoms = u.select_atoms("protein and name CA")
    not_protein = u.select_atoms("not protein")
    
    print(f"\nProtein atoms: {len(protein)}")
    print(f"Backbone atoms: {len(backbone)}")
    print(f"CA atoms (for RMSF): {len(ca_atoms)}")
    
    # =========================================================================
    # Apply PBC Corrections using Transformations
    # =========================================================================
    print("\n--- Applying PBC Corrections ---")
    
    # Create transformation workflow:
    # 1. Unwrap - make molecules whole across PBC
    # 2. Center - center protein in box
    # 3. Wrap - wrap everything back into box
    
    workflow = [
        transformations.unwrap(protein),  # Make protein whole
        transformations.center_in_box(protein, center='geometry'),  # Center protein
        transformations.wrap(not_protein, compound='atoms'),  # Wrap solvent back
    ]
    
    u.trajectory.add_transformations(*workflow)
    print("✓ PBC transformations applied: unwrap → center → wrap")
    
    # =========================================================================
    # Create Reference Structure (first frame, aligned)
    # =========================================================================
    print("\n--- Creating Reference Structure ---")
    
    # Go to first frame
    u.trajectory[0]
    
    # Create reference universe
    ref = mda.Universe(str(TPR_FILE), xtc_files)
    ref.trajectory.add_transformations(*[
        transformations.unwrap(ref.select_atoms("protein")),
        transformations.center_in_box(ref.select_atoms("protein"), center='geometry'),
    ])
    ref.trajectory[0]
    
    print("✓ Reference structure created from frame 0")
    
    # =========================================================================
    # RMSD Analysis with Alignment
    # =========================================================================
    print("\n--- RMSD Analysis (with PBC fix) ---")
    
    # Align trajectory to reference first
    print("  Aligning trajectory to reference...")
    aligner = align.AlignTraj(u, ref, select="backbone", in_memory=True)
    aligner.run()
    print("  ✓ Alignment complete")
    
    # Calculate RMSD
    R = RMSD(u, ref, select="backbone", ref_frame=0)
    R.run()
    
    rmsd_time = R.results.rmsd[:, 1] / 1000  # Convert ps to ns
    rmsd_values = R.results.rmsd[:, 2] / 10  # Convert Å to nm
    
    # Running average
    window = 100
    rmsd_smooth = running_average(rmsd_values, window)
    rmsd_time_smooth = rmsd_time[window//2:len(rmsd_smooth)+window//2]
    
    # Save RMSD data
    np.savetxt(OUTPUT_DIR / "rmsd_backbone.dat", 
               np.column_stack([rmsd_time, rmsd_values]),
               header="Time(ns) RMSD(nm)", comments="")
    
    print(f"  RMSD saved: rmsd_backbone.dat")
    print(f"  Average: {rmsd_values.mean():.3f} ± {rmsd_values.std():.3f} nm")
    print(f"  Min: {rmsd_values.min():.3f} nm, Max: {rmsd_values.max():.3f} nm")
    
    # =========================================================================
    # RMSF Analysis (already aligned)
    # =========================================================================
    print("\n--- RMSF Analysis ---")
    
    rmsf_calc = RMSF(ca_atoms).run()
    rmsf_values = rmsf_calc.results.rmsf / 10  # Convert Å to nm
    residue_ids = ca_atoms.resids
    
    # Save RMSF data
    np.savetxt(OUTPUT_DIR / "rmsf.dat",
               np.column_stack([residue_ids, rmsf_values]),
               header="Residue RMSF(nm)", comments="")
    
    print(f"  RMSF saved: rmsf.dat")
    print(f"  Average: {rmsf_values.mean():.3f} nm")
    
    # =========================================================================
    # Radius of Gyration
    # =========================================================================
    print("\n--- Radius of Gyration ---")
    
    rg_time = []
    rg_values = []
    
    for ts in u.trajectory:
        rg_time.append(ts.time / 1000)  # ps to ns
        rg_values.append(protein.radius_of_gyration() / 10)  # Å to nm
    
    rg_time = np.array(rg_time)
    rg_values = np.array(rg_values)
    
    # Running average
    rg_smooth = running_average(rg_values, window)
    rg_time_smooth = rg_time[window//2:len(rg_smooth)+window//2]
    
    # Save Rg data
    np.savetxt(OUTPUT_DIR / "gyrate.dat",
               np.column_stack([rg_time, rg_values]),
               header="Time(ns) Rg(nm)", comments="")
    
    print(f"  Rg saved: gyrate.dat")
    print(f"  Average: {rg_values.mean():.3f} ± {rg_values.std():.3f} nm")
    
    # =========================================================================
    # Plotting
    # =========================================================================
    print("\n--- Generating Plots ---")
    
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    fig.suptitle("264THM_PPARG - MD Analysis (50 ns, PBC Fixed)", fontsize=16, fontweight='bold')
    
    # Colors
    raw_alpha = 0.3
    smooth_color = '#2E86AB'
    rg_color = '#F18F01'
    
    # RMSD Plot
    axes[0].plot(rmsd_time, rmsd_values, color=smooth_color, linewidth=0.5, alpha=raw_alpha, label='Raw')
    axes[0].plot(rmsd_time_smooth, rmsd_smooth, color=smooth_color, linewidth=2, label=f'Running Avg')
    axes[0].axhline(y=rmsd_values.mean(), color='red', linestyle='--', alpha=0.5)
    axes[0].set_xlabel('Time (ns)', fontsize=12)
    axes[0].set_ylabel('RMSD (nm)', fontsize=12)
    axes[0].set_title(f'Backbone RMSD\n(Avg: {rmsd_values.mean():.3f} ± {rmsd_values.std():.3f} nm)', fontsize=12)
    axes[0].set_xlim(0, 50)
    axes[0].set_ylim(0, None)
    axes[0].legend(loc='lower right')
    axes[0].grid(True, alpha=0.3)
    
    # RMSF Plot
    axes[1].fill_between(residue_ids, 0, rmsf_values, color=smooth_color, alpha=0.3)
    axes[1].plot(residue_ids, rmsf_values, color=smooth_color, linewidth=1)
    axes[1].axhline(y=rmsf_values.mean(), color='red', linestyle='--', alpha=0.5)
    axes[1].set_xlabel('Residue Number', fontsize=12)
    axes[1].set_ylabel('RMSF (nm)', fontsize=12)
    axes[1].set_title(f'RMSF per Residue\n(Avg: {rmsf_values.mean():.3f} nm)', fontsize=12)
    axes[1].grid(True, alpha=0.3)
    
    # Radius of Gyration Plot
    axes[2].plot(rg_time, rg_values, color=rg_color, linewidth=0.5, alpha=raw_alpha, label='Raw')
    axes[2].plot(rg_time_smooth, rg_smooth, color=rg_color, linewidth=2, label=f'Running Avg')
    axes[2].axhline(y=rg_values.mean(), color='red', linestyle='--', alpha=0.5)
    axes[2].set_xlabel('Time (ns)', fontsize=12)
    axes[2].set_ylabel('Rg (nm)', fontsize=12)
    axes[2].set_title(f'Radius of Gyration\n(Avg: {rg_values.mean():.3f} ± {rg_values.std():.3f} nm)', fontsize=12)
    axes[2].set_xlim(0, 50)
    axes[2].legend(loc='upper right')
    axes[2].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "264THM_PPARG_analysis_50ns_fixed.png", dpi=300, bbox_inches='tight')
    print(f"  Plot saved: 264THM_PPARG_analysis_50ns_fixed.png")
    
    plt.show()
    
    # =========================================================================
    # Summary
    # =========================================================================
    print("\n" + "=" * 60)
    print("ANALYSIS SUMMARY (PBC Fixed)")
    print("=" * 60)
    print(f"Total time: {rmsd_time[-1]:.1f} ns")
    print(f"Total frames: {len(rmsd_time)}")
    
    print(f"\n--- Full Trajectory ---")
    print(f"RMSD:  {rmsd_values.mean():.3f} ± {rmsd_values.std():.3f} nm")
    print(f"Rg:    {rg_values.mean():.3f} ± {rg_values.std():.3f} nm")
    print(f"RMSF:  {rmsf_values.mean():.3f} nm (average)")
    
    # Equilibration check (last 20 ns)
    last_20ns = rmsd_time >= 30
    if np.any(last_20ns):
        print(f"\n--- Equilibrated (30-50 ns) ---")
        print(f"RMSD:  {rmsd_values[last_20ns].mean():.3f} ± {rmsd_values[last_20ns].std():.3f} nm")
        print(f"Rg:    {rg_values[last_20ns].mean():.3f} ± {rg_values[last_20ns].std():.3f} nm")
    
    print(f"\nOutput: {OUTPUT_DIR}")
    print("=" * 60)

if __name__ == "__main__":
    main()
