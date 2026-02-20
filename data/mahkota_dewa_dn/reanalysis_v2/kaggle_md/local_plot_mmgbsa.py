"""
Generate Publication-Quality MM-GBSA Plots from CSV Data.

This script reads 'mmgbsa_per_frame.csv' and 'mmgbsa_summary.csv'
to generate the 2x2 figure panel:
(A) MM-GBSA Distribution Histogram
(B) Binding Free Energy Time Series
(C) Energy Decomposition Bar Chart
(D) Energy Components over Time

Run this script LOCALLY after downloading output from Kaggle.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# --- Configuration ---
CSV_DIR = Path(r"e:\Project\herbal-network-pharmacology\data\mahkota_dewa_dn\reanalysis_v2\kaggle_md\mmgbsa_output")
OUTPUT_PLOT = CSV_DIR / "mmgbsa_plot_local.png"

# --- Load Data ---
def load_data():
    try:
        df_frame = pd.read_csv(CSV_DIR / "mmgbsa_per_frame.csv")
        df_summary = pd.read_csv(CSV_DIR / "mmgbsa_summary.csv")
        print(f"✅ Loaded {len(df_frame)} frames of MM-GBSA data")
        return df_frame, df_summary.iloc[0]
    except FileNotFoundError as e:
        print(f"❌ Error: Could not find CSV files in {CSV_DIR}")
        print(f"   Ensure you downloaded the output from Kaggle correctly.")
        print(f"   Missing: {e.filename}")
        exit(1)

# --- Plotting Function ---
def plot_mmgbsa(df, summary):
    print("📊 Generating 2x2 plot panel...")
    
    # Extract data
    dG_total = df["dG_total"].values
    dE_gas = df["dE_gas"].values
    dG_solv = df["dG_solv"].values
    t_ns = df["time_ns"].values
    
    # Extract stats from summary
    dG_m = summary["dG_total_mean"]
    dG_s = summary["dG_total_sem"]
    dG_sd = summary["dG_total_std"]
    
    # Setup plot style
    plt.rcParams.update({
        "font.family": "sans-serif", "font.size": 10,
        "axes.labelsize": 12, "axes.titlesize": 13,
        "figure.dpi": 300, "axes.linewidth": 1.2,
    })
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # (A) Histogram of total ΔG
    ax = axes[0, 0]
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
    means = [summary["dE_gas_mean"], summary["dG_solv_mean"], summary["dG_total_mean"]]
    sems = [summary["dE_gas_sem"], summary["dG_solv_sem"], summary["dG_total_sem"]]
    colors = ["#2196F3" if m < 0 else "#FF5722" for m in means]
    # For error bars, ensure positive values (matplotlib requires positive errors)
    sems = [abs(s) for s in sems] 
    
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
    
    fig.suptitle(f"MM-GBSA Binding Free Energy: {summary['Ligand']} – {summary['Receptor']}\n(RELA/NF-κB p65) | Solute Dielectric ε=4.0",
                 fontsize=14, fontweight="bold")
    plt.tight_layout()
    plt.savefig(OUTPUT_PLOT, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"✅ Plot saved to: {OUTPUT_PLOT}")

# --- Main ---
if __name__ == "__main__":
    if not CSV_DIR.exists():
        print(f"Creating directory: {CSV_DIR}")
        CSV_DIR.mkdir(parents=True, exist_ok=True)
        print("⚠️ Place the downloaded CSV files here and run again!")
    else:
        df, summary = load_data()
        plot_mmgbsa(df, summary)
