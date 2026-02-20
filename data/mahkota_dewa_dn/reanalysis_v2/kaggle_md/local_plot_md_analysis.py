"""
Generate Publication-Quality MD Analysis Plots from CSV Data.

This script reads analysis CSVs from the 'Mangiferin_RELA' directory
and generates high-quality plots for:
1. RMSD (Protein vs Ligand)
2. Core RMSD (Backbone vs Core vs Binding Site)
3. RMSF (with binding site highlighting)
4. Radius of Gyration (Rg)
5. Solvent Accessible Surface Area (SASA)
6. Hydrogen Bonds (Protein-Ligand)
7. Contact Frequency (Key Reserves)
8. Secondary Structure (DSSP) over time

Run this script LOCALLY.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# --- Configuration ---
# Directory containing the CSV files
DATA_DIR = Path(r"e:\Project\herbal-network-pharmacology\data\mahkota_dewa_dn\reanalysis_v2\kaggle_md\Mangiferin_RELA")
OUTPUT_DIR = DATA_DIR / "plots"

# Create output directory
OUTPUT_DIR.mkdir(exist_ok=True)

# Plot styling
plt.style.use('seaborn-v0_8-paper')
plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 10,
    "axes.labelsize": 12, "axes.titlesize": 14,
    "figure.dpi": 300, "axes.linewidth": 1.2,
    "lines.linewidth": 1.5
})

def load_csv(filename):
    path = DATA_DIR / filename
    if not path.exists():
        print(f"⚠️ Warning: {filename} not found. Skipping related plots.")
        return None
    return pd.read_csv(path)

# --- 1. RMSD Plot (Protein & Ligand) ---
def plot_rmsd():
    df = load_csv("timeseries.csv")
    if df is None: return

    plt.figure(figsize=(10, 6))
    plt.plot(df["Time_ns"], df["RMSD_Protein_A"], label="Protein Backbone", color="#1f77b4", alpha=0.8)
    plt.plot(df["Time_ns"], df["RMSD_Ligand_A"], label="Ligand (Fit on Prot)", color="#d62728", alpha=0.8)
    
    # Moving average
    if len(df) > 20:
        ma_prot = df["RMSD_Protein_A"].rolling(window=10).mean()
        ma_lig = df["RMSD_Ligand_A"].rolling(window=10).mean()
        plt.plot(df["Time_ns"], ma_prot, color="#0b4a78", linewidth=2)
        plt.plot(df["Time_ns"], ma_lig, color="#8c1b1b", linewidth=2)

    plt.xlabel("Time (ns)")
    plt.ylabel("RMSD (Å)")
    plt.title("RMSD Evolution: Protein Backbone vs Ligand")
    plt.legend()
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "RMSD_Protein_Ligand.png")
    plt.close()
    print("✅ Generated RMSD_Protein_Ligand.png")

# --- 2. Core RMSD Plot ---
def plot_core_rmsd():
    df = load_csv("core_rmsd_timeseries.csv")
    if df is None: return

    plt.figure(figsize=(10, 6))
    plt.plot(df["Time_ns"], df["RMSD_Backbone_A"], label="Full Backbone", color="gray", alpha=0.5)
    plt.plot(df["Time_ns"], df["RMSD_Core_HE_A"], label="Core (Helix+Sheet)", color="#2ca02c", linewidth=2)
    plt.plot(df["Time_ns"], df["RMSD_BindingSite_A"], label="Binding Site", color="#ff7f0e", linewidth=2)
    
    plt.xlabel("Time (ns)")
    plt.ylabel("RMSD (Å)")
    plt.title("Stability Analysis: Core & Binding Site")
    plt.legend()
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "RMSD_Stability_Analysis.png")
    plt.close()
    print("✅ Generated RMSD_Stability_Analysis.png")

# --- 3. RMSF Plot ---
def plot_rmsf():
    df_rmsf = load_csv("rmsf_per_residue.csv")
    df_binding = load_csv("binding_site_residues.csv")
    
    if df_rmsf is None: return

    plt.figure(figsize=(12, 6))
    
    # Extract residue numbers
    # Assuming 'Residue' column is in format 'ALA123'
    import re
    def get_resnum(s):
        m = re.search(r'\d+', str(s))
        return int(m.group()) if m else 0

    df_rmsf['ResSeq'] = df_rmsf['Residue'].apply(get_resnum)
    
    plt.plot(df_rmsf['ResSeq'], df_rmsf['RMSF_A'], color="#1f77b4", linewidth=1.5, label="RMSF")
    plt.fill_between(df_rmsf['ResSeq'], df_rmsf['RMSF_A'], color="#1f77b4", alpha=0.1)

    # Highlight binding site
    if df_binding is not None:
        binding_indices = []
        for res in df_binding['Residue']:
             res_idx = get_resnum(res)
             if res_idx > 0:
                 binding_indices.append(res_idx)
                 # Highlight with vertical bars or dots
                 plt.axvline(x=res_idx, color='orange', alpha=0.3, linewidth=1)
        
        # Add legend entry for binding site
        plt.plot([], [], color='orange', alpha=0.5, linewidth=1, label='Binding Site Residues')

    plt.xlabel("Residue Index")
    plt.ylabel("RMSF (Å)")
    plt.title("Root Mean Square Fluctuation (RMSF)")
    plt.legend()
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "RMSF_Residue_Fluctuation.png")
    plt.close()
    print("✅ Generated RMSF_Residue_Fluctuation.png")

# --- 4. Radius of Gyration (Rg) ---
def plot_rg():
    df = load_csv("timeseries.csv")
    if df is None: return

    plt.figure(figsize=(10, 6))
    plt.plot(df["Time_ns"], df["Rg_A"], color="#9467bd", label="Rg")
    
    # Overall mean
    mean_rg = df["Rg_A"].mean()
    plt.axhline(mean_rg, color="black", linestyle="--", label=f"Mean: {mean_rg:.2f} Å")
    
    plt.xlabel("Time (ns)")
    plt.ylabel("Radius of Gyration (Å)")
    plt.title("Protein Compactness (Radius of Gyration)")
    plt.legend()
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "Radius_of_Gyration.png")
    plt.close()
    print("✅ Generated Radius_of_Gyration.png")

# --- 5. SASA Plot ---
def plot_sasa():
    df = load_csv("timeseries.csv")
    if df is None: return
    
    plt.figure(figsize=(10, 6))
    plt.plot(df["Time_ns"], df["SASA_A2"], color="#17becf")
    plt.xlabel("Time (ns)")
    plt.ylabel("SASA ($Å^2$)")
    plt.title("Solvent Accessible Surface Area (SASA)")
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "SASA.png")
    plt.close()
    print("✅ Generated SASA.png")

# --- 6. H-Bonds Plot ---
def plot_hbonds():
    df = load_csv("timeseries.csv")
    if df is None: return
    
    plt.figure(figsize=(10, 6))
    plt.plot(df["Time_ns"], df["Hbonds"], color="#e377c2", alpha=0.6, label="Instantaneous")
    
    # Rolling average
    if len(df) > 20:
        ma = df["Hbonds"].rolling(window=20).mean()
        plt.plot(df["Time_ns"], ma, color="#c51b7d", linewidth=2, label="Rolling Avg (20 frames)")
        
    plt.xlabel("Time (ns)")
    plt.ylabel("Number of Hydrogen Bonds")
    plt.title("Protein-Ligand Hydrogen Bonds")
    plt.legend()
    plt.grid(True, axis='y', linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "Hydrogen_Bonds.png")
    plt.close()
    print("✅ Generated Hydrogen_Bonds.png")

# --- 7. Contact Frequency ---
def plot_contacts():
    df = load_csv("contact_frequency.csv")
    if df is None: return
    
    # Take top 15 residues
    top_df = df.sort_values("Contact_Frequency", ascending=False).head(15)
    
    plt.figure(figsize=(12, 6))
    sns.barplot(x="Contact_Frequency", y="Residue", data=top_df, hue="Residue", palette="viridis")
    plt.xlabel("Contact Frequency (Fraction of Frames)")
    plt.ylabel("Residue")
    plt.title("Top Protein-Ligand Contacts")
    plt.xlim(0, 1.0)
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "Contact_Frequency.png")
    plt.close()
    print("✅ Generated Contact_Frequency.png")

# --- 8. Secondary Structure (DSSP) ---
def plot_dssp():
    df = load_csv("dssp_over_time.csv")
    if df is None: return
    
    plt.figure(figsize=(10, 6))
    plt.stackplot(df["Time_ns"], 
                  df["N_Helix"], df["N_Sheet"], df["N_Coil"],
                  labels=["Helix", "Sheet", "Coil"],
                  colors=["#ff7f0e", "#1f77b4", "#2ca02c"], alpha=0.7)
    
    plt.xlabel("Time (ns)")
    plt.ylabel("Number of Residues")
    plt.title("Secondary Structure Evolution (DSSP)")
    plt.legend(loc="upper left")
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "Secondary_Structure_DSSP.png")
    plt.close()
    print("✅ Generated Secondary_Structure_DSSP.png")

# --- Main Execution ---
if __name__ == "__main__":
    print(f"📂 Reading data from: {DATA_DIR}")
    print(f"📂 Saving plots to: {OUTPUT_DIR}")
    
    plot_rmsd()
    plot_core_rmsd()
    plot_rmsf()
    plot_rg()
    plot_sasa()
    plot_hbonds()
    plot_contacts()
    plot_dssp()
    
    print("\n🎉 All plots generated successfully!")
