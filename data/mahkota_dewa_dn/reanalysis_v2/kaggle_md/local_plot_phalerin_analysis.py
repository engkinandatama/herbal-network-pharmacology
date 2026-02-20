"""
Generate MD Analysis Plots for Phalerin-AGTR1.

This script reads analysis CSVs from 'analysis_output/analysis/Phalerin_AGTR1'
and generates plots for available data:
1. RMSD (Protein vs Ligand)
2. RMSF
3. Rg
4. SASA
5. H-Bonds
6. Contacts

Note: Core RMSD and specific binding site residues might be missing for this system,
so those plots will be skipped if data is unavailable.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# --- Configuration ---
# Directory containing the CSV files
DATA_DIR = Path(r"e:\Project\herbal-network-pharmacology\data\mahkota_dewa_dn\reanalysis_v2\kaggle_md\analysis_output\analysis\Phalerin_AGTR1")
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
    if "RMSD_Protein_A" in df.columns:
        plt.plot(df["Time_ns"], df["RMSD_Protein_A"], label="Protein Backbone", color="#1f77b4", alpha=0.8)
    if "RMSD_Ligand_A" in df.columns:
        plt.plot(df["Time_ns"], df["RMSD_Ligand_A"], label="Ligand (Fit on Prot)", color="#d62728", alpha=0.8)
    
    # Moving average
    if len(df) > 20 and "RMSD_Protein_A" in df.columns:
        ma_prot = df["RMSD_Protein_A"].rolling(window=10).mean()
        plt.plot(df["Time_ns"], ma_prot, color="#0b4a78", linewidth=2)
    if len(df) > 20 and "RMSD_Ligand_A" in df.columns:
        ma_lig = df["RMSD_Ligand_A"].rolling(window=10).mean()
        plt.plot(df["Time_ns"], ma_lig, color="#8c1b1b", linewidth=2)

    plt.xlabel("Time (ns)")
    plt.ylabel("RMSD (Å)")
    plt.title("RMSD Evolution: Protein Backbone vs Ligand (Phalerin-AGTR1)")
    plt.legend()
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "RMSD_Protein_Ligand.png")
    plt.close()
    print("✅ Generated RMSD_Protein_Ligand.png")

# --- 2. RMSF Plot ---
def plot_rmsf():
    df_rmsf = load_csv("rmsf_per_residue.csv")
    # Binding site might not be defined for Phalerin_AGTR1 in the same way
    # Check if we have a binding site file, otherwise skip highlighting
    df_binding = load_csv("binding_site_residues.csv") 
    
    if df_rmsf is None: return

    plt.figure(figsize=(12, 6))
    
    # Extract residue numbers
    import re
    def get_resnum(s):
        m = re.search(r'\d+', str(s))
        return int(m.group()) if m else 0

    df_rmsf['ResSeq'] = df_rmsf['Residue'].apply(get_resnum)
    
    plt.plot(df_rmsf['ResSeq'], df_rmsf['RMSF_A'], color="#1f77b4", linewidth=1.5, label="RMSF")
    plt.fill_between(df_rmsf['ResSeq'], df_rmsf['RMSF_A'], color="#1f77b4", alpha=0.1)

    # Highlight binding site if available
    if df_binding is not None:
        binding_indices = []
        for res in df_binding['Residue']:
             res_idx = get_resnum(res)
             if res_idx > 0:
                 binding_indices.append(res_idx)
                 plt.axvline(x=res_idx, color='orange', alpha=0.3, linewidth=1)
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

# --- 3. Radius of Gyration (Rg) ---
def plot_rg():
    df = load_csv("timeseries.csv")
    if df is None: return

    if "Rg_A" not in df.columns: return

    plt.figure(figsize=(10, 6))
    plt.plot(df["Time_ns"], df["Rg_A"], color="#9467bd", label="Rg")
    
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

# --- 4. SASA Plot ---
def plot_sasa():
    df = load_csv("timeseries.csv")
    if df is None: return
    
    if "SASA_A2" not in df.columns: return

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

# --- 5. H-Bonds Plot ---
def plot_hbonds():
    df = load_csv("timeseries.csv")
    if df is None: return
    
    if "Hbonds" not in df.columns: return

    plt.figure(figsize=(10, 6))
    plt.plot(df["Time_ns"], df["Hbonds"], color="#e377c2", alpha=0.6, label="Instantaneous")
    
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

# --- 6. Contact Frequency ---
def plot_contacts():
    df = load_csv("contact_frequency.csv")
    if df is None: return
    
    # Check column name
    col_name = "Contact_Frequency"
    if col_name not in df.columns:
        if "Contact_Freq" in df.columns:
            col_name = "Contact_Freq"
        else:
            print(f"⚠️ Warning: Contact frequency column not found in {df.columns}")
            return

    # Take top 15 residues
    top_df = df.sort_values(col_name, ascending=False).head(15)
    
    plt.figure(figsize=(12, 6))
    sns.barplot(x=col_name, y="Residue", data=top_df, hue="Residue", palette="viridis")
    plt.xlabel("Contact Frequency (Fraction of Frames)")
    plt.ylabel("Residue")
    plt.title("Top Protein-Ligand Contacts")
    plt.xlim(0, 1.0)
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "Contact_Frequency.png")
    plt.close()
    print("✅ Generated Contact_Frequency.png")

# --- Main Execution ---
if __name__ == "__main__":
    print(f"📂 Reading data from: {DATA_DIR}")
    print(f"📂 Saving plots to: {OUTPUT_DIR}")
    
    plot_rmsd()
    plot_rmsf()
    plot_rg()
    plot_sasa()
    plot_hbonds()
    plot_contacts()
    
    print("\n🎉 Phalerin analysis plots generated!")
