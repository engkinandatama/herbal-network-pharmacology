"""
Local plotting script for MD analysis results.
Reads CSV files downloaded from Kaggle and generates publication-quality figures.
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from pathlib import Path

# ==============================================================
# Configuration
# ==============================================================
DATA_DIR = Path(__file__).parent / "analysis_output" / "analysis" / "Phalerin_AGTR1"
OUTPUT_DIR = Path(__file__).parent / "figures"
OUTPUT_DIR.mkdir(exist_ok=True)

LIGAND_NAME = "Phalerin"
RECEPTOR_PDB_ID = "7YMH"
SYSTEM_NAME = f"{LIGAND_NAME}_{RECEPTOR_PDB_ID}"

# ==============================================================
# Load CSV data
# ==============================================================
print("📂 Loading CSV data...")

ts = pd.read_csv(DATA_DIR / "timeseries.csv")
time_ns      = ts["Time_ns"].values
rmsd_protein = ts["RMSD_Protein_A"].values
rmsd_ligand  = ts["RMSD_Ligand_A"].values
rg           = ts["Rg_A"].values
hbond_counts = ts["Hbonds"].values
sasa_full    = ts["SASA_A2"].values
print(f"  Timeseries: {len(ts)} frames, {time_ns[-1]:.1f} ns")

rmsf_df = pd.read_csv(DATA_DIR / "rmsf_per_residue.csv")
residue_ids = rmsf_df["Residue"].values
rmsf        = rmsf_df["RMSF_A"].values
print(f"  RMSF: {len(rmsf_df)} residues")

cf = pd.read_csv(DATA_DIR / "contact_frequency.csv")
cf = cf.dropna(subset=["Contact_Frequency"])
cf = cf.sort_values("Contact_Frequency", ascending=False)
print(f"  Contacts: {len(cf)} residues")

kd = pd.read_csv(DATA_DIR / "key_distances.csv")
dist_cols = [c for c in kd.columns if c != "Time_ns"]
print(f"  Key distances: {len(dist_cols)} residues ({', '.join(dist_cols)})")

# ==============================================================
# Plot settings
# ==============================================================
plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 10,
    "axes.labelsize": 12, "axes.titlesize": 13,
    "xtick.labelsize": 9, "ytick.labelsize": 9,
    "legend.fontsize": 9, "figure.dpi": 300, "axes.linewidth": 1.2,
})

PRODUCTION_NS = time_ns[-1]
window = max(1, len(rmsd_protein) // 100)

# ==============================================================
# Figure 1: 4-Panel (RMSD Protein, RMSD Ligand, RMSF, Rg)
# ==============================================================
print("\n🎨 Generating Figure 1: 4-panel...")
fig = plt.figure(figsize=(14, 10))
gs = GridSpec(2, 2, hspace=0.35, wspace=0.3)

# (A) Protein RMSD
ax1 = fig.add_subplot(gs[0, 0])
ax1.plot(time_ns, rmsd_protein, color="#2196F3", linewidth=0.3, alpha=0.4)
rmsd_smooth = np.convolve(rmsd_protein, np.ones(window)/window, mode="valid")
time_smooth = time_ns[:len(rmsd_smooth)]
ax1.plot(time_smooth, rmsd_smooth, color="#1565C0", linewidth=1.5, label="Running avg")
ax1.set_xlabel("Time (ns)"); ax1.set_ylabel("RMSD (Å)")
ax1.set_title(f"(A) Protein Backbone RMSD\nmean={np.mean(rmsd_protein):.2f}±{np.std(rmsd_protein):.2f} Å")
ax1.legend(); ax1.set_xlim(0, time_ns[-1])

# (B) Ligand RMSD
ax2 = fig.add_subplot(gs[0, 1])
ax2.plot(time_ns, rmsd_ligand, color="#FF9800", linewidth=0.3, alpha=0.4)
rmsd_lig_smooth = np.convolve(rmsd_ligand, np.ones(window)/window, mode="valid")
ax2.plot(time_smooth, rmsd_lig_smooth, color="#E65100", linewidth=1.5, label="Running avg")
ax2.set_xlabel("Time (ns)"); ax2.set_ylabel("RMSD (Å)")
ax2.set_title(f"(B) Ligand RMSD\nmean={np.mean(rmsd_ligand):.2f}±{np.std(rmsd_ligand):.2f} Å")
ax2.legend(); ax2.set_xlim(0, time_ns[-1])

# (C) RMSF per Residue
ax3 = fig.add_subplot(gs[1, 0])
ax3.plot(residue_ids, rmsf, color="#4CAF50", linewidth=1.0)
ax3.fill_between(residue_ids, 0, rmsf, color="#4CAF50", alpha=0.15)
ax3.set_xlabel("Residue Number"); ax3.set_ylabel("RMSF (Å)")
ax3.set_title("(C) Cα RMSF per Residue")
threshold = np.mean(rmsf) + 2 * np.std(rmsf)
ax3.axhline(threshold, color="red", linestyle="--", alpha=0.5, label=f"Mean+2σ ({threshold:.1f} Å)")
ax3.legend()

# (D) Radius of Gyration
ax4 = fig.add_subplot(gs[1, 1])
ax4.plot(time_ns, rg, color="#9C27B0", linewidth=0.3, alpha=0.4)
rg_smooth = np.convolve(rg, np.ones(window)/window, mode="valid")
ax4.plot(time_smooth, rg_smooth, color="#6A1B9A", linewidth=1.5, label="Running avg")
ax4.set_xlabel("Time (ns)"); ax4.set_ylabel("Rg (Å)")
ax4.set_title(f"(D) Radius of Gyration\nmean={np.mean(rg):.2f}±{np.std(rg):.2f} Å")
ax4.legend(); ax4.set_xlim(0, time_ns[-1])

fig.suptitle(f"MD Simulation: {LIGAND_NAME} – {RECEPTOR_PDB_ID} (AGTR1) | {PRODUCTION_NS:.0f} ns",
             fontsize=14, fontweight="bold", y=0.98)
out1 = OUTPUT_DIR / "md_analysis_4panel.png"
plt.savefig(out1, dpi=300, bbox_inches="tight", facecolor="white")
plt.close()
print(f"  ✅ Saved: {out1}")

# ==============================================================
# Figure 2: H-bond + SASA
# ==============================================================
print("🎨 Generating Figure 2: H-bond + SASA...")
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
out2 = OUTPUT_DIR / "md_hbond_sasa.png"
plt.savefig(out2, dpi=300, bbox_inches="tight", facecolor="white")
plt.close()
print(f"  ✅ Saved: {out2}")

# ==============================================================
# Figure 3: Contacts + Key Distances
# ==============================================================
print("🎨 Generating Figure 3: Contacts + Distances...")
fig3 = plt.figure(figsize=(14, 5))
gs3 = GridSpec(1, 2, wspace=0.35)

# (G) Contact Frequency Bar Chart
ax7 = fig3.add_subplot(gs3[0, 0])
top_n = min(15, len(cf))
if top_n > 0:
    labels_c = cf["Residue"].values[:top_n]
    freqs_c  = cf["Contact_Frequency"].values[:top_n] * 100
    colors_c = plt.cm.YlOrRd(np.linspace(0.3, 0.9, top_n))
    ax7.barh(range(top_n), freqs_c, color=colors_c, edgecolor="black", linewidth=0.5)
    ax7.set_yticks(range(top_n)); ax7.set_yticklabels(labels_c)
    ax7.set_xlabel("Contact Frequency (%)"); ax7.set_title("(G) Protein-Ligand Contacts")
    ax7.invert_yaxis()

# (H) Key Residue–Ligand Distances
ax8 = fig3.add_subplot(gs3[0, 1])
dist_colors = ["#E91E63", "#3F51B5", "#009688"]
for i, col in enumerate(dist_cols):
    c = dist_colors[i % len(dist_colors)]
    dists = kd[col].values
    ax8.plot(time_ns, dists, color=c, linewidth=0.3, alpha=0.3)
    d_smooth = np.convolve(dists, np.ones(window)/window, mode="valid")
    ax8.plot(time_smooth, d_smooth, color=c, linewidth=1.5, label=col)
ax8.axhline(4.5, color="gray", linestyle=":", alpha=0.5, label="Contact cutoff")
ax8.set_xlabel("Time (ns)"); ax8.set_ylabel("Min Distance (Å)")
ax8.set_title("(H) Key Residue–Ligand Distances")
ax8.legend(fontsize=8); ax8.set_xlim(0, time_ns[-1])

fig3.suptitle(f"Contact Analysis: {LIGAND_NAME} – {RECEPTOR_PDB_ID}", fontsize=14, fontweight="bold")
plt.tight_layout()
out3 = OUTPUT_DIR / "md_contacts_distances.png"
plt.savefig(out3, dpi=300, bbox_inches="tight", facecolor="white")
plt.close()
print(f"  ✅ Saved: {out3}")

# ==============================================================
# Summary statistics
# ==============================================================
print(f"\n{'='*60}")
print(f"📋 SUMMARY: {SYSTEM_NAME}")
print(f"{'='*60}")
print(f"""
System: {LIGAND_NAME} bound to {RECEPTOR_PDB_ID} (AGTR1)
Simulation: {PRODUCTION_NS:.0f} ns ({len(ts):,} frames)

Stability:
  Protein RMSD: {np.mean(rmsd_protein):.2f} ± {np.std(rmsd_protein):.2f} Å
  Ligand RMSD:  {np.mean(rmsd_ligand):.2f} ± {np.std(rmsd_ligand):.2f} Å
  Rg:           {np.mean(rg):.2f} ± {np.std(rg):.2f} Å
  SASA:         {np.mean(sasa_full):.0f} ± {np.std(sasa_full):.0f} Å²

Interactions:
  H-bonds:      {np.mean(hbond_counts):.1f} ± {np.std(hbond_counts):.1f}
  Top contacts:  {', '.join(cf['Residue'].values[:5])}

Figures saved to: {OUTPUT_DIR.resolve()}
""")
print("🎉 Local plotting complete!")
