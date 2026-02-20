
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from pathlib import Path

# --- Configuration ---
BASE_DIR = Path(r"e:\Project\herbal-network-pharmacology\data\mahkota_dewa_dn\reanalysis_v2\kaggle_md")

# Input Directories
MANGIFERIN_DIR = BASE_DIR / "Mangiferin_RELA/plots"
PHALERIN_DIR = BASE_DIR / "analysis_output/analysis/Phalerin_AGTR1/plots"

# Output Directory
OUTPUT_DIR = BASE_DIR / "combined_plots"
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

# Pairs to combine (Mangiferin on top, Phalerin on bottom)
# Filenames must match what is in the directories
pairs = [
    ("RMSD", "Mangiferin_RMSD_Protein_Ligand.png", "Phalerin_RMSD_Protein_Ligand.png"),
    ("RMSF", "Mangiferin_RMSF_Residue_Fluctuation.png", "Phalerin_RMSF_Residue_Fluctuation.png"),
    ("Rg", "Mangiferin_Radius_of_Gyration.png", "Phalerin_Radius_of_Gyration.png"),
    ("H-Bonds", "Mangiferin_Hydrogen_Bonds.png", "Phalerin_Hydrogen_Bonds.png"),
    ("Contacts", "Mangiferin_Mangiferin_Contact_Frequency.png", "Phalerin_Phalerin_Contact_Frequency.png"),
    ("SASA", "Mangiferin_SASA.png", "Phalerin_SASA.png"),
]

def combine_images_vertical(name, img1_filename, img2_filename):
    p1 = MANGIFERIN_DIR / img1_filename
    p2 = PHALERIN_DIR / img2_filename
    
    # Check existence
    if not p1.exists():
        print(f"⚠️ Missing Mangiferin plot: {p1}")
        return
    if not p2.exists():
        print(f"⚠️ Missing Phalerin plot: {p2}")
        return

    try:
        img1 = mpimg.imread(str(p1))
        img2 = mpimg.imread(str(p2))

        # Create figure with 2 subplots vertically
        fig, axes = plt.subplots(2, 1, figsize=(10, 12))
        
        axes[0].imshow(img1)
        axes[0].axis('off')
        axes[0].set_title("Mangiferin-RELA (Top)", fontsize=14, pad=10)
        
        axes[1].imshow(img2)
        axes[1].axis('off')
        axes[1].set_title("Phalerin-AGTR1 (Bottom)", fontsize=14, pad=10)
        
        plt.tight_layout()
        output_path = OUTPUT_DIR / f"Combined_{name}.png"
        plt.savefig(output_path, bbox_inches='tight', dpi=150)
        plt.close()
        print(f"✅ Created {output_path.name}")
    except Exception as e:
        print(f"❌ Error combining {name}: {e}")

if __name__ == "__main__":
    print(f"Processing plots...")
    print(f"  Mangiferin: {MANGIFERIN_DIR}")
    print(f"  Phalerin:   {PHALERIN_DIR}")
    print(f"  Output:     {OUTPUT_DIR}")
    
    for name, f1, f2 in pairs:
        combine_images_vertical(name, f1, f2)
        
    print("\nDone! Combined plots are in the output directory.")
