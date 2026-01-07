"""
Prepare Input Files for MD Simulation

This script extracts and prepares the input files needed for MD simulation
from the docking results.

Output:
- Receptor PDB (clean, without ligand)
- Ligand PDB (from best docking pose)
"""

import os
import shutil
from pathlib import Path

# ============================================
# CONFIGURATION
# ============================================

# Complex 1: 264THM + PPARG
COMPLEX_1 = {
    "name": "264THM_PPARG",
    "compound": "264-trihydroxy-4-methoxybenzophenone",
    "target": "PPARG",
    "pdb_id": "6MS7",
}

# Complex 2: Luteolin + PDE5A
COMPLEX_2 = {
    "name": "Luteolin_PDE5A",
    "compound": "Luteolin",
    "target": "PDE5A",
    "pdb_id": "1TBF",
}

# ============================================
# PATHS
# ============================================

BASE_DIR = Path(__file__).parent.parent
DATA_DIR = BASE_DIR / "data" / "mahkota_dewa_dn"
DOCKING_DIR = DATA_DIR / "results" / "docking"
RECEPTOR_DIR = DOCKING_DIR / "receptors" / "prepared"
LIGAND_DIR = DOCKING_DIR / "ligands" / "pdb"
DOCKED_DIR = DOCKING_DIR / "mahkota_dewa_dn_docking_results"

OUTPUT_DIR = BASE_DIR / "data" / "mahkota_dewa_dn" / "md_input"


def extract_first_pose_from_pdbqt(pdbqt_file: Path, output_pdb: Path):
    """
    Extract the first (best) pose from a docked PDBQT file and convert to PDB.
    """
    with open(pdbqt_file, 'r') as f:
        lines = f.readlines()
    
    pdb_lines = []
    in_model = False
    model_count = 0
    
    for line in lines:
        if line.startswith("MODEL"):
            model_count += 1
            if model_count == 1:
                in_model = True
            continue
        
        if line.startswith("ENDMDL"):
            if model_count == 1:
                break  # Only take first model
            continue
        
        if in_model or model_count == 0:  # Handle files without MODEL markers
            if line.startswith(("ATOM", "HETATM")):
                # Convert PDBQT to PDB (remove charge column)
                pdb_line = line[:66] + "\n" if len(line) > 66 else line
                # Replace HETATM with proper residue name
                if "LIG" not in pdb_line and "UNL" not in pdb_line:
                    pdb_line = pdb_line[:17] + "LIG" + pdb_line[20:]
                pdb_lines.append(pdb_line)
    
    with open(output_pdb, 'w') as f:
        f.writelines(pdb_lines)
        f.write("END\n")
    
    print(f"  ✅ Extracted best pose: {output_pdb.name}")


def prepare_complex(config: dict):
    """Prepare input files for one complex."""
    print(f"\n{'='*60}")
    print(f"Preparing: {config['name']}")
    print(f"{'='*60}")
    
    complex_dir = OUTPUT_DIR / config['name']
    complex_dir.mkdir(parents=True, exist_ok=True)
    
    # 1. Copy receptor PDB
    receptor_pdb = RECEPTOR_DIR / f"{config['pdb_id']}_clean.pdb"
    if receptor_pdb.exists():
        output_receptor = complex_dir / f"{config['target']}_{config['pdb_id']}.pdb"
        shutil.copy(receptor_pdb, output_receptor)
        print(f"  ✅ Receptor: {output_receptor.name}")
    else:
        print(f"  ❌ Receptor not found: {receptor_pdb}")
    
    # 2. Extract ligand from docked PDBQT
    docked_pdbqt = DOCKED_DIR / f"{config['target']}_{config['compound']}_docked.pdbqt"
    if docked_pdbqt.exists():
        output_ligand = complex_dir / f"{config['compound']}_docked.pdb"
        extract_first_pose_from_pdbqt(docked_pdbqt, output_ligand)
    else:
        # Try original ligand PDB
        original_ligand = LIGAND_DIR / f"{config['compound']}.pdb"
        if original_ligand.exists():
            output_ligand = complex_dir / f"{config['compound']}.pdb"
            shutil.copy(original_ligand, output_ligand)
            print(f"  ⚠️  Using original ligand (not docked pose): {output_ligand.name}")
        else:
            print(f"  ❌ Ligand not found: {docked_pdbqt}")
    
    # 3. Create README
    readme = f"""# MD Simulation Input Files - {config['name']}

## Complex
- **Compound**: {config['compound']}
- **Target**: {config['target']}
- **PDB ID**: {config['pdb_id']}

## Files
- `{config['target']}_{config['pdb_id']}.pdb` - Receptor structure (clean)
- `{config['compound']}_docked.pdb` - Ligand (best docking pose)

## Usage
1. Upload this folder as a Kaggle dataset
2. Update paths in the MD simulation notebook
3. Run the notebook

## Notes
- Ligand is extracted from the best docking pose (-9.53 kcal/mol for PPARG)
- Receptor has been cleaned (remove waters, ions, heteroatoms)
"""
    
    with open(complex_dir / "README.md", 'w') as f:
        f.write(readme)
    
    print(f"  ✅ README created")
    print(f"\n📁 Output: {complex_dir}")


def main():
    print("=" * 60)
    print("MD Simulation Input Preparation")
    print("=" * 60)
    
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    # Prepare both complexes
    prepare_complex(COMPLEX_1)
    prepare_complex(COMPLEX_2)
    
    print("\n" + "=" * 60)
    print("✅ All input files prepared!")
    print("=" * 60)
    print(f"\nUpload these folders to Kaggle as datasets:")
    print(f"  1. {OUTPUT_DIR / COMPLEX_1['name']}")
    print(f"  2. {OUTPUT_DIR / COMPLEX_2['name']}")


if __name__ == "__main__":
    main()
