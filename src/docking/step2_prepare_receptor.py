"""
Step 2: Prepare Receptor Structures
====================================

Cleans PDB files by removing ligands, water, and selecting chains.
Converts to PDBQT format for docking.

Input: raw PDB files from Step 1
Output: clean PDBQT files in receptors/prepared/

Usage:
    python -m src.docking.step2_prepare_receptor --input-dir receptors/raw
"""

import os
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import click


def parse_pdb_components(pdb_file: Path) -> Dict:
    """
    Parse PDB file and identify components.
    
    Returns:
        Dictionary with chains, ligands, waters, metals info
    """
    chains = set()
    ligands = []
    waters = 0
    metals = []
    protein_atoms = 0
    
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                chains.add(line[21])
                protein_atoms += 1
            elif line.startswith('HETATM'):
                res_name = line[17:20].strip()
                chain = line[21]
                if res_name == 'HOH':
                    waters += 1
                elif res_name in ['ZN', 'MG', 'CA', 'FE', 'MN', 'CU', 'NA', 'K', 'CL']:
                    metals.append(res_name)
                else:
                    ligands.append(res_name)
    
    return {
        'chains': sorted(chains),
        'ligands': list(set(ligands)),
        'waters': waters,
        'metals': list(set(metals)),
        'protein_atoms': protein_atoms
    }


def extract_ligand_coordinates(pdb_file: Path, ligand_name: str) -> Optional[Tuple[float, float, float]]:
    """
    Extract center coordinates of a ligand for binding site definition.
    
    Returns:
        (center_x, center_y, center_z) or None
    """
    coords = []
    
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('HETATM'):
                res_name = line[17:20].strip()
                if res_name == ligand_name:
                    try:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        coords.append((x, y, z))
                    except ValueError:
                        continue
    
    if not coords:
        return None
    
    center_x = sum(c[0] for c in coords) / len(coords)
    center_y = sum(c[1] for c in coords) / len(coords)
    center_z = sum(c[2] for c in coords) / len(coords)
    
    return (center_x, center_y, center_z)


def clean_pdb(
    input_file: Path,
    output_file: Path,
    keep_chain: str = 'A',
    remove_hetatm: bool = True,
    keep_metals: bool = False
) -> Path:
    """
    Clean PDB file by removing unwanted components.
    
    Args:
        input_file: Input PDB path
        output_file: Output PDB path
        keep_chain: Chain ID to keep (default: 'A')
        remove_hetatm: Remove all HETATM records
        keep_metals: Keep metal ions
        
    Returns:
        Path to cleaned PDB
    """
    output_file.parent.mkdir(parents=True, exist_ok=True)
    
    metal_residues = {'ZN', 'MG', 'CA', 'FE', 'MN', 'CU', 'NA', 'K', 'CL'}
    
    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        for line in f_in:
            # Keep ATOM records for specified chain
            if line.startswith('ATOM'):
                chain = line[21]
                if chain == keep_chain or keep_chain == 'ALL':
                    f_out.write(line)
            
            # Handle HETATM
            elif line.startswith('HETATM'):
                if not remove_hetatm:
                    f_out.write(line)
                elif keep_metals:
                    res_name = line[17:20].strip()
                    if res_name in metal_residues:
                        f_out.write(line)
            
            # Keep TER and END
            elif line.startswith(('TER', 'END')):
                f_out.write(line)
    
    return output_file


def pdb_to_pdbqt(pdb_file: Path, pdbqt_file: Path) -> Optional[Path]:
    """
    Convert PDB to PDBQT format.
    
    Uses Open Babel if available, otherwise simple conversion.
    """
    import subprocess
    
    pdbqt_file.parent.mkdir(parents=True, exist_ok=True)
    
    # Try Open Babel first
    try:
        result = subprocess.run(
            ['obabel', str(pdb_file), '-O', str(pdbqt_file), '-xr'],
            capture_output=True,
            text=True
        )
        if result.returncode == 0 and pdbqt_file.exists():
            return pdbqt_file
    except FileNotFoundError:
        pass
    
    # Fallback: simple conversion
    print("  ⚠️  Open Babel not found, using simple conversion...")
    
    with open(pdb_file, 'r') as f_in, open(pdbqt_file, 'w') as f_out:
        for line in f_in:
            if line.startswith(('ATOM', 'HETATM')):
                # Get atom type from atom name
                atom_name = line[12:16].strip()
                atom_type = atom_name[0] if atom_name else 'C'
                
                # Pad line and add charge + type
                padded = line.rstrip().ljust(70)
                f_out.write(f"{padded}  0.000 {atom_type:>2}\n")
            elif line.startswith(('TER', 'END')):
                f_out.write(line)
    
    return pdbqt_file


@click.command()
@click.option('--input-dir', '-i', required=True, help='Directory with raw PDB files')
@click.option('--output-dir', '-o', default=None, help='Output directory for prepared files')
@click.option('--chain', '-c', default='A', help='Chain to keep (default: A, use ALL for all)')
@click.option('--keep-metals', is_flag=True, help='Keep metal ions')
def main(input_dir: str, output_dir: str, chain: str, keep_metals: bool):
    """Step 2: Prepare receptor structures for docking."""
    print("=" * 50)
    print("STEP 2: Prepare Receptor Structures")
    print("=" * 50)
    
    input_path = Path(input_dir)
    if output_dir is None:
        output_path = input_path.parent / "prepared"
    else:
        output_path = Path(output_dir)
    
    pdb_files = list(input_path.glob("*.pdb"))
    print(f"\nFound {len(pdb_files)} PDB files\n")
    
    results = []
    binding_sites = {}
    
    for pdb_file in pdb_files:
        pdb_id = pdb_file.stem
        print(f"[{pdb_id}]")
        
        # Analyze components
        components = parse_pdb_components(pdb_file)
        print(f"  Chains: {components['chains']}")
        print(f"  Ligands: {components['ligands']}")
        print(f"  Waters: {components['waters']}")
        print(f"  Metals: {components['metals']}")
        
        # Extract binding site from first ligand
        if components['ligands']:
            first_ligand = components['ligands'][0]
            center = extract_ligand_coordinates(pdb_file, first_ligand)
            if center:
                binding_sites[pdb_id] = {
                    'ligand': first_ligand,
                    'center_x': round(center[0], 2),
                    'center_y': round(center[1], 2),
                    'center_z': round(center[2], 2)
                }
                print(f"  📍 Binding site from {first_ligand}: ({center[0]:.2f}, {center[1]:.2f}, {center[2]:.2f})")
        
        # Clean PDB
        clean_file = output_path / f"{pdb_id}_clean.pdb"
        clean_pdb(pdb_file, clean_file, keep_chain=chain, keep_metals=keep_metals)
        print(f"  ✅ Cleaned: {clean_file.name}")
        
        # Convert to PDBQT
        pdbqt_file = output_path / f"{pdb_id}.pdbqt"
        pdb_to_pdbqt(clean_file, pdbqt_file)
        print(f"  ✅ PDBQT: {pdbqt_file.name}")
        
        results.append((pdb_id, clean_file, pdbqt_file))
    
    # Save binding sites
    if binding_sites:
        import json
        sites_file = output_path / "binding_sites.json"
        with open(sites_file, 'w') as f:
            json.dump(binding_sites, f, indent=2)
        print(f"\n📍 Binding sites saved to: {sites_file}")
    
    # Save manifest
    manifest_file = output_path / "receptor_manifest.txt"
    with open(manifest_file, 'w') as f:
        f.write("# Step 2: Receptor Preparation Manifest\n\n")
        for pdb_id, clean_file, pdbqt_file in results:
            f.write(f"{pdb_id}\t{pdbqt_file}\n")
    
    print("\n" + "=" * 50)
    print("SUMMARY")
    print("=" * 50)
    print(f"✅ Prepared: {len(results)} receptors")
    print(f"📁 Output: {output_path}")
    print(f"\n➡️  Next: Run step3_prepare_ligand.py")


if __name__ == "__main__":
    main()
