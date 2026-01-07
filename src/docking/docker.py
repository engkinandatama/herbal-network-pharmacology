"""
Molecular Docking Module
========================

Performs molecular docking using AutoDock Vina.
Can be run locally on laptop or in Colab.

Requirements:
- AutoDock Vina installed and in PATH
- RDKit for ligand preparation
- Open Babel (optional) for file conversion
"""

import os
import subprocess
import pandas as pd
from pathlib import Path
from typing import List, Dict, Optional, Tuple
from dataclasses import dataclass
import tempfile
import shutil

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False
    print("Warning: RDKit not installed. Ligand preparation will be limited.")


@dataclass
class DockingResult:
    """Container for docking results."""
    compound_name: str
    target_name: str
    binding_affinity: float  # kcal/mol
    pose_file: Optional[str] = None
    rmsd_lb: float = 0.0
    rmsd_ub: float = 0.0


class MolecularDocker:
    """
    Performs molecular docking using AutoDock Vina.
    
    Usage:
        docker = MolecularDocker(config)
        docker.prepare_ligands()
        docker.prepare_receptors()
        results = docker.run_docking()
        docker.save_results()
    """
    
    def __init__(self, config):
        """
        Initialize the docking module.
        
        Args:
            config: Configuration object with paths and settings
        """
        self.config = config
        self.data_dir = Path(config.data_dir)
        self.results_dir = self.data_dir / "results" / "docking"
        self.results_dir.mkdir(parents=True, exist_ok=True)
        
        # Directories
        self.ligand_dir = self.results_dir / "ligands"
        self.receptor_dir = self.results_dir / "receptors"
        self.output_dir = self.results_dir / "outputs"
        
        for d in [self.ligand_dir, self.receptor_dir, self.output_dir]:
            d.mkdir(parents=True, exist_ok=True)
        
        # Settings
        self.exhaustiveness = config.get("docking.exhaustiveness", 8)
        self.num_modes = config.get("docking.num_modes", 9)
        self.energy_range = config.get("docking.energy_range", 3)
        
        # Results storage
        self.docking_results: List[DockingResult] = []
        self.prepared_ligands: Dict[str, Path] = {}
        self.prepared_receptors: Dict[str, Path] = {}
    
    def check_vina_installed(self) -> bool:
        """Check if AutoDock Vina is installed and accessible."""
        try:
            result = subprocess.run(
                ["vina", "--version"],
                capture_output=True,
                text=True
            )
            return result.returncode == 0
        except FileNotFoundError:
            return False
    
    def smiles_to_pdbqt(self, smiles: str, name: str, output_dir: Optional[Path] = None) -> Optional[Path]:
        """
        Convert SMILES to PDBQT format for docking.
        
        Args:
            smiles: SMILES string
            name: Compound name
            output_dir: Output directory
            
        Returns:
            Path to PDBQT file or None if failed
        """
        if not HAS_RDKIT:
            print("RDKit required for SMILES to PDBQT conversion")
            return None
        
        output_dir = output_dir or self.ligand_dir
        
        try:
            # Create molecule from SMILES
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                print(f"Failed to parse SMILES for {name}")
                return None
            
            # Add hydrogens
            mol = Chem.AddHs(mol)
            
            # Generate 3D coordinates
            result = AllChem.EmbedMolecule(mol, randomSeed=42)
            if result == -1:
                print(f"Failed to generate 3D coordinates for {name}")
                return None
            
            # Optimize geometry
            AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
            
            # Save as PDB first
            pdb_file = output_dir / f"{name}.pdb"
            Chem.MolToPDBFile(mol, str(pdb_file))
            
            # Convert to PDBQT using Open Babel or manual conversion
            pdbqt_file = output_dir / f"{name}.pdbqt"
            
            # Try using Open Babel
            try:
                subprocess.run(
                    ["obabel", str(pdb_file), "-O", str(pdbqt_file), "-xr"],
                    capture_output=True,
                    check=True
                )
            except (FileNotFoundError, subprocess.CalledProcessError):
                # Fallback: create basic PDBQT manually
                self._pdb_to_pdbqt_simple(pdb_file, pdbqt_file)
            
            return pdbqt_file
            
        except Exception as e:
            print(f"Error converting {name}: {e}")
            return None
    
    def _pdb_to_pdbqt_simple(self, pdb_file: Path, pdbqt_file: Path):
        """Simple PDB to PDBQT conversion (basic, no Gasteiger charges)."""
        with open(pdb_file, 'r') as f:
            lines = f.readlines()
        
        with open(pdbqt_file, 'w') as f:
            for line in lines:
                if line.startswith(('ATOM', 'HETATM')):
                    # Add dummy charge and atom type
                    atom_name = line[12:16].strip()
                    atom_type = atom_name[0]  # Simple: first letter
                    # Pad line to 78 chars and add charge + type
                    padded = line.rstrip().ljust(70)
                    f.write(f"{padded}  0.000 {atom_type:>2}\n")
                elif line.startswith('END'):
                    f.write(line)
            f.write("TORSDOF 0\n")
    
    def prepare_ligands(self, compounds_file: Optional[Path] = None) -> Dict[str, Path]:
        """
        Prepare all ligands for docking.
        
        Args:
            compounds_file: Path to compounds CSV with name and smiles columns
            
        Returns:
            Dictionary mapping compound names to PDBQT paths
        """
        if compounds_file is None:
            compounds_file = self.data_dir / "raw" / "compounds.csv"
        
        df = pd.read_csv(compounds_file)
        
        print(f"Preparing {len(df)} ligands...")
        
        for _, row in df.iterrows():
            name = row['name'].replace(' ', '_').replace("'", "")
            smiles = row['smiles']
            
            pdbqt_path = self.smiles_to_pdbqt(smiles, name)
            if pdbqt_path:
                self.prepared_ligands[name] = pdbqt_path
                print(f"  ✓ {name}")
            else:
                print(f"  ✗ {name} (failed)")
        
        print(f"Prepared {len(self.prepared_ligands)} ligands")
        return self.prepared_ligands
    
    def download_receptor(self, pdb_id: str, chain: str = "A") -> Optional[Path]:
        """
        Download receptor structure from PDB.
        
        Args:
            pdb_id: PDB ID (e.g., "1FM6")
            chain: Chain ID to extract
            
        Returns:
            Path to cleaned PDB file
        """
        import urllib.request
        
        pdb_file = self.receptor_dir / f"{pdb_id}.pdb"
        
        if not pdb_file.exists():
            url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
            try:
                urllib.request.urlretrieve(url, pdb_file)
                print(f"Downloaded {pdb_id} from PDB")
            except Exception as e:
                print(f"Failed to download {pdb_id}: {e}")
                return None
        
        return pdb_file
    
    def prepare_receptor(self, pdb_file: Path, name: str) -> Optional[Path]:
        """
        Prepare receptor PDBQT file for docking.
        
        Args:
            pdb_file: Path to PDB file
            name: Receptor name for output
            
        Returns:
            Path to PDBQT file
        """
        pdbqt_file = self.receptor_dir / f"{name}.pdbqt"
        
        # Try using Open Babel or MGLTools
        try:
            # Try prepare_receptor4.py from MGLTools
            subprocess.run(
                ["pythonsh", "prepare_receptor4.py", "-r", str(pdb_file), "-o", str(pdbqt_file)],
                capture_output=True,
                check=True
            )
        except (FileNotFoundError, subprocess.CalledProcessError):
            # Fallback: try Open Babel
            try:
                subprocess.run(
                    ["obabel", str(pdb_file), "-O", str(pdbqt_file), "-xr"],
                    capture_output=True,
                    check=True
                )
            except (FileNotFoundError, subprocess.CalledProcessError):
                # Last resort: simple conversion
                self._pdb_to_pdbqt_simple(pdb_file, pdbqt_file)
        
        if pdbqt_file.exists():
            self.prepared_receptors[name] = pdbqt_file
            return pdbqt_file
        
        return None
    
    def get_binding_site(self, pdb_file: Path, ligand_name: str = None) -> Dict[str, float]:
        """
        Estimate binding site center and dimensions.
        
        If ligand_name is provided, use ligand coordinates.
        Otherwise, use geometric center of protein.
        
        Returns:
            Dictionary with center_x, center_y, center_z, size_x, size_y, size_z
        """
        coords = []
        
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith(('ATOM', 'HETATM')):
                    try:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        coords.append((x, y, z))
                    except ValueError:
                        continue
        
        if not coords:
            return {"center_x": 0, "center_y": 0, "center_z": 0,
                    "size_x": 30, "size_y": 30, "size_z": 30}
        
        # Calculate center
        center_x = sum(c[0] for c in coords) / len(coords)
        center_y = sum(c[1] for c in coords) / len(coords)
        center_z = sum(c[2] for c in coords) / len(coords)
        
        # Default box size (can be refined based on known binding site)
        return {
            "center_x": center_x,
            "center_y": center_y,
            "center_z": center_z,
            "size_x": 25,
            "size_y": 25,
            "size_z": 25
        }
    
    def run_vina(
        self,
        receptor_pdbqt: Path,
        ligand_pdbqt: Path,
        output_pdbqt: Path,
        binding_site: Dict[str, float]
    ) -> Optional[List[DockingResult]]:
        """
        Run AutoDock Vina docking.
        
        Args:
            receptor_pdbqt: Path to receptor PDBQT
            ligand_pdbqt: Path to ligand PDBQT
            output_pdbqt: Path for output poses
            binding_site: Dictionary with center and size
            
        Returns:
            List of DockingResult objects
        """
        log_file = output_pdbqt.with_suffix('.log')
        
        cmd = [
            "vina",
            "--receptor", str(receptor_pdbqt),
            "--ligand", str(ligand_pdbqt),
            "--out", str(output_pdbqt),
            "--log", str(log_file),
            "--center_x", str(binding_site["center_x"]),
            "--center_y", str(binding_site["center_y"]),
            "--center_z", str(binding_site["center_z"]),
            "--size_x", str(binding_site["size_x"]),
            "--size_y", str(binding_site["size_y"]),
            "--size_z", str(binding_site["size_z"]),
            "--exhaustiveness", str(self.exhaustiveness),
            "--num_modes", str(self.num_modes),
            "--energy_range", str(self.energy_range)
        ]
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
            
            if result.returncode != 0:
                print(f"Vina error: {result.stderr}")
                return None
            
            # Parse log file for results
            return self._parse_vina_log(log_file, ligand_pdbqt.stem, receptor_pdbqt.stem)
            
        except subprocess.TimeoutExpired:
            print(f"Docking timeout for {ligand_pdbqt.stem}")
            return None
        except FileNotFoundError:
            print("AutoDock Vina not found. Please install it.")
            return None
    
    def _parse_vina_log(self, log_file: Path, ligand_name: str, receptor_name: str) -> List[DockingResult]:
        """Parse Vina log file for binding affinities."""
        results = []
        
        if not log_file.exists():
            return results
        
        with open(log_file, 'r') as f:
            in_results = False
            for line in f:
                if "-----+------------" in line:
                    in_results = True
                    continue
                
                if in_results and line.strip():
                    parts = line.split()
                    if len(parts) >= 4:
                        try:
                            mode = int(parts[0])
                            affinity = float(parts[1])
                            rmsd_lb = float(parts[2])
                            rmsd_ub = float(parts[3])
                            
                            results.append(DockingResult(
                                compound_name=ligand_name,
                                target_name=receptor_name,
                                binding_affinity=affinity,
                                rmsd_lb=rmsd_lb,
                                rmsd_ub=rmsd_ub
                            ))
                        except (ValueError, IndexError):
                            continue
        
        return results
    
    def run_docking(
        self,
        target_genes: Optional[List[str]] = None,
        top_compounds: Optional[List[str]] = None
    ) -> List[DockingResult]:
        """
        Run docking for selected targets and compounds.
        
        Args:
            target_genes: List of gene symbols to dock against
            top_compounds: List of compound names to dock
            
        Returns:
            List of all docking results
        """
        # Load hub genes if not specified
        if target_genes is None:
            hub_file = self.data_dir / "results" / "hub_genes.csv"
            if hub_file.exists():
                df = pd.read_csv(hub_file)
                target_genes = df['gene'].head(5).tolist()  # Top 5 hub genes
            else:
                print("No hub genes file found")
                return []
        
        # Load drug-like compounds if not specified
        if top_compounds is None:
            admet_file = self.data_dir / "results" / "admet_predictions.csv"
            if admet_file.exists():
                df = pd.read_csv(admet_file)
                df = df[df['drug_like'] == True]
                top_compounds = df['name'].head(10).tolist()  # Top 10 drug-like
            else:
                top_compounds = list(self.prepared_ligands.keys())[:10]
        
        print(f"Running docking: {len(top_compounds)} compounds × {len(target_genes)} targets")
        
        all_results = []
        
        for target in target_genes:
            # Get PDB structure for target (would need mapping)
            pdb_id = self._get_pdb_for_gene(target)
            if not pdb_id:
                print(f"No PDB structure found for {target}")
                continue
            
            # Download and prepare receptor
            pdb_file = self.download_receptor(pdb_id)
            if not pdb_file:
                continue
            
            receptor_pdbqt = self.prepare_receptor(pdb_file, target)
            if not receptor_pdbqt:
                continue
            
            # Get binding site
            binding_site = self.get_binding_site(pdb_file)
            
            # Dock each compound
            for compound in top_compounds:
                if compound not in self.prepared_ligands:
                    continue
                
                ligand_pdbqt = self.prepared_ligands[compound]
                output_file = self.output_dir / f"{compound}_{target}.pdbqt"
                
                print(f"  Docking {compound} → {target}...")
                results = self.run_vina(receptor_pdbqt, ligand_pdbqt, output_file, binding_site)
                
                if results:
                    all_results.extend(results)
                    print(f"    Best affinity: {results[0].binding_affinity:.2f} kcal/mol")
        
        self.docking_results = all_results
        return all_results
    
    def _get_pdb_for_gene(self, gene_symbol: str) -> Optional[str]:
        """
        Get PDB ID for a gene symbol.
        
        This is a simplified mapping. In production, use UniProt or PDB API.
        """
        # Manual mapping for common DN targets
        pdb_mapping = {
            "PPARG": "2PRG",  # PPAR-gamma
            "HMGCR": "1HWK",  # HMG-CoA reductase
            "AGTR1": "4YAY",  # Angiotensin II receptor type 1
            "RELA": "1NFI",   # NF-kB p65
            "KDR": "4AGD",    # VEGFR2
            "PDE5A": "1UDT",  # Phosphodiesterase 5A
            "ADORA2A": "3EML",# Adenosine A2A receptor
            "ADORA2B": "5MZP",# Adenosine A2B receptor
            "SERPINE1": "1DVN",# PAI-1
            "SLC5A2": "7VSI", # SGLT2
            "VDR": "1DB1",    # Vitamin D receptor
            "AXL": "5U6B",    # AXL kinase
        }
        return pdb_mapping.get(gene_symbol)
    
    def save_results(self, output_file: Optional[Path] = None) -> Path:
        """Save docking results to CSV."""
        if output_file is None:
            output_file = self.results_dir / "docking_results.csv"
        
        data = []
        for r in self.docking_results:
            data.append({
                "compound": r.compound_name,
                "target": r.target_name,
                "binding_affinity_kcal_mol": r.binding_affinity,
                "rmsd_lb": r.rmsd_lb,
                "rmsd_ub": r.rmsd_ub
            })
        
        df = pd.DataFrame(data)
        
        # Sort by binding affinity (more negative = better)
        df = df.sort_values("binding_affinity_kcal_mol")
        
        df.to_csv(output_file, index=False)
        print(f"Results saved to: {output_file}")
        
        return output_file
    
    def get_best_interactions(self, top_n: int = 20) -> pd.DataFrame:
        """Get top binding interactions."""
        if not self.docking_results:
            return pd.DataFrame()
        
        df = pd.DataFrame([
            {
                "compound": r.compound_name,
                "target": r.target_name,
                "binding_affinity": r.binding_affinity
            }
            for r in self.docking_results
        ])
        
        # Get best pose per compound-target pair
        df = df.sort_values("binding_affinity").groupby(["compound", "target"]).first().reset_index()
        
        return df.head(top_n)
