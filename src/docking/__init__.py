"""
Docking Module
==============

Modular molecular docking workflow using AutoDock Vina.

Steps:
    1. step1_download     - Download PDB structures
    2. step2_prepare_receptor - Clean and prepare receptors
    3. step3_prepare_ligand   - Convert SMILES to 3D PDBQT
    4. step4_run_docking      - Run AutoDock Vina
    5. step5_analyze          - Analyze results

Usage:
    python -m src.docking.step1_download --config config/docking_config.yaml
    python -m src.docking.step2_prepare_receptor --input-dir receptors/raw
    python -m src.docking.step3_prepare_ligand --compounds compounds.csv
    python -m src.docking.step4_run_docking --receptors prepared --ligands pdbqt
    python -m src.docking.step5_analyze --results docking_results.json
"""

from .docker import MolecularDocker, DockingResult

__all__ = [
    "MolecularDocker", 
    "DockingResult",
]
