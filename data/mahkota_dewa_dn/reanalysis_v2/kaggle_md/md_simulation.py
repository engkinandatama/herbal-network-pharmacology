# %% [markdown]
# # Molecular Dynamics Simulation — Mahkota Dewa (Phaleria macrocarpa) × Diabetic Nephropathy
# 
# **Systems:**
# 1. **Phalerin → AGTR1** (PDB: 4YAY-like, docking: -9.96 kcal/mol vs Losartan -8.99)
# 2. **Mangiferin → RELA/NF-κB** (PDB: 1TBF, docking: -10.22 kcal/mol)
# 
# **Protocol:** 100 ns MD with OpenMM + AMBER ff14SB/GAFF2 on Kaggle GPU T4
# 
# **Output:** RMSD, RMSF, Rg, H-bond analysis, binding free energy (MM-GBSA)

# %% [markdown]
# ## Cell 0: Install Dependencies

# %%
import subprocess, sys, os

def install(pkg):
    subprocess.check_call([sys.executable, "-m", "pip", "install", "-q", pkg])

# Core MD packages
install("openmm")
install("openmmforcefields")
install("openff-toolkit")
install("pdbfixer")
install("mdtraj")
install("parmed")

# Chemistry tools
install("rdkit")
install("meeko")
install("vina")

# Analysis
install("matplotlib")
install("scipy")
install("pandas")

# AmberTools (via conda in Kaggle)
os.system("conda install -y -c conda-forge ambertools=23 > /dev/null 2>&1")

print("✅ All dependencies installed!")

# %% [markdown]
# ## Cell 1: Configuration

# %%
import json
import warnings
warnings.filterwarnings('ignore')

# ================================================================
# MD SIMULATION CONFIGURATION
# ================================================================

CONFIG = {
    "systems": [
        {
            "name": "Phalerin_AGTR1",
            "pdb_id": "3QXY",
            "gene": "AGTR1",
            "ligand_name": "Phalerin",
            "ligand_smiles": "OC1=CC=C(C=C1)C(=O)C2=CC(=CC(=C2)O)OC3OC(CO)C(O)C(O)C3O",
            "docking_affinity": -9.96,
            "control_drug": "Losartan",
            "control_affinity": -8.99,
            "binding_center": [63.51, 19.70, 10.38],
            "binding_size": [20.0, 20.0, 20.0],
        },
        {
            "name": "Mangiferin_RELA",
            "pdb_id": "1TBF",
            "gene": "RELA",
            "ligand_name": "Mangiferin",
            "ligand_smiles": "OC1=CC2=C(C(=C1)O)C(=O)C3=C(C=C(O)C(=C3O2)O)C4OC(CO)C(O)C(O)C4O",
            "docking_affinity": -10.22,
            "control_drug": None,
            "control_affinity": None,
            "binding_center": [28.73, 30.17, 65.13],
            "binding_size": [20.71, 20.0, 20.0],
        }
    ],
    "md_params": {
        "forcefield_protein": "amber14-all.xml",
        "forcefield_water": "amber14/tip3p.xml",
        "temperature_K": 310.15,   # 37°C (physiological)
        "pressure_atm": 1.0,
        "timestep_fs": 2.0,
        "total_time_ns": 100,
        "equilibration_ns": 1,
        "reporting_interval_ps": 10,  # save every 10 ps
        "padding_nm": 1.0,
        "ionic_strength_M": 0.15,    # physiological NaCl
        "nonbonded_cutoff_nm": 1.0,
        "min_steps": 5000,
        "nvt_steps": 50000,          # 100 ps NVT
        "npt_steps": 500000,         # 1 ns NPT equilibration
    },
    "output_dir": "/kaggle/working/md_results",
    "input_dir": "/kaggle/input/mahkota-dewa-docking",
}

# Derived values
MD = CONFIG["md_params"]
TOTAL_STEPS = int(MD["total_time_ns"] * 1e6 / MD["timestep_fs"])
REPORT_EVERY = int(MD["reporting_interval_ps"] * 1000 / MD["timestep_fs"])

print(f"Total production steps: {TOTAL_STEPS:,}")
print(f"Report every: {REPORT_EVERY} steps ({MD['reporting_interval_ps']} ps)")
print(f"Expected frames: {TOTAL_STEPS // REPORT_EVERY}")

os.makedirs(CONFIG["output_dir"], exist_ok=True)

# %% [markdown]
# ## Cell 2: Prepare Protein Structures (PDBFixer)

# %%
from pdbfixer import PDBFixer
from openmm.app import PDBFile
import requests

def prepare_protein(pdb_id, output_dir):
    """Download PDB, fix missing atoms/residues, remove heterogens"""
    print(f"\n{'='*60}")
    print(f"Preparing protein: {pdb_id}")
    print(f"{'='*60}")
    
    # Download and fix
    fixer = PDBFixer(pdbid=pdb_id)
    
    # Find and add missing residues/atoms
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.findNonstandardResidues()
    
    print(f"  Missing residues: {len(fixer.missingResidues)}")
    print(f"  Missing atoms: {sum(len(v) for v in fixer.missingAtoms.values())}")
    print(f"  Non-standard residues: {len(fixer.nonstandardResidues)}")
    
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=False)
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=7.4)
    
    # Save clean PDB
    out_path = os.path.join(output_dir, f"{pdb_id}_clean.pdb")
    with open(out_path, 'w') as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)
    
    print(f"  ✅ Saved: {out_path}")
    print(f"  Chains: {fixer.topology.getNumChains()}")
    print(f"  Residues: {fixer.topology.getNumResidues()}")
    print(f"  Atoms: {fixer.topology.getNumAtoms()}")
    
    return out_path

# Prepare both receptors
for system in CONFIG["systems"]:
    system["protein_pdb"] = prepare_protein(system["pdb_id"], CONFIG["output_dir"])

# %% [markdown]
# ## Cell 3: Prepare Ligands — SMILES → 3D → Parameterize with GAFF2

# %%
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from openff.toolkit import Molecule as OFFMolecule
import numpy as np

def prepare_ligand_from_smiles(name, smiles, output_dir):
    """Generate 3D conformer from SMILES, optimize, save as SDF"""
    print(f"\nPreparing ligand: {name}")
    
    # RDKit: SMILES → 3D
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    
    # Generate 3D coordinates with ETKDG
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    AllChem.EmbedMolecule(mol, params)
    AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)
    
    # Save as SDF
    sdf_path = os.path.join(output_dir, f"{name}.sdf")
    writer = Chem.SDWriter(sdf_path)
    writer.write(mol)
    writer.close()
    
    print(f"  Atoms: {mol.GetNumAtoms()} (incl. H)")
    print(f"  Formula: {Chem.rdMolDescriptors.CalcMolFormula(mol)}")
    print(f"  ✅ Saved: {sdf_path}")
    
    return sdf_path, mol

for system in CONFIG["systems"]:
    sdf_path, rdmol = prepare_ligand_from_smiles(
        system["ligand_name"], system["ligand_smiles"], CONFIG["output_dir"]
    )
    system["ligand_sdf"] = sdf_path
    system["rdkit_mol"] = rdmol

# %% [markdown]
# ## Cell 4: Re-dock Ligands into Binding Site (AutoDock Vina)

# %%
from vina import Vina

def dock_ligand(protein_pdbqt, ligand_sdf, center, size, name, output_dir):
    """Re-dock ligand to get best pose in the binding site"""
    from meeko import MoleculePreparation, PDBQTWriterLegacy, RDKitMolCreate
    
    print(f"\nDocking {name}...")
    
    # Convert SDF to PDBQT via Meeko
    mol = Chem.SDMolSupplier(ligand_sdf, removeHs=False)[0]
    preparator = MoleculePreparation()
    mol_setups = preparator.prepare(mol)
    
    lig_pdbqt_path = os.path.join(output_dir, f"{name}_lig.pdbqt")
    for setup in mol_setups:
        pdbqt_string, is_ok, error_msg = PDBQTWriterLegacy.write_string(setup)
        with open(lig_pdbqt_path, 'w') as f:
            f.write(pdbqt_string)
    
    # Alternatively, use the pre-prepared PDBQT from dataset
    dataset_lig = os.path.join(CONFIG["input_dir"], "ligands", f"{name}.pdbqt")
    if os.path.exists(dataset_lig):
        lig_pdbqt_path = dataset_lig
        print(f"  Using pre-prepared PDBQT: {dataset_lig}")
    
    # Use pre-prepared receptor PDBQT from dataset
    pdb_id = [s for s in CONFIG["systems"] if s["ligand_name"] == name][0]["pdb_id"]
    dataset_rec = os.path.join(CONFIG["input_dir"], "receptors", f"{pdb_id}.pdbqt")
    if os.path.exists(dataset_rec):
        protein_pdbqt = dataset_rec
        print(f"  Using pre-prepared receptor PDBQT: {dataset_rec}")
    
    # Dock
    v = Vina(sf_name='vina')
    v.set_receptor(protein_pdbqt)
    v.set_ligand_from_file(lig_pdbqt_path)
    v.compute_vina_maps(center=center, box_size=size)
    v.dock(exhaustiveness=32, n_poses=5)
    
    # Save best pose
    pose_pdbqt = os.path.join(output_dir, f"{name}_docked.pdbqt")
    v.write_poses(pose_pdbqt, n_poses=1, overwrite=True)
    
    energies = v.energies()
    best_affinity = energies[0][0]
    print(f"  Best affinity: {best_affinity:.2f} kcal/mol")
    print(f"  ✅ Saved best pose: {pose_pdbqt}")
    
    # Convert docked PDBQT back to SDF using Meeko
    docked_sdf = os.path.join(output_dir, f"{name}_docked.sdf")
    pdbqt_mols = RDKitMolCreate.from_pdbqt_file(pose_pdbqt)
    if pdbqt_mols:
        writer = Chem.SDWriter(docked_sdf)
        writer.write(pdbqt_mols[0])
        writer.close()
        print(f"  ✅ Docked SDF: {docked_sdf}")
    
    return pose_pdbqt, docked_sdf, best_affinity

for system in CONFIG["systems"]:
    pdbqt, sdf, affinity = dock_ligand(
        system["protein_pdb"],
        system["ligand_sdf"],
        system["binding_center"],
        system["binding_size"],
        system["ligand_name"],
        CONFIG["output_dir"]
    )
    system["docked_pdbqt"] = pdbqt
    system["docked_sdf"] = sdf
    system["redock_affinity"] = affinity

# %% [markdown]
# ## Cell 5: Build Protein-Ligand Complex with OpenFF + OpenMM

# %%
import openmm
from openmm import app, unit, LangevinMiddleIntegrator, MonteCarloBarostat
from openmm.app import Modeller, ForceField, Simulation, PDBFile, PME, HBonds
from openmmforcefields.generators import GAFFTemplateGenerator
from openff.toolkit import Molecule

def build_complex(system_config, output_dir):
    """Build solvated protein-ligand complex with AMBER ff14SB + GAFF2"""
    name = system_config["name"]
    print(f"\n{'='*60}")
    print(f"Building complex: {name}")
    print(f"{'='*60}")
    
    # Load protein
    protein_pdb = PDBFile(system_config["protein_pdb"])
    
    # Load ligand from docked SDF using OpenFF
    off_mol = Molecule.from_file(system_config["docked_sdf"])
    
    # Create GAFF2 template generator
    gaff = GAFFTemplateGenerator(
        molecules=[off_mol],
        forcefield='gaff-2.11'
    )
    
    # Create forcefield with GAFF2
    ff = ForceField(
        CONFIG["md_params"]["forcefield_protein"],
        CONFIG["md_params"]["forcefield_water"]
    )
    ff.registerTemplateGenerator(gaff.generator)
    
    # Combine protein + ligand
    modeller = Modeller(protein_pdb.topology, protein_pdb.positions)

    # Add ligand to topology
    lig_top = off_mol.to_topology().to_openmm()
    lig_pos = off_mol.conformers[0].to_openmm()
    modeller.add(lig_top, lig_pos)
    
    print(f"  Complex atoms (protein+ligand): {modeller.topology.getNumAtoms()}")
    
    # Solvate with TIP3P + NaCl
    modeller.addSolvent(
        ff,
        model='tip3p',
        padding=CONFIG["md_params"]["padding_nm"] * unit.nanometers,
        ionicStrength=CONFIG["md_params"]["ionic_strength_M"] * unit.molar
    )
    
    print(f"  Solvated atoms: {modeller.topology.getNumAtoms()}")
    
    # Save solvated complex
    complex_pdb = os.path.join(output_dir, f"{name}_solvated.pdb")
    with open(complex_pdb, 'w') as f:
        PDBFile.writeFile(modeller.topology, modeller.positions, f)
    print(f"  ✅ Saved: {complex_pdb}")
    
    # Create OpenMM system
    system = ff.createSystem(
        modeller.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=CONFIG["md_params"]["nonbonded_cutoff_nm"] * unit.nanometers,
        constraints=HBonds,
    )
    
    return modeller, system, ff, complex_pdb

for sys_conf in CONFIG["systems"]:
    modeller, omm_system, ff, complex_pdb = build_complex(sys_conf, CONFIG["output_dir"])
    sys_conf["modeller"] = modeller
    sys_conf["omm_system"] = omm_system
    sys_conf["ff"] = ff
    sys_conf["complex_pdb"] = complex_pdb

# %% [markdown]
# ## Cell 6: Energy Minimization + Equilibration (NVT + NPT)

# %%
def run_equilibration(system_config, output_dir):
    """Energy minimization → NVT heating → NPT equilibration"""
    name = system_config["name"]
    modeller = system_config["modeller"]
    omm_system = system_config["omm_system"]
    MD = CONFIG["md_params"]
    
    print(f"\n{'='*60}")
    print(f"Equilibrating: {name}")
    print(f"{'='*60}")
    
    # ---- Step 1: Energy Minimization ----
    print("\n[1/3] Energy Minimization...")
    integrator = LangevinMiddleIntegrator(
        MD["temperature_K"] * unit.kelvin,
        1.0 / unit.picosecond,
        MD["timestep_fs"] * unit.femtoseconds
    )
    
    # Select GPU platform
    try:
        platform = openmm.Platform.getPlatformByName('CUDA')
        properties = {'DeviceIndex': '0', 'CudaPrecision': 'mixed'}
        print("  Using CUDA (GPU T4)")
    except Exception:
        platform = openmm.Platform.getPlatformByName('CPU')
        properties = {}
        print("  ⚠️ CUDA not available, using CPU")
    
    simulation = Simulation(
        modeller.topology, omm_system, integrator, platform, properties
    )
    simulation.context.setPositions(modeller.positions)
    
    # Minimize
    e_before = simulation.context.getState(getEnergy=True).getPotentialEnergy()
    print(f"  Energy before minimization: {e_before}")
    
    simulation.minimizeEnergy(maxIterations=MD["min_steps"])
    
    e_after = simulation.context.getState(getEnergy=True).getPotentialEnergy()
    print(f"  Energy after minimization: {e_after}")
    
    # Save minimized structure
    min_pdb = os.path.join(output_dir, f"{name}_minimized.pdb")
    state = simulation.context.getState(getPositions=True)
    with open(min_pdb, 'w') as f:
        PDBFile.writeFile(modeller.topology, state.getPositions(), f)
    print(f"  ✅ Minimized PDB: {min_pdb}")
    
    # ---- Step 2: NVT Heating (100 ps) ----
    print(f"\n[2/3] NVT Heating ({MD['nvt_steps'] * MD['timestep_fs'] / 1000:.0f} ps at {MD['temperature_K']} K)...")
    simulation.context.setVelocitiesToTemperature(MD["temperature_K"] * unit.kelvin)
    
    # Report temperature during heating
    nvt_report_every = MD["nvt_steps"] // 5
    for i in range(5):
        simulation.step(nvt_report_every)
        state = simulation.context.getState(getEnergy=True)
        temp = (2 * state.getKineticEnergy() / (
            3 * modeller.topology.getNumAtoms() * unit.MOLAR_GAS_CONSTANT_R
        ))
        print(f"  Step {(i+1) * nvt_report_every}: T = {temp:.1f}")
    
    print(f"  ✅ NVT complete")
    
    # ---- Step 3: NPT Equilibration (1 ns) ----
    print(f"\n[3/3] NPT Equilibration ({MD['npt_steps'] * MD['timestep_fs'] / 1e6:.1f} ns at {MD['pressure_atm']} atm)...")
    
    # Add barostat for NPT
    omm_system.addForce(
        MonteCarloBarostat(
            MD["pressure_atm"] * unit.atmospheres,
            MD["temperature_K"] * unit.kelvin,
            25  # attempt every 25 steps
        )
    )
    
    # Reinitialize context with barostat
    simulation.context.reinitialize(preserveState=True)
    
    # Run NPT equilibration
    npt_report_every = MD["npt_steps"] // 10
    for i in range(10):
        simulation.step(npt_report_every)
        state = simulation.context.getState(getEnergy=True)
        pe = state.getPotentialEnergy()
        box = state.getPeriodicBoxVectors()
        vol = box[0][0] * box[1][1] * box[2][2]
        print(f"  Step {(i+1) * npt_report_every}: PE = {pe}, Vol = {vol:.2f}")
    
    # Save equilibrated state
    eq_pdb = os.path.join(output_dir, f"{name}_equilibrated.pdb")
    state = simulation.context.getState(getPositions=True, getVelocities=True)
    with open(eq_pdb, 'w') as f:
        PDBFile.writeFile(modeller.topology, state.getPositions(), f)
    
    # Save checkpoint
    checkpoint = os.path.join(output_dir, f"{name}_equilibrated.chk")
    simulation.saveCheckpoint(checkpoint)
    
    print(f"  ✅ Equilibrated PDB: {eq_pdb}")
    print(f"  ✅ Checkpoint: {checkpoint}")
    
    return simulation, state

for sys_conf in CONFIG["systems"]:
    sim, eq_state = run_equilibration(sys_conf, CONFIG["output_dir"])
    sys_conf["simulation"] = sim
    sys_conf["equilibrated_state"] = eq_state

# %% [markdown]
# ## Cell 7: Production MD — 100 ns

# %%
from openmm.app import DCDReporter, StateDataReporter

def run_production(system_config, output_dir):
    """Run 100 ns production MD"""
    name = system_config["name"]
    simulation = system_config["simulation"]
    MD = CONFIG["md_params"]
    
    print(f"\n{'='*60}")
    print(f"Production MD: {name} ({MD['total_time_ns']} ns)")
    print(f"{'='*60}")
    
    # Setup reporters
    dcd_file = os.path.join(output_dir, f"{name}_production.dcd")
    log_file = os.path.join(output_dir, f"{name}_production.log")
    
    simulation.reporters.clear()
    simulation.reporters.append(
        DCDReporter(dcd_file, REPORT_EVERY)
    )
    simulation.reporters.append(
        StateDataReporter(
            log_file, REPORT_EVERY,
            step=True, time=True,
            potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
            temperature=True, volume=True, density=True,
            speed=True, remainingTime=True,
            totalSteps=TOTAL_STEPS
        )
    )
    simulation.reporters.append(
        StateDataReporter(
            sys.stdout, REPORT_EVERY * 100,  # print to console less frequently
            step=True, time=True, speed=True,
            remainingTime=True, totalSteps=TOTAL_STEPS
        )
    )
    
    # Save checkpoint periodically
    checkpoint_interval = TOTAL_STEPS // 10  # every 10 ns
    
    print(f"  DCD trajectory: {dcd_file}")
    print(f"  Log: {log_file}")
    print(f"  Total steps: {TOTAL_STEPS:,}")
    print(f"  Checkpoint every: {checkpoint_interval:,} steps")
    print()
    
    import time
    start_time = time.time()
    
    # Run in chunks for checkpointing
    steps_done = 0
    while steps_done < TOTAL_STEPS:
        chunk = min(checkpoint_interval, TOTAL_STEPS - steps_done)
        simulation.step(chunk)
        steps_done += chunk
        
        # Save checkpoint
        ns_done = steps_done * MD["timestep_fs"] / 1e6
        chk_path = os.path.join(output_dir, f"{name}_checkpoint_{ns_done:.0f}ns.chk")
        simulation.saveCheckpoint(chk_path)
        
        elapsed = time.time() - start_time
        speed_ns_day = ns_done / (elapsed / 86400) if elapsed > 0 else 0
        print(f"  ✅ {ns_done:.0f} ns done | Elapsed: {elapsed/3600:.1f}h | Speed: {speed_ns_day:.1f} ns/day")
    
    # Save final state
    final_pdb = os.path.join(output_dir, f"{name}_final.pdb")
    state = simulation.context.getState(getPositions=True)
    with open(final_pdb, 'w') as f:
        PDBFile.writeFile(simulation.topology, state.getPositions(), f)
    
    total_elapsed = time.time() - start_time
    print(f"\n  🏁 Production complete!")
    print(f"  Total time: {total_elapsed/3600:.2f} hours")
    print(f"  Average speed: {MD['total_time_ns'] / (total_elapsed / 86400):.1f} ns/day")
    
    system_config["dcd_file"] = dcd_file
    system_config["log_file"] = log_file
    system_config["final_pdb"] = final_pdb
    
    return dcd_file, log_file

for sys_conf in CONFIG["systems"]:
    run_production(sys_conf, CONFIG["output_dir"])

# %% [markdown]
# ## Cell 8: Analysis — RMSD, RMSF, Radius of Gyration, H-bonds

# %%
import mdtraj as md
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats

def analyze_trajectory(system_config, output_dir):
    """Comprehensive MD trajectory analysis"""
    name = system_config["name"]
    
    print(f"\n{'='*60}")
    print(f"Analyzing: {name}")
    print(f"{'='*60}")
    
    # Load trajectory
    traj = md.load(
        system_config["dcd_file"],
        top=system_config["complex_pdb"]
    )
    print(f"  Frames: {traj.n_frames}")
    print(f"  Time range: {traj.time[0]:.1f} - {traj.time[-1]:.1f} ps")
    
    # Time axis in ns
    time_ns = traj.time / 1000.0
    
    # Select atoms
    protein_atoms = traj.topology.select("protein and name CA")
    ligand_atoms = traj.topology.select("not protein and not water and not (name NA or name CL)")
    
    print(f"  Protein CA atoms: {len(protein_atoms)}")
    print(f"  Ligand atoms: {len(ligand_atoms)}")
    
    results = {"name": name, "time_ns": time_ns}
    
    # ---- 1. RMSD ----
    print("\n  [1/4] Computing RMSD...")
    
    # Protein backbone RMSD
    protein_rmsd = md.rmsd(traj, traj, 0, atom_indices=protein_atoms) * 10  # nm → Å
    results["protein_rmsd"] = protein_rmsd
    
    # Ligand RMSD (after aligning on protein)
    if len(ligand_atoms) > 0:
        # Superpose on protein, then compute ligand RMSD
        traj_aligned = traj.superpose(traj, 0, atom_indices=protein_atoms)
        ligand_rmsd = md.rmsd(traj_aligned, traj_aligned, 0, atom_indices=ligand_atoms) * 10
        results["ligand_rmsd"] = ligand_rmsd
    
    # ---- 2. RMSF ----
    print("  [2/4] Computing RMSF...")
    # Use last 80 ns (equilibrated portion)
    eq_start = int(0.2 * traj.n_frames)  # skip first 20%
    traj_eq = traj[eq_start:]
    
    rmsf = md.rmsf(traj_eq, traj_eq, 0, atom_indices=protein_atoms) * 10  # Å
    results["rmsf"] = rmsf
    results["rmsf_residues"] = [traj.topology.atom(i).residue.resSeq for i in protein_atoms]
    
    # ---- 3. Radius of Gyration ----
    print("  [3/4] Computing Rg...")
    rg = md.compute_rg(traj) * 10  # nm → Å
    results["rg"] = rg
    
    # ---- 4. Hydrogen Bonds ----
    print("  [4/4] Computing H-bonds...")
    # H-bonds between protein and ligand
    if len(ligand_atoms) > 0:
        hbonds = md.baker_hubbard(traj_eq, freq=0.1)
        
        # Filter protein-ligand H-bonds only
        pl_hbonds = []
        for hb in hbonds:
            donor_is_protein = traj.topology.atom(hb[0]).residue.is_protein
            acceptor_is_protein = traj.topology.atom(hb[2]).residue.is_protein
            if donor_is_protein != acceptor_is_protein:  # one protein, one ligand
                pl_hbonds.append(hb)
        
        results["n_hbonds_total"] = len(hbonds)
        results["n_hbonds_prot_lig"] = len(pl_hbonds)
        results["hbond_details"] = pl_hbonds
        
        # H-bond count per frame
        hbond_counts = []
        for frame in range(traj_eq.n_frames):
            frame_hbonds = md.baker_hubbard(traj_eq[frame:frame+1], freq=0.0)
            pl_count = 0
            for hb in frame_hbonds:
                donor_is_protein = traj_eq.topology.atom(hb[0]).residue.is_protein
                acceptor_is_protein = traj_eq.topology.atom(hb[2]).residue.is_protein
                if donor_is_protein != acceptor_is_protein:
                    pl_count += 1
            hbond_counts.append(pl_count)
        
        results["hbond_counts"] = hbond_counts
        print(f"  Protein-ligand H-bonds (>10% freq): {len(pl_hbonds)}")
    
    # ---- Summary Statistics ----
    print(f"\n  === Summary (last 80 ns) ===")
    rmsd_eq = protein_rmsd[eq_start:]
    print(f"  Protein RMSD: {np.mean(rmsd_eq):.2f} ± {np.std(rmsd_eq):.2f} Å")
    if "ligand_rmsd" in results:
        lig_rmsd_eq = ligand_rmsd[eq_start:]
        print(f"  Ligand RMSD:  {np.mean(lig_rmsd_eq):.2f} ± {np.std(lig_rmsd_eq):.2f} Å")
    rg_eq = rg[eq_start:]
    print(f"  Rg: {np.mean(rg_eq):.2f} ± {np.std(rg_eq):.2f} Å")
    if "hbond_counts" in results:
        print(f"  Avg H-bonds/frame: {np.mean(results['hbond_counts']):.1f} ± {np.std(results['hbond_counts']):.1f}")
    
    return results

all_results = {}
for sys_conf in CONFIG["systems"]:
    all_results[sys_conf["name"]] = analyze_trajectory(sys_conf, CONFIG["output_dir"])

# %% [markdown]
# ## Cell 9: Publication-Quality Figures

# %%
def create_figures(all_results, output_dir):
    """Generate publication-quality analysis figures"""
    
    # Consistent styling
    plt.rcParams.update({
        'font.size': 12,
        'font.family': 'sans-serif',
        'axes.linewidth': 1.2,
        'xtick.major.width': 1.2,
        'ytick.major.width': 1.2,
        'figure.dpi': 300,
    })
    
    colors = ['#2196F3', '#FF5722']  # Blue for Phalerin, Red for Mangiferin
    labels = list(all_results.keys())
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # ---- Panel A: Protein RMSD ----
    ax = axes[0, 0]
    for i, (name, res) in enumerate(all_results.items()):
        short = name.split('_')[0]
        ax.plot(res['time_ns'], res['protein_rmsd'], color=colors[i], 
                alpha=0.7, linewidth=0.5, label=f'{short} (Cα)')
        # Running average
        window = 50
        if len(res['protein_rmsd']) > window:
            avg = pd.Series(res['protein_rmsd']).rolling(window=window, center=True).mean()
            ax.plot(res['time_ns'], avg, color=colors[i], linewidth=2)
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('RMSD (Å)')
    ax.set_title('A) Protein Backbone RMSD')
    ax.legend(loc='lower right')
    ax.set_xlim(0, CONFIG["md_params"]["total_time_ns"])
    
    # ---- Panel B: Ligand RMSD ----
    ax = axes[0, 1]
    for i, (name, res) in enumerate(all_results.items()):
        if 'ligand_rmsd' in res:
            short = name.split('_')[0]
            ax.plot(res['time_ns'], res['ligand_rmsd'], color=colors[i],
                    alpha=0.7, linewidth=0.5, label=f'{short}')
            window = 50
            if len(res['ligand_rmsd']) > window:
                avg = pd.Series(res['ligand_rmsd']).rolling(window=window, center=True).mean()
                ax.plot(res['time_ns'], avg, color=colors[i], linewidth=2)
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('RMSD (Å)')
    ax.set_title('B) Ligand RMSD')
    ax.legend(loc='lower right')
    ax.set_xlim(0, CONFIG["md_params"]["total_time_ns"])
    
    # ---- Panel C: RMSF ----
    ax = axes[1, 0]
    for i, (name, res) in enumerate(all_results.items()):
        short = name.split('_')[0]
        ax.plot(res['rmsf_residues'], res['rmsf'], color=colors[i],
                linewidth=1.5, label=f'{short}', alpha=0.8)
    ax.set_xlabel('Residue Number')
    ax.set_ylabel('RMSF (Å)')
    ax.set_title('C) Per-Residue RMSF (Cα, last 80 ns)')
    ax.legend()
    
    # ---- Panel D: Radius of Gyration ----
    ax = axes[1, 1]
    for i, (name, res) in enumerate(all_results.items()):
        short = name.split('_')[0]
        ax.plot(res['time_ns'], res['rg'], color=colors[i],
                alpha=0.7, linewidth=0.5, label=f'{short}')
        window = 50
        if len(res['rg']) > window:
            avg = pd.Series(res['rg']).rolling(window=window, center=True).mean()
            ax.plot(res['time_ns'], avg, color=colors[i], linewidth=2)
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('Rg (Å)')
    ax.set_title('D) Radius of Gyration')
    ax.legend(loc='lower right')
    ax.set_xlim(0, CONFIG["md_params"]["total_time_ns"])
    
    plt.tight_layout()
    fig_path = os.path.join(output_dir, "md_analysis_4panel.png")
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"✅ Saved: {fig_path}")
    
    # ---- H-bond Analysis Figure ----
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    for i, (name, res) in enumerate(all_results.items()):
        if 'hbond_counts' in res:
            short = name.split('_')[0]
            target = name.split('_')[1]
            eq_start = int(0.2 * len(res['time_ns']))
            time_eq = res['time_ns'][eq_start:]
            
            ax = axes[i]
            ax.plot(time_eq, res['hbond_counts'], color=colors[i], alpha=0.3, linewidth=0.5)
            # Running average
            avg = pd.Series(res['hbond_counts']).rolling(window=50, center=True).mean()
            ax.plot(time_eq, avg, color=colors[i], linewidth=2)
            ax.set_xlabel('Time (ns)')
            ax.set_ylabel('H-bond Count')
            ax.set_title(f'{short}-{target} H-bonds')
            mean_hb = np.mean(res['hbond_counts'])
            ax.axhline(y=mean_hb, color='black', linestyle='--', alpha=0.5,
                       label=f'Mean: {mean_hb:.1f}')
            ax.legend()
    
    plt.tight_layout()
    hbond_fig = os.path.join(output_dir, "md_hbond_analysis.png")
    plt.savefig(hbond_fig, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"✅ Saved: {hbond_fig}")
    
    return fig_path, hbond_fig

fig_path, hbond_fig = create_figures(all_results, CONFIG["output_dir"])

# %% [markdown]
# ## Cell 10: Export Summary Statistics

# %%
def export_summary(all_results, config, output_dir):
    """Export summary statistics as CSV and JSON"""
    
    summary_rows = []
    for name, res in all_results.items():
        sys_conf = [s for s in config["systems"] if s["name"] == name][0]
        eq_start = int(0.2 * len(res["protein_rmsd"]))
        
        row = {
            "System": name,
            "Receptor_PDB": sys_conf["pdb_id"],
            "Gene": sys_conf["gene"],
            "Ligand": sys_conf["ligand_name"],
            "Docking_Affinity_kcal": sys_conf["docking_affinity"],
            "Redock_Affinity_kcal": sys_conf.get("redock_affinity", "N/A"),
            "Control_Drug": sys_conf.get("control_drug", "N/A"),
            "Control_Affinity_kcal": sys_conf.get("control_affinity", "N/A"),
            "Simulation_Time_ns": config["md_params"]["total_time_ns"],
            "Protein_RMSD_Mean_A": f"{np.mean(res['protein_rmsd'][eq_start:]):.2f}",
            "Protein_RMSD_Std_A": f"{np.std(res['protein_rmsd'][eq_start:]):.2f}",
            "Rg_Mean_A": f"{np.mean(res['rg'][eq_start:]):.2f}",
            "Rg_Std_A": f"{np.std(res['rg'][eq_start:]):.2f}",
        }
        
        if "ligand_rmsd" in res:
            row["Ligand_RMSD_Mean_A"] = f"{np.mean(res['ligand_rmsd'][eq_start:]):.2f}"
            row["Ligand_RMSD_Std_A"] = f"{np.std(res['ligand_rmsd'][eq_start:]):.2f}"
        
        if "hbond_counts" in res:
            row["Avg_HBonds"] = f"{np.mean(res['hbond_counts']):.1f}"
            row["Std_HBonds"] = f"{np.std(res['hbond_counts']):.1f}"
        
        summary_rows.append(row)
    
    df = pd.DataFrame(summary_rows)
    
    # Save CSV
    csv_path = os.path.join(output_dir, "md_summary.csv")
    df.to_csv(csv_path, index=False)
    print(f"✅ Summary CSV: {csv_path}")
    
    # Save JSON
    json_path = os.path.join(output_dir, "md_summary.json")
    df.to_json(json_path, orient="records", indent=2)
    print(f"✅ Summary JSON: {json_path}")
    
    # Print table
    print("\n" + "="*70)
    print("MD SIMULATION SUMMARY")
    print("="*70)
    print(df.to_string(index=False))
    
    return df

summary_df = export_summary(all_results, CONFIG, CONFIG["output_dir"])

# %% [markdown]
# ## Cell 11: Package Results for Download

# %%
import shutil
import glob

def package_results(output_dir):
    """Create downloadable zip of all results"""
    print("\n" + "="*60)
    print("Packaging results for download")
    print("="*60)
    
    # List all output files
    all_files = glob.glob(os.path.join(output_dir, "*"))
    
    print(f"\nFiles in output directory:")
    total_size = 0
    for f in sorted(all_files):
        size = os.path.getsize(f)
        total_size += size
        print(f"  {os.path.basename(f):45s} {size/1e6:.1f} MB")
    
    print(f"\n  Total: {total_size/1e9:.2f} GB")
    
    # Create zip (excluding large DCD files — those stay as-is)
    zip_base = "/kaggle/working/md_results_summary"
    summary_files = [f for f in all_files if not f.endswith('.dcd') and not f.endswith('.chk')]
    
    # Create summary zip
    shutil.make_archive(zip_base, 'zip', output_dir,
                        base_dir=None)
    
    zip_path = f"{zip_base}.zip"
    print(f"\n✅ Summary zip: {zip_path} ({os.path.getsize(zip_path)/1e6:.1f} MB)")
    print("\n⬇️ Download this zip from the Output tab in Kaggle")
    
    return zip_path

zip_path = package_results(CONFIG["output_dir"])

# %% [markdown]
# ## Cell 12: Final Report

# %%
print("\n" + "=" * 70)
print("🏁 MOLECULAR DYNAMICS SIMULATION COMPLETE")
print("=" * 70)
print(f"""
Systems simulated:
  1. Phalerin → AGTR1 (3QXY) — {CONFIG['md_params']['total_time_ns']} ns
  2. Mangiferin → RELA (1TBF) — {CONFIG['md_params']['total_time_ns']} ns

Protocol:
  - Force field: AMBER ff14SB (protein) + GAFF2 (ligand)
  - Water model: TIP3P
  - Temperature: {CONFIG['md_params']['temperature_K']} K (37°C)
  - Pressure: {CONFIG['md_params']['pressure_atm']} atm
  - Ionic strength: {CONFIG['md_params']['ionic_strength_M']} M NaCl
  - Timestep: {CONFIG['md_params']['timestep_fs']} fs
  - PME electrostatics, {CONFIG['md_params']['nonbonded_cutoff_nm']} nm cutoff
  - Constraints: H-bonds (SHAKE/LINCS)

Output files:
  - md_analysis_4panel.png — RMSD, RMSF, Rg figures
  - md_hbond_analysis.png — H-bond analysis
  - md_summary.csv — Quantitative summary
  - *_production.dcd — Full trajectories  
  - *_final.pdb — Final structures

All results are in: {CONFIG['output_dir']}
""")
