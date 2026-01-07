"""
Step 4: Run Molecular Docking
==============================

Runs AutoDock Vina docking for all receptor-ligand pairs.

Input: Prepared receptors and ligands from Steps 2 & 3
Output: Docking poses and scores

Usage:
    python -m src.docking.step4_run_docking --receptors receptors/prepared --ligands ligands/pdbqt
"""

import os
import subprocess
import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import click
from dataclasses import dataclass


@dataclass
class DockingJob:
    """A single docking job."""
    receptor_name: str
    receptor_file: Path
    ligand_name: str
    ligand_file: Path
    output_file: Path
    log_file: Path
    binding_site: Dict[str, float]


def check_vina():
    """Check if AutoDock Vina is available."""
    try:
        result = subprocess.run(['vina', '--version'], capture_output=True, text=True)
        return result.returncode == 0
    except FileNotFoundError:
        return False


def run_vina(job: DockingJob, exhaustiveness: int = 8, num_modes: int = 9) -> Optional[Dict]:
    """
    Run AutoDock Vina for a single job.
    
    Returns:
        Dictionary with results or None if failed
    """
    job.output_file.parent.mkdir(parents=True, exist_ok=True)
    
    cmd = [
        'vina',
        '--receptor', str(job.receptor_file),
        '--ligand', str(job.ligand_file),
        '--out', str(job.output_file),
        '--log', str(job.log_file),
        '--center_x', str(job.binding_site.get('center_x', 0)),
        '--center_y', str(job.binding_site.get('center_y', 0)),
        '--center_z', str(job.binding_site.get('center_z', 0)),
        '--size_x', str(job.binding_site.get('size_x', 25)),
        '--size_y', str(job.binding_site.get('size_y', 25)),
        '--size_z', str(job.binding_site.get('size_z', 25)),
        '--exhaustiveness', str(exhaustiveness),
        '--num_modes', str(num_modes),
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        
        if result.returncode != 0:
            return {'error': result.stderr}
        
        # Parse log for results
        return parse_vina_log(job.log_file)
        
    except subprocess.TimeoutExpired:
        return {'error': 'Timeout (10 min)'}
    except FileNotFoundError:
        return {'error': 'Vina not found'}


def parse_vina_log(log_file: Path) -> Dict:
    """Parse Vina log file for binding affinities."""
    results = {
        'modes': [],
        'best_affinity': None
    }
    
    if not log_file.exists():
        return {'error': 'Log file not found'}
    
    with open(log_file, 'r') as f:
        in_results = False
        for line in f:
            if '-----+------------' in line:
                in_results = True
                continue
            
            if in_results and line.strip():
                parts = line.split()
                if len(parts) >= 4 and parts[0].isdigit():
                    try:
                        mode = int(parts[0])
                        affinity = float(parts[1])
                        rmsd_lb = float(parts[2])
                        rmsd_ub = float(parts[3])
                        
                        results['modes'].append({
                            'mode': mode,
                            'affinity': affinity,
                            'rmsd_lb': rmsd_lb,
                            'rmsd_ub': rmsd_ub
                        })
                        
                        if results['best_affinity'] is None or affinity < results['best_affinity']:
                            results['best_affinity'] = affinity
                    except (ValueError, IndexError):
                        continue
    
    return results


@click.command()
@click.option('--receptors', '-r', required=True, help='Directory with receptor PDBQT files')
@click.option('--ligands', '-l', required=True, help='Directory with ligand PDBQT files')
@click.option('--output-dir', '-o', default=None, help='Output directory')
@click.option('--binding-sites', '-b', default=None, help='JSON file with binding site coordinates')
@click.option('--exhaustiveness', '-e', default=8, help='Vina exhaustiveness')
@click.option('--dry-run', is_flag=True, help='Show jobs without running')
def main(receptors: str, ligands: str, output_dir: str, binding_sites: str, 
         exhaustiveness: int, dry_run: bool):
    """Step 4: Run molecular docking with AutoDock Vina."""
    print("=" * 50)
    print("STEP 4: Run Molecular Docking")
    print("=" * 50)
    
    # Check Vina
    if not check_vina():
        print("\n❌ ERROR: AutoDock Vina not found in PATH!")
        print("Download from: https://vina.scripps.edu/")
        return
    print("✅ AutoDock Vina found")
    
    receptor_dir = Path(receptors)
    ligand_dir = Path(ligands)
    
    if output_dir is None:
        output_dir = Path("data/mahkota_dewa_dn/results/docking/outputs")
    else:
        output_dir = Path(output_dir)
    
    # Load binding sites
    sites = {}
    if binding_sites:
        with open(binding_sites) as f:
            sites = json.load(f)
    else:
        # Look for binding_sites.json in receptor dir parent
        sites_file = receptor_dir / "binding_sites.json"
        if sites_file.exists():
            with open(sites_file) as f:
                sites = json.load(f)
    
    # Get files
    receptor_files = list(receptor_dir.glob("*.pdbqt"))
    ligand_files = list(ligand_dir.glob("*.pdbqt"))
    
    print(f"\n📦 Receptors: {len(receptor_files)}")
    print(f"💊 Ligands: {len(ligand_files)}")
    print(f"🔢 Total jobs: {len(receptor_files) * len(ligand_files)}")
    
    # Create jobs
    jobs = []
    for receptor_file in receptor_files:
        receptor_name = receptor_file.stem.replace('_clean', '')
        
        # Get binding site
        site = sites.get(receptor_name, {
            'center_x': 0, 'center_y': 0, 'center_z': 0,
            'size_x': 25, 'size_y': 25, 'size_z': 25
        })
        
        for ligand_file in ligand_files:
            ligand_name = ligand_file.stem
            
            job = DockingJob(
                receptor_name=receptor_name,
                receptor_file=receptor_file,
                ligand_name=ligand_name,
                ligand_file=ligand_file,
                output_file=output_dir / f"{ligand_name}_{receptor_name}.pdbqt",
                log_file=output_dir / f"{ligand_name}_{receptor_name}.log",
                binding_site=site
            )
            jobs.append(job)
    
    if dry_run:
        print("\n🔍 DRY RUN - Jobs that would run:")
        for job in jobs[:10]:
            print(f"  {job.ligand_name} → {job.receptor_name}")
        if len(jobs) > 10:
            print(f"  ... and {len(jobs) - 10} more")
        return
    
    # Run docking
    results = []
    failed = []
    
    print(f"\n🚀 Running {len(jobs)} docking jobs...\n")
    
    for i, job in enumerate(jobs, 1):
        print(f"[{i}/{len(jobs)}] {job.ligand_name} → {job.receptor_name}", end=" ")
        
        result = run_vina(job, exhaustiveness=exhaustiveness)
        
        if 'error' in result:
            print(f"❌ {result['error']}")
            failed.append((job, result['error']))
        else:
            best = result.get('best_affinity', 'N/A')
            print(f"✅ {best} kcal/mol")
            results.append({
                'receptor': job.receptor_name,
                'ligand': job.ligand_name,
                'best_affinity': best,
                'modes': result.get('modes', [])
            })
    
    # Save results
    results_file = output_dir / "docking_results.json"
    with open(results_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    # Save CSV summary
    csv_file = output_dir / "docking_summary.csv"
    with open(csv_file, 'w') as f:
        f.write("ligand,receptor,best_affinity_kcal_mol\n")
        for r in sorted(results, key=lambda x: x['best_affinity'] or 0):
            f.write(f"{r['ligand']},{r['receptor']},{r['best_affinity']}\n")
    
    print("\n" + "=" * 50)
    print("SUMMARY")
    print("=" * 50)
    print(f"✅ Completed: {len(results)}")
    print(f"❌ Failed: {len(failed)}")
    print(f"📁 Results: {results_file}")
    print(f"📊 Summary: {csv_file}")
    print(f"\n➡️  Next: Run step5_analyze.py")


if __name__ == "__main__":
    main()
