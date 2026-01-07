"""
Step 5: Analyze Docking Results
================================

Analyzes docking results and generates summary reports.

Input: Docking results from Step 4
Output: Analysis reports, rankings, visualizations

Usage:
    python -m src.docking.step5_analyze --results docking_results.json
"""

import json
from pathlib import Path
from typing import Dict, List
import click

try:
    import pandas as pd
    HAS_PANDAS = True
except ImportError:
    HAS_PANDAS = False


def load_results(results_file: Path) -> List[Dict]:
    """Load docking results from JSON."""
    with open(results_file) as f:
        return json.load(f)


def create_affinity_matrix(results: List[Dict]) -> Dict:
    """Create ligand × receptor affinity matrix."""
    matrix = {}
    
    for r in results:
        ligand = r['ligand']
        receptor = r['receptor']
        affinity = r.get('best_affinity')
        
        if ligand not in matrix:
            matrix[ligand] = {}
        matrix[ligand][receptor] = affinity
    
    return matrix


def rank_by_affinity(results: List[Dict]) -> List[Dict]:
    """Rank all results by binding affinity."""
    valid = [r for r in results if r.get('best_affinity') is not None]
    return sorted(valid, key=lambda x: x['best_affinity'])


def rank_by_target(results: List[Dict]) -> Dict[str, List]:
    """Rank ligands for each target."""
    by_target = {}
    
    for r in results:
        receptor = r['receptor']
        if receptor not in by_target:
            by_target[receptor] = []
        by_target[receptor].append(r)
    
    # Sort each target's results
    for receptor in by_target:
        by_target[receptor] = sorted(
            by_target[receptor],
            key=lambda x: x.get('best_affinity') or 0
        )
    
    return by_target


def compare_with_controls(results: List[Dict], control_names: List[str]) -> Dict:
    """Compare test compounds with control drugs."""
    controls = {}
    test_compounds = {}
    
    for r in results:
        ligand = r['ligand']
        if ligand in control_names:
            controls[ligand] = r
        else:
            test_compounds[ligand] = r
    
    # Find compounds that beat controls
    better_than_control = []
    for ligand, result in test_compounds.items():
        affinity = result.get('best_affinity')
        control_affinity = None
        
        # Find control for this target
        for ctrl_name, ctrl_result in controls.items():
            if ctrl_result['receptor'] == result['receptor']:
                control_affinity = ctrl_result.get('best_affinity')
                break
        
        if affinity and control_affinity and affinity < control_affinity:
            better_than_control.append({
                'ligand': ligand,
                'receptor': result['receptor'],
                'affinity': affinity,
                'control': ctrl_name,
                'control_affinity': control_affinity,
                'improvement': control_affinity - affinity
            })
    
    return {
        'controls': controls,
        'test_compounds': test_compounds,
        'better_than_control': better_than_control
    }


@click.command()
@click.option('--results', '-r', required=True, help='Docking results JSON file')
@click.option('--output-dir', '-o', default=None, help='Output directory')
@click.option('--controls', '-c', multiple=True, default=['Pioglitazone', 'Atorvastatin', 'Losartan', 'Sildenafil'],
              help='Control compound names')
def main(results: str, output_dir: str, controls: tuple):
    """Step 5: Analyze docking results."""
    print("=" * 50)
    print("STEP 5: Analyze Docking Results")
    print("=" * 50)
    
    results_file = Path(results)
    data = load_results(results_file)
    
    if output_dir is None:
        output_dir = results_file.parent / "analysis"
    else:
        output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"\nLoaded {len(data)} docking results")
    
    # 1. Overall ranking
    print("\n📊 Top 10 Overall Binding Affinities:")
    print("-" * 50)
    ranked = rank_by_affinity(data)
    for i, r in enumerate(ranked[:10], 1):
        print(f"  {i}. {r['ligand']} → {r['receptor']}: {r['best_affinity']:.2f} kcal/mol")
    
    # 2. Ranking by target
    print("\n📊 Best Ligand per Target:")
    print("-" * 50)
    by_target = rank_by_target(data)
    for receptor, ligands in by_target.items():
        if ligands:
            best = ligands[0]
            print(f"  {receptor}: {best['ligand']} ({best['best_affinity']:.2f} kcal/mol)")
    
    # 3. Compare with controls
    print("\n📊 Comparison with Control Drugs:")
    print("-" * 50)
    comparison = compare_with_controls(data, list(controls))
    
    better = comparison['better_than_control']
    if better:
        print(f"  ✅ {len(better)} compounds outperform controls:")
        for item in better:
            print(f"     {item['ligand']} ({item['affinity']:.2f}) > {item['control']} ({item['control_affinity']:.2f})")
    else:
        print("  ⚠️  No compounds significantly outperform controls")
    
    # 4. Save reports
    
    # JSON full results
    analysis_file = output_dir / "analysis_full.json"
    analysis_data = {
        'total_results': len(data),
        'top_10_overall': ranked[:10],
        'best_per_target': {k: v[0] if v else None for k, v in by_target.items()},
        'comparison': comparison
    }
    with open(analysis_file, 'w') as f:
        json.dump(analysis_data, f, indent=2, default=str)
    
    # CSV ranking
    ranking_file = output_dir / "ranking_overall.csv"
    with open(ranking_file, 'w') as f:
        f.write("rank,ligand,receptor,affinity_kcal_mol\n")
        for i, r in enumerate(ranked, 1):
            f.write(f"{i},{r['ligand']},{r['receptor']},{r['best_affinity']}\n")
    
    # Affinity matrix
    matrix = create_affinity_matrix(data)
    matrix_file = output_dir / "affinity_matrix.csv"
    if HAS_PANDAS:
        df = pd.DataFrame(matrix).T
        df.to_csv(matrix_file)
    else:
        # Simple CSV
        receptors = list(set(r['receptor'] for r in data))
        with open(matrix_file, 'w') as f:
            f.write("ligand," + ",".join(receptors) + "\n")
            for ligand, values in matrix.items():
                row = [str(values.get(r, '')) for r in receptors]
                f.write(f"{ligand}," + ",".join(row) + "\n")
    
    print("\n" + "=" * 50)
    print("OUTPUTS")
    print("=" * 50)
    print(f"📄 Full analysis: {analysis_file}")
    print(f"📊 Overall ranking: {ranking_file}")
    print(f"📈 Affinity matrix: {matrix_file}")
    print(f"\n🎉 Docking analysis complete!")


if __name__ == "__main__":
    main()
