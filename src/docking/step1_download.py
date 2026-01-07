"""
Step 1: Download PDB Structures
================================

Downloads protein structures from RCSB PDB database.
Output: raw PDB files in receptors/raw/

Usage:
    python -m src.docking.step1_download --config config/docking_config.yaml
"""

import os
import urllib.request
from pathlib import Path
from typing import Dict, List, Optional
import yaml
import click


def download_pdb(pdb_id: str, output_dir: Path) -> Optional[Path]:
    """
    Download PDB file from RCSB.
    
    Args:
        pdb_id: PDB ID (e.g., "6MS7")
        output_dir: Directory to save file
        
    Returns:
        Path to downloaded file or None if failed
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    pdb_file = output_dir / f"{pdb_id}.pdb"
    
    if pdb_file.exists():
        print(f"  ⏭️  {pdb_id}.pdb already exists, skipping...")
        return pdb_file
    
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    
    try:
        print(f"  📥 Downloading {pdb_id} from RCSB...")
        urllib.request.urlretrieve(url, pdb_file)
        print(f"  ✅ Saved to {pdb_file}")
        return pdb_file
    except Exception as e:
        print(f"  ❌ Failed to download {pdb_id}: {e}")
        return None


def load_config(config_path: str) -> Dict:
    """Load docking configuration."""
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)


@click.command()
@click.option('--config', '-c', required=True, help='Path to docking_config.yaml')
@click.option('--output-dir', '-o', default=None, help='Output directory')
def main(config: str, output_dir: str):
    """Step 1: Download PDB structures for docking."""
    print("=" * 50)
    print("STEP 1: Download PDB Structures")
    print("=" * 50)
    
    cfg = load_config(config)
    targets = cfg.get('targets', {})
    
    if output_dir is None:
        output_dir = Path("data/mahkota_dewa_dn/results/docking/receptors/raw")
    else:
        output_dir = Path(output_dir)
    
    downloaded = []
    failed = []
    
    for target_name, target_info in targets.items():
        pdb_id = target_info.get('pdb_id')
        print(f"\n[{target_name}] PDB: {pdb_id}")
        
        result = download_pdb(pdb_id, output_dir)
        if result:
            downloaded.append((target_name, pdb_id, result))
        else:
            failed.append((target_name, pdb_id))
    
    # Summary
    print("\n" + "=" * 50)
    print("SUMMARY")
    print("=" * 50)
    print(f"✅ Downloaded: {len(downloaded)}")
    print(f"❌ Failed: {len(failed)}")
    print(f"📁 Output: {output_dir}")
    
    # Save manifest
    manifest_file = output_dir / "download_manifest.txt"
    with open(manifest_file, 'w') as f:
        f.write("# Step 1: Download PDB Manifest\n")
        f.write(f"# Total: {len(downloaded)} downloaded, {len(failed)} failed\n\n")
        for name, pdb_id, path in downloaded:
            f.write(f"{name}\t{pdb_id}\t{path}\n")
    
    print(f"\n📄 Manifest saved to: {manifest_file}")
    print("\n➡️  Next: Run step2_prepare_receptor.py")


if __name__ == "__main__":
    main()
