"""
Kaggle API Helper for MD Simulation

Commands:
    python kaggle_helper.py push      - Push notebook to Kaggle
    python kaggle_helper.py status    - Check kernel status  
    python kaggle_helper.py output    - Get kernel output/errors
    python kaggle_helper.py download  - Download results

Requirements:
    1. pip install kaggle python-dotenv
    2. Set KAGGLE_API_TOKEN in .env file
    3. Get token from: https://www.kaggle.com/settings -> API -> Create New Token
"""

import argparse
import json
import os
import sys
from pathlib import Path

# Load environment variables from .env
try:
    from dotenv import load_dotenv
    env_path = Path(__file__).parent.parent / ".env"
    load_dotenv(env_path)
except ImportError:
    pass  # dotenv not installed, use system env vars

# Set Kaggle API token from environment
if os.getenv("KAGGLE_API_TOKEN"):
    os.environ["KAGGLE_KEY"] = os.getenv("KAGGLE_API_TOKEN")

try:
    from kaggle.api.kaggle_api_extended import KaggleApi
except ImportError:
    print("❌ Kaggle not installed. Run: pip install kaggle")
    sys.exit(1)

# ============================================
# CONFIGURATION
# ============================================

CONFIG = {
    # Kaggle username (from your profile URL)
    "username": "engkinandatama",
    
    # Kernel/Notebook slug (will be created from notebook name)
    "kernel_slug": "md-simulation-264thm-pparg",
    
    # Notebook path
    "notebook_path": Path(__file__).parent.parent / "notebooks" / "kaggle" / "md_simulation_264THM_PPARG.ipynb",
    
    # Dataset dependencies
    "dataset_sources": [
        "engkinandatama/md-simulation-264thm-pparg-input"
    ],
    
    # GPU settings
    "enable_gpu": True,
    "enable_internet": True,
}

# ============================================
# API FUNCTIONS
# ============================================

def get_api():
    """Get authenticated Kaggle API."""
    api = KaggleApi()
    api.authenticate()
    return api


def push_notebook():
    """Push notebook to Kaggle."""
    api = get_api()
    
    notebook_path = CONFIG["notebook_path"]
    if not notebook_path.exists():
        print(f"❌ Notebook not found: {notebook_path}")
        return False
    
    # Create kernel metadata
    kernel_metadata = {
        "id": f"{CONFIG['username']}/{CONFIG['kernel_slug']}",
        "title": "MD Simulation - 264THM PPARG",
        "code_file": str(notebook_path),
        "language": "python",
        "kernel_type": "notebook",
        "is_private": True,
        "enable_gpu": CONFIG["enable_gpu"],
        "enable_internet": CONFIG["enable_internet"],
        "dataset_sources": CONFIG["dataset_sources"],
        "competition_sources": [],
        "kernel_sources": [],
    }
    
    # Save metadata
    metadata_path = notebook_path.parent / "kernel-metadata.json"
    with open(metadata_path, "w", encoding="utf-8") as f:
        json.dump(kernel_metadata, f, indent=2)
    
    print(f"📤 Pushing notebook: {CONFIG['kernel_slug']}")
    
    try:
        api.kernels_push(str(notebook_path.parent))
        print("✅ Notebook pushed successfully!")
        print(f"🔗 URL: https://www.kaggle.com/code/{CONFIG['username']}/{CONFIG['kernel_slug']}")
        return True
    except Exception as e:
        print(f"❌ Push failed: {e}")
        return False


def check_status():
    """Check kernel status."""
    api = get_api()
    
    kernel_ref = f"{CONFIG['username']}/{CONFIG['kernel_slug']}"
    print(f"📊 Checking status: {kernel_ref}")
    
    try:
        status = api.kernel_status(kernel_ref)
        print(f"\n{'='*50}")
        print(f"Status: {status.get('status', 'Unknown')}")
        print(f"Failed: {status.get('failureMessage', 'None')}")
        print(f"{'='*50}")
        return status
    except Exception as e:
        print(f"❌ Status check failed: {e}")
        return None


def get_output():
    """Get kernel output and errors."""
    api = get_api()
    
    kernel_ref = f"{CONFIG['username']}/{CONFIG['kernel_slug']}"
    print(f"📋 Getting output: {kernel_ref}")
    
    try:
        output = api.kernel_output(kernel_ref, path="./kaggle_output")
        
        # Try to get log
        log_path = Path("./kaggle_output") / f"{CONFIG['kernel_slug']}.log"
        if log_path.exists():
            print(f"\n{'='*50}")
            print("📜 KERNEL LOG:")
            print(f"{'='*50}")
            with open(log_path) as f:
                # Print last 100 lines
                lines = f.readlines()
                for line in lines[-100:]:
                    print(line, end="")
        
        print(f"\n✅ Output saved to: ./kaggle_output/")
        return True
    except Exception as e:
        print(f"❌ Output retrieval failed: {e}")
        return False


def download_results():
    """Download kernel results."""
    api = get_api()
    
    kernel_ref = f"{CONFIG['username']}/{CONFIG['kernel_slug']}"
    output_dir = Path("./kaggle_results")
    output_dir.mkdir(exist_ok=True)
    
    print(f"📥 Downloading results: {kernel_ref}")
    
    try:
        api.kernel_output(kernel_ref, path=str(output_dir))
        
        print(f"✅ Results downloaded to: {output_dir}")
        
        # List files
        print("\n📁 Downloaded files:")
        for f in output_dir.rglob("*"):
            if f.is_file():
                size_kb = f.stat().st_size / 1024
                print(f"  - {f.relative_to(output_dir)} ({size_kb:.1f} KB)")
        
        return True
    except Exception as e:
        print(f"❌ Download failed: {e}")
        return False


def list_kernels():
    """List all kernels."""
    api = get_api()
    
    print(f"📋 Listing kernels for: {CONFIG['username']}")
    
    try:
        kernels = api.kernels_list(user=CONFIG['username'])
        
        print(f"\n{'='*60}")
        for k in kernels:
            print(f"📓 {k.ref}")
            # Handle different API versions
            if hasattr(k, 'status'):
                print(f"   Status: {k.status}")
            if hasattr(k, 'lastRunTime'):
                print(f"   Last run: {k.lastRunTime}")
            print()
        print(f"{'='*60}")
        print(f"✅ Found {len(kernels)} kernel(s)")
        
        return kernels
    except Exception as e:
        print(f"❌ List failed: {e}")
        return None


def create_dataset():
    """Create dataset from input files."""
    api = get_api()
    
    input_dir = Path(__file__).parent.parent / "data" / "mahkota_dewa_dn" / "md_input" / "264THM_PPARG"
    
    if not input_dir.exists():
        print(f"❌ Input directory not found: {input_dir}")
        return False
    
    # Create dataset metadata
    dataset_metadata = {
        "title": "MD Simulation Input - 264THM PPARG",
        "id": f"{CONFIG['username']}/md-simulation-264thm-pparg-input",
        "licenses": [{"name": "CC0-1.0"}],
    }
    
    metadata_path = input_dir / "dataset-metadata.json"
    with open(metadata_path, "w") as f:
        json.dump(dataset_metadata, f, indent=2)
    
    print(f"📤 Creating dataset from: {input_dir}")
    
    try:
        api.dataset_create_new(str(input_dir), dir_mode="zip")
        print("✅ Dataset created successfully!")
        print(f"🔗 URL: https://www.kaggle.com/datasets/{CONFIG['username']}/md-simulation-264thm-pparg-input")
        return True
    except Exception as e:
        print(f"❌ Dataset creation failed: {e}")
        return False


# ============================================
# MAIN
# ============================================

def main():
    parser = argparse.ArgumentParser(description="Kaggle API Helper")
    parser.add_argument(
        "command",
        choices=["push", "status", "output", "download", "list", "create-dataset"],
        help="Command to run"
    )
    
    args = parser.parse_args()
    
    print(f"🚀 Kaggle Helper - {args.command}")
    print(f"{'='*50}\n")
    
    if args.command == "push":
        push_notebook()
    elif args.command == "status":
        check_status()
    elif args.command == "output":
        get_output()
    elif args.command == "download":
        download_results()
    elif args.command == "list":
        list_kernels()
    elif args.command == "create-dataset":
        create_dataset()


if __name__ == "__main__":
    main()
