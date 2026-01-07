"""
PyMOL Docking Pose Visualization Script
Mahkota Dewa - Diabetic Nephropathy Study

This script generates publication-quality images of docking poses.
Run with: pymol -cq generate_poses.py (headless mode)
Or: pymol generate_poses.py (with GUI)
"""

import os
from pymol import cmd, util

# ============================================
# CONFIGURATION
# ============================================

# Base directories (WSL paths)
BASE_DIR = "/mnt/e/Project/herbal-network-pharmacology/data/mahkota_dewa_dn/results/docking"
RECEPTOR_DIR = f"{BASE_DIR}/receptors/prepared"
DOCKED_DIR = f"{BASE_DIR}/mahkota_dewa_dn_docking_results"
OUTPUT_DIR = f"{BASE_DIR}/figures"

# Top 3 docking results (from docking analysis)
TOP_POSES = [
    {
        "name": "264-trihydroxy-4-methoxybenzophenone_PPARG",
        "receptor": "6MS7",
        "ligand": "PPARG_264-trihydroxy-4-methoxybenzophenone_docked",
        "gene": "PPARG",
        "compound": "264-trihydroxy-4-methoxybenzophenone",
        "score": -9.53
    },
    {
        "name": "Phalerin_RELA",
        "receptor": "3QXY",
        "ligand": "RELA_Phalerin_docked",
        "gene": "RELA",
        "compound": "Phalerin",
        "score": -9.26
    },
    {
        "name": "264-trihydroxy-4-methoxybenzophenone_RELA",
        "receptor": "3QXY",
        "ligand": "RELA_264-trihydroxy-4-methoxybenzophenone_docked",
        "gene": "RELA",
        "compound": "264-trihydroxy-4-methoxybenzophenone",
        "score": -9.31
    },
]

# Image settings
IMAGE_WIDTH = 1280
IMAGE_HEIGHT = 720
RAY_TRACE = False  # Disable for faster preview (enable for publication)

# ============================================
# STYLE FUNCTIONS
# ============================================

def setup_style():
    """Configure PyMOL display settings."""
    cmd.bg_color("white")
    cmd.set("ray_trace_mode", 1)
    cmd.set("ray_shadows", 1)
    cmd.set("antialias", 2)
    cmd.set("orthoscopic", 1)
    cmd.set("depth_cue", 0)
    cmd.set("spec_reflect", 0.3)
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("cartoon_side_chain_helper", 1)
    cmd.set("stick_radius", 0.2)

def style_protein(selection="protein"):
    """Apply publication-style to protein."""
    cmd.hide("everything", selection)
    cmd.show("cartoon", selection)
    cmd.color("lightblue", selection)
    cmd.set("cartoon_transparency", 0.3, selection)

def style_ligand(selection="ligand"):
    """Apply publication-style to ligand."""
    cmd.hide("everything", selection)
    cmd.show("sticks", selection)
    cmd.util.cbag(selection)  # Color by atom, green carbons
    cmd.set("stick_radius", 0.25, selection)

def style_binding_site(protein_sel="protein", ligand_sel="ligand", distance=5.0):
    """Show and style residues within binding site."""
    binding_site = f"byres ({protein_sel}) within {distance} of ({ligand_sel})"
    cmd.select("binding_site", binding_site)
    cmd.show("sticks", "binding_site")
    cmd.util.cbay("binding_site")  # Color by atom, yellow carbons
    cmd.set("stick_radius", 0.15, "binding_site")
    cmd.label("binding_site and name CA", "resn+resi")

def add_surface_around_ligand(protein_sel="protein", ligand_sel="ligand", distance=6.0):
    """Add surface representation around binding site."""
    pocket = f"byres ({protein_sel}) within {distance} of ({ligand_sel})"
    cmd.select("pocket_surface", pocket)
    cmd.show("surface", "pocket_surface")
    cmd.set("surface_color", "white", "pocket_surface")
    cmd.set("transparency", 0.7, "pocket_surface")

# ============================================
# MAIN RENDERING FUNCTION
# ============================================

def render_pose(pose_info, output_dir):
    """Render a single docking pose."""
    print(f"\n{'='*60}")
    print(f"Rendering: {pose_info['name']}")
    print(f"{'='*60}")
    
    # Clear previous
    cmd.delete("all")
    
    # Load files
    receptor_file = f"{RECEPTOR_DIR}/{pose_info['receptor']}.pdbqt"
    ligand_file = f"{DOCKED_DIR}/{pose_info['ligand']}.pdbqt"
    
    print(f"Loading receptor: {receptor_file}")
    print(f"Loading ligand: {ligand_file}")
    
    if not os.path.exists(receptor_file):
        print(f"ERROR: Receptor not found: {receptor_file}")
        return
    if not os.path.exists(ligand_file):
        print(f"ERROR: Ligand not found: {ligand_file}")
        return
    
    cmd.load(receptor_file, "protein")
    cmd.load(ligand_file, "ligand")
    
    # Apply styles
    setup_style()
    style_protein("protein")
    style_ligand("ligand")
    style_binding_site("protein", "ligand", distance=5.0)
    # add_surface_around_ligand("protein", "ligand")  # Uncomment for surface
    
    # Position camera
    cmd.orient("ligand")
    cmd.zoom("ligand", 8)
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Generate image
    output_file = f"{output_dir}/{pose_info['name']}.png"
    
    if RAY_TRACE:
        cmd.ray(IMAGE_WIDTH, IMAGE_HEIGHT)
    
    cmd.png(output_file, width=IMAGE_WIDTH, height=IMAGE_HEIGHT, dpi=300)
    print(f"Saved: {output_file}")
    
    # Also save a session file for later editing
    session_file = f"{output_dir}/{pose_info['name']}.pse"
    cmd.save(session_file)
    print(f"Session saved: {session_file}")

# ============================================
# MAIN EXECUTION
# ============================================

def main():
    print("\n" + "="*60)
    print("DOCKING POSE VISUALIZATION")
    print("="*60)
    print(f"\nGenerating {len(TOP_POSES)} images...")
    print(f"Output directory: {OUTPUT_DIR}")
    
    for pose in TOP_POSES:
        render_pose(pose, OUTPUT_DIR)
    
    print("\n" + "="*60)
    print("COMPLETE!")
    print("="*60)
    print(f"\nImages saved to: {OUTPUT_DIR}")
    print("\nTip: Open .pse session files in PyMOL GUI to fine-tune views")

if __name__ == "__main__" or True:  # Always run in PyMOL
    main()
