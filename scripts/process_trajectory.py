import MDAnalysis as mda
from MDAnalysis.transformations import fit_rot_trans
from MDAnalysis.analysis import align
import warnings
import os

# Suppress warnings
warnings.filterwarnings('ignore')

def process_trajectory():
    print("=== MD Trajectory Pre-processor (Python/MDAnalysis) ===")
    
    # Path Configuration
    top_file = r"kaggle-output-sim\264THM_PPARG\md\md.gro"
    traj_file = r"kaggle-output-sim\264THM_PPARG\md\full_trajectory.xtc"
    output_file = r"visual_ready.xtc"
    
    # Parameters
    SKIP_FRAMES = 10  # Downsample factor
    
    print(f"Loading topology: {top_file}")
    print(f"Loading trajectory: {traj_file}")
    
    if not os.path.exists(top_file) or not os.path.exists(traj_file):
        print("ERROR: Input files not found!")
        return

    try:
        u = mda.Universe(top_file, traj_file)
        print(f"Total frames: {len(u.trajectory)}")
        print(f"Processing with skip={SKIP_FRAMES} (Expected output: {len(u.trajectory)//SKIP_FRAMES} frames)")
        
        # Define atoms
        protein = u.select_atoms("protein")
        all_atoms = u.select_atoms("all")
        
        # Setup workflow to write new trajectory
        # 1. Align protein to ref (removes rotation/translation)
        # We use the first frame as reference
        ref = mda.Universe(top_file, top_file)
        
        print("Starting processing... (This may take a few minutes)")
        
        # We will iterate and write
        with mda.Writer(output_file, all_atoms.n_atoms) as W:
            for ts in u.trajectory[::SKIP_FRAMES]:
                # On-the-fly alignment
                # Align current frame protein to reference protein
                # This centers the protein and removes global rotation
                align.alignto(u, ref, select="protein", weights="mass")
                
                # Write frame
                W.write(all_atoms)
                
                if ts.frame % (SKIP_FRAMES * 50) == 0:
                    print(f"Processed frame {ts.frame}/{len(u.trajectory)}")
                    
        print(f"Done! Saved to: {output_file}")
        print("You can now open 'md.gro' and 'visual_ready.xtc' in ChimeraX.")
        
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    process_trajectory()
