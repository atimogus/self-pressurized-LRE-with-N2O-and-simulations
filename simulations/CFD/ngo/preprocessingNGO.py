#!/usr/bin/env python3

import os
import shutil
import subprocess
import glob

def move_and_generate_mesh():
    # Get directory where main python file is run
    current_dir = os.path.dirname(os.path.abspath(__file__))
    print(f"Current directory: {current_dir}")
    
    # Read all files in the directory
    files = os.listdir(current_dir)
    print(f"Files in directory: {files}")
    
    # Find .msh2 files
    msh2_files = [f for f in files if f.endswith('.msh2')]
    
    if not msh2_files:
        print("No .msh2 files found in the current directory!")
        return
    
    # Use the first .msh2 file found
    msh_file = msh2_files[0]
    print(f"Found .msh2 file: {msh_file}")
    
    # Define the target directory (raketa is the main directory)
    raketa_dir = os.path.dirname(current_dir)  # Go up to raketa directory
    target_dir = os.path.join(raketa_dir, "ngo", "sonicFoam")
    
    if not os.path.exists(target_dir):
        print(f"Target directory does not exist: {target_dir}")
        return
    
    # Copy the .msh2 file to the target directory (overwrite if exists)
    source_path = os.path.join(current_dir, msh_file)
    target_path = os.path.join(target_dir, msh_file)
    
    try:
        shutil.copy2(source_path, target_path)
        print(f"Successfully copied {msh_file} to {target_dir}")
    except Exception as e:
        print(f"Error copying file: {e}")
        return
    
    # Change to the target directory
    os.chdir(target_dir)
    print(f"Changed directory to: {target_dir}")
    
    # # Source OpenFOAM environment and run gmshToFoam
    # try:
    #     # First, source OpenFOAM and run gmshToFoam
    #     cmd = f"source /usr/lib/openfoam/openfoam2412/etc/bashrc && gmshToFoam {msh_file}"
    #     result = subprocess.run(cmd, shell=True, capture_output=True, text=True, executable='/bin/bash')
        
    #     if result.returncode == 0:
    #         print("gmshToFoam conversion successful!")
    #         print(result.stdout)
    #     else:
    #         print("gmshToFoam conversion failed!")
    #         print(result.stderr)
    #         return
            
    # except Exception as e:
    #     print(f"Error running gmshToFoam: {e}")
    #     return
    
    # # Run Allrun script
    # try:
    #     allrun_path = os.path.join(target_dir, "openfoamslurm.sh")
    #     if os.path.exists(allrun_path):
    #         # Make sure Allrun is executable
    #         os.chmod(allrun_path, 0o755)
            
    #         # Run Allrun script with OpenFOAM environment
    #         cmd = f"sbatch {allrun_path}"
    #         print("Starting simulation with Allrun...")
            
    #         # Run in background since this might take a long time
    #         process = subprocess.Popen(cmd, shell=True, executable='/bin/bash')
    #         print(f"Simulation started with PID: {process.pid}")
    #         print("The simulation is running in the background.")
            
    #     else:
    #         print(f"Allrun script not found in {target_dir}")
            
    # except Exception as e:
    #     print(f"Error running Allrun: {e}")

if __name__ == "__main__":
    move_and_generate_mesh()
    print("Preprocessing NGO script executed successfully.")
