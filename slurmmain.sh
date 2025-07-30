#!/bin/bash
#SBATCH --job-name=openfoam_ngo_simulation
#SBATCH --partition=computes_thin
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --time=12:00:00
#SBATCH --mem=48G
#SBATCH --output=openfoam_%j.log
#SBATCH --error=openfoam_%j.log   # same file for error

echo "=========================================="
echo "Starting SLURM job at: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURMD_NODENAME"
echo "Working directory: $(pwd)"
echo "Python version: $(python --version)"
echo "Python path: $(which python)"
echo "Python path in venv: $(which python3)"
echo "Conda environment: $CONDA_DEFAULT_ENV"
echo "=========================================="

module load OpenFOAM/OpenFOAM_2312-IntelMPI_2021.5 

echo "==========================================" # rm -r processor* log.*
echo "Starting main.py execution..."
echo "=========================================="

export OMP_NUM_THREADS=24  # or the number of cores you want to use

# Run the main Python script that will handle everything
python3 -u main.py


echo "=========================================="
echo "✓ Pipeline completed successfully!"
echo "Job finished at: $(date)"
echo "=========================================="


cd ./simulations/CFD/ngo/sonicFoam/

rm -r processor* log.* *00 

echo "Workcurrent changed directory: $(pwd)"

. ${WM_PROJECT_DIR:?}/bin/tools/RunFunctions 

runApplication gmshToFoam structured_nozzle_3D.msh2

runApplication changeDictionary

runApplication decomposePar 

runParallel $(getApplication)

runApplication reconstructPar


echo "=========================================="
echo "✓ code completed successfully!"
echo "Job finished at: $(date)"
echo "=========================================="