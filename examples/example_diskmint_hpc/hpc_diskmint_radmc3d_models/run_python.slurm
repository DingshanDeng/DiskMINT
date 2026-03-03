#!/bin/bash
#SBATCH --job-name=py_helper
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --output=logs/py_helper_%j.out
#SBATCH --error=logs/py_helper_%j.err
#SBATCH --partition=standard
#SBATCH --account=pascucci

# Load your conda env
# source activate your_env  # or: module load anaconda; conda activate your_env
load anaconda/2024.06 

python "$@"