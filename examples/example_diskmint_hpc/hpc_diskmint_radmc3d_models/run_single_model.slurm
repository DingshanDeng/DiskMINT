#!/bin/bash
#SBATCH --job-name=diskmint_model
#SBATCH --output=logs/diskmint_%A_%a.out
#SBATCH --error=logs/diskmint_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=94
#SBATCH --partition=standard
#SBATCH --account=pascucci

# Load modules
# Activate aconda environment
# source ~/.bashrc
# conda activate diskmint_stable
module load anaconda/2024.06  # or your env

# Read parameters from CSV file based on array task ID
PARAM_FILE=$1 # this is the grid csv
MSTAR=$2
MODEL_DATE=$3
MODEL_MARK=$4
MODEL_PARAM_FILE=$5 # this is the base parameter file this model

# Get the specific row for this array task (skip header)
ROW=$((SLURM_ARRAY_TASK_ID + 2)) # +2 to account for 0-indexing and header row
PARAMS=$(sed -n "${ROW}p" ${PARAM_FILE})

# Parse parameters
RATIO_MDISK=$(echo $PARAMS | cut -d',' -f1)
RATIO_G2D=$(echo $PARAMS   | cut -d',' -f2)
RTAP=$(echo $PARAMS        | cut -d',' -f3)
MODEL_MARK_SHORT=$(echo $PARAMS | cut -d',' -f4)
MODEL_DIR=$(echo $PARAMS   | cut -d',' -f5)
echo "=================================================="
echo "SLURM Job ID: ${SLURM_JOB_ID}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Running on: $(hostname)"
echo "Parameters: MSTAR=${MSTAR}, Mdisk/Mstar=${RATIO_MDISK}, g2d=${RATIO_G2D}, Rtap=${RTAP}"
echo "Model mark: ${MODEL_MARK_SHORT}"
echo "Model date: ${MODEL_DATE}"
echo "Model parameter file: ${MODEL_PARAM_FILE}"
echo "Model directory: ${MODEL_DIR}"
echo "=================================================="

# Run the model
python run_model.py \
    --mstar ${MSTAR} \
    --ratio_mdisk ${RATIO_MDISK} \
    --ratio_g2d ${RATIO_G2D} \
    --rtap ${RTAP} \
    --model_date ${MODEL_DATE} \
    --model_mark ${MODEL_MARK} \
    --model_id ${SLURM_ARRAY_TASK_ID} \
    --param_file ${MODEL_PARAM_FILE} \
    --model_dir ${MODEL_DIR} \

echo "Job completed at: $(date)"
