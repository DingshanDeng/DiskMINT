#!/bin/bash
# Master script to submit DiskMINT model grids for different stellar masses

set -e  # Exit on error

# Configuration
MODEL_DATE="20260222_0" # Model date identifier for organizing outputs
# MODEL_DIR="/groups/pascucci/dingshandeng/diskmint_model_grids_0p3ms_test/radmc3d_model_grids_test/" # Base
# Now the model dir is set up directly by running the generate_grid.py on the local machine,
# And it is stored inside the 4th col of $GRID_FILE (grids/grid_x.csv file).
# so we don't need to specify it here. --- IGNORE ---
WORK_DIR=$(pwd)
MSTAR=0.7 # Stellar mass for this grid (in Msun)
GRID_TYPE="supplementary" # Type of grid to generate (e.g., "default", "extended", "test", "supplementary")
MODEL_MARK="hpc_grid0p7msolar" # Model mark/label for this grid (used in output naming)
MODEL_PARAM_FILE="model_parameters_0p7Msolar.csv" # Base parameter CSV file for this grid (should be in the same directory as this script or specify path)
# GRID_FILE="grids/grid_${MSTAR}msolar_${GRID_TYPE}.csv" # Output grid file path
GRID_FILE="grids/grid_0p7ms_rerun_20260222.csv"

# """
# before uploading the scripts to hpc, 
# we run the generate_grid.py on the local machine to generate the grid file (grids/grid_x.csv) for each grid, 
# and then we can directly use the generated grid file to submit the jobs on hpc.
# This way we can avoid the time-consuming grid generation step on hpc and also have more control over the grid parameters and the generated grid file.

# To generate the grid file, you can run the following command on your local machine (make sure to adjust the parameters as needed):
# python generate_grid.py --output ${GRID_FILE} --grid ${GRID_TYPE} --mstar ${MSTAR} --model_dir ${MODEL_DIR}

# and past executions
# python generate_grid.py --output grids/grid_0p3ms_20260219.csv --grid extended --mstar 0.3 --model_dir /groups/pascucci/dingshandeng/diskmint_model_grids_0p3ms/radmc3d_model_grids/

# python generate_grid.py --output grids/grid_0p3ms_rerun_20260220.csv --grid supplementary --mstar 0.3 --model_dir /groups/pascucci/dingshandeng/diskmint_model_grids_0p3ms_rerun/radmc3d_model_grids/ 

# python generate_grid.py --output grids/grid_0p7ms_20260220.csv --grid extended --mstar 0.7 --model_dir /groups/pascucci/dingshandeng/diskmint_model_grids_0p7ms/radmc3d_model_grids/ 

# python generate_grid.py --output grids/grid_0p7ms_rerun_20260222.csv --grid supplementary --mstar 0.7 --model_dir /groups/pascucci/dingshandeng/diskmint_model_grids_0p7ms_rerun/radmc3d_model_grids/

# """

# Create logs directory
mkdir -p logs

# Function to submit a grid
submit_grid() {
    local MSTAR=$1
    local GRID_TYPE=$2
    local MODEL_MARK=$3
    local MODEL_PARAM_FILE=$4

    echo "=================================================="
    echo "Setting up grid for Mstar = ${MSTAR} Msun"
    echo "Grid type: ${GRID_TYPE}"
    echo "Model mark: ${MODEL_MARK}"
    echo "Model parameter file: ${MODEL_PARAM_FILE}"
    echo "=================================================="
    
    # NOTE: Now the parameter grid is generated before uploading the scripts to hpc
    # # Generate parameter grid
    # GRID_FILE="grids/grid_${MSTAR}msolar_${GRID_TYPE}.csv"
    # mkdir -p grids
    # python generate_grid.py --output ${GRID_FILE} --grid ${GRID_TYPE} --mstar ${MSTAR} --model_dir ${MODEL_DIR}

    # sbatch --wait run_python.slurm generate_grid.py --output ${GRID_FILE} --grid ${GRID_TYPE} --mstar ${MSTAR} --model_dir ${MODEL_DIR}
    
    # Count number of models
    N_MODELS=$(($(wc -l < ${GRID_FILE}) - 1))  # Subtract header

    # echo "Generated grid file inside: ${GRID_FILE}"

    echo "Submitting ${N_MODELS} models for Mstar=${MSTAR}"
    
    # Submit array job
    JOB_ID=$(sbatch --parsable \
        --array=0-$((N_MODELS-1)) \
        --job-name=diskmint_${MSTAR}ms \
        run_single_model.slurm \
        ${GRID_FILE} ${MSTAR} ${MODEL_DATE} ${MODEL_MARK} ${MODEL_PARAM_FILE})
    
    echo "Submitted job array: ${JOB_ID}"
    echo "Array size: ${N_MODELS}"
    echo ""
    
    # Save job info
    echo "${MSTAR},${GRID_TYPE},${MODEL_MARK},${JOB_ID},${N_MODELS}" >> job_submissions.log
}

# Initialize log file
echo "Mstar,GridType,ModelMark,JobID,NumModels" > job_submissions.log

# =============================================================================
# Submit different grids - CUSTOMIZE THIS SECTION
# =============================================================================

# USAGE:
submit_grid $MSTAR $GRID_TYPE $MODEL_MARK $MODEL_PARAM_FILE

# Example 1: 2.0 Msolar with BT-Settl model (your current setup)
# submit_grid 2.0 "default" "hpc_grid2p0msolar" "model_parameters_2p0Msolar.csv"

# Example 2: 1.0 Msolar grid
# submit_grid 1.0 "default" "hpc_grid1p0msolar" "model_parameters_1p0Msolar.csv"

# Example 3: 0.5 Msolar grid
# submit_grid 0.5 "extended" "hpc_grid0p5msolar" "model_parameters_0p5Msolar.csv"

# Example 4: Test grid for quick verification
# submit_grid 2.0 "test" "test_grid2p0msolar" "model_parameters_2p0Msolar.csv"

# =============================================================================

echo "=================================================="
echo "All jobs submitted!"
echo "Check status with: squeue -u $USER"
echo "Monitor logs in: ${WORK_DIR}/logs/"
echo "Job summary saved to: job_submissions.log"
echo "=================================================="
