#!/bin/bash

#SBATCH --array=1#1-9999%500
#SBATCH --mem=16G         # Memory per task
#SBATCH -c 4              # Number of CPU cores per task
#SBATCH --time=48:00:00   # Maximum runtime
#SBATCH --error=log/bacteroides_pul/isabl1_mags_signalp_%A_%a.err
#SBATCH --output=log/bacteroides_pul/isabl1_mags_signalp_%A_%a.out
#SBATCH --partition=cpu

# Activate conda
source $(conda info --base)/etc/profile.d/conda.sh
conda activate signalp

# Get the MAG name for the current task
MAG=$(cat analyses/bacteroides_pul/isabl1_nc_mags.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)

# Set paths
PUL_DIR="analyses/bacteroides_pul/mag_pul_prediction/isabl1/${MAG}"
PROTEINS="${PUL_DIR}/CGC.faa"
OUT_DIR="${PUL_DIR}/signalp"

# Run SignalP
signalp6 --fastafile ${PROTEINS} \
    --organism other \
    --output_dir ${OUT_DIR} \
    --format txt \
    --mode fast \
    --torch_num_threads 4