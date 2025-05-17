#!/bin/bash

#SBATCH --array=1-307%20  # Process 307 MAGs, with up to 20 jobs running concurrently
#SBATCH --mem=16G         # Memory per task
#SBATCH -c 4              # Number of CPU cores per task
#SBATCH --time=24:00:00   # Maximum runtime
#SBATCH --error=log/bacteroides_variation/bacteroides_dbcan_%A_%a.err
#SBATCH --output=log/bacteroides_variation/bacteroides_dbcan_%A_%a.out
#SBATCH --partition=cpu

# Activate conda
source $(conda info --base)/etc/profile.d/conda.sh
conda activate dbcan

# Set paths
MAG_DIR="analyses/bacteroides_variation/genomes/mags"
OUTPUT_DIR="analyses/bacteroides_variation/cgc/dbcan"
DB_DIR="/data1/xavierj/carlos/dbs/dbcan"

# Get the MAG file for the current task
MAG=$(ls ${MAG_DIR} | sed -n ${SLURM_ARRAY_TASK_ID}p)
MAG_PATH="${MAG_DIR}/${MAG}"

# Create output directory for the current MAG
MAG_OUTPUT_DIR="${OUTPUT_DIR}/${MAG%%.fna}"
mkdir -p ${MAG_OUTPUT_DIR}

# Run dbCAN
run_dbcan easy_substrate \
  --input_raw_data ${MAG_PATH} \
  --mode prok \
  --output_dir ${MAG_OUTPUT_DIR} \
  --db_dir ${DB_DIR} \
  --input_gff ${MAG_OUTPUT_DIR}/uniInput.gff \
  --gff_type prodigal