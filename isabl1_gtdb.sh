#!/bin/bash

#SBATCH --output=log/enterococcus_diversity/isabl1_gtdb.out
#SBATCH --error=log/enterococcus_diversity/isabl1_gtdb.err
#SBATCH --time=168:00:00  # Adjust the time as needed
#SBATCH --cpus-per-task=64
#SBATCH --mem=256G
#SBATCH --partition=cpu

# Activate the GTDB environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate gtdb

# Set paths
BIN_DIR="analyses/enterococcus_diversity/binning/short/isabl1/vamb/multi_bins"
OUTPUT_DIR="analyses/enterococcus_diversity/binning/short/isabl1/gtdb"
MASH_DB="/data1/xavierj/carlos/dbs/gtdbtk-2.4.1/db/gtdb_ref_sketch.msh"

# Run the download-db.sh script
gtdbtk classify_wf --genome_dir ${BIN_DIR} \
    --out_dir ${OUTPUT_DIR} \
    --cpus 64 \
    --mash_db ${MASH_DB}