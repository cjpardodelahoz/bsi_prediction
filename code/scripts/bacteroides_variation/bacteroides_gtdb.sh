#!/bin/bash

#SBATCH --output=log/bacteroides_variation/bacteroides_gtdb.out
#SBATCH --error=log/bacteroides_variation/bacteroides_gtdb.err
#SBATCH --time=96:00:00  # Adjust the time as needed
#SBATCH --cpus-per-task=12
#SBATCH --mem=128G
#SBATCH --partition=cpu

# Activate the GTDB environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate gtdb

# Set paths
BIN_DIR="analyses/bacteroides_variation/binning/compiled"
OUTPUT_DIR="analyses/bacteroides_variation/binning/gtdb"
MASH_DB="/data1/xavierj/carlos/dbs/gtdbtk-2.4.1/db/gtdb_ref_sketch.msh"

# Run the download-db.sh script
gtdbtk classify_wf --genome_dir ${BIN_DIR} \
    --out_dir ${OUTPUT_DIR} \
    --cpus 12 \
    --mash_db ${MASH_DB}