#!/bin/bash

#SBATCH --output=log/gtdb_mash_db.out
#SBATCH --error=log/gtdb_mash_db.err
#SBATCH --time=24:00:00  # Adjust the time as needed
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --partition=cpu

# Activate the GTDB environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate gtdb

# Run the download-db.sh script
gtdbtk classify_wf --genome_dir test/ \
    --out_dir test_out \
    --cpus 2 \
    --mash_db /data1/xavierj/carlos/dbs/gtdbtk-2.4.1/db