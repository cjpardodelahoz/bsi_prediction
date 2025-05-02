#!/bin/bash

#SBATCH --mem-per-cpu=8G
#SBATCH -c 8
#SBATCH --error=log/bacteroides_variation/build_silva_db.err
#SBATCH --output=log/bacteroides_variation/build_silva_db.out
#SBATCH --partition=cpu
#SBATCH --time=96:00:00

# Load dependencies
source $(conda info --base)/etc/profile.d/conda.sh
conda activate kraken2

# Set paths
DBPATH="/data1/xavierj/carlos/dbs/kraken/silva"

# Build silva database
k2 build --db $DBPATH \
    --threads 8 \
    --special silva