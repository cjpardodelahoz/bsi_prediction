#!/bin/bash

#SBATCH --mem-per-cpu=4G
#SBATCH -c 16
#SBATCH --error=log/bacteroides_variation/build_kraken_db.err
#SBATCH --output=log/bacteroides_variation/build_kraken_db.err
#SBATCH --partition=cpu
#SBATCH --time=96:00:00

# Load dependencies
source $(conda info --base)/etc/profile.d/conda.sh
conda activate kraken2

# Kraken databases
mkdir -p /data1/xavierj/carlos/dbs/kraken/18042025
DBNAME="/data1/xavierj/carlos/dbs/kraken/18042025"
kraken2-build --download-taxonomy --db ${DBNAME}
kraken2-build --download-library bacteria --db ${DBNAME} --threads 16
kraken2-build --download-library archaea --db ${DBNAME} --threads 16
kraken2-build --download-library plasmid --db ${DBNAME} --threads 16
kraken2-build --download-library viral --db ${DBNAME} --threads 16
kraken2-build --download-library fungi --db ${DBNAME} --threads 16