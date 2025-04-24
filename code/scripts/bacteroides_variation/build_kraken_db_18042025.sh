#!/bin/bash

#SBATCH --mem-per-cpu=8G
#SBATCH -c 8
#SBATCH --error=log/bacteroides_variation/build_kraken_db.err
#SBATCH --output=log/bacteroides_variation/build_kraken_db.out
#SBATCH --partition=cpu
#SBATCH --time=96:00:00

# Load dependencies
source $(conda info --base)/etc/profile.d/conda.sh
conda activate kraken2

# Set paths
DBNAME="/data1/xavierj/carlos/dbs/kraken/18042025"
CUSTOM_GENOMES="analyses/bacteroides_variation/genomes/ncbi/ncbi_dataset/kraken_header_fastas/*.fasta"

# Download Kraken databases
k2 download-taxonomy --db ${DBNAME}
k2 download-library --library bacteria --db ${DBNAME} --threads 8
k2 download-library --library archaea --db ${DBNAME} --threads 8
k2 download-library --library plasmid --db ${DBNAME} --threads 8
k2 download-library --library viral --db ${DBNAME} --threads 8
k2 download-library --library fungi --db ${DBNAME} --threads 8
k2 download-library --library UniVec --db ${DBNAME} --threads 8

# Adding Bacteroides genomes
k2 add-to-library --file ${CUSTOM_GENOMES} --db ${DBNAME} --threads 8

# Make a copy of the db library for future updates
mkdir -p /data1/xavierj/carlos/dbs/kraken/base
cp -r ${DBNAME}/* /data1/xavierj/carlos/dbs/kraken/base

# Build the database
k2 build --db ${DBNAME} --minimizer-spaces 4 --kmer-len 49 --threads 8