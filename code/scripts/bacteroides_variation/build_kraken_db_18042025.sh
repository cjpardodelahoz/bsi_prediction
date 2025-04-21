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
#mkdir -p /data1/xavierj/carlos/dbs/kraken/18042025
DBNAME="/data1/xavierj/carlos/dbs/kraken/18042025"
kraken2-build --download-taxonomy --db ${DBNAME}
kraken2-build --download-library bacteria --db ${DBNAME} --threads 16
kraken2-build --download-library archaea --db ${DBNAME} --threads 16
kraken2-build --download-library plasmid --db ${DBNAME} --threads 16
kraken2-build --download-library viral --db ${DBNAME} --threads 16
kraken2-build --download-library fungi --db ${DBNAME} --threads 16

# Adding Bacteroides genomes
ls analyses/bacteroides_variation/genomes/ncbi/ncbi_dataset/kraken_header_fastas/*.fasta | while read line; do
	file=$(echo ${line});
	kraken2-build --add-to-library ${file} --db ${DBNAME} --threads 16;
done

# Make a copy of the db library for future updates
mkdir -p /data1/xavierj/carlos/dbs/kraken/base
cp -r ${DBNAME}/* /data1/xavierj/carlos/dbs/kraken/base

# Build the database
kraken2-build --build --db ${DBNAME} --minimizer-spaces 4 --kmer-len 49 --threads 16