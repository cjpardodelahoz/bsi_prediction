#!/bin/bash

#SBATCH --mem=16G
#SBATCH -c 4
#SBATCH --time=12:00:00
#SBATCH --error=log/bacteroides_variation/classify_bacteroides_16S.err
#SBATCH --output=log/bacteroides_variation/classify_bacteroides_16S.out
#SBATCH --partition=cpu

# Activate conda
source $(conda info --base)/etc/profile.d/conda.sh
conda activate kraken2

# Set paths
SILVA_DB="/data1/xavierj/carlos/dbs/kraken/silva"  # Path to the SILVA Kraken2 database
INPUT_FASTA="analyses/bacteroides_variation/16S/compiled/compiled_16S.fasta"
KRAKEN_OUTPUT="analyses/bacteroides_variation/16S/compiled/kraken_output.txt"
KRAKEN_REPORT="analyses/bacteroides_variation/16S/compiled/kraken_report.txt"
BACTEROIDES_FASTA="analyses/bacteroides_variation/16S/compiled/bacteroides_16S.fasta"

# Run Kraken2 classification
kraken2 --db ${SILVA_DB} \
        --threads 4 \
        --report ${KRAKEN_REPORT} \
        --output ${KRAKEN_OUTPUT} \
        ${INPUT_FASTA}

# Extract sequences classified as Bacteroidaceae
extract_kraken_reads.py -k ${KRAKEN_OUTPUT} \
                        -r ${KRAKEN_REPORT} \
                        -s ${INPUT_FASTA} \
                        -o ${BACTEROIDES_FASTA} \
                        --include-children --taxid 43873 # TaxID for Bacteroidaceae