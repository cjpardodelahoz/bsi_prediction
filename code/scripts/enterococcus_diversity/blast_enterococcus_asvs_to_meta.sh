#!/bin/bash

#SBATCH --mem-per-cpu=8G
#SBATCH -c 4
#SBATCH --time=12:00:00
#SBATCH --error=log/enterococcus_diversity/blast_asvs_to_meta.err
#SBATCH --output=log/enterococcus_diversity/blast_asvs_to_meta.out
#SBATCH --partition=cpu

# Activate conda
source $(conda info --base)/etc/profile.d/conda.sh
conda activate blast

# Set paths
QUERY_FASTA="analyses/enterococcus_diversity/16S/compiled/enterococcus_asvs.fasta"
SUBJECT_FASTA="analyses/enterococcus_diversity/16S/compiled/enterococcus_16S.fasta"
OUTPUT_TABLE="analyses/enterococcus_diversity/16S/compiled/enterococcus_asv_to_meta_blast.txt"
BLAST_DB_DIR="analyses/enterococcus_diversity/16S/compiled/meta_blast_db"

# Create a BLAST database from the subject FASTA
mkdir -p ${BLAST_DB_DIR}
makeblastdb -in ${SUBJECT_FASTA} -dbtype nucl -parse_seqids -out ${BLAST_DB_DIR}/enterococcus_meta_db

# Run BLASTn search
blastn -query ${QUERY_FASTA} \
       -db ${BLAST_DB_DIR}/enterococcus_meta_db \
       -out ${OUTPUT_TABLE} \
       -outfmt "6 qseqid sseqid pident qlen slen length mismatch gaps" \
       -num_threads 4