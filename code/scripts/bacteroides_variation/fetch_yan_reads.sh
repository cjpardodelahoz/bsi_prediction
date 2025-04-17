#!/bin/bash

#SBATCH --array=1-395%20
#SBATCH --mem-per-cpu=16G
#SBATCH -c 1
#SBATCH --error=log/bacteroides_variation/fetch_yan_reads_%A_%a.err
#SBATCH --output=log/bacteroides_variation/fetch_yan_reads_%A_%a.err
#SBATCH --partition=cpu
#SBATCH --time=06:00:00


# Set path for SRA toolkits
export PATH=/usersoftware/xavierj/sratoolkit.3.2.0-centos_linux64/bin:$PATH

# Variable with SRA accession
sra_accession=$(sed -n ${SLURM_ARRAY_TASK_ID}p misc_files/bacteroides_variation/yan_sd_2022_accessions.txt | cut -f 1)
# Variable with sample name
sample_name=$(sed -n ${SLURM_ARRAY_TASK_ID}p misc_files/bacteroides_variation/yan_sd_2022_accessions.txt | cut -f 2)

# Download raw reads from SRA
prefetch ${sra_accession}
fasterq-dump -O data/reads/yan_sd_2022 ${sra_accession}

# Rename and sort reads
mkdir -p data/reads/yan_sd_2022/${sampqle_name}
mv data/reads/yan_sd_2022/${sra_accession}_1.fastq data/reads/yan_sd_2022/${sample_name}/${sample_name}_R1_all.fastq
mv data/reads/yan_sd_2022/${sra_accession}_2.fastq data/reads/yan_sd_2022/${sample_name}/${sample_name}_R2_all.fastq
# Gzip the files
gzip data/reads/yan_sd_2022/${sample_name}/${sample_name}_R1_all.fastq
gzip data/reads/yan_sd_2022/${sample_name}/${sample_name}_R2_all.fastq