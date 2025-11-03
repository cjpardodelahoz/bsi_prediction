#!/bin/bash

#SBATCH --array=1-565%50
#SBATCH --output=log/enterococcus_diversity/vamb_multi_isabl1_%A_%a.out
#SBATCH --error=log/enterococcus_diversity/vamb_multi_isabl1_%A_%a.err
#SBATCH --time=120:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=48G
#SBATCH --partition=cpu

# Activate conda
source $(conda info --base)/etc/profile.d/conda.sh

# Define paths
ASSEMBLY_DIR="analyses/enterococcus_diversity/metagenomes/assembly/metaspades"
BINNING_DIR="analyses/enterococcus_diversity/binning/short/isabl1"
PATIENT_SAMPLES_DIR="${BINNING_DIR}/patient_samples"
VAMB_OUT_DIR="${BINNING_DIR}/vamb/multi"
READS_PATHS="misc_files/enterococcus_diversity/all_target_samples_paths.txt"

# Get the patient ID for the current task
patient_file=$(ls ${PATIENT_SAMPLES_DIR} | sed -n ${SLURM_ARRAY_TASK_ID}p)
patient_id=${patient_file%%.txt}

# Create output directories
mkdir -p ${VAMB_OUT_DIR}/${patient_id}
AEMB_DIR="${VAMB_OUT_DIR}/${patient_id}/aemb"
mkdir -p ${AEMB_DIR}

# Step 1: Filter contigs >= 1000 bp for all samples in the patient
conda activate hmmerseqkit
PATIENT_SAMPLES=$(cat ${PATIENT_SAMPLES_DIR}/${patient_file})
for sample in ${PATIENT_SAMPLES}; do
    INPUT_FILE="${ASSEMBLY_DIR}/${sample}/contigs.fasta"
    OUTPUT_FILE="${ASSEMBLY_DIR}/${sample}/contigs_1000.fasta"
    if [[ -f "${INPUT_FILE}" ]]; then
        seqkit seq -m 1000 ${INPUT_FILE} -o ${OUTPUT_FILE}
    fi
done
conda deactivate

# Step 2: Concatenate assemblies (contigs >1000 bp) for all samples in the patient
CONCATENATED_FASTA="${VAMB_OUT_DIR}/${patient_id}/contigs.fna.gz"
software/custom/concatenate_assemblies_for_vamb.sh \
    ${CONCATENATED_FASTA} $(for sample in ${PATIENT_SAMPLES}; do echo ${ASSEMBLY_DIR}/${sample}/contigs_1000.fasta; done)

# Step 3: Map reads from all samples to the concatenated assembly using strobealign
conda activate strobealign
for sample in ${PATIENT_SAMPLES}; do
    # Get R1 and R2 paths from the paths table (first record for each sample)
    record=$(awk -v s="['${sample}']" '$2==s{print; exit}' ${READS_PATHS})
    R1=$(echo "$record" | cut -f3 | sed 's|/data/brinkvd/|/data1/collab004/|')
    R2=$(echo "$record" | cut -f4 | sed 's|/data/brinkvd/|/data1/collab004/|')
    strobealign -t 16 --aemb ${CONCATENATED_FASTA} ${R1} ${R2} > ${AEMB_DIR}/${sample}.tsv
done
conda deactivate

# Step 4: Generate the abundance TSV file
conda activate python3
ABUNDANCE_TSV="${VAMB_OUT_DIR}/${patient_id}/abundance.tsv"
python software/vamb/merge_aemb.py ${AEMB_DIR} ${ABUNDANCE_TSV}
conda deactivate

# Step 4: Run VAMB
vamb bin default \
    --outdir ${VAMB_OUT_DIR}/${patient_id}/vambout \
    --fasta ${CONCATENATED_FASTA} \
    --abundance_tsv ${ABUNDANCE_TSV} \
    --minfasta 250000 \
    -m 1000 \
    -p 24 \
    -o _C_