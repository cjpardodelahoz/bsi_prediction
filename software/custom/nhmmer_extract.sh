#!/bin/bash

# Command-line application to search for sequences from a custom HMM profile and extract matching regions

# Function to display the help menu
function show_help() {
    echo "Usage: $0 -p <HMM_PROFILE> -e <EVALUE_THRESHOLD> -i <INPUT_FASTA> -o <OUTPUT_FASTA> [-t <OUTPUT_HITS>] [-c <CPUS>]"
    echo
    echo "Options:"
    echo "  -p, --profile         Path to the HMM profile file (e.g., 16S_rRNA.hmm)"
    echo "  -e, --evalue          E-value threshold for filtering hits (e.g., 0)"
    echo "  -i, --input           Input FASTA file to search"
    echo "  -o, --output          Output FASTA file to save extracted sequences"
    echo "  -t, --table           Path to save the HMMER output table (optional, default: temporary file)"
    echo "  -c, --cpus            Number of CPU cores to use for nhmmer (optional, default: 4)"
    echo "  -h, --help            Display this help menu"
    echo
    echo "This script requires 'nhmmer' and 'seqkit' to be installed and available in the PATH."
}

# Parse command-line arguments
CPUS=4  # Default number of CPU cores
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -p|--profile) HMM_PROFILE="$2"; shift ;;
        -e|--evalue) EVALUE_THRESHOLD="$2"; shift ;;
        -i|--input) INPUT_FASTA="$2"; shift ;;
        -o|--output) OUTPUT_FASTA="$2"; shift ;;
        -t|--table) OUTPUT_HITS="$2"; shift ;;
        -c|--cpus) CPUS="$2"; shift ;;
        -h|--help) show_help; exit 0 ;;
        *) echo "Unknown option: $1"; show_help; exit 1 ;;
    esac
    shift
done

# Check if all required arguments are provided
if [[ -z "$HMM_PROFILE" || -z "$EVALUE_THRESHOLD" || -z "$INPUT_FASTA" || -z "$OUTPUT_FASTA" ]]; then
    echo "Error: Missing required arguments."
    show_help
    exit 1
fi

# Check if required binaries are in the PATH
if ! command -v nhmmer &> /dev/null; then
    echo "Error: 'nhmmer' is not installed or not in the PATH."
    exit 1
fi

if ! command -v seqkit &> /dev/null; then
    echo "Error: 'seqkit' is not installed or not in the PATH."
    exit 1
fi

# Set default for OUTPUT_HITS if not provided
if [[ -z "$OUTPUT_HITS" ]]; then
    OUTPUT_HITS=$(mktemp)
fi

# Step 1: Search for sequences using nhmmer
nhmmer --tblout ${OUTPUT_HITS} --cpu ${CPUS} ${HMM_PROFILE} ${INPUT_FASTA}

# Step 2: Extract matching regions based on coordinates
grep -v "^#" ${OUTPUT_HITS} | awk -v evalue="$EVALUE_THRESHOLD" '$13 <= evalue {print $1, $9, $10, $12}' | while read -r target envfrom envto strand; do
    # Handle reverse complement if strand is '-'
    if [[ "$strand" == "-" ]]; then
        seqkit subseq --chr "$target" -r "$envto:$envfrom" ${INPUT_FASTA} | seqkit seq -r -p >> ${OUTPUT_FASTA}
    else
        seqkit subseq --chr "$target" -r "$envfrom:$envto" ${INPUT_FASTA} >> ${OUTPUT_FASTA}
    fi
done

# Check if any sequences were extracted
if [[ ! -s ${OUTPUT_FASTA} ]]; then
    echo "No sequences extracted from ${INPUT_FASTA} with E-value <= ${EVALUE_THRESHOLD}."
    [[ -z "$OUTPUT_HITS" ]] && rm -f ${OUTPUT_HITS}  # Remove temporary file if used
    exit 1
fi

# Clean up temporary files if a temporary file was used
[[ -z "$OUTPUT_HITS" ]] && rm -f ${OUTPUT_HITS}

echo "Extracted sequences saved to ${OUTPUT_FASTA}."