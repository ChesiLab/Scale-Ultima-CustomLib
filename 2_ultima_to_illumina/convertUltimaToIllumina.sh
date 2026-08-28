#!/usr/bin/env bash
set -euo pipefail # could add back 'x' for debugging

# =========================================================
# Usage
# =========================================================

# ./convertUltimaToIllumina.sh <INPUT_FASTQ> <WORKDIR> [SAMPLE_BC]

# Example:
# ./convertUltimaToIllumina.sh \
#     /home/sjpl/github/scratchData/in/Z0001_02p/Z0001_02p_rna.fastq.gz \
#     /home/sjpl/github/scratchData/work \
#     "GCATCGTATG"

# Produces in WORKDIR:
#   ScaleRNA_R1_001.fastq.gz
#   ScaleRNA_R2_001.fastq.gz
#   ScaleRNA_I1_001.fastq.gz
#   ScaleRNA_I2_001.fastq.gz

# =========================================================
# Arg Definition
# =========================================================

INPUT_FASTQ="$1"             # Input fastq that came from cram
WORKDIR="$2"                 # local SSD scratch for intermediates + finals
BC="$3"                      # Scale barcode for I2 file
# BC="GCATCGTATG"            # Stand-in Scale barcode

# Input check:
if [ ! -f "$INPUT_FASTQ" ]; then
    echo "ERROR: Input FASTQ not found: $INPUT_FASTQ" >&2
    exit 1
fi

# Create string for the I2 Quality Score. Will use D because that's the most common.
QUAL=$(printf 'D%.0s' $(seq ${#BC}))

# ---- Intermediate R1 and R2 fastqs (stay on local SSD) ----
OUTPUT_R1_LONG="${WORKDIR}/R1_LONG.fastq.gz"
OUTPUT_R2_LONG="${WORKDIR}/R2_LONG.fastq.gz"
OUTPUT_R2_REV="${WORKDIR}/R2_REV.fastq.gz"

# ---- Final outputs (written locally, then uploaded to S3) ----
OUTPUT_I1="${WORKDIR}/ScaleRNA_I1_001.fastq.gz"
OUTPUT_I2="${WORKDIR}/ScaleRNA_I2_001.fastq.gz"
OUTPUT_R1="${WORKDIR}/ScaleRNA_R1_001.fastq.gz"
OUTPUT_R2="${WORKDIR}/ScaleRNA_R2_001.fastq.gz"

# >>>>>> TODO: the 001 etc needs to match the Z0001 format!! <<<<<

# =========================================================
# R1 & R2: trim adapters and filter for pairs reaching min lengths (i.e. really just R2 >= 76bp)
# looking for the partial TruSeq adapter and the Nextera adapter
# =========================================================

args1=(
    # -j 0 # use all available CPU cores
    --discard-untrimmed # remove any read w/o an adapter
    --pair-filter any # using cram twice, treated as a "pair". Remove read if either in "pair" fail adapter check
    # Define Read 1 5'+3' adapters; define error_rate and min_overlap for each adapter; 5' adapter is required.
    -a "CTACACGACGCTCTTCCGATCT;max_error_rate=0.2;min_overlap=10;required...CTGTCTCTTATACACATCTC;max_error_rate=0.2;min_overlap=6"
    # Define Read 2 5'+3' "adapters"; define error_rate and min_overlap for each adapter; 5' adapter is required.
    -A "CTACACGACGCTCTTCCGATCT;max_error_rate=0.2;min_overlap=10;required...CTGTCTCTTATACACATCTC;max_error_rate=0.2;min_overlap=6"
    -o "${OUTPUT_R1_LONG}" # Output FASTQ file for R1
    -p "${OUTPUT_R2_LONG}" # Output FASTQ file for R2
    --minimum-length 34:76 # Discard pairs where R1 < 34 *or* R2 < 76
    "${INPUT_FASTQ}" "${INPUT_FASTQ}" # Ultima input corresponding to R1 and R2 from Illumina
)

# R1 & R2: Run trimming and filtering
cutadapt "${args1[@]}"

# =========================================================
# R1: trim to standard Illumina R1 length: 34bp for 34 cycles
# =========================================================

cutadapt \
    --minimum-length 34 \
    --maximum-length 34 \
    --length 34 \
    -o "${OUTPUT_R1}" \
    "${OUTPUT_R1_LONG}"

# =========================================================
# R2: get reverse complement, then trim to 76bp for 76 cycles
# =========================================================

# Get reverse complement
zcat "${OUTPUT_R2_LONG}" | awk '{print $0 }' | seqkit seq -p -r -t DNA | gzip > "${OUTPUT_R2_REV}"

# Trim to standard Illumina R2 length: 76bp for 76 cycles
cutadapt \
    --minimum-length 76 \
    --maximum-length 76 \
    --length 76 \
    -o "${OUTPUT_R2}" \
    "${OUTPUT_R2_REV}"

# =========================================================
# I2: Create fastq file for I2 that matches headers of R1 (and R2)
# =========================================================

zcat "${OUTPUT_R1}" | awk 'NR%4==1' | \
awk -v bc="$BC" -v qual="$QUAL" '{print $0 "\n" bc "\n+\n" qual}' | \
gzip > "${OUTPUT_I2}"

# =========================================================
# I1: Create fastq for I1 that matches the expected random i7 indices
# =========================================================

zcat "${OUTPUT_R1}" | \
awk '
    BEGIN {
        srand(42)                     
        bc[0]="ATCTGCAGTC"
        bc[1]="GCTCTCGCCT"
        bc[2]="TGAGCATAGA"
        bc[3]="CAGAAGCTAG"
        qual="DDDDDDDDDD"             
    }
    NR%4==1 {print}                   # header from R1, unchanged
    NR%4==2 {print bc[int(rand()*4)]} # random barcode as sequence
    NR%4==3 {print "+"}               # plus line
    NR%4==0 {print qual}              # constant quality
' | \
gzip > "${OUTPUT_I1}"

# =========================================================
# Clean up intermediates only (finals left for wrapper to upload)
# =========================================================
rm -f "${OUTPUT_R1_LONG}" "${OUTPUT_R2_LONG}" "${OUTPUT_R2_REV}"

echo "Done: Illumina FASTQs from ${INPUT_FASTQ} in ${WORKDIR}"

