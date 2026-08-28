# Input fastq from cram
INPUT_FASTQ="/home/sjpl/scale_preprocess/cutadapt/cram/Z0001_02p_rna.fastq.gz"

# Setting sample barcode for I2 fastq file from this cram/fastq
# BC="CAGCTCGAATGCGAT" # Ultima barcode
BC="GCATCGTATG" # Stand-in Scale barcode
# Create string for the I2 Quality Score. Will use D because that's the most common.
QUAL=$(printf 'D%.0s' $(seq ${#BC}))

# Intermediate R1 and R2 fastqs
OUTPUT_R1_LONG="/home/sjpl/scale_preprocess/cutadapt/fastq/R1_LONG.fastq.gz"
OUTPUT_R2_LONG="/home/sjpl/scale_preprocess/cutadapt/fastq/R2_LONG.fastq.gz"
OUTPUT_R2_REV="/home/sjpl/scale_preprocess/cutadapt/fastq/R2_REV.fastq.gz"

# Outputs: R1, R2, and I2
OUTPUT_R1="/home/sjpl/scale_preprocess/cutadapt/fastq/ScaleRNA_R1_001.fastq.gz"
OUTPUT_R2="/home/sjpl/scale_preprocess/cutadapt/fastq/ScaleRNA_R2_001.fastq.gz"

OUTPUT_I1="/home/sjpl/scale_preprocess/cutadapt/fastq/ScaleRNA_I1_001.fastq.gz"
OUTPUT_I2="/home/sjpl/scale_preprocess/cutadapt/fastq/ScaleRNA_I2_001.fastq.gz"

# R1 & R2: Set up to trim adapters and filter for R1 and R2 pairs that reach minimum lengths. (so really R2>=76bp)
# looking for the partial TruSeq adapter and the Nextera adapter
args1=(
    # -j 0 # use all available CPU cores
    --discard-untrimmed # remove any read w/o an adapter
    --pair-filter any # using cram twice, treated as a "pair". Remove read if either in "pair" fail adapter check
    # Define Read 1 5'+3' adapters; define error_rate and min_overlap for each adapter; 5' adapter is required.
    -a "CTACACGACGCTCTTCCGATCT;max_error_rate=0.2;min_overlap=10;required...CTGTCTCTTATACACATCTC;max_error_rate=0.2;min_overlap=6"
    # -U 50 # Unconditionally remove the first 50bp from R2
    # -q 25 # Quality-trim 30bp on R1 from the 3' end, using Phred score cutoff of 30
    # Define Read 2 5'+3' "adapters"; define error_rate and min_overlap for each adapter; 5' adapter is required.
    -A "CTACACGACGCTCTTCCGATCT;max_error_rate=0.2;min_overlap=10;required...CTGTCTCTTATACACATCTC;max_error_rate=0.2;min_overlap=6"
    -o "${OUTPUT_R1_LONG}" # Output FASTQ file for R1
    -p "${OUTPUT_R2_LONG}" # Output FASTQ file for R2
    --minimum-length 34:76 # Discard pairs where R1 < 34 *or* R2 < 76
    "${INPUT_FASTQ}" "${INPUT_FASTQ}" # Ultima input corresponding to R1 and R2 from Illumina
)

# R1 & R2: Run trimming and filtering
cutadapt "${args1[@]}"

# R1: Trimming R1 down to standard Illumina R1 length: 34bp for 34 cycles
cutadapt \
    --minimum-length 34 \
    --maximum-length 34 \
    --length 34 \
    -o "${OUTPUT_R1}" \
    "${OUTPUT_R1_LONG}"

# R2: Getting reverse complement of R2
zcat "${OUTPUT_R2_LONG}" | awk '{print $0 }' | seqkit seq -p -r -t DNA | gzip > "${OUTPUT_R2_REV}"

# R2: Trimming R2 to standard Illumina R2 length: 76bp for 76 cycles
cutadapt \
    --minimum-length 76 \
    --maximum-length 76 \
    --length 76 \
    -o "${OUTPUT_R2}" \
    "${OUTPUT_R2_REV}"


# I2: Create fastq file for I2 that matches headers of R1 (and R2)
zcat "${OUTPUT_R1}" | awk 'NR%4==1' | \
awk -v bc="$BC" -v qual="$QUAL" '{print $0 "\n" bc "\n+\n" qual}' | \
gzip > "${OUTPUT_I2}"

# I1
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



