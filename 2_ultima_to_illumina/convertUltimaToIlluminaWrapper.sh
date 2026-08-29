# Needs to take the S3 location, output save location, and work directory as input arguments

# TODO: Need to have a barcode conversion file!! i.e Ultima to Scale. Probably a .csv file.
# First column for the Ultima barcodes, 2nd column for the corresponding Scale barcodes.
# Probably just keep this in the github repo

# Then download the file from S3

# after downloading, then run convertUltimaToIllumina
# will need the path for the downloaded INPUT_FASTQ file, the work directory
# will need to get the PCR barcode from the conversion file.

# after running, then upload the fastq files back to the S3 bucket

# Then clean up the fastq files on the SSD

#!/usr/bin/env bash
set -uo pipefail
# Note: intentionally NOT using `set -e` — we want one failed sample to be
# logged and skipped, not kill the whole batch run.

# =========================================================
# Example to run wrapper
# =========================================================
# [script] [S3 input folder] [S3 output folder] [EC2 staging dir] [optional barcode CSV]
./convertUltimaToIlluminaWrapper.sh \
    s3://chesilab-testbucket/ultimatestout \
    s3://chesilab-testbucket/scalein \
    /home/ec2-user/scratch \
    ./PCR_Barcodes.csv

# =========================================================
# User config / args
# =========================================================
BUCKET1=$1                                   # S3 folder with source files, one _rna.fastq.gz per subfolder
BUCKET2=$2                                   # S3 folder to receive Illumina fastq.gz outputs
WORKDIR=$3                                   # local SSD scratch space
# Barcode conversion CSV: col1 = Ultima BC, col2 = Scale BC. Defaults to repo copy.
BARCODE_CSV="${4:-$(dirname "$(realpath "${BASH_SOURCE[0]}")")/barcodeMap.csv}"

# Path to the per-file conversion script (assumed to sit next to this wrapper)
SCRIPT="$(dirname "$(realpath "${BASH_SOURCE[0]}")")/convertUltimaToIllumina.sh"

# mkdir -p "$WORKDIR"
LOGFILE="$WORKDIR/convert_log_$(date +%Y%m%d_%H%M%S).txt"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOGFILE"
}

# =========================================================
# Pre-flight checks
# =========================================================
if [ ! -f "$SCRIPT" ]; then
    log "ERROR: Conversion script not found: $SCRIPT"
    exit 1
fi
if [ ! -f "$BARCODE_CSV" ]; then
    log "ERROR: Barcode conversion CSV not found: $BARCODE_CSV"
    exit 1
fi

log "=== Starting batch run ==="
log "Bucket1 (input):  $BUCKET1"
log "Bucket2 (output): $BUCKET2"
log "Workdir:          $WORKDIR"
log "Barcode CSV:      $BARCODE_CSV"
log "Conversion script: $SCRIPT"

# Helper: look up Scale barcode for a given Ultima barcode.
# Matches by column *header name*, so column order/count can change safely.
lookup_scale_bc() {
    local ultima="$1"
    awk -F',' -v u="$ultima" '
        NR==1 {
            for (i=1; i<=NF; i++) {
                h=$i; gsub(/\r/,"",h); gsub(/^[[:space:]]+|[[:space:]]+$/,"",h)
                if (tolower(h)=="ultima_bc_name") seq_col=i
                if (tolower(h)=="scale_bc_index2") idx_col=i
            }
            if (!seq_col || !idx_col) {
                print "LOOKUP_ERROR: missing required columns" > "/dev/stderr"
                exit 2
            }
            next
        }
        {
            s=$seq_col; gsub(/[[:space:]\r]/,"",s)
            v=$idx_col; gsub(/[[:space:]\r]/,"",v)
            if (tolower(s)==tolower(u)) { print v; exit }
        }
    ' "$BARCODE_CSV"
}

# =========================================================
# Get list of top-level "folder" names in Bucket1
# =========================================================
folders=$(aws s3 ls "${BUCKET1}/" | awk '/PRE/ {print $2}' | sed 's#/$##')

if [ -z "$folders" ]; then
    log "ERROR: No folders found in $BUCKET1 — check the bucket path/permissions."
    exit 1
fi

total=$(echo "$folders" | wc -l)
count=0

for folder in $folders; do
    count=$((count + 1))
    log "--- [$count/$total] Processing folder: $folder ---"

    out_prefix_local="${WORKDIR}/${folder}"     # local output dir for this sample
    out_s3_prefix="${BUCKET2}/${folder}"

    # ---- 1. Find the input _rna.fastq.gz in this folder ----
    # (filename is NOT the folder name, so we list and match the pattern)
    fastq_name=$(aws s3 ls "${BUCKET1}/${folder}/" 2>/dev/null \
                 | awk '{print $4}' | grep -E '_rna\.fastq\.gz$' | head -n1 || true)

    if [ -z "$fastq_name" ]; then
        log "WARNING: No *_rna.fastq.gz found in ${BUCKET1}/${folder}/ — skipping."
        continue
    fi

    fastq_key="${BUCKET1}/${folder}/${fastq_name}"
    local_fastq="${out_prefix_local}/${fastq_name}"

 # ---- 2. Extract the full Ultima barcode NAME (Z-code) from the filename ----
    # e.g. 430978-ChesiLab1-Z0001-CAGCTCGAATGCGAT_rna.fastq.gz -> Z0001
    #   grep -oP 'Z[0-9]+'  keeps the leading 'Z' (no \K), so it matches the CSV.
    ultima_bc=$(echo "$fastq_name" | grep -oP 'Z[0-9]+' | sed -n '1p')
    if [ -z "$ultima_bc" ]; then
        log "WARNING: Could not parse Z-code from '$fastq_name' — skipping."
        continue
    fi

    # ---- 3. Map Ultima -> Scale barcode via the conversion CSV ----
    scale_bc=$(lookup_scale_bc "$ultima_bc")
    if [ -z "$scale_bc" ]; then
        log "WARNING: No Scale barcode mapping for Ultima BC '$ultima_bc' (folder $folder) — skipping."
        continue
    fi
    log "Ultima BC $ultima_bc -> Scale BC $scale_bc"

    # ---- 4. Derive numeric sample number from the same Z-code (no second grep) ----
    # Z0001 -> 0001 -> force base-10 -> zero-pad to 3 -> 001
    raw_num="${ultima_bc#Z}"                       # strip leading 'Z' -> 0001
    if ! [[ "$raw_num" =~ ^[0-9]+$ ]]; then
        log "WARNING: Z-code '$ultima_bc' has no numeric part — skipping."
        continue
    fi
    sample_num=$(printf '%03d' $((10#$raw_num)))   # 10# avoids octal on leading zeros

    # ---- 5. Download the input FASTQ to local SSD ----
    mkdir -p "$out_prefix_local"
    log "Copying $fastq_key -> $local_fastq"
    if ! aws s3 cp "$fastq_key" "$local_fastq" --only-show-errors; then
        log "ERROR: Failed to copy $fastq_key — skipping."
        rm -rf "$out_prefix_local"
        continue
    fi

    # ---- 6. Run the conversion script ----
    log "Running conversion on $local_fastq (sample $sample_num)"
    if ! bash "$SCRIPT" "$local_fastq" "$out_prefix_local" "$scale_bc" >> "$LOGFILE" 2>&1; then
        log "ERROR: Conversion failed for $folder — leaving local files for inspection, skipping upload."
        continue
    fi

    # ---- 7. Upload the four Illumina FASTQs to Bucket2 ----
    upload_ok=true
    for read in R1 R2 I1 I2; do
        out_file="${out_prefix_local}/ScaleRNA_${read}_${sample_num}.fastq.gz"
        if [ -f "$out_file" ]; then
            aws s3 cp "$out_file" "${out_s3_prefix}/ScaleRNA_${read}_${sample_num}.fastq.gz" \
                --only-show-errors || upload_ok=false
        else
            log "WARNING: Expected output missing: $out_file"
            upload_ok=false
        fi
    done

    if [ "$upload_ok" = false ]; then
        log "ERROR: Upload incomplete for $folder — local files kept for inspection."
        continue
    fi
    log "Upload confirmed for $folder."

    # ---- 8. Clean up local SSD ----
    rm -rf "$out_prefix_local"
    log "Cleaned up local files for $folder."

done

log "=== Batch run complete: $count/$total folders processed ==="