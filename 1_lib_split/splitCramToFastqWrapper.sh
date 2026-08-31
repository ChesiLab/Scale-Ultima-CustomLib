#!/usr/bin/env bash
set -uo pipefail
# Note: intentionally NOT using `set -e` — we want one failed CRAM to be
# logged and skipped, not kill the whole batch run.

# =========================================================
# Example to run wrapper
# =========================================================

# conda activate libsplit # or whatever env you're using, i.e. ngstools
# [script] [AWS S3 folder with source files] [AWS S3 folder to put output] [EC2 folder for staging intermediate files]
# nohup ./splitCramToFastqWrapper.sh s3://upenn-chesi-lab-broad-data-transfer s3://chesilab-testbucket/ultima_fastq /home/ec2-user/scratch &

# =========================================================
# Troubleshooting: Folders to skip (by Z000# code) — add more as issues arise
# =========================================================
SKIP_CODES=(
    "Z0002"
    # "Z0007"
    # "Z0013"
)

# =========================================================
# Argument config
# =========================================================
BUCKET1=$1 # i.e. "s3://chesilab-testbucket/ultimatest"          # folders live here, one .cram per folder
BUCKET2=$2 # i.e. "s3://chesilab-testpail/ultimatestout"         # fastq.gz outputs go here
WORKDIR=$3 # i.e. "/home/ec2-user/scratch"             # local SSD scratch space
SCRIPT="$(dirname "$(realpath "${BASH_SOURCE[0]}")")/splitCramToFastq.py" #i.e. /home/sjpl/github/Scale-Ultima-CustomLib/1_lib_split/splitCramToFastq.py

mkdir -p "$WORKDIR"
LOGFILE="$WORKDIR/process_log_$(date +%Y%m%d_%H%M%S).txt"

# Helper function for logging.
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOGFILE"
}

# Helper function for skipping specific folders by SKIP_CODE
should_skip() {
    local folder="$1"
    for code in "${SKIP_CODES[@]}"; do
        if [[ "$folder" == *"$code"* ]]; then
            return 0   # true, should skip
        fi
    done
    return 1   # false, don't skip
}

log "=== Starting batch run ==="
log "Bucket1 (input):  $BUCKET1"
log "Bucket2 (output): $BUCKET2"
log "Workdir:          $WORKDIR"

# =========================================================
# Get list of top-level "folder" names in Bucket1
# =========================================================
# aws s3 ls with a trailing slash lists common prefixes (folders) one level deep.
folders=$(aws s3 ls "${BUCKET1}/" | awk '/PRE/ {print $2}' | sed 's#/$##') # if on AWS
# folders=$(ls "${BUCKET1}/" | awk '/PRE/ {print $2}' | sed 's#/$##') # if local


if [ -z "$folders" ]; then
    log "ERROR: No folders found in $BUCKET1 — check the bucket path/permissions."
    exit 1
fi

total=$(echo "$folders" | wc -l)
count=0

for folder in $folders; do
    count=$((count + 1))
    log "--- [$count/$total] Processing folder: $folder ---"

    # ---- Skip folders that don't contain "ChesiLab" ----
    if [[ "$folder" != *ChesiLab* ]]; then
        log "SKIP: ${folder} does not contain 'ChesiLab'."
        continue
    fi

    # ---- Skip known "bad" folders by Z000# code ----
    if should_skip "$folder"; then
        log "SKIP: ${folder} matches a known Z000# code in SKIP_CODES."
        continue
    fi

    cram_key="${BUCKET1}/${folder}/${folder}.cram"
    out_prefix_local="${WORKDIR}/${folder}"          # local output dir for this sample
    local_cram="${WORKDIR}/${folder}.cram"
    out_s3_prefix="${BUCKET2}/${folder}"

    # ---- Skip if already processed (idempotency / safe to re-run) ----
    already_done=$(aws s3 ls "${out_s3_prefix}/${folder}_crispr.fastq.gz" 2>/dev/null || true)
    if [ -n "$already_done" ]; then
        log "SKIP: ${folder} already has output in Bucket2."
        continue
    fi

    # ---- 1. Check the source CRAM exists ----
    exists=$(aws s3 ls "$cram_key" 2>/dev/null || true)
    if [ -z "$exists" ]; then
        log "WARNING: No CRAM found at $cram_key — skipping."
        continue
    fi

    # ---- 2. Copy CRAM down to local SSD ----
    log "Copying $cram_key -> $local_cram"
    if ! aws s3 cp "$cram_key" "$local_cram" --only-show-errors; then
        log "ERROR: Failed to copy $cram_key — skipping."
        rm -f "$local_cram"
        continue
    fi

    # ---- 3. Run processing script ----
    mkdir -p "$out_prefix_local"
    log "Running processing script on $local_cram"
    if ! python3 "$SCRIPT" "$local_cram" "$out_prefix_local" >> "$LOGFILE" 2>&1; then
        log "ERROR: Processing failed for $folder — leaving local files for inspection, skipping upload."
        continue
    fi

    # ---- 4. Upload results to Bucket2 ----
    crispr_fastq="${out_prefix_local}/${folder}_crispr.fastq.gz"
    rna_fastq="${out_prefix_local}/${folder}_rna.fastq.gz"

    upload_ok=true
    if [ -f "$crispr_fastq" ]; then
        aws s3 cp "$crispr_fastq" "${out_s3_prefix}/${folder}_crispr.fastq.gz" --only-show-errors || upload_ok=false
    else
        log "WARNING: Expected output missing: $crispr_fastq"
        upload_ok=false
    fi

    if [ -f "$rna_fastq" ]; then
        aws s3 cp "$rna_fastq" "${out_s3_prefix}/${folder}_rna.fastq.gz" --only-show-errors || upload_ok=false
    else
        log "WARNING: Expected output missing: $rna_fastq"
        upload_ok=false
    fi

    if [ "$upload_ok" = false ]; then
        log "ERROR: Upload incomplete for $folder — local files kept for inspection."
        continue
    fi

    log "Upload confirmed for $folder."

    # ---- 5. Clean up local SSD (cram + fastqs + folder) ----
    rm -f "$local_cram"
    rm -rf "$out_prefix_local"
    log "Cleaned up local files for $folder."

done

log "=== Batch run complete: $count/$total folders processed ==="