#!/usr/bin/env bash
#
# count_cram_reads.sh
# Count reads in every .cram file (one per folder) in an S3 bucket,
# skipping loose files at the bucket root.
#
# Requirements on the EC2 instance:
#   - awscli   (aws --version)
#   - samtools (samtools --version)  # v1.10+ recommended
#   - an IAM role / credentials with s3:ListBucket and s3:GetObject
#
# Usage:
#   ./count_cram_reads.sh s3://my-bucket-name [optional/prefix/]

set -euo pipefail

BUCKET_URI="${1:?Usage: $0 s3://bucket-name [prefix/]}"
PREFIX="${2:-}"
OUTFILE="cram_read_counts.csv"

# Strip the s3:// scheme to get the bare bucket name for `aws s3api`
BUCKET="${BUCKET_URI#s3://}"
BUCKET="${BUCKET%%/*}"

echo "folder,cram_key,read_count" > "$OUTFILE"

echo "Listing .cram objects in s3://${BUCKET}/${PREFIX} ..." >&2

# Enumerate every object key, keep only .cram files that live inside a folder
aws s3api list-objects-v2 \
    --bucket "$BUCKET" \
    --prefix "$PREFIX" \
    --query 'Contents[].Key' \
    --output text \
  | tr '\t' '\n' \
  | grep -E '\.cram$' \
  | grep -E '/' \
  | while read -r KEY; do

        # The folder is everything up to the last slash
        FOLDER="${KEY%/*}"

        echo "Counting: ${KEY}" >&2

        # Stream the CRAM from S3 and count records.
        # -@ 4 gives samtools 4 decode threads; bump/lower to match your instance.
        COUNT=$(aws s3 cp "s3://${BUCKET}/${KEY}" - \
                  | samtools view -c -@ 4 - 2>/dev/null) || COUNT="ERROR"

        echo "${FOLDER},${KEY},${COUNT}" >> "$OUTFILE"
        echo "  -> ${COUNT} reads" >&2
    done

echo "Done. Results written to ${OUTFILE}" >&2