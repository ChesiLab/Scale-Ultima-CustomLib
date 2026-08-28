import pysam # pip install pysam; python=3.12; conda install stalls out..
import edlib # pip install edlib
import gzip
import argparse 
from pathlib import Path
# -----------------------
# Description
# -----------------------
# Takes folder name from Ultima as input.
# Processes contained .cram file.

# Set CRISPR sequence that you want to use to extract CRISPR reads.
# Set the path.

# Outputs fastq.gz for downstream processing.

# >>>> Running script
# ./processCram.sh s3://chesilab-testbucket/ultimatest s3://chesilab-testbucket/ultimatestout /home/ec2-user/scratch
# >>>>> If you want to run just ONE file, i.e. for testing: <<<<<
# python3 splitCramToFastq.py ~/scale_preprocess/cutadapt/cram/Z0001_02p.cram Z0001_02p

# -----------------------
# User parameters
# -----------------------

# Comment this out if there are problems! Suppressess "index" error being thrown by AlignmentFile
pysam.set_verbosity(0) 


# Sequence specific to CRISPR lib. 2nd seq is the reverse complement.
crispr_marker = ["TTCCAGCATAGCTCTTAAAC","GTTTAAGAGCTATGCTGGAA"]
max_mismatch = 3 # in bp

# root_dir = Path("home/sjpl/scale_preprocess/lib_split/s3bucket")

# -----------------------
# Helper Functions
# -----------------------

def matches(seq, motifs, max_mismatches=3):
    for m in motifs:
        result = edlib.align(
            m,
            seq,
            mode="HW",  # in-fix mode - i.e. finding a word in a sentence
            k=max_mismatches
        )
        if result["editDistance"] != -1:
            return True
    return False

def write_fastq(fh, read):
    fh.write(
        f"@{read.query_name}\n"
        f"{read.query_sequence}\n"
        "+\n"
        f"{read.qual}\n"
    )

def process_cram(cram_path, out_prefix):
    crispr_out_path = f"{out_prefix}_crispr.fastq.gz"
    rna_out_path = f"{out_prefix}_rna.fastq.gz"

    with pysam.AlignmentFile(
        cram_path,
        "rc", # to open a cram file
        reference_filename= None,
        require_index = False, # don't really need index, we can go sequentially
        check_sq = False # might not need this, usually for unaligned reads.
    ) as cram, \
         gzip.open(crispr_out_path, "wt") as crispr_out, \
         gzip.open(rna_out_path, "wt") as rna_out:
        
        # until_eof = False lets us get away w/ not indexing the .cram file.
        for i, read in enumerate(cram.fetch(until_eof=True)):
            # Skips read with no sequence. (Edge-case)
            if read.query_sequence is None:
                continue
            
            # Read the sequence
            seq = read.query_sequence

            # Check for CRISPR motif, then write sequence to fastq's
            if matches(seq, crispr_marker, max_mismatch):
                write_fastq(crispr_out, read)
            else:
                write_fastq(rna_out, read)

            # Output progress every 100,000 processed reads.
            if i % 100000 == 0 and i > 0:
                print(f"{cram_path.name}: {i:,} reads processed")

# -----------------------
# Set up args
# -----------------------

parser = argparse.ArgumentParser()
parser.add_argument("cram")
parser.add_argument("outputDir")

args = parser.parse_args()

cramPath = Path(args.cram)
outputDir = Path(args.outputDir)

# Make and set output directory.
outputDir.mkdir(parents=True, exist_ok=True)
outPrefix = outputDir / cramPath.stem

# -----------------------
# Main
# -----------------------

# Process the .cram
process_cram(cramPath, outPrefix)