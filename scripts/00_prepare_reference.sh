#!/bin/bash
# Reference preparation: CanFam3.1 (GCF_000002285.3), subset to chr6 + chr28,
# bwa index. Not re-run by this repo (reference/ is gitignored, ~2.4GB full
# genome + indexes); documented here for reproducibility.
set -euo pipefail

source ~/miniforge3/etc/profile.d/conda.sh
conda activate fennec-genomics

mkdir -p reference
cd reference

datasets download genome accession GCF_000002285.3 --include genome --filename GCF_000002285.3.zip
unzip -o GCF_000002285.3.zip -d GCF_000002285.3_genome
cp GCF_000002285.3_genome/ncbi_dataset/data/GCF_000002285.3/*_genomic.fna reference.fna
samtools faidx reference.fna

# Subset to the two chromosomes used in this project:
#   NC_006588.3 = chr6  (GTF2I, GTF2IRD1)
#   NC_006610.3 = chr28 (SORCS1)
samtools faidx reference.fna NC_006588.3 NC_006610.3 > reference_subset.fna
samtools faidx reference_subset.fna
bwa index reference_subset.fna
