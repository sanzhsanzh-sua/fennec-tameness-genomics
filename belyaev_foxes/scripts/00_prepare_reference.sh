#!/bin/bash
# Reference preparation: VulVul3 (GCF_048418805.1, PRJNA1206188), a 2025
# PacBio HiFi/Hifiasm chromosome-level assembly.
#
# This SUPERSEDES the initial choice of VulVul2.2 (GCF_003160815.1,
# PRJNA378561): VulVul2.2 is a 2018 Illumina/SOAPdenovo assembly (contig
# N50 = 55kb, 82,423 scaffolds) that RefSeq now marks "suppressed" in
# favor of VulVul3 (contig N50 = 55.7Mb, 289 scaffolds, 17 full
# chromosomes). The switch resolves three concrete problems the
# VulVul2.2-based analysis had already documented: GTF2I had no gene
# symbol under VulVul2.2's annotation (Gnomon left it as LOC112935683,
# confirmed by transcript product + synteny) but IS named correctly under
# VulVul3; GTF2I sat only 4.4kb from the edge of its (4Mb) VulVul2.2
# scaffold vs 72Mb of margin on a 152Mb VulVul3 chromosome; and the
# VulVul2.2 SORCS1 transcript model was flagged partial=true
# (assembly-gap truncated) while its VulVul3 counterpart
# (XM_025990732.2, an updated version of the very same transcript ID) is
# not.
#
# Subset to the two chromosomes carrying our candidate genes:
#   NC_132782.1 (152,625,573 bp) -- GTF2I, GTF2IRD1
#   NC_132794.1 (118,630,945 bp) -- SORCS1
set -euo pipefail

source ~/miniforge3/etc/profile.d/conda.sh
conda activate fennec-genomics

mkdir -p reference_fox3
cd reference_fox3

datasets download genome accession GCF_048418805.1 --include genome,gff3 --filename GCF_048418805.1.zip
unzip -o GCF_048418805.1.zip -d GCF_048418805.1_genome

FNA=GCF_048418805.1_genome/ncbi_dataset/data/GCF_048418805.1/GCF_048418805.1_VulVul3_genomic.fna
GFF=GCF_048418805.1_genome/ncbi_dataset/data/GCF_048418805.1/genomic.gff

samtools faidx "$FNA"
samtools faidx "$FNA" NC_132782.1 NC_132794.1 > reference_subset.fna
samtools faidx reference_subset.fna
bwa index reference_subset.fna

# Keep the GFF3 (needed later for intron-boundary lookup); the full 2.4GB
# genome fasta and its zip are not needed once the subset is built and
# are NOT kept (gitignored / deleted) -- see .gitignore and README.
mv "$GFF" GCF_048418805.1_genomic.gff
rm -rf GCF_048418805.1_genome GCF_048418805.1.zip
