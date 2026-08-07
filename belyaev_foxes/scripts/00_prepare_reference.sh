#!/bin/bash
# Reference preparation: VulVul2.2 (GCF_003160815.1, PRJNA378561 -- the fox
# reference *assembly* bioproject, not to be confused with PRJNA376561, the
# WGS resequencing bioproject the reads come from). RefSeq flags this
# assembly "suppressed" (superseded by the 2025 VulVul3 chromosome-level
# assembly) but it is still fully downloadable and has a complete
# Gnomon/RefSeq annotation (Annotation Release 100) -- the only annotation
# whose coordinates match this scaffold numbering.
#
# Subset to the two scaffolds carrying our candidate genes:
#   NW_020356435.1 (55,683,013 bp) -- SORCS1
#   NW_020356599.1 ( 4,086,530 bp) -- GTF2IRD1 and GTF2I
# GTF2I is NOT annotated under its gene symbol in this assembly -- Gnomon
# left it as LOC112935683 ("general transcription factor II-I-like"),
# almost certainly because GTF2I/GTF2IRD1/GTF2IRD2 sit in a segmental
# duplication that automated annotation is conservative about naming.
# Confirmed by transcript product description and synteny (ELN-LIMK1-
# EIF4H-LAT2-RFC2-CLIP2-GTF2IRD1-GTF2I, matching the canonical WBS region
# gene order).
set -euo pipefail

source ~/miniforge3/etc/profile.d/conda.sh
conda activate fennec-genomics

mkdir -p reference_fox
cd reference_fox

datasets download genome accession GCF_003160815.1 --include genome,gff3 --filename GCF_003160815.1.zip
unzip -o GCF_003160815.1.zip -d GCF_003160815.1_genome

FNA=GCF_003160815.1_genome/ncbi_dataset/data/GCF_003160815.1/GCF_003160815.1_VulVul2.2_genomic.fna
GFF=GCF_003160815.1_genome/ncbi_dataset/data/GCF_003160815.1/genomic.gff

samtools faidx "$FNA"
samtools faidx "$FNA" NW_020356435.1 NW_020356599.1 > reference_subset.fna
samtools faidx reference_subset.fna
bwa index reference_subset.fna

# Keep the GFF3 (needed later for intron-boundary lookup); the full 2.4GB
# genome fasta and its zip are not needed once the subset is built and
# are NOT kept (gitignored / deleted) -- see .gitignore and README.
mv "$GFF" GCF_003160815.1_genomic.gff
rm -rf GCF_003160815.1_genome GCF_003160815.1.zip
