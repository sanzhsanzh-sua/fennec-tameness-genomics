#!/bin/bash
# Build the candidate-region BED file (SORCS1 chr28, GTF2I/GTF2IRD1 chr6,
# CanFam3.1 / GCF_000002285.3 coordinates, +/-5000bp padding) and extract
# those regions from every sorted BAM in aligned/.
#
# Gene coordinates were taken from the official GFF3 annotation for the
# exact assembly version used as reference (GCF_000002285.3) -- NOT from
# NCBI Gene's default summary, which is pinned to a newer, coordinate-
# incompatible chromosome build (NC_006588.4 / NC_006610.4).
set -euo pipefail

source ~/miniforge3/etc/profile.d/conda.sh
conda activate fennec-genomics

mkdir -p variants variants/region_bams

cat > variants/candidate_regions.bed << 'EOF'
NC_006610.3	18614015	18938785	SORCS1_chr28
NC_006588.3	5721652	5845740	GTF2I_chr6
NC_006588.3	5882757	6002627	GTF2IRD1_chr6
EOF

BED=variants/candidate_regions.bed
mapfile -t ACCESSIONS < SRR_Acc_List.txt

for ACC in "${ACCESSIONS[@]}"; do
  echo "=== [$ACC] extracting candidate regions ==="
  samtools view -b -h -L "$BED" "aligned/${ACC}.sorted.bam" \
    > "variants/region_bams/${ACC}.candidate_regions.bam"
  samtools index "variants/region_bams/${ACC}.candidate_regions.bam"
  samtools quickcheck -v "variants/region_bams/${ACC}.candidate_regions.bam"
done

echo "=== all region extractions complete ==="
