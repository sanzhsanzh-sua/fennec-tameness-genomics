#!/bin/bash
# Compute per-base, multi-sample depth (samtools depth) and jointly call
# variants across all 10 samples (bcftools mpileup + call) in the three
# candidate regions. Applies a soft MQ<40 filter (many "variants" in these
# regions turn out to be multi-mapping artifacts from LINE/SINE repeats --
# see report/ for details) and writes a PASS-only VCF as well.
set -euo pipefail

source ~/miniforge3/etc/profile.d/conda.sh
conda activate fennec-genomics

cd "$(dirname "$0")/.."

BED=variants/candidate_regions.bed
mapfile -t ACCESSIONS < SRR_Acc_List.txt
BAMS=()
for ACC in "${ACCESSIONS[@]}"; do
  BAMS+=("variants/region_bams/${ACC}.candidate_regions.bam")
done

# --- depth ---
{
  printf "chrom\tpos"
  for ACC in "${ACCESSIONS[@]}"; do printf "\t%s" "$ACC"; done
  printf "\n"
} > variants/candidate_regions.depth.txt

samtools depth -a -b "$BED" "${BAMS[@]}" >> variants/candidate_regions.depth.txt

# --- joint variant calling ---
bcftools mpileup -f reference/reference_subset.fna \
  -R "$BED" \
  -a AD,DP \
  -Ou "${BAMS[@]}" \
  | bcftools call -mv -Ov -o variants/fennec_candidate_regions.raw.vcf

cp variants/fennec_candidate_regions.raw.vcf variants/fennec_candidate_regions.unfiltered.vcf

# --- soft MQ<40 filter (sites kept, FILTER column marked LowMQ vs PASS) ---
bcftools filter \
  -e 'INFO/MQ<40' \
  -s LowMQ \
  -m + \
  -Ov \
  -o variants/fennec_candidate_regions.vcf \
  variants/fennec_candidate_regions.unfiltered.vcf

bcftools view -f PASS variants/fennec_candidate_regions.vcf \
  -Ov -o variants/fennec_candidate_regions.PASS.vcf

rm -f variants/fennec_candidate_regions.raw.vcf

echo "Done. See variants/fennec_candidate_regions.vcf (all sites, FILTER=PASS/LowMQ)"
echo "and variants/fennec_candidate_regions.PASS.vcf (PASS only)."
