#!/bin/bash
# Same approach as the fennec stage (scripts/03_depth_and_call_variants.sh):
# per-base depth (samtools depth) + joint variant calling (bcftools
# mpileup/call) across the 9 fox samples in the expanded pilot (5 tame,
# 4 aggressive -- see metadata/pilot5_samples.tsv; asymmetric because two
# ~7.5Gbp aggressive samples, AGGR03 and AGGR05, both hit a reproducible
# fasterq-dump disk-exhaustion failure at ~6.1-6.2Gbp processed,
# consistent with shrinking available disk headroom as more BAMs
# accumulated over the course of the run), restricted to the 3
# candidate-gene regions, with the same soft MQ<40 -> FILTER=LowMQ rule.
# Reference: VulVul3 (reference_fox3/reference_subset.fna).
set -euo pipefail

source ~/miniforge3/etc/profile.d/conda.sh
conda activate fennec-genomics

cd "$(dirname "$0")/../.."

REF=reference_fox3/reference_subset.fna
BED=variants_fox/candidate_regions.bed
mapfile -t RUNS < <(tail -n +2 belyaev_foxes/metadata/pilot5_samples.tsv | cut -f1)

BAMS=()
for RUN in "${RUNS[@]}"; do
  BAMS+=("aligned_fox/${RUN}.sorted.bam")
done

# --- depth ---
{
  printf "chrom\tpos"
  for RUN in "${RUNS[@]}"; do printf "\t%s" "$RUN"; done
  printf "\n"
} > variants_fox/candidate_regions.depth.txt

samtools depth -a -b "$BED" "${BAMS[@]}" >> variants_fox/candidate_regions.depth.txt

# --- joint variant calling ---
bcftools mpileup -f "$REF" \
  -R "$BED" \
  -a AD,DP \
  -Ou "${BAMS[@]}" \
  | bcftools call -mv -Ov -o variants_fox/fox_candidate_regions.raw.vcf

cp variants_fox/fox_candidate_regions.raw.vcf variants_fox/fox_candidate_regions.unfiltered.vcf

# --- soft MQ<40 filter, same rule as the fennec stage ---
bcftools filter \
  -e 'INFO/MQ<40' \
  -s LowMQ \
  -m + \
  -Ov \
  -o variants_fox/fox_candidate_regions.vcf \
  variants_fox/fox_candidate_regions.unfiltered.vcf

bcftools view -f PASS variants_fox/fox_candidate_regions.vcf \
  -Ov -o variants_fox/fox_candidate_regions.PASS.vcf

rm -f variants_fox/fox_candidate_regions.raw.vcf

echo "Done. See variants_fox/fox_candidate_regions.vcf (all sites, FILTER=PASS/LowMQ)"
echo "and variants_fox/fox_candidate_regions.PASS.vcf (PASS only)."
