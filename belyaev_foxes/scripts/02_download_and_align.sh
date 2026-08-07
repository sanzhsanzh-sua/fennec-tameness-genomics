#!/bin/bash
# Download + align the 6 pilot Belyaev-fox samples (3 tame, 3 aggressive;
# see metadata/pilot_samples.tsv) WITHOUT ever writing a full-genome fastq
# to disk. Each sample's SRA object is whole-genome (5.7-12.8 Gbp of
# reads), which as plain fastq would be tens of GB -- far more than we
# need for 3 candidate genes and more than fits in the available disk
# budget for 6 samples at once.
#
# Instead: prefetch the .sra (~4-8GB), stream it through
# `fasterq-dump --stdout --split-spot` (interleaved paired fastq on
# stdout, nothing written to disk) straight into `bwa mem -p` against the
# 2-scaffold subset reference, drop unmapped reads immediately
# (samtools view -F 4), sort, index, then delete the .sra. Peak extra
# disk per sample is just the one .sra file; >99% of each sample's reads
# don't match either scaffold and never touch disk as fastq.
set -euo pipefail

source ~/miniforge3/etc/profile.d/conda.sh
conda activate fennec-genomics

cd "$(dirname "$0")/../.."

REF=reference_fox/reference_subset.fna
SRA_TMP=sra_tmp
mkdir -p aligned_fox "$SRA_TMP"
MASTER_LOG=aligned_fox/align.master.log

tail -n +2 belyaev_foxes/metadata/pilot_samples.tsv | while IFS=$'\t' read -r RUN SAMPLE POP SEX BIOSAMPLE SIZE_MB BASES; do
  BAM="aligned_fox/${RUN}.sorted.bam"

  if [ -f "$BAM" ] && samtools quickcheck -v "$BAM" 2>/dev/null; then
    echo "=== [$RUN/$SAMPLE] already present and verified, skipping ===" >> "$MASTER_LOG"
    continue
  fi

  echo "=== [$RUN/$SAMPLE/$POP] prefetch: $(date) ===" >> "$MASTER_LOG"
  rm -rf "${SRA_TMP:?}/${RUN}"
  prefetch "$RUN" -O "$SRA_TMP" >> "$MASTER_LOG" 2>&1

  echo "=== [$RUN/$SAMPLE/$POP] stream fasterq-dump | bwa mem | samtools: $(date) ===" >> "$MASTER_LOG"
  fasterq-dump --stdout --split-spot --skip-technical -e 4 "${SRA_TMP}/${RUN}/${RUN}.sra" 2>> "aligned_fox/${RUN}.fqd.log" \
    | bwa mem -p -t 4 -R "@RG\tID:${RUN}\tSM:${RUN}\tPL:ILLUMINA" "$REF" - 2> "aligned_fox/${RUN}.bwa.log" \
    | samtools view -b -F 4 - \
    | samtools sort -o "$BAM" -

  samtools index "$BAM"

  if samtools quickcheck -v "$BAM" >> "$MASTER_LOG" 2>&1; then
    echo "=== [$RUN/$SAMPLE/$POP] DONE and verified: $(date) ===" >> "$MASTER_LOG"
  else
    echo "=== [$RUN/$SAMPLE/$POP] FAILED quickcheck: $(date) ===" >> "$MASTER_LOG"
  fi

  # reclaim disk immediately -- the .sra is not needed once the region-
  # filtered BAM exists
  rm -rf "${SRA_TMP:?}/${RUN}"
done

echo "=== Fox alignment run finished: $(date) ===" >> "$MASTER_LOG"
