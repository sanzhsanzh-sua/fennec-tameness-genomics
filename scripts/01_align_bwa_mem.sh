#!/bin/bash
# Align all accessions in SRR_Acc_List.txt to reference/reference_subset.fna
# with bwa mem, coordinate-sort with samtools, index, and verify with
# samtools quickcheck. Safe to re-run: already-verified BAMs are skipped,
# so an interrupted run can simply be restarted.
set -uo pipefail

source ~/miniforge3/etc/profile.d/conda.sh
conda activate fennec-genomics

mkdir -p aligned
MASTER_LOG=aligned/align.master.log
mapfile -t ACCESSIONS < SRR_Acc_List.txt

echo "=== Alignment run started: $(date) ===" >> "$MASTER_LOG"

for ACC in "${ACCESSIONS[@]}"; do
  BAM="aligned/${ACC}.sorted.bam"

  if [ -f "$BAM" ] && samtools quickcheck -v "$BAM" 2>/dev/null; then
    echo "=== [$ACC] already present and verified, skipping ===" >> "$MASTER_LOG"
    continue
  fi

  # clean up any partial/broken output from a previous interrupted run
  rm -f "$BAM" "${BAM}.bai" "${BAM}".tmp.*.bam "aligned/${ACC}.bwa.log"

  echo "=== [$ACC] starting: $(date) ===" >> "$MASTER_LOG"

  bwa mem -t 4 \
    -R "@RG\tID:${ACC}\tSM:${ACC}\tPL:ILLUMINA" \
    reference/reference_subset.fna \
    "${ACC}/${ACC}_1.fastq.gz" "${ACC}/${ACC}_2.fastq.gz" \
    2> "aligned/${ACC}.bwa.log" \
    | samtools sort -o "$BAM" -

  if [ $? -ne 0 ]; then
    echo "=== [$ACC] FAILED (bwa/sort pipeline): $(date) ===" >> "$MASTER_LOG"
    continue
  fi

  samtools index "$BAM"

  if samtools quickcheck -v "$BAM" >> "$MASTER_LOG" 2>&1; then
    echo "=== [$ACC] DONE and verified: $(date) ===" >> "$MASTER_LOG"
  else
    echo "=== [$ACC] FAILED quickcheck after completion: $(date) ===" >> "$MASTER_LOG"
  fi
done

echo "=== Alignment run finished: $(date) ===" >> "$MASTER_LOG"
