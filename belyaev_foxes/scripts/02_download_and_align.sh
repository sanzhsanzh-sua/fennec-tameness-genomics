#!/bin/bash
# Download + align Belyaev-fox samples against the VulVul3 subset
# reference, WITHOUT ever writing a full-genome fastq to disk. Sample
# list defaults to metadata/full_samples.tsv (all 30: 10 tame, 10
# aggressive, 10 conventional) or pass an alternate TSV path as $1, e.g.
# metadata/pilot5_samples.tsv (5 tame + 5 aggressive, no conventional --
# used for the expanded-pilot stage to keep runtime to ~2 days instead
# of ~6-7 on unreliable power).
#
# Each sample's SRA object is whole-genome (5.7-12.8 Gbp of reads), which
# as plain fastq would be tens of GB. Instead: prefetch the .sra
# (~3.6-7.8GB), stream it through `fasterq-dump --stdout --split-spot`
# (interleaved paired fastq on stdout, nothing written to disk as fastq)
# straight into `bwa mem -p` against the 2-chromosome subset reference,
# drop unmapped reads immediately (samtools view -F 4), sort, index, then
# delete the .sra. Peak extra disk per sample is one .sra file.
#
# Network in this environment is slow (~0.5-0.8 MB/s observed) and prone
# to transient failures, so every step below retries:
#   - prefetch: up to 3 attempts, clearing a stale .sra.lock between tries
#     (sra-tools itself will often self-heal a mid-download timeout on
#     retry, resuming from its own checkpoint rather than starting over)
#   - fasterq-dump: run with --size-check off, since its default
#     disk-limit heuristic throws a false-positive "disk-limit exeeded!"
#     even with tens of GB genuinely free, silently producing an empty
#     (header-only) BAM that `samtools quickcheck` does NOT catch
#   - post-alignment: verified via `samtools flagstat` (mapped reads > 0),
#     not just quickcheck (which only checks structural validity) --  if
#     a sample comes out with 0 mapped reads, the whole sample (prefetch
#     through alignment) is retried up to 2 more times
set -uo pipefail

source ~/miniforge3/etc/profile.d/conda.sh
conda activate fennec-genomics

cd "$(dirname "$0")/../.."

SAMPLE_SHEET="${1:-belyaev_foxes/metadata/full_samples.tsv}"
REF=reference_fox3/reference_subset.fna
SRA_TMP=sra_tmp
FQD_TMP=fqd_tmp
mkdir -p aligned_fox "$SRA_TMP" "$FQD_TMP"
MASTER_LOG=aligned_fox/align.master.log

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$MASTER_LOG"; }

align_one_sample() {
  local RUN=$1 SAMPLE=$2 POP=$3
  local BAM="aligned_fox/${RUN}.sorted.bam"

  # Deliberately NOT wiping ${SRA_TMP}/${RUN} here: prefetch supports
  # resuming a partial .sra.tmp via HTTP range requests (confirmed
  # working in practice -- e.g. "Continue download of 'X' from
  # <byte-offset>"), so a prior interrupted/failed attempt's partial
  # download is left in place on purpose. Only the lock file (which goes
  # stale if a previous prefetch process was killed rather than exiting
  # cleanly) is cleared before each attempt.
  local ATTEMPT
  for ATTEMPT in 1 2 3; do
    log "[$RUN/$SAMPLE/$POP] prefetch attempt $ATTEMPT"
    rm -f "${SRA_TMP}/${RUN}/${RUN}.sra.lock"
    prefetch "$RUN" -O "$SRA_TMP" >> "$MASTER_LOG" 2>&1 && break
    log "[$RUN/$SAMPLE/$POP] prefetch attempt $ATTEMPT failed"
    if [ "$ATTEMPT" = 3 ]; then
      log "[$RUN/$SAMPLE/$POP] GIVING UP on prefetch after 3 attempts"
      return 1
    fi
  done

  if [ ! -f "${SRA_TMP}/${RUN}/${RUN}.sra" ]; then
    log "[$RUN/$SAMPLE/$POP] no .sra file after prefetch, skipping"
    return 1
  fi

  log "[$RUN/$SAMPLE/$POP] stream fasterq-dump | bwa mem | samtools sort"
  fasterq-dump --stdout --split-spot --skip-technical --size-check off \
      -t "$FQD_TMP" -e 4 "${SRA_TMP}/${RUN}/${RUN}.sra" \
      2> "aligned_fox/${RUN}.fqd.log" \
    | bwa mem -p -t 4 -R "@RG\tID:${RUN}\tSM:${RUN}\tPL:ILLUMINA" "$REF" - \
      2> "aligned_fox/${RUN}.bwa.log" \
    | samtools view -b -F 4 - \
    | samtools sort -o "$BAM" -

  # fqd_tmp is fasterq-dump's own scratch space, never resumable/reusable
  # -- always safe to clear. The .sra, by contrast, is NOT cleared here:
  # if fasterq-dump/bwa/sort failed downstream of a successful download
  # (e.g. the recurring "storage exhausted" fasterq-dump failure), the
  # .sra is left in place so a retry skips re-downloading it entirely.
  # It's only removed once this sample is fully DONE, below.
  rm -rf "${FQD_TMP:?}"/*

  if ! samtools quickcheck -v "$BAM" 2>> "$MASTER_LOG"; then
    log "[$RUN/$SAMPLE/$POP] FAILED quickcheck"
    rm -f "$BAM"
    return 1
  fi

  local MAPPED
  MAPPED=$(samtools flagstat "$BAM" | head -1 | awk '{print $1}')
  if [ -z "$MAPPED" ] || [ "$MAPPED" -eq 0 ]; then
    log "[$RUN/$SAMPLE/$POP] FAILED: 0 mapped reads in BAM (likely fasterq-dump/bwa pipeline error)"
    rm -f "$BAM"
    return 1
  fi

  samtools index "$BAM"
  rm -rf "${SRA_TMP:?}/${RUN}"
  log "[$RUN/$SAMPLE/$POP] DONE: $MAPPED mapped reads"
  return 0
}

tail -n +2 "$SAMPLE_SHEET" | while IFS=$'\t' read -r RUN SAMPLE POP SEX BIOSAMPLE SIZE_MB BASES; do
  BAM="aligned_fox/${RUN}.sorted.bam"

  if [ -f "$BAM" ] && samtools quickcheck -v "$BAM" 2>/dev/null; then
    MAPPED=$(samtools flagstat "$BAM" 2>/dev/null | head -1 | awk '{print $1}')
    if [ -n "$MAPPED" ] && [ "$MAPPED" -gt 0 ]; then
      log "[$RUN/$SAMPLE/$POP] already present and verified ($MAPPED mapped reads), skipping"
      continue
    fi
  fi

  SAMPLE_ATTEMPT_OK=0
  for SAMPLE_ATTEMPT in 1 2; do
    if align_one_sample "$RUN" "$SAMPLE" "$POP"; then
      SAMPLE_ATTEMPT_OK=1
      break
    fi
    log "[$RUN/$SAMPLE/$POP] sample-level attempt $SAMPLE_ATTEMPT failed, retrying whole sample"
  done

  if [ "$SAMPLE_ATTEMPT_OK" -eq 0 ]; then
    log "[$RUN/$SAMPLE/$POP] *** GIVING UP after repeated failures -- needs manual attention ***"
  fi
done

log "=== Fox alignment run finished (sample sheet: $SAMPLE_SHEET) ==="
