#!/usr/bin/env python3
"""Compute +/-5000bp candidate regions for SORCS1/GTF2I/GTF2IRD1 in
VulVul3 (GCF_048418805.1), clamped to each chromosome's actual length so
padding never runs off the end. Writes variants_fox/candidate_regions.bed.

Superseded VulVul2.2 (GCF_003160815.1, PRJNA378561) as the reference for
this stage: VulVul2.2 is a 2018 Illumina/SOAPdenovo scaffold assembly
(contig N50 = 55kb, 82,423 scaffolds) and is RefSeq-suppressed in favor
of VulVul3 (GCF_048418805.1, PRJNA1206188), a 2025 PacBio HiFi/Hifiasm
chromosome-level assembly (contig N50 = 55.7Mb, 17 chromosomes). The
switch resolves three concrete problems documented in the VulVul2.2-stage
findings: GTF2I had no gene symbol in VulVul2.2's annotation (Gnomon left
it as LOC112935683) but is named correctly in VulVul3; GTF2I sat only
4.4kb from its (4Mb) scaffold's edge in VulVul2.2 vs 72Mb of margin on a
152Mb chromosome in VulVul3; and the VulVul2.2 SORCS1 transcript model was
flagged partial=true (assembly-gap truncated) while its VulVul3
counterpart (XM_025990732.2, an updated version of the very same
transcript ID) is not.

Gene coordinates below were read directly out of
reference_fox3/GCF_048418805.1_genomic.gff (`gene` features), not from
NCBI Gene's default summary.
"""
GENES = [
    # name,      scaffold,      start,    end,      strand, scaffold_len
    ("SORCS1",   "NC_132794.1", 96013244, 96508980, "-",    118630945),
    ("GTF2I",    "NC_132782.1", 80253639, 80368744, "-",    152625573),
    ("GTF2IRD1", "NC_132782.1", 80415590, 80526170, "-",    152625573),
]
PAD = 5000

def main():
    rows = []
    for name, scaf, start, end, strand, slen in GENES:
        p_start = max(1, start - PAD)
        p_end = min(slen, end + PAD)
        clamped = (start - PAD) < 1 or (end + PAD) > slen
        margin_left = start - PAD - 1
        margin_right = slen - (end + PAD)
        rows.append((name, scaf, start, end, strand, slen, p_start, p_end))
        print(f"{name:10s} {scaf} gene={start}-{end} ({strand}) scaffold_len={slen:,}")
        print(f"    padded = {p_start}-{p_end}  clamped={clamped}  "
              f"margin_left={margin_left:,}bp  margin_right={margin_right:,}bp")

    with open("variants_fox/candidate_regions.bed", "w") as f:
        for name, scaf, start, end, strand, slen, p_start, p_end in rows:
            # BED is 0-based half-open
            f.write(f"{scaf}\t{p_start - 1}\t{p_end}\t{name}_{scaf}\n")
    print("\nWrote variants_fox/candidate_regions.bed")

if __name__ == "__main__":
    main()
