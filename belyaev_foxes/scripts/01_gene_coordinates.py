#!/usr/bin/env python3
"""Compute +/-5000bp candidate regions for SORCS1/GTF2I/GTF2IRD1 in
VulVul2.2, clamped to each scaffold's actual length so padding never runs
off the end of a short scaffold. Writes variants_fox/candidate_regions.bed.

Gene coordinates below were read directly out of
reference_fox/GCF_003160815.1_genomic.gff (`gene` features), not from
NCBI Gene's default summary -- Gene's default view points at the newer,
coordinate-incompatible VulVul3 assembly.
"""
GENES = [
    # name,      scaffold,          start,    end,      strand, scaffold_len
    ("SORCS1",   "NW_020356435.1",  41874895, 42378408, "-",    55683013),
    ("GTF2IRD1", "NW_020356599.1",  3841685,  3952986,  "+",    4086530),
    # Annotated as LOC112935683 ("general transcription factor II-I-like")
    ("GTF2I",    "NW_020356599.1",  3999936,  4077133,  "+",    4086530),
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
