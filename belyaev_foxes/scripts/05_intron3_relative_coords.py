#!/usr/bin/env python3
"""
Cross-species comparison in RELATIVE coordinates, per the task requirement
not to compare raw bp positions across assemblies/species directly.

SORCS1 intron 3 boundaries (from GFF3 mRNA exon models, both single-
isoform, both 27 exons -- see metadata/sorcs1_intron3_boundaries.tsv):

  dog (CanFam3.1, NC_006610.3, XM_005637877.2, - strand):
    intron 3 = 18,797,250-18,848,312 (51,063 bp)
    5' boundary (adjacent to exon3) = 18,848,312

  fox (VulVul2.2, NW_020356435.1, XM_025990732.1, - strand):
    intron 3 = 42,053,244-42,102,689 (49,446 bp)
    5' boundary (adjacent to exon3) = 42,102,689

Both genes are minus-strand, so transcript 5'->3' direction runs from
high to low genomic coordinate. Relative offset = distance from the
intron's 5' boundary (exon3 side), i.e.:

    rel(pos) = intron_5prime_boundary - pos

which is always in [0, intron_length] for a position inside the intron,
and is directly comparable between species regardless of assembly/strand.
"""
import subprocess
import sys
from pathlib import Path

DOG_INTRON3_5P = 18_848_312
DOG_INTRON3_3P = 18_797_250
DOG_LEN = DOG_INTRON3_5P - DOG_INTRON3_3P

FOX_INTRON3_5P = 42_102_689
FOX_INTRON3_3P = 42_053_244
FOX_LEN = FOX_INTRON3_5P - FOX_INTRON3_3P

DOG_VCF = Path("variants/fennec_candidate_regions.PASS.vcf")
FOX_VCF = Path("variants_fox/fox_candidate_regions.PASS.vcf")
DOG_CHROM = "NC_006610.3"
FOX_CHROM = "NW_020356435.1"


def rel_offset(pos, boundary5p):
    return boundary5p - pos


def query_vcf_sites(vcf, chrom, lo, hi):
    out = subprocess.run(
        ["bcftools", "query", "-t", f"{chrom}:{lo}-{hi}",
         "-f", "%POS\t%REF\t%ALT\t%INFO/DP[\t%GT]\n", str(vcf)],
        capture_output=True, text=True, check=True,
    ).stdout.strip()
    return [line.split("\t") for line in out.split("\n")] if out else []


def main():
    print("=== SORCS1 intron 3, relative coordinates (distance from exon3-adjacent boundary) ===")
    print(f"dog:  {DOG_INTRON3_3P}-{DOG_INTRON3_5P}  length={DOG_LEN}bp")
    print(f"fox:  {FOX_INTRON3_3P}-{FOX_INTRON3_5P}  length={FOX_LEN}bp")
    print(f"length ratio fox/dog = {FOX_LEN / DOG_LEN:.3f}")
    print()

    # known dog SNP rs23402730, liftover position from findings_summary.md
    dog_known_snp = 18_801_488
    print(f"dog rs23402730 (liftover CanFam3.1 pos {dog_known_snp}): "
          f"rel={rel_offset(dog_known_snp, DOG_INTRON3_5P)}bp "
          f"({100 * rel_offset(dog_known_snp, DOG_INTRON3_5P) / DOG_LEN:.1f}% into intron)")
    print()

    if not DOG_VCF.exists() or not FOX_VCF.exists():
        print("[!] One or both PASS VCFs not found yet -- run the variant-calling "
              "steps for both stages before the final comparison.")
        print(f"    dog VCF present: {DOG_VCF.exists()}  ({DOG_VCF})")
        print(f"    fox VCF present: {FOX_VCF.exists()}  ({FOX_VCF})")
        return 1

    print("=== dog PASS sites inside SORCS1 intron 3 ===")
    for pos, ref, alt, dp, *gts in query_vcf_sites(DOG_VCF, DOG_CHROM, DOG_INTRON3_3P, DOG_INTRON3_5P):
        rel = rel_offset(int(pos), DOG_INTRON3_5P)
        pct = 100 * rel / DOG_LEN
        print(f"  {pos} {ref}>{alt}  rel={rel}bp ({pct:.1f}%)  DP={dp}")

    print()
    print("=== fox PASS sites inside SORCS1 intron 3 ===")
    for pos, ref, alt, dp, *gts in query_vcf_sites(FOX_VCF, FOX_CHROM, FOX_INTRON3_3P, FOX_INTRON3_5P):
        rel = rel_offset(int(pos), FOX_INTRON3_5P)
        pct = 100 * rel / FOX_LEN
        print(f"  {pos} {ref}>{alt}  rel={rel}bp ({pct:.1f}%)  DP={dp}")


if __name__ == "__main__":
    sys.exit(main())
