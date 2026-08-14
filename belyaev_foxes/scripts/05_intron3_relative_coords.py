#!/usr/bin/env python3
"""
Cross-species comparison in RELATIVE coordinates, per the task requirement
not to compare raw bp positions across assemblies/species directly.

SORCS1 intron 3 boundaries (from GFF3 mRNA exon models, both single-
isoform, both 27 exons -- see metadata/sorcs1_intron3_boundaries.tsv):

  dog (CanFam3.1, NC_006610.3, XM_005637877.2, - strand):
    intron 3 = 18,797,250-18,848,312 (51,063 bp)
    5' boundary (adjacent to exon3) = 18,848,312

  fox (VulVul3, NC_132794.1, XM_025990732.2, - strand):
    intron 3 = 96,192,951-96,241,834 (48,884 bp)
    5' boundary (adjacent to exon3) = 96,241,834

(VulVul2.2's corresponding intron -- 42,053,244-42,102,689,
NW_020356435.1 -- is kept for the record in
metadata/sorcs1_intron3_boundaries.tsv but superseded once the fox stage
switched reference to VulVul3.)

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
DOG_N_SAMPLES = 10

FOX_INTRON3_5P = 96_241_834
FOX_INTRON3_3P = 96_192_951
FOX_LEN = FOX_INTRON3_5P - FOX_INTRON3_3P
FOX_N_SAMPLES = 9  # 5 tame + 4 aggressive (expanded pilot, no conventional)

DOG_VCF = Path("variants/fennec_candidate_regions.PASS.vcf")
FOX_VCF = Path("variants_fox/fox_candidate_regions.PASS.vcf")
DOG_CHROM = "NC_006610.3"
FOX_CHROM = "NC_132794.1"


def rel_offset(pos, boundary5p):
    return boundary5p - pos


def carriers(gt_list):
    return sum(1 for gt in gt_list if any(a not in (".", "0") for a in gt.replace("|", "/").split("/")))


def query_vcf_sites(vcf, chrom, lo, hi):
    out = subprocess.run(
        ["bcftools", "query", "-t", f"{chrom}:{lo}-{hi}",
         "-f", "%POS[\t%GT]\n", str(vcf)],
        capture_output=True, text=True, check=True,
    ).stdout.strip()
    result = []
    for line in out.split("\n") if out else []:
        parts = line.split("\t")
        result.append((int(parts[0]), carriers(parts[1:])))
    return result


def main():
    print("=== SORCS1 intron 3, relative coordinates (distance from exon3-adjacent boundary) ===")
    print(f"dog:  {DOG_INTRON3_3P}-{DOG_INTRON3_5P}  length={DOG_LEN}bp")
    print(f"fox:  {FOX_INTRON3_3P}-{FOX_INTRON3_5P}  length={FOX_LEN}bp")
    print(f"length ratio fox/dog = {FOX_LEN / DOG_LEN:.3f}")
    print()

    dog_known_snp = 18_801_488
    known_snp_pct = 100 * rel_offset(dog_known_snp, DOG_INTRON3_5P) / DOG_LEN
    print(f"dog rs23402730 (liftover CanFam3.1 pos {dog_known_snp}): "
          f"rel={rel_offset(dog_known_snp, DOG_INTRON3_5P)}bp ({known_snp_pct:.1f}% into intron)")
    print()

    if not DOG_VCF.exists() or not FOX_VCF.exists():
        print("[!] One or both PASS VCFs not found yet -- run the variant-calling "
              "steps for both stages before the final comparison.")
        print(f"    dog VCF present: {DOG_VCF.exists()}  ({DOG_VCF})")
        print(f"    fox VCF present: {FOX_VCF.exists()}  ({FOX_VCF})")
        return 1

    dog_sites = [(pos, rel_offset(pos, DOG_INTRON3_5P), 100 * rel_offset(pos, DOG_INTRON3_5P) / DOG_LEN, n)
                 for pos, n in query_vcf_sites(DOG_VCF, DOG_CHROM, DOG_INTRON3_3P, DOG_INTRON3_5P)]
    fox_sites = [(pos, rel_offset(pos, FOX_INTRON3_5P), 100 * rel_offset(pos, FOX_INTRON3_5P) / FOX_LEN, n)
                 for pos, n in query_vcf_sites(FOX_VCF, FOX_CHROM, FOX_INTRON3_3P, FOX_INTRON3_5P)]

    dog_core = [s for s in dog_sites if s[3] == DOG_N_SAMPLES]
    fox_core = [s for s in fox_sites if s[3] == FOX_N_SAMPLES]
    print(f"dog intron3 PASS sites: {len(dog_sites)}  (core={DOG_N_SAMPLES}/{DOG_N_SAMPLES}: {len(dog_core)})")
    print(f"fox intron3 PASS sites: {len(fox_sites)}  (core={FOX_N_SAMPLES}/{FOX_N_SAMPLES}: {len(fox_core)})")
    print()

    print(f"=== known rs23402730 position: {known_snp_pct:.1f}% into intron 3 ===")
    print("nearest dog PASS sites:")
    for pos, rel, pct, n in sorted(dog_sites, key=lambda s: abs(s[2] - known_snp_pct))[:5]:
        print(f"  dog {pos} ({pct:5.1f}%, {n}/{DOG_N_SAMPLES})")
    print("nearest fox PASS sites:")
    for pos, rel, pct, n in sorted(fox_sites, key=lambda s: abs(s[2] - known_snp_pct))[:5]:
        print(f"  fox {pos} ({pct:5.1f}%, {n}/{FOX_N_SAMPLES})")
    print()

    print(f"=== dog CORE sites ({DOG_N_SAMPLES}/{DOG_N_SAMPLES}) vs nearest fox site by %% into intron ===")
    for pos, rel, pct, n in dog_core:
        if fox_sites:
            fpos, frel, fpct, fn = min(fox_sites, key=lambda f: abs(f[2] - pct))
            print(f"  dog {pos} ({pct:5.1f}%, {n}/{DOG_N_SAMPLES})  <->  "
                  f"fox {fpos} ({fpct:5.1f}%, {fn}/{FOX_N_SAMPLES})  delta%={abs(pct - fpct):.2f}")


if __name__ == "__main__":
    sys.exit(main())
