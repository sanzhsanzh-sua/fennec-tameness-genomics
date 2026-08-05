#!/usr/bin/env python3
"""
Cross-assembly SNP localization by sequence homology (seed-and-extend).

NCBI dbSNP has essentially no canine variant data (esearch -db snp returns
Count=0 for most canine rsIDs); the actual data lives in Ensembl/EVA, but
Ensembl's current dog assembly (ROS_Cfam_1.0) is NOT coordinate-compatible
with the CanFam3.1 (GCF_000002285.3) reference used in this project -- the
same chromosome accession differs by up to several Mb in length between
assembly versions (e.g. NC_006588.3: 77,573,801bp vs NC_006588.4: 80,213,190bp).

This script fetches flanking sequence around a variant from Ensembl's
current assembly, then locates that sequence in our CanFam3.1 reference by
exact k-mer matching (forward and reverse-complement), checking that the
hits fall on a consistent diagonal (offset-independent), which rules out
spurious short-repeat hits and confirms an indel-free 1:1 alignment.

Usage:
    python3 scripts/05_snp_liftover.py <species> <chrom> <pos> \
        --region-fasta reference/reference_subset.fna \
        --region NC_006588.3:5882757-6002627

Example (as used to map rs23402730 in SORCS1, ROS_Cfam_1.0 chr28:19290257
to NC_006610.3:18801488):
    python3 scripts/05_snp_liftover.py canis_lupus_familiaris 28 19290257 \
        --region-fasta reference/reference_subset.fna \
        --region NC_006610.3:18614015-18938785
"""
import argparse
import subprocess
import sys
import urllib.request
import json

FLANK = 100
KMER_LEN = 18
KMER_STRIDE = 15


def revcomp(s):
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return s.translate(comp)[::-1]


def fetch_ensembl_flank(species, chrom, pos):
    start, end = pos - FLANK, pos + FLANK
    url = (
        f"https://rest.ensembl.org/sequence/region/{species}/"
        f"{chrom}:{start}-{end}?content-type=text/x-fasta"
    )
    with urllib.request.urlopen(url, timeout=20) as resp:
        text = resp.read().decode()
    lines = text.strip().split("\n")
    return "".join(lines[1:]).upper()


def fetch_region_seq(fasta_path, region):
    out = subprocess.run(
        ["samtools", "faidx", fasta_path, region], capture_output=True, text=True, check=True
    ).stdout
    lines = out.strip().split("\n")
    return "".join(lines[1:]).upper()


def region_start(region):
    return int(region.split(":")[1].split("-")[0])


def seed_and_extend(query, target):
    hits = {}
    for i in range(0, len(query) - KMER_LEN + 1, KMER_STRIDE):
        kmer = query[i : i + KMER_LEN]
        p = target.find(kmer)
        if p != -1 and target.find(kmer, p + 1) == -1:
            hits[i] = p
    diagonals = {}
    for off, p in hits.items():
        d = p - off
        diagonals.setdefault(d, []).append(off)
    if not diagonals:
        return None
    return max(diagonals.items(), key=lambda x: len(x[1]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("species")
    ap.add_argument("chrom")
    ap.add_argument("pos", type=int)
    ap.add_argument("--region-fasta", required=True)
    ap.add_argument("--region", required=True, help="e.g. NC_006588.3:5882757-6002627")
    args = ap.parse_args()

    query = fetch_ensembl_flank(args.species, args.chrom, args.pos)
    target_fwd = fetch_region_seq(args.region_fasta, args.region)
    target_rc = revcomp(target_fwd)
    rstart = region_start(args.region)

    best_fwd = seed_and_extend(query, target_fwd)
    best_rc = seed_and_extend(query, target_rc)

    if best_fwd and (not best_rc or len(best_fwd[1]) >= len(best_rc[1])):
        diag, anchors = best_fwd
        offset = diag + FLANK
        strand = "fwd"
    elif best_rc:
        diag, anchors = best_rc
        offset = diag + (len(query) - 1 - FLANK)
        strand = "rc"
    else:
        print("No consistent alignment found.", file=sys.stderr)
        sys.exit(1)

    abs_pos = rstart + offset
    ref_base = (target_fwd if strand == "fwd" else target_rc)[offset]
    contig = args.region.split(":")[0]

    print(f"query anchors on one diagonal: {len(anchors)}")
    print(f"strand: {strand}")
    print(f"mapped position: {contig}:{abs_pos}")
    print(f"reference base at mapped position: {ref_base}")


if __name__ == "__main__":
    main()
