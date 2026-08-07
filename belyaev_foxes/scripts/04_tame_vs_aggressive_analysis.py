#!/usr/bin/env python3
"""
Classify PASS variants in variants_fox/fox_candidate_regions.PASS.vcf as
private/shared/core (same logic as the fennec stage's
scripts/04_private_shared_analysis.py), plus split by line (tame vs
aggressive, from metadata/pilot_samples.tsv) to flag sites where ALT
carriage is concordant within a line but differs between lines --
candidate line-associated sites for this pilot's 3 tame + 3 aggressive
samples. NOTE: n=3 per line, so this is a pilot-scale hypothesis-
generating screen, not a population-genetics test.
"""
import csv
import subprocess
import sys
from pathlib import Path

VCF = Path("variants_fox/fox_candidate_regions.PASS.vcf")
SAMPLE_SHEET = Path("belyaev_foxes/metadata/pilot_samples.tsv")


def run(cmd):
    return subprocess.run(cmd, capture_output=True, text=True, check=True).stdout


def carries_alt(gt):
    return any(a not in (".", "0") for a in gt.replace("|", "/").split("/"))


def main():
    pop_of = {}
    with open(SAMPLE_SHEET) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            pop_of[row["Run"]] = row["Population"]

    samples = run(["bcftools", "query", "-l", str(VCF)]).strip().split("\n")
    gt_matrix = run(
        ["bcftools", "query", "-f", "%CHROM\t%POS[\t%GT]\n", str(VCF)]
    ).strip().split("\n")

    tame = [s for s in samples if pop_of.get(s) == "tame"]
    aggr = [s for s in samples if pop_of.get(s) == "aggressive"]

    private_count = {s: 0 for s in samples}
    shared_count = {s: 0 for s in samples}
    site_carrier_count = {}
    core_sites = []
    line_concordant_sites = []  # all tame same genotype-class, all aggressive the other

    for line in gt_matrix:
        parts = line.split("\t")
        chrom, pos = parts[0], parts[1]
        gts = dict(zip(samples, parts[2:]))
        carriers = [s for s, gt in gts.items() if carries_alt(gt)]
        n = len(carriers)
        site_carrier_count[n] = site_carrier_count.get(n, 0) + 1
        if n == 1:
            private_count[carriers[0]] += 1
        elif n > 1:
            for s in carriers:
                shared_count[s] += 1
        if n == len(samples):
            core_sites.append((chrom, pos))

        tame_carriers = sum(1 for s in tame if carries_alt(gts[s]))
        aggr_carriers = sum(1 for s in aggr if carries_alt(gts[s]))
        if (tame_carriers, aggr_carriers) in ((len(tame), 0), (0, len(aggr))):
            line_concordant_sites.append((chrom, pos, tame_carriers, aggr_carriers))

    total = sum(site_carrier_count.values())
    print(f"Total PASS sites: {total}")
    print(f"  private (exactly 1 sample):  {site_carrier_count.get(1, 0)}")
    print(f"  shared  (2+ samples):        {total - site_carrier_count.get(0, 0) - site_carrier_count.get(1, 0)}")
    print(f"  core    (all {len(samples)} samples):    {len(core_sites)}")
    print()
    print(f"{'Sample':<14}{'Line':<12}{'Private':>10}{'Shared':>10}{'Total':>8}")
    for s in samples:
        print(f"{s:<14}{pop_of.get(s, '?'):<12}{private_count[s]:>10}{shared_count[s]:>10}{private_count[s] + shared_count[s]:>8}")

    print()
    print("Sites by number of carrier samples:")
    for n in sorted(site_carrier_count):
        print(f"  {n:>2} samples: {site_carrier_count[n]} sites")

    print()
    print(f"Line-concordant sites (ALT fixed in one line, absent in the other, n=3+3 pilot): {len(line_concordant_sites)}")
    for chrom, pos, tc, ac in line_concordant_sites:
        which = "tame-fixed" if tc == len(tame) else "aggressive-fixed"
        print(f"  {chrom}:{pos}  {which}  (tame={tc}/{len(tame)}, aggressive={ac}/{len(aggr)})")


if __name__ == "__main__":
    sys.exit(main())
