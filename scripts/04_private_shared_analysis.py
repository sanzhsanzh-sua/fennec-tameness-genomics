#!/usr/bin/env python3
"""
Classify PASS variants in variants/fennec_candidate_regions.PASS.vcf as
private (exactly one sample carries an ALT allele) vs shared (2+ samples),
and report per-sample private/shared counts plus the "core" set of sites
fixed (ALT present) across all 10 samples.
"""
import subprocess
import sys
from pathlib import Path

VCF = Path("variants/fennec_candidate_regions.PASS.vcf")


def run(cmd):
    return subprocess.run(cmd, capture_output=True, text=True, check=True).stdout


def main():
    samples = run(["bcftools", "query", "-l", str(VCF)]).strip().split("\n")
    gt_matrix = run(
        ["bcftools", "query", "-f", "%CHROM\t%POS[\t%GT]\n", str(VCF)]
    ).strip().split("\n")

    private_count = {s: 0 for s in samples}
    shared_count = {s: 0 for s in samples}
    site_carrier_count = {}
    core_sites = []

    for line in gt_matrix:
        parts = line.split("\t")
        chrom, pos = parts[0], parts[1]
        gts = parts[2:]
        carriers = [
            s
            for s, gt in zip(samples, gts)
            if any(a not in (".", "0") for a in gt.replace("|", "/").split("/"))
        ]
        n = len(carriers)
        site_carrier_count[n] = site_carrier_count.get(n, 0) + 1
        if n == 1:
            private_count[carriers[0]] += 1
        elif n > 1:
            for s in carriers:
                shared_count[s] += 1
        if n == len(samples):
            core_sites.append((chrom, pos))

    total = sum(site_carrier_count.values())
    print(f"Total PASS sites: {total}")
    print(f"  private (exactly 1 sample):  {site_carrier_count.get(1, 0)}")
    print(f"  shared  (2+ samples):        {total - site_carrier_count.get(0, 0) - site_carrier_count.get(1, 0)}")
    print(f"  core    (all {len(samples)} samples):    {len(core_sites)}")
    print()
    print(f"{'Sample':<16}{'Private':>10}{'Shared':>10}{'Total carried':>16}")
    for s in samples:
        print(f"{s:<16}{private_count[s]:>10}{shared_count[s]:>10}{private_count[s] + shared_count[s]:>16}")

    print()
    print("Sites by number of carrier samples:")
    for n in sorted(site_carrier_count):
        print(f"  {n:>2} samples: {site_carrier_count[n]} sites")


if __name__ == "__main__":
    sys.exit(main())
