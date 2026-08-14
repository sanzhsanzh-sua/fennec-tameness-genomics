#!/usr/bin/env python3
"""
Classify PASS variants in variants_fox/fox_candidate_regions.PASS.vcf as
private/shared/core (same logic as the fennec stage's
scripts/04_private_shared_analysis.py) across the 9-sample expanded
pilot: 5 tame + 4 aggressive, no conventional line (see
metadata/pilot5_samples.tsv). Asymmetric n because two ~7.5Gbp
aggressive samples (AGGR03, AGGR05) both hit a reproducible
fasterq-dump disk-exhaustion failure and were dropped rather than
retried indefinitely -- see findings_summary.md. Conventional line was
descoped for this stage (network/time constraints); the "line-fixed"
check below degrades gracefully to a tame-vs-aggressive-only comparison
when the conventional group is empty.
"""
import csv
import subprocess
import sys
from pathlib import Path

VCF = Path("variants_fox/fox_candidate_regions.PASS.vcf")
SAMPLE_SHEET = Path("belyaev_foxes/metadata/pilot5_samples.tsv")


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
    conv = [s for s in samples if pop_of.get(s) == "conventional"]

    private_count = {s: 0 for s in samples}
    shared_count = {s: 0 for s in samples}
    site_carrier_count = {}
    core_sites = []
    line_fixed_sites = []  # fixed in exactly one line, absent from both others

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

        tc = sum(1 for s in tame if carries_alt(gts[s]))
        ac = sum(1 for s in aggr if carries_alt(gts[s]))
        cc = sum(1 for s in conv if carries_alt(gts[s]))
        counts = {"tame": (tc, len(tame)), "aggressive": (ac, len(aggr)), "conventional": (cc, len(conv))}
        fixed_lines = [name for name, (k, n_) in counts.items() if k == n_ and n_ > 0]
        absent_lines = [name for name, (k, n_) in counts.items() if k == 0]
        if len(fixed_lines) == 1 and len(absent_lines) == 2:
            line_fixed_sites.append((chrom, pos, fixed_lines[0], tc, ac, cc))

    total = sum(site_carrier_count.values())
    print(f"Total PASS sites: {total}  (samples: {len(samples)} = {len(tame)} tame + {len(aggr)} aggressive + {len(conv)} conventional)")
    print(f"  private (exactly 1 sample):  {site_carrier_count.get(1, 0)}")
    print(f"  shared  (2+ samples):        {total - site_carrier_count.get(0, 0) - site_carrier_count.get(1, 0)}")
    print(f"  core    (all {len(samples)} samples):    {len(core_sites)}")
    print()
    print(f"{'Sample':<14}{'Line':<14}{'Private':>10}{'Shared':>10}{'Total':>8}")
    for s in samples:
        print(f"{s:<14}{pop_of.get(s, '?'):<14}{private_count[s]:>10}{shared_count[s]:>10}{private_count[s] + shared_count[s]:>8}")

    print()
    print("Sites by number of carrier samples:")
    for n in sorted(site_carrier_count):
        print(f"  {n:>2} samples: {site_carrier_count[n]} sites")

    print()
    print(f"Line-fixed sites (ALT fixed in exactly one of tame/aggressive/conventional, absent from both others, n={len(tame)}+{len(aggr)}+{len(conv)}): {len(line_fixed_sites)}")
    for chrom, pos, which, tc, ac, cc in line_fixed_sites:
        print(f"  {chrom}:{pos}  {which}-fixed  (tame={tc}/{len(tame)}, aggressive={ac}/{len(aggr)}, conventional={cc}/{len(conv)})")


if __name__ == "__main__":
    sys.exit(main())
