# fennec-tameness-genomics

Small variant-calling pipeline over 10 fennec fox (*Vulpes zerda*) WGS
samples, restricted to three candidate genes implicated in tameness/social
behavior in canids: **SORCS1** (chr28) and **GTF2I / GTF2IRD1** (chr6, the
Williams-Beuren syndrome region linked to hypersociability in dogs by
[vonHoldt et al. 2017](https://www.science.org/doi/10.1126/sciadv.1700398)).

## Data source

- Reads: [BioProject PRJNA951250](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA951250) (SRA), 10 runs listed in [`SRR_Acc_List.txt`](SRR_Acc_List.txt)
- Reference: CanFam3.1, [GCF_000002285.3](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000002285.3/) (dog), subset to chr6 (NC_006588.3) + chr28 (NC_006610.3) — fennec has no chromosome-level assembly of its own, so alignment used the closest well-annotated canid reference

Raw sequencing data, the reference FASTA/bwa index, and intermediate BAMs
are **not** included in this repository (tens of GB) — see
[Reproducing](#reproducing) below. What *is* included are the pipeline
scripts and the compact, final analysis outputs (VCFs, depth table,
findings summary, report).

## Pipeline

```
scripts/00_prepare_reference.sh        download CanFam3.1, subset to chr6+chr28, bwa index
scripts/01_align_bwa_mem.sh            bwa mem | samtools sort, per accession, with quickcheck verification
scripts/02_extract_candidate_regions.sh   extract SORCS1 / GTF2I / GTF2IRD1 (+/-5000bp) from each BAM
scripts/03_depth_and_call_variants.sh  samtools depth + bcftools mpileup/call (joint, 10 samples) + MQ<40 soft filter
scripts/04_private_shared_analysis.py  classify PASS variants as private / shared / core across samples
scripts/05_snp_liftover.py             map a known dog SNP (Ensembl/EVA) onto the CanFam3.1 reference by sequence homology
```

Environment: `environment.yaml` (conda), plus `bcftools`, `bedtools`, and
`entrez-direct` installed into the same env during analysis.

Gene coordinates were pulled from the GFF3 annotation matching the exact
reference assembly version (GCF_000002285.3) — not from NCBI Gene's
default summary, which by default returns coordinates for a newer,
coordinate-incompatible chromosome build.

## Key findings

Full writeup: [`variants/findings_summary.md`](variants/findings_summary.md) · [`report/fennec_candidate_regions_report.docx`](report/fennec_candidate_regions_report.docx)

- **17,370** joint-called variant sites across the three regions in 10 samples; **18.7%** flagged `LowMQ` (MQ<40) and excluded from the PASS set (14,129 sites).
- Two clusters of superficially high-QUAL variants turned out to be **repeat-driven mapping artifacts**: a LINE/L1MC3 element upstream of GTF2I (confirmed via UCSC RepeatMasker + bedtools intersect) and the SINEC_Cf transposon in GTF2I intron 17 described by vonHoldt et al. — both show 3-5x inflated depth and MQ collapsed to 9-38.
- **1,284** sites are fixed (ALT allele) across all 10 fennec samples — likely fennec-vs-dog-reference divergence rather than within-species polymorphism.
- NCBI dbSNP has essentially no canine variant data; known dog SNPs (e.g. rs23402730 in SORCS1, GTF2IRD1 missense variants) were located via Ensembl/EVA and mapped onto CanFam3.1 by direct sequence-homology liftover (`scripts/05_snp_liftover.py`), since the assemblies are not coordinate-compatible.

## Repository layout

```
environment.yaml                 conda environment (fennec-genomics)
SRR_Acc_List.txt                 10 SRA run accessions used
scripts/                         pipeline scripts (see above)
variants/
  candidate_regions.bed          the 3 regions analyzed (+/-5000bp padding)
  candidate_regions.depth.txt    samtools depth, all 10 samples, per base
  fennec_candidate_regions.vcf              joint VCF, all sites, FILTER=PASS/LowMQ
  fennec_candidate_regions.unfiltered.vcf   same, before the MQ<40 filter
  fennec_candidate_regions.PASS.vcf         PASS-only subset
  annotation/                    UCSC RepeatMasker annotation for the repeat cluster
  findings_summary.md            full findings writeup
report/
  fennec_candidate_regions_report.docx      same writeup, as a document
```

Not included (see `.gitignore`): raw reads (`.fastq.gz`/`.sra`), the
reference FASTA and its bwa index, and per-sample BAM/BAI files.

## Reproducing

```bash
conda env create -f environment.yaml
conda activate fennec-genomics
conda install -c bioconda -c conda-forge bcftools bedtools entrez-direct

# fetch reads for the accessions in SRR_Acc_List.txt with sra-tools (prefetch + fasterq-dump),
# then:
bash scripts/00_prepare_reference.sh
bash scripts/01_align_bwa_mem.sh
bash scripts/02_extract_candidate_regions.sh
bash scripts/03_depth_and_call_variants.sh
python3 scripts/04_private_shared_analysis.py
```
