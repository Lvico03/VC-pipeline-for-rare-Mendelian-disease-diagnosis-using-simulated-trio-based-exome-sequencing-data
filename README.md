# A Variant Calling Pipeline for Rare Mendelian Disease Diagnosis Using Simulated Trio-Based Exome Sequencing Data

Bioinformatics for Computational Genomics (BCG) course project — MSc in Bioinformatics for Computational Genomics, University of Milan – Politecnico di Milano.

- **Authors:** Francesco Matassa, Ludovico Salmin (both authors contributed equally)

## Overview

Trio-based exome sequencing — sequencing a proband (child) together with both parents — enables identification of disease-causing variants through inheritance-based filtering, and is an effective diagnostic strategy for rare Mendelian diseases. This project implements a complete bioinformatics pipeline for diagnosing simulated trio-based exome sequencing data, restricted to chromosome 20 and aligned to GRCh38.

Ten trios derived from 1000 Genomes Project samples were analyzed (five trios per student), each with disease-causing variants simulated in silico under a specific mode of inheritance:

- **Autosomal recessive (AR):** trios 1, 4, 5
- **Autosomal dominant, inherited (AD inh., mother affected):** trio 2
- **Autosomal dominant, de novo (AD dn):** trio 3

## Pipeline

1. **Read alignment** — Bowtie2 (paired-end, default parameters), with read groups assigned per sample
2. **BAM processing** — SAMtools (conversion, coordinate sorting, indexing)
3. **Quality control** — FastQC per BAM, Qualimap per sample against the exome target BED file, aggregated with MultiQC
4. **Variant calling** — joint calling per trio with FreeBayes (`-m 20 -C 5 -Q 10 --min-coverage 10`), output compressed with bgzip and indexed with bcftools
5. **Inheritance-based filtering** — bcftools, genotype filters tailored to each trio's mode of inheritance, restricted to exome target regions, `QUAL > 20`
6. **Annotation** — Ensembl VEP (GRCh38, MANE Select transcripts, one consequence per variant), with gnomAD exomes and 1000 Genomes allele frequencies
7. **Final filtering** — `filter_vep`, retaining HIGH/MODERATE impact variants with `MAX_AF < 0.0001` or not reported
8. **Interpretation** — candidate genes cross-referenced against the simulated disorder list, OMIM, and ClinVar
9. **Coverage visualization** — per-sample BEDGraph tracks via `bedtools genomecov`, visualized on the UCSC Genome Browser

## Results

All 30 samples (3 per trio × 10 trios) aligned at ~100%, with mean coverage of ~38–42x across all samples. FreeBayes identified ~9,000–9,100 raw variant sites per trio on chromosome 20; inheritance and quality filtering reduced this to ~100–260 candidates per trio, and VEP impact/frequency filtering narrowed this further to 1–2 high-confidence candidates in affected trios.

Diagnostic yield was 8/10 trios (80%), consistent with the ~20% of cases simulated without a disease-causing variant.

| Student / Trio | Mode | Gene | Variant type | Diagnosis |
|---|---|---|---|---|
| Matassa, Trio 1 | AR | – | – | Unaffected |
| Matassa, Trio 2 | AD inh. | *KCNQ2* | Frameshift insertion | Seizures Benign Familial Neonatal 1 / DEE 7 |
| Matassa, Trio 3 | AD dn | *ASXL1* | Frameshift insertion | Bohring-Opitz Syndrome |
| Matassa, Trio 4 | AR | *SLC2A10* | Frameshift deletion | Arterial Tortuosity Syndrome |
| Matassa, Trio 5 | AR | *GSS* | Stop gained (p.Q308*) | Glutathione Synthetase Deficiency |
| Salmin, Trio 1 | AR | – | – | Unaffected |
| Salmin, Trio 2 | AD inh. | *KCNQ2* | Frameshift deletion | Seizures Benign Familial Neonatal 1 / DEE 7 |
| Salmin, Trio 3 | AD dn | *SLC12A5* | Frameshift deletion | Developmental and Epileptic Encephalopathy 34 |
| Salmin, Trio 4 | AR | *SLC52A3* | Stop gained (p.E249*) | Brown-Vialetto-Van Laere Syndrome 1 |
| Salmin, Trio 5 | AR | *SLC52A3* | Stop gained (p.E249*) | Brown-Vialetto-Van Laere Syndrome 1 |

## Repository structure

```
.
├── cands/              # Candidate variants after inheritance-based genotype + quality filtering
├── extra_info/         # Supplementary reference/annotation files (e.g. target BED, sample list)
├── final_report/        # Final written report (PDF)
├── script/              # Full pipeline script (alignment through coverage tracks)
├── trio_cov/             # Per-sample BEDGraph coverage tracks
├── vcf_filtered/         # Inheritance- and quality-filtered VCFs
├── vepcmd_filtered/       # VEP-annotated and filter_vep-filtered variants (command-line VEP)
└── vepweb_filtered/       # VEP web-interface annotation/verification output
```

## Tools used

Bowtie2 · SAMtools · FastQC · Qualimap · MultiQC · FreeBayes · bcftools · Ensembl VEP · filter_vep · bedtools · UCSC Genome Browser

## Data availability

Raw reads were simulated from real 1000 Genomes Project genotypes and are not included in this repository. Intermediate files (BAMs, VCFs, annotated VCFs, QC reports, coverage tracks) were generated and stored on the course's `leon` compute cluster during the exam.

## References

Key methods references: Langmead & Salzberg (Bowtie2), Li et al. (SAMtools), García-Alcalde et al. (Qualimap), Ewels et al. (MultiQC), Garrison & Marth (FreeBayes), McLaren et al. (Ensembl VEP), Quinlan & Hall (BEDTools). Full citation list in the final report.
