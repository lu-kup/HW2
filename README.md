# Home work assignment II - WGBS data analysis

## Project Description

This project contains the analysis workflow and results for a Whole Genome Bisulfite Sequencing (WGBS) homework assignment focused on identifying DNA methylation differences between esophageal squamous cell carcinoma (ESCC) and matched normal esophageal tissue samples. The project performs a complete methylation analysis workflow from raw FASTQ files to functional enrichment analysis and biological interpretation.

The analysis results can be reproduced by running the following scripts in order: download_script.sh, script.sh, dataQC.R, collect_data.R, dmc.R.

---

# Dataset Description

## Samples

| Sample ID | Description |
|---|---|
| GSM4505857 | Normal esophagus tissue (R2) |
| GSM4505859 | Normal esophagus tissue (R7) |
| GSM4505867 | Tumor tissue (R2) |
| GSM4505869 | Tumor tissue (R7) |

## Sequencing Information

- Sequencing type: Whole Genome Bisulfite Sequencing (WGBS)
- Platform: Illumina HiSeq2500
- Library type: Paired-end
- Read length: 2 × 125 bp

## Public Data Source

The raw sequencing data was obtained from European Nucleotide Archive at https://www.ebi.ac.uk/.

Original study:
- Cao W et al., "Multi-faceted epigenetic dysregulation of gene expression promotes esophageal squamous cell carcinoma.", Nat Commun, 2020 Jul 22;11(1):3675

---

# Project Structure

```text
HW2/
├── raw_data/              # Raw FASTQ files (not committed to repository)
├── references/            # Reference genome files (not committed to repository)
├── reads_trimmed/         # Intermediate trimmed FASTQ files (not committed to repository)
├── qc/                    # Outputs related to the quality control steps
├── qc_trimmed/            # Outputs related to the quality control of intermediate trimmed files
├── outputs/               # All non-graphical outputs of the workflow
├── diagrams/              # Graphical outputs - generated plots and figures
├── download_script.sh     # Commands used in 'Getting datasets' part of assignment
├── script.sh              # Commands used in 'Raw data QC and methylation calling' part of assignment
├── dataQC.R               # Code used in 'Data QC and statistical analysis' part of assignment
├── collect_data.R         # Code used to prepare intermediate files methylation_dt, coverage_dt, sample_info (based on examples provided online)
└── dmc.R                  # Code used in 'DMC calling', 'DMR calling', "Biological interpretation' parts of the analysis
```

---

# Software and Packages Used

## Command-Line Tools

| Tool | Purpose |
|---|---|
| FastQC | Raw read quality control |
| Trim Galore | Adapter and quality trimming |
| BWA-Meth | Bisulfite-aware read alignment |
| Samtools | BAM/SAM processing |
| Picard MarkDuplicates | Duplication analysis |
| MethylDackel | Methylation calling |

## R Packages

| Package | Purpose |
|---|---|
| DSS | DMC and DMR calling |
| methylKit | Methylation, coverage analysis |
| ggplot2 | Data visualization |
| pheatmap | Heatmaps |
| VennDiagram | Venn diagrams |
| clusterProfiler | Functional enrichment |
| org.Hs.eg.db | Human gene annotation |
| GenomicRanges | Genomic interval handling |
| annotatr | Genomic annotation |
| ComplexHeatmap | Advanced heatmaps |
| dplyr | Data manipulation |

## Reference Genome

- Human genome assembly: GRCh38

