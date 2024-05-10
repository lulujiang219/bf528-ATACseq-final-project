# RNA-seq Analysis Documentation

This repository contains materials and methods for analyzing RNA-seq data comparing wild-type (WT) and knockout (KO) samples. The analysis pipeline uses Snakemake to automate tasks from quality control to differential expression analysis, followed by detailed visualization and enrichment analysis using R.

## Overview

The project analyzes six RNA-seq samples, three WT and three KO, derived from a human source. The primary objective is to identify differentially expressed genes between the WT and KO conditions and interpret these results through various bioinformatics tools and techniques.

## Pipeline and Analysis Workflow

### 1. Snakemake Pipeline

The Snakemake pipeline manages the workflow automation for the following tasks:
- **Quality Control (QC):** FastQC and MultiQC are used for sequencing reads and alignment QC.
- **Read Alignment:** Reads are aligned to the human reference genome using the STAR aligner.
- **Count Matrix Generation:** FeatureCounts is used to generate a matrix of read counts per gene.

### 2. R-based Analysis and Visualization

Post-pipeline analysis is conducted in R, focusing on:
- **Differential Expression Analysis:** Utilizing DESeq2 to compare expression levels between WT and KO samples.
- **Visualization:**
  - **PCA Biplot** (`PCA_biplot.png`): Visualizes sample relationships.
  - **Histogram** (`histogram.png`): Shows distribution of log2 fold changes for significantly expressed genes.
  - **Volcano Plot** (`volcano.png`): Distinguishes significant from non-significant DE genes, with emphasis on the top ten most significant genes.
- **Gene Set Enrichment Analysis (GSEA):**
  - Performed using `fgseaMultilevel` and Enrichr to understand the biological implications of differentially expressed genes.

## Deliverables

The deliverables from this analysis include:
1. **QC Reports:** Quality control reports for sequencing reads and alignments.
2. **PCA Biplot:** A visual representation of sample relationships.
3. **DE Analysis Results:** Detailed results provided in `results.csv`.
4. **Histogram:** Distribution of log2FoldChanges for DE genes.
5. **Volcano Plot:** Highlights significant DE genes.
6. **FGSEA Results:** Results from the Functional Gene Set Enrichment Analysis.
7. **Gene Enrichment Analysis:** Outcomes of the gene enrichment analysis using DE genes.

## Installation and Usage

### Prerequisites

Installation of necessary software via Conda:
```bash
conda install -c conda-forge -c bioconda verse
```
```
verse -S -a gencode.vM33.primary_assembly.annotation.gtf -o KOrep2 KOrep2.Aligned.out.bam
```
```
python concat_df.py -i CTLrep1.exon.txt CTLrep2.exon.txt CTLrep3.exon.txt KOrep3.exon.txt KOrep1.exon.txt KOrep2.exon.txt -o output.txt
```
```
python parse_gtf.py -i gencode.vM33.primary_assembly.annotation.gtf -o id2gene.txt
```
