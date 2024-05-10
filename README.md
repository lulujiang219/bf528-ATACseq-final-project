# RNAseq Analysis Documentation

This repository contains materials related to an RNA-seq experiment analysis comparing wild-type (WT) and knockout (KO) samples. Below are the deliverables and instructions for replicating the analysis.

## Overview
This project involves the analysis of six RNAseq samples: three wild-type (WT) and three knockout (KO) derived from a human source. The primary objective is to perform a basic differential expression analysis comparing the WT and KO conditions. Quality control (QC) will be conducted at both the read and alignment levels, followed by the generation of a counts matrix, PCA or sample-to-sample distance plot, and differential expression (DE) analysis. Gene set enrichment analysis (GSEA) will be performed using both a ranked list of all genes and a list of statistically significant DE genes.


## Methods

The analysis of the RNAseq data will begin with quality control of the raw sequencing reads (FASTQ files) for the 6 samples (3 WT and 3 KO) using FastQC v0.11.9 and MultiQC v1.11. Reads will be aligned to the human reference genome (GRCh38, release 104) using STAR v2.7.10a, guided by the GENCODE v42 annotation file. Alignment quality will be assessed using RSeQC v4.0.0 and MultiQC.

A counts matrix will be generated using STAR's gene-level counts. The matrix will be analyzed in R v4.2.2 using DESeq2 v1.38.3. Variance stabilizing transformation will be applied, and sample-to-sample distance plot and PCA will be generated to assess the overall structure of the experiment and identify potential outliers or batch effects.

Differential expression (DE) analysis will be performed using DESeq2, comparing KO to WT samples, with an FDR threshold of 0.05. Gene set enrichment analysis will be conducted using FGSEA v1.24.0 on the ranked gene list based on log2FoldChange values and separately on the list of significant DE genes using enrichR v3.0. The MSigDB v7.5.1 C2 canonical pathways gene sets will be used for FGSEA. Results from both analyses will be compared and discussed, inferring potential biological functions of the factor of interest.

## Deliverables
1. **QC Reports**: Quality control reports for sequencing reads and alignments.
2. **PCA Biplot**: Visualizes sample relationships (`PCA_biplot.png`).
3. **DE Analysis Results**: Comprehensively listed in `results.csv`.
4. **Histogram**: Shows the distribution of log2FoldChanges (`histogram.png`).
5. **Volcano Plot**: Distinguishes significant from non-significant DE genes, highlighting the top ten most significant (`volcano.png`).
6. **FGSEA Results**: Functional Gene Set Enrichment Analysis results presented in `fgsea.png`.
7. **Gene Enrichment Analysis**: Results of the gene enrichment analysis using DE genes (`david.png`).

## Installation and Usage

### Prerequisites
Ensure the necessary tools are installed:
```bash
conda install -c conda-forge -c bioconda verse

verse -S -a gencode.vM33.primary_assembly.annotation.gtf -o KOrep2 KOrep2.Aligned.out.bam
python concat_df.py -i CTLrep1.exon.txt CTLrep2.exon.txt CTLrep3.exon.txt KOrep3.exon.txt KOrep1.exon.txt KOrep2.exon.txt -o output.txt
python parse_gtf.py -i gencode.vM33.primary_assembly.annotation.gtf -o id2gene.txt

