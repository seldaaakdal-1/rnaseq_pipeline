# RNA-seq Analysis Pipeline

This repository contains an end-to-end RNA-seq analysis workflow including preprocessing, alignment, gene quantification, differential expression analysis, and functional enrichment.

## Bash Pipeline
- SRA data download (prefetch, fasterq-dump)
- Quality control (FastQC, MultiQC)
- Read trimming (Trimmomatic)
- Genome alignment (STAR)
- Gene quantification (featureCounts)
- Processing of featureCounts output to generate a clean count matrix (`gene_count_cutted.txt`) for downstream analysis

## R Analysis
- Differential expression analysis (DESeq2)
- Data visualization (PCA, MA plot, volcano plot, heatmap)
- Functional enrichment analysis (GO and KEGG)

## Tools
- Bash, STAR, FastQC, MultiQC, Trimmomatic
- R (DESeq2, clusterProfiler, pheatmap, ggplot2)
