# Bioinformatics NGS Pipeline

## Project Overview
This repository contains a series of bioinformatics pipeline scripts for next-generation sequencing (NGS) data analysis, which I personally developed and validated through extensive practical application in my research group during my graduate studies. These scripts cover multiple aspects from raw data processing to downstream analysis, including Chromatin Immunoprecipitation Sequencing (ChIP-seq), Bisulfite Sequencing (Bisulfite-seq), RNA Sequencing (RNA-seq), and Transcription Start Site (TSS) enrichment analysis.

## Included Pipelines

### TSS Enrichment Heatmap Analysis
The `TSS_enrich_heatMap` directory contains scripts for generating heatmaps of enrichment levels around transcription start sites (TSS) and transcription end sites (TES). This is particularly useful for evaluating the quality of ChIP-seq or ATAC-seq data and visualizing protein binding patterns in promoter regions.

### Bisulfite Sequencing Analysis
The `bisulfite-seq` directory provides scripts for processing and analyzing bisulfite sequencing data for DNA methylation studies.

### ChIP-seq Histone Analysis
Scripts in the `chip-seq` directory are specifically designed for processing ChIP-seq data, with a focus on histone modification analysis.

### tDNA Insertion Resequencing
The `resequencing_tDNA_insert` directory contains scripts for analyzing transgenic DNA (tDNA) insertion sites using paired-end (PE) sequencing data.

### RNA-seq Analysis
The `rna-seq` directory provides a complete pipeline for RNA sequencing data analysis for gene expression studies.

### Single-cell Sequencing Analysis
The `singlecell` directory contains analysis pipelines for single-cell sequencing data.

### Metagenomics Analysis
The `meta-genomics` directory provides resources for analyzing metagenomic sequencing data.

## Directory Structure Note
**Important**: This repository is transitioning to use consistent lowercase naming for all directories. Please use these current folders:
- Use `singlecell` instead of `Singlecell`
- Use `meta-genomics` instead of `Meta-genomics`

## Usage Instructions

### Environment Requirements
* Bash shell environment (I use Linux; other systems have not been tested)
* Bioinformatics tools: All software mentioned in the pipelines can be downloaded via bioconda
* R statistical environment (for data visualization)
* Python 3

### How to Use
The scripts in this repository are of two types:
1. **Executable programs**: These scripts can be run directly, using customized functions for quick large-scale raw data processing
2. **Code notes**: Collections of code snippets similar to work records, convenient for copy-paste usage

Please check the comments and header descriptions in each script for specific usage instructions. Most scripts include detailed parameter descriptions and usage examples.

## File Descriptions

### Core Functionality
Brief description of the main scripts in each directory:

* **TSS_enrich_heatMap/**
   * Generate signal heatmaps for gene TSS/TES regions ±x kb
* **bisulfite-seq/**
   * Complete pipeline from raw FASTQ files to methylation site analysis
* **chip-seq/**
   * Scripts specifically for histone modification ChIP-seq analysis
* **resequencing_tDNA_insert/**
   * Precise localization and identification of transgenic insertion sites
* **rna-seq/**
   * Complete pipeline from sequencing data to differential expression analysis

## Frequently Asked Questions

* **Q: Which parameters need to be modified in the scripts?** A: Most scripts require modification of input/output paths and reference genome paths. I was lazy and didn't set these paths as variables at the beginning, but direct modification isn't too troublesome.
* **Q: How do I install the dependent tools?** A: It's recommended to use conda/bioconda for installation, for example: `conda install -c bioconda bwa samtools bedtools`

## Contributions and Feedback
Questions or suggestions are welcome through Issues. If you're interested in improving these scripts, Pull Requests are also welcome.

## License
This project is licensed under the MIT License - see the LICENSE file for details.

## Citation
If you use these scripts in your research, please consider citing this repository: FROST-web3. (2025). Bioinformatics_NGS_pipeline: A collection of validated NGS analysis pipelines. GitHub.
