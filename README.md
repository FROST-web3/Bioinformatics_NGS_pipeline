# Bioinformatics NGS Pipeline

## Project Overview
This repository contains a series of bioinformatics pipeline scripts for next-generation sequencing (NGS) data analysis, which I personally developed and validated through extensive practical application in my research group during my graduate studies. These scripts cover multiple aspects from raw data processing to downstream analysis, including Chromatin Immunoprecipitation Sequencing (ChIP-seq), Bisulfite Sequencing (Bisulfite-seq), RNA Sequencing (RNA-seq), and Transcription Start Site (TSS) enrichment analysis.

## Included Pipelines

### TSS Enrichment Heatmap Analysis
The `TSS_enrich_heatMap` directory contains scripts for generating heatmaps of enrichment levels around transcription start sites (TSS) and transcription end sites (TES) using deepTools. It includes:
- `bs_TSS_enrich.sh`: For bisulfite sequencing data visualization
- `chip_TSS_enrich.sh`: For ChIP-seq data visualization

Both scripts generate heatmaps with a 1000bp gene body region, 5000bp upstream region and 3000bp downstream region using the RdBu_r color scheme.

### Bisulfite Sequencing Analysis
The `bisulfite-seq` directory provides comprehensive pipelines for processing and analyzing bisulfite sequencing data for DNA methylation studies:
- `bisulfite-seq_PE.sh`: Pipeline for paired-end bisulfite sequencing data
- `bisulfite-seq_SE.sh`: Pipeline for single-end bisulfite sequencing data

Both scripts handle the complete workflow from raw SRA files to methylation analysis, including:
1. SRA conversion to FASTQ
2. Quality control with fastp
3. Alignment with Bismark 
4. Deduplication
5. Methylation extraction and report generation
6. Conversion to bedGraph and bigWig for visualization
7. Context-specific methylation level calculation (CpG, CHG, CHH)
8. Hypermethylated region identification using MethPipe
9. CSV output generation for further analysis

### ChIP-seq Analysis
Scripts in the `chip-seq` directory process both transcription factor and histone modification ChIP-seq data:
- `chip-seq_PE_histone.sh`: Processing paired-end histone ChIP-seq data
- `chip-seq_PE_TF.sh`: Processing paired-end transcription factor ChIP-seq data  
- `chip-seq_SE_histone.sh`: Processing single-end histone ChIP-seq data
- `chip-seq_SE_TF.sh`: Processing single-end transcription factor ChIP-seq data
- `chip-seq_annotation.sh`: Annotation of peak regions
- `chip-seq_change_name.sh`: Utility for batch renaming files

The pipeline includes:
1. SRA to FASTQ conversion
2. Bowtie2 alignment
3. SAM to BAM conversion with sorting
4. Paired-read processing with fixmate
5. Deduplication
6. Peak calling with MACS2 (with histone-specific or TF-specific parameters)
7. BigWig file generation for visualization

### tDNA Insertion Resequencing
The `resequencing_tDNA_insert` directory contains specialized scripts for analyzing transgenic DNA (tDNA) insertion sites using paired-end (PE) sequencing data:
- `resequencing_PE.sh`: Identifies and localizes tDNA insertion sites in genome

The workflow:
1. Alignment of reads to tDNA reference
2. Identification of reads with soft-clipped regions (indicating genomic junctions)
3. Extraction of soft-clipped sequence portions
4. BLAST alignment of these sequences to the genome
5. Filtering for high-quality matches
6. Conversion to BED format for visualization
7. Statistical analysis of insertion sites

### RNA-seq Analysis
The `rna-seq` directory provides pipelines for RNA sequencing data analysis:
- `rna-seq_PE.sh`: Complete paired-end RNA-seq workflow
- `rna-seq_SE.sh`: Complete single-end RNA-seq workflow
- `quick_rna-seq_PE.sh`: Streamlined pipeline for paired-end RNA-seq

The pipelines include:
1. Quality control with FastQC and MultiQC
2. Adapter and quality trimming with fastp
3. Alignment with HISAT2
4. Sorting and indexing BAM files
5. Read counting with featureCounts
6. BigWig file generation for visualization

### Single-cell Sequencing Analysis
The `singlecell` directory contains notes and code examples for single-cell sequencing data analysis. **Please note that these are notes rather than directly executable pipelines**. They include references for quality control, preprocessing, dimensionality reduction, clustering, and cell type annotation.

### Metagenomics Analysis
The `meta-genomics` directory provides notes and code examples for metagenomic sequencing data analysis. **Please note that these are notes rather than directly executable pipelines**. They include references for taxonomic classification, functional annotation, and comparative analysis tools.

## Usage Instructions

### Environment Requirements
* Bash shell environment (Linux required; other systems have not been tested)
* Bioinformatics tools: All software mentioned in the pipelines can be downloaded via bioconda
  * Core tools: bowtie2, hisat2, samtools, bedtools, deepTools, MACS2, featureCounts, etc.
  * For bisulfite-seq: Bismark, MethPipe
* R statistical environment (for data visualization and downstream analysis)
* Python 3 (required for some analysis scripts)

### How to Use
The scripts in this repository are of two types:
1. **Executable pipelines**: These scripts can be run directly, using customized functions for quick large-scale raw data processing. Most follow a modular design with main() functions that can process multiple directories in batch.
2. **Code notes**: Collections of code snippets similar to work records, convenient for copy-paste usage.

Most scripts require some customization before use:
1. Reference genome paths need to be updated to your local paths
2. Parameters may need adjustment based on your specific experiment design
3. Parallel processing settings can be adjusted to match your computational resources

### Path Configuration
Before running any script, you'll need to modify the reference paths. Look for paths marked with `/path/to/reference/` in the scripts and replace them with your actual paths.

## File Descriptions

### Core Functionality
Detailed description of the main scripts in each directory:

* **TSS_enrich_heatMap/**
   * `bs_TSS_enrich.sh`: Takes bedGraph files from methylation analysis and creates TSS enrichment heatmaps, specifically normalized for bisulfite data
   * `chip_TSS_enrich.sh`: Takes bigWig files from ChIP-seq analysis and creates TSS enrichment heatmaps
   
* **bisulfite-seq/**
   * `bisulfite-seq_PE.sh`: Complete pipeline for paired-end bisulfite sequencing data
   * `bisulfite-seq_SE.sh`: Complete pipeline for single-end bisulfite sequencing data

* **chip-seq/**
   * `chip-seq_PE_histone.sh`: Processing paired-end data with parameters optimized for histone marks (broad peaks)
   * `chip-seq_PE_TF.sh`: Processing paired-end data with parameters for transcription factors (narrow peaks)
   * `chip-seq_SE_histone.sh`: Single-end version for histone marks
   * `chip-seq_SE_TF.sh`: Single-end version for transcription factors
   * `chip-seq_annotation.sh`: Annotates called peaks with genomic features
   * `chip-seq_change_name.sh`: Utility script for batch renaming files

* **resequencing_tDNA_insert/**
   * `resequencing_PE.sh`: Processes PE sequencing data to identify genomic tDNA insertion sites using soft-clipped reads

* **rna-seq/**
   * `rna-seq_PE.sh`: Full pipeline for paired-end RNA-seq with extensive QC and visualization
   * `rna-seq_SE.sh`: Full pipeline for single-end RNA-seq
   * `quick_rna-seq_PE.sh`: Streamlined pipeline focusing on essential steps for faster processing

* **singlecell/**
   * Notes and code examples for single-cell RNA-seq analysis (not executable pipelines)

* **meta-genomics/**
   * Notes and code examples for metagenomic data analysis (not executable pipelines)

## Frequently Asked Questions

* **Q: Which parameters need to be modified in the scripts?** A: Most scripts require modification of input/output paths and reference genome paths. All pathways marked with `/path/to/reference/` should be updated to your actual reference paths.

* **Q: How do I install the dependent tools?** A: It's recommended to use conda/bioconda for installation. Example:
  ```
  conda create -n ngs-analysis
  conda activate ngs-analysis
  conda install -c bioconda bowtie2 hisat2 samtools bedtools macs2 subread
  conda install -c bioconda bismark methpipe deeptools
  ```

* **Q: Can I run these scripts on a high-performance computing cluster?** A: Yes, the scripts use parallel processing and can be adapted for cluster environments by adjusting the number of cores and memory parameters.

* **Q: What reference genomes are these scripts designed for?** A: The scripts were originally designed for plant genome analysis (specifically tomato), but can be adapted to any organism by changing the reference genome paths and adjusting parameters like genome size.

## Contributions and Feedback
Questions or suggestions are welcome through Issues. If you're interested in improving these scripts, Pull Requests are also welcome. Areas for improvement include:
- Adding more detailed documentation
- Parameterizing reference paths
- Optimizing parallel processing
- Adding additional downstream analysis options

## License
This project is licensed under the MIT License - see the LICENSE file for details.

## Citation
If you use these scripts in your research, please consider citing this repository: FROST-web3. (2025). Bioinformatics_NGS_pipeline: A collection of validated NGS analysis pipelines. GitHub.
