# WGBS Processing Pipeline with Bismark (Nextflow)

This Nextflow pipeline automates the processing of paired-end Whole Genome Bisulfite Sequencing (WGBS) data using the Bismark suite. It handles steps from raw FASTQ files to methylation calling.

## Table of Contents

1.  [Overview](#overview)
2.  [Features](#features)
3.  [Requirements](#requirements)
4.  [Setup](#setup)
    *   [Pipeline Code](#pipeline-code)
    *   [Input Data](#input-data)
    *   [Bismark Genome Preparation](#bismark-genome-preparation)
5.  [Usage](#usage)
    *   [Input Samplesheet](#input-samplesheet)
    *   [Running the Pipeline](#running-the-pipeline)
6.  [Parameters](#parameters)
    *   [Mandatory Parameters](#mandatory-parameters)
    *   [Optional Parameters](#optional-parameters)
    *   [Trim Galore Parameters](#trim-galore-parameters)
    *   [Bismark Parameters](#bismark-parameters)
    *   [Methylation Extractor Parameters](#methylation-extractor-parameters)
    *   [QC Parameters](#qc-parameters)
7.  [Output Directory Structure](#output-directory-structure)
8.  [Dependency Management](#dependency-management)
9.  [Troubleshooting](#troubleshooting)

## Overview

The pipeline performs the following key steps for paired-end WGBS data:
1.  **FASTQ Merging:** Concatenates FASTQ files from different lanes for each sample (R1 and R2).
2.  **Quality Control (Optional Raw):** Runs FastQC on merged raw FASTQ files.
3.  **Adapter & Quality Trimming:** Uses Trim Galore for adapter removal, quality trimming, and optional 5'/3' end clipping. Trim Galore also runs FastQC on trimmed reads.
4.  **Alignment:** Aligns trimmed reads to a bisulfite-converted reference genome using Bismark (with Bowtie2 or HISAT2 as the aligner).
5.  **Deduplication:** Removes PCR duplicates using `deduplicate_bismark`.
6.  **Methylation Extraction:** Extracts methylation calls from the deduplicated BAM files using `bismark_methylation_extractor`, generating various report formats including BedGraph and coverage files.
7.  **Aggregate QC:** Generates a MultiQC report summarizing QC metrics from various steps.

## Features

*   **Automated Workflow:** Streamlines the entire WGBS analysis from raw reads to methylation calls.
*   **Paired-End Support:** Specifically designed for paired-end WGBS data.
*   **Parallelization:** Leverages Nextflow for efficient parallel processing of multiple samples and steps.
*   **Reproducibility:** Uses Conda for managing software dependencies, ensuring consistent environments.
*   **Configurable:** Offers various parameters to customize trimming, alignment, and methylation extraction steps.
*   **Comprehensive QC:** Integrates FastQC and MultiQC for thorough quality assessment.

## Requirements

*   **Nextflow:** (version 21.10.x or later recommended).
*   **Conda:** To manage software dependencies.
*   **Bismark-Prepared Genome:** The reference genome must be prepared using `bismark_genome_preparation` **before** running this pipeline (see [Bismark Genome Preparation](#bismark-genome-preparation)).
*   **Required Tools (installed via Conda by the pipeline):**
    *   `fastqc`
    *   `multiqc`
    *   `trim-galore` (and its dependency `cutadapt`)
    *   `bismark` (which requires `samtools` and an aligner: `bowtie2` or `hisat2`)

## Setup

### 1. Pipeline Code
Clone or download this repository/pipeline script (`wgbs_pipeline.nf`).

```bash
# If it's in a git repository:
# git clone <repository_url>
# cd <repository_directory>
```

### 2. Input Data
*   Place your gzipped paired-end FASTQ files (R1 and R2) in a directory. The path to this directory will be specified by `params.input_dir`.
*   The pipeline uses a samplesheet (see [Input Samplesheet](#input-samplesheet)) to identify these files.

### 3. Bismark Genome Preparation
**This is a critical prerequisite.** Before running the pipeline, you must prepare your reference genome using Bismark.

```bash
# Example for preparing a genome with Bowtie2 (default aligner for this pipeline)
bismark_genome_preparation --bowtie2 /path/to/your/reference_fasta_directory/

# Example for HISAT2
# bismark_genome_preparation --hisat2 /path/to/your/reference_fasta_directory/
```
This command will create a new directory (often named `Bisulfite_Genome` within your reference FASTA directory). The path to this *Bismark-prepared genome directory* is what you need for the `params.bismark_genome_folder` parameter.

## Usage

### 1. Input Samplesheet
Create a CSV file (e.g., `samples_wgbs.csv`) with the following header and structure:

```csv
sample,fastq_r1_glob,fastq_r2_glob
Tag-2,Tag-2_S1_L00*_R1_001.fastq.gz,Tag-2_S1_L00*_R2_001.fastq.gz
Tag-3,Tag-3_S2_L00*_R1_001.fastq.gz,Tag-3_S2_L00*_R2_001.fastq.gz
Tag-6,Tag-6_S3_L00*_R1_001.fastq.gz,Tag-6_S3_L00*_R2_001.fastq.gz
# ... and so on for all your samples
```
*   `sample`: A unique identifier for your sample. This will be used for naming output files.
*   `fastq_r1_glob`: A glob pattern relative to `params.input_dir` that matches the R1 FASTQ file(s) for that sample.
*   `fastq_r2_glob`: A glob pattern relative to `params.input_dir` that matches the R2 FASTQ file(s) for that sample.

### 2. Running the Pipeline
Execute the pipeline using the `nextflow run` command:

```bash
nextflow run wgbs_pipeline.nf \
    --samplesheet samples_wgbs.csv \
    --input_dir /path/to/your/fastq/directory \
    --bismark_genome_folder /path/to/your/Bismark_Prepared_Genome_Folder \
    --outdir ./results_wgbs_runX \
    # Add other parameters as needed (see Parameters section)
    # e.g., --bismark_aligner_cores 8
```
You can also use a `nextflow.config` file to specify parameters and process configurations (e.g., executor, CPU, memory).

## Parameters

### Mandatory Parameters
*   `--samplesheet`: Path to the input samplesheet CSV file (e.g., `samples_wgbs.csv`).
*   `--input_dir`: Path to the directory containing raw FASTQ files (e.g., `data/fastq_wgbs`).
*   `--bismark_genome_folder`: Path to the Bismark-prepared genome directory (e.g., `genome/Bisulfite_Genome`).

### Optional Parameters
*   `--outdir`: Directory where results will be published (Default: `results_wgbs`).

### Trim Galore Parameters
*   `--clip_r1`: Number of bases to clip from the 5' end of R1 reads (Default: `15`). Set to `0` to disable.
*   `--clip_r2`: Number of bases to clip from the 5' end of R2 reads (Default: `15`). Set to `0` to disable.
*   `--trim_galore_cores`: Number of cores for Trim Galore's internal Cutadapt calls (Default: `2`).

### Bismark Parameters
*   `--bismark_non_directional`: Specify if the library is non-directional (Default: `true`). Set to `false` for directional libraries.
*   `--bismark_aligner`: Aligner to be used by Bismark (`bowtie2` or `hisat2`) (Default: `bowtie2`).
*   `--bismark_aligner_cores`: Number of cores allocated for the aligner (e.g., Bowtie2/HISAT2 threads). Bismark uses `--parallel` and aligner's `--threads` based on this. (Default: `4`).

### Methylation Extractor Parameters
*   `--meth_extractor_no_overlap`: For paired-end reads, ignore overlapping parts of alignments (Default: `true`).
*   `--meth_extractor_paired_end`: Specify that reads are paired-end (Default: `true`). Should generally remain true for this WGBS pipeline.
*   `--meth_extractor_cytosine_report`: Output genome-wide cytosine methylation report (Default: `true`).
*   `--meth_extractor_bedgraph`: Output BedGraph file (Default: `true`).
*   `--meth_extractor_zero_based`: Use 0-based coordinates for BedGraph output (Default: `true`).
*   `--meth_extractor_multicore`: Number of cores for `bismark_methylation_extractor` sub-processes (Default: `2`).

### QC Parameters
*   `--skip_fastqc_raw`: Skip FastQC on raw merged FASTQ files (Default: `false`).
*   `--skip_multiqc`: Skip MultiQC report generation (Default: `false`).

## Output Directory Structure

The pipeline will create an output directory (specified by `params.outdir`, default `results_wgbs/`) with the following structure:

```
results_wgbs/
├── 01_merged_fastq/            # Merged R1 and R2 FASTQ files per sample
├── 02_fastqc_raw/              # (If not skipped) FastQC reports for raw merged FASTQs
├── 03_trimmed_fastq/           # Trimmed R1 and R2 FASTQ files
│   └── reports/                # TrimGalore reports and FastQC reports for trimmed reads
├── 04_bismark_alignment/       # Bismark alignment BAM files and reports
├── 05_deduplicated_bam/        # Deduplicated BAM files and deduplication reports
├── 06_methylation_calls/       # Methylation extraction results (BedGraph, coverage, M-bias, etc.)
└── multiqc/                    # MultiQC report aggregating QC results
    ├── multiqc_report.html
    └── multiqc_data/
```

## Dependency Management

Software dependencies are managed using Conda. Each process in the Nextflow script specifies its required Conda environment. Nextflow will automatically create (or use cached) these environments.
Ensure Conda is installed and accessible in your `PATH`.

## Troubleshooting

*   **Bismark Genome Not Found/Incorrect:** Ensure `params.bismark_genome_folder` points to the directory created by `bismark_genome_preparation` and that this preparation step was successful.
*   **Conda Environment Issues:** Check your Conda setup and internet connection. You can try creating the environments manually first if issues persist.
*   **File/Path Issues:** Double-check all input file paths and glob patterns in your samplesheet. Nextflow is case-sensitive.
*   **Resource Limits:** WGBS alignment, especially with Bismark, can be resource-intensive (CPU and memory). Adjust process configurations in a `nextflow.config` file if jobs fail due to resource limitations. For `BISMARK_ALIGN`, consider the memory requirements of the chosen aligner for your reference genome size.
*   **MultiQC Report Empty/Missing Sections:** Ensure that the upstream processes generating reports (FastQC, TrimGalore, Bismark) completed successfully and produced the expected output files.
```
