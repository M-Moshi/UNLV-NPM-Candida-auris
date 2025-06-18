# UNLV Lab of Neuogenetics and Precision Medicine Amplicon Pipeline

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A robust and parallelized bioinformatics pipeline for processing paired-end amplicon sequencing data (e.g., for SARS-CoV-2 or other targeted sequencing).

This pipeline takes raw FASTQ files and performs trimming, alignment, primer sequence removal, variant calling, and generates comprehensive quality control reports.

## Pipeline Workflow

The script automates the following sequence of analysis steps:

```
Input FASTQ files (R1/R2)
│
└───▶ 1. fastp: Trim adapters and low-quality bases
    │
    └───▶ 2. bwa-mem: Align trimmed reads to a reference genome
        │
        └───▶ 3. samtools: Convert SAM to sorted, indexed BAM
            │
            └───▶ 4. fgbio: Trim primer sequences from alignments
                │
                └───▶ 5. ivar: Call variants and generate consensus
                    │
                    └───▶ 6. samtools/qualimap: Generate QC & coverage metrics
                        │
                        └───▶ 7. Python/Pandas: Aggregate statistics
                            │
                            └───▶ Final Reports (TSV, Qualimap, Logs)
```

## Installation

The recommended method for installation is to use **Conda**, which will handle all software dependencies automatically within a self-contained environment.

**Step 1: Clone the Repository**

First, clone this repository to your local machine.

```bash
git clone https://github.com/moshi321/UNLV-NPM-Candida-auris.git
cd UNLV-NPM-Candida-auris
```

**Step 2: Create and Activate the Conda Environment**

Use the provided `environment.yml` file to create the Conda environment. This command installs all the required tools (`bwa`, `samtools`, `pandas`, etc.).

```bash
# Create the environment (this may take several minutes)
conda env create -f environment.yml

# Activate the environment to use the pipeline
conda activate unlv-amplicon-pipeline
```

> **Note:** You must activate the `unlv-amplicon-pipeline` environment every time you open a new terminal session to run the pipeline.

**Step 3: Make the Pipeline Executable**

Make the main script executable.

```bash
chmod +x Amplicon_pipeline.sh
```

## Usage

Once the environment is activated, you can run the pipeline using the following command structure.

```bash
./Amplicon_pipeline.sh <input_fastq_folder> <output_folder> <reference.fasta> <primers.tab>
```

### Arguments

* `<input_fastq_folder>`: Path to the directory containing paired-end `_R1_001.fastq.gz` and `_R2_001.fastq.gz` files.
* `<output_folder>`: Path to the directory where all results will be saved. It will be created if it doesn't exist.
* `<reference.fasta>`: Path to the reference genome in FASTA format.
* `<primers.tab>`: Path to a tab-separated file defining the primer scheme, required by `fgbio TrimPrimers`.

### Example Command

```bash
# Ensure the conda environment is active
conda activate unlv-amplicon-pipeline

# Run the pipeline
./Amplicon_pipeline.sh \
  ./data/raw_fastqs/ \
  ./results/run01/ \
  ./reference/MN908947.3.fasta \
  ./reference/artic_v3_primers.tsv
```

## Output Structure

The pipeline will generate a structured output directory containing all results:

```
<output_folder>/
├── bam/                # Sorted BAM files after alignment
├── filtered_bam/       # Final BAM files after primer trimming
├── cov/                # Per-sample coverage reports from samtools
├── depth/              # Per-base depth files
├── filtered_tsv/       # Per-sample variant calls from iVar
├── logs/               # Log files for the pipeline run and individual tools
├── qualimap/           # Detailed HTML reports from Qualimap
├── stats/              # Per-sample summary statistics
├── trimmed_fastq/      # FASTQ files after fastp trimming
├── all_filtered.tsv    # Merged TSV of all variants from all samples
└── combined_stats.tsv  # A summary report of key metrics for all samples
```
