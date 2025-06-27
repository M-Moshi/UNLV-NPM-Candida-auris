# UNLV Lab of Neuogenetics and Precision Medicine Amplicon Pipeline


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
                        └───▶ 7. Annotation: Coversion of iVar to VCF and annotation using snpEff and snpSift
                            │
                            └───▶ 8. Python: Aggregate statistics
                                │
                                └───▶ Final Reports (TSV, Qualimap, Logs)
```

## Installation

The recommended method for installation is to use **Conda**, which will handle all software dependencies automatically within a self-contained environment.

**Step 1: Clone the Repository**

First, clone this repository to your local machine.

```bash
git clone https://github.com/m-moshi/UNLV-NPM-Candida-auris.git
cd UNLV-NPM-Candida-auris
```

**Step 1: Make the Pipeline Executable**

Make the main script executable.

```bash
chmod +x /amplicon_sequencing/Amplicon_pipeline.sh
```

**Step 2: Run Amplicon Pipeline**

## Usage

Once the environment is activated, you can run the pipeline using the following command structure.

```bash
./amplicon_sequencing/Amplicon_pipeline.sh <input_fastq_folder> <output_folder> <reference.fasta> <primers.tab> <genes.tsv>
```

### Arguments

* `<input_fastq_folder>`: Path to the directory containing paired-end `_R1_001.fastq.gz` and `_R2_001.fastq.gz` files.
* `<output_folder>`: Path to the directory where all results will be saved. It will be created if it doesn't exist.
* `<reference.fasta>`: Path to the reference genome in FASTA format.
* `<primers.tab>`: Path to a tab-separated file defining the primer scheme, required by `fgbio TrimPrimers`.
* `<genes.tsv>`: Path to a tab-separated file defining the coordinates of the genes within the panel (ex. `C_auris_panel.tsv`).

### Example Command

```bash
# Ensure the conda environment is active
conda activate npm-candida

# Run the pipeline
./amplicon_sequencing/Amplicon_pipeline.sh \
  ./raw_fastqs/ \
  ./results/ \
  ./reference_genomes/C_auris_panel/C_auris_panel.fasta \
  ./reference_genomes/C_auris_panel/C_auris_panel.primers.tab \
  ./reference_genomes/C_auris_panel/C_auris_panel.tsv
```

**Step 3: Annotate Amplicon Variantse**

### Usage: ./snpeff.sh <pipeline_results_dir> <genome_name> <config_path>

```bash
./amplicon_sequencing/snpeff.sh results/ C_auris_panel reference_genomes/snpEff/snpEff.config
```




## Output Structure

The pipeline will generate a structured output directory containing all results:

```
<output_folder>/
├── annotation/         # Annotated VCFs 
├── bam/                # Sorted BAM files after alignment
├── cov/                # Per-sample coverage reports from samtools
├── coverage_analysis/  # Gene coverage report sheet
├── depth/              # Per-base depth files
├── filtered_bam/       # Final BAM files after primer trimming
├── filtered_tsv/       # Per-sample variant calls from iVar
├── histogram/          # Histogram of coverage across reference fasta
├── logs/               # Log files for the pipeline run and individual tools
├── qualimap/           # Detailed HTML reports from Qualimap
├── rawtsv/             # Unfiltered raw TSV variant calls from iVar
├── sam/                # Raw SAM file from BWA-MEM
├── stats/              # Per-sample summary statistics
├── trimmed_fastq/      # FASTQ files after fastp trimming
├── tsv/                # Final, combined reports for each sample
└── vcf/                # VCF files created during annotation
```
