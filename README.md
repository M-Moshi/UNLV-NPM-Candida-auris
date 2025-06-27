# UNLV NPM Candida auris Bioinformatics Workflows

This repository hosts a collection of robust and parallelized bioinformatics pipelines developed at the UNLV Laboratory of Neuogenetics and Precision Medicine for the analysis of *Candida auris* sequencing data. It includes distinct workflows for amplicon, whole-genome (WGS), and RNA-Seq data.

## Project Overview

The primary goal of this repository is to provide a standardized, reproducible, and easy-to-use set of tools for processing various types of Candida auris sequencing data. Each workflow is designed to be self-contained and comes with its own set of instructions and reference files. 
These pipelines were used to generate the data for the NCBI BioProjects listed in the **[Associated Datasets](#associated-datasets)** section at the bottom of this document.

### Directory Structure
```
├── amplicon_sequencing/     # Scripts and README for the amplicon panel pipeline
├── rna_seq/                 # Scripts and README for RNA-Seq & differential abundance analysis
├── whole_genome_sequencing/ # Scripts and README for the WGS MycoSNP-nf pipeline
├── reference_genomes/       # All reference genomes, panel files, and the SnpEff build script
├── environment.yml          # Conda environment file with all software dependencies
└── README.md                # This file
```

---

## Analysis Workflows

This repository contains three distinct analysis workflows. Please refer to the specific `README.md` file within each directory for detailed instructions.

### 1. [Amplicon Sequencing Analysis](./amplicon_sequencing/README.md)

A parallelized pipeline for processing paired-end amplicon sequencing data. This workflow takes raw FASTQ files and performs trimming, alignment, primer trimming, variant calling (iVar), and automated variant annotation using a custom SnpEff database.

### 2. [WGS Analysis using `CDCgov/mycosnp-nf`](./whole_genome_sequencing/README.md)

Instructions for performing whole-genome SNP analysis using the `CDCgov/mycosnp-nf` Nextflow pipeline. This workflow is ideal for phylogenetic analysis and outbreak investigation, and it leverages our custom-built SnpEff databases for *C. auris*.

### 3. [RNA-Seq and Differential Abundance Analysis](./rna_seq/README.md)

A multi-stage workflow to process RNA-Seq data, from raw reads to differential gene expression analysis. This process uses the `nf-core/rnaseq` and `nf-core/differentialabundance` pipelines with specially prepared reference annotation files.

---

## [Reference Genomes & SnpEff Database Management](./reference_genomes/README.md)

This project relies on custom-built reference genomes and SnpEff databases. The `reference_genomes` directory contains all the necessary source files and an automated script (`build_custom_snpeff_dbs.sh`) to build the databases required by the annotation steps in the other pipelines.

### Available References

* **`candida_auris_b11205` & `candida_auris_b11221`**: Primary references for **Whole Genome Sequencing (WGS)** analysis.
* **`candida_auris_b8441`**: The reference genome used for **RNA-Seq analysis**.
* **`C_auris_panel`**: A custom-built reference panel for targeted **amplicon sequencing**.

Before running the WGS or Amplicon pipelines, you must first build these databases. Please see the **[detailed instructions in the `reference_genomes` directory](./reference_genomes/README.md)**.

---

## Installation and Setup

The recommended method for installation is to use **Conda**, which will handle all software dependencies automatically within a self-contained environment.

**Step 1: Clone the Repository**

First, clone this repository to your local machine.

```bash
git clone [https://github.com/m-moshi/UNLV-NPM-Candida-auris.git](https://github.com/m-moshi/UNLV-NPM-Candida-auris.git)
cd UNLV-NPM-Candida-auris
```

**Step 2: Create and Activate the Conda Environment**

Use the provided environment.yml file to create the Conda environment. This command installs all the required tools (bwa, samtools, snpEff, nextflow, pandas, etc.).

```bash
# Create the environment (this may take several minutes)
conda env create -f environment.yml

# Activate the environment to use the pipelines
conda activate npm-candida
```


---

## Associated Datasets

The workflows in this repository were used for the analysis of the following publicly available datasets:

* **Whole Genome Sequencing (WGS) of *Candida auris* isolated from wastewater**
    * [NCBI BioProject: PRJNA1279182](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1279182)

* **RNA sequencing of *Candida auris* isolated from wastewater**
    * [NCBI BioProject: PRJNA1279245](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1279245)

* **Targeted amplicon sequencing of *Candida auris* antifungal resistance genes from wastewater**
    * [NCBI BioProject: PRJNA1279255](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1279255)
