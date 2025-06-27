# Reference Genomes

This directory contains the reference genome and annotation files required for the various analysis pipelines in this project. Each subdirectory corresponds to a specific reference genome or panel used for a distinct purpose.

---

### `candida_auris_b11205` and `candida_auris_b11221`

-   **Use Case**: These two directories contain the primary reference genomes for **Whole Genome Sequencing (WGS)** analysis.
-   **Key Files**: Each contains a genome sequence (**`.fna`**) and a corresponding rich annotation file in GenBank format (**`.gbk`**). These are used for variant calling and annotation in the WGS pipeline.

---

### `candida_auris_b8441`

-   **Use Case**: This reference genome is specifically used for **RNA-Seq analysis**, such as differential gene expression.
-   **Key Files**: Contains the genome sequence (**`.fasta`**) and gene models in GFF/GTF format (**`.gtf`**). The GTF format is essential for correctly quantifying transcript abundance.

---

### `C_auris_panel`

-   **Use Case**: This is a custom-built reference used for a targeted **amplicon panel analysis**. Instead of whole chromosomes, it contains specific gene or marker sequences of interest.
-   **Key Files**: Contains the panel sequences (**`sequences.fa`**) and a corresponding simple annotation in GFF format (**`genes.gff`**) that defines the regions of interest within the panel.


# Automated SnpEff Database Management

This document explains how to use the automated `build_custom_snpeff_dbs.sh` script to prepare all necessary SnpEff databases for this project. This script replaces the need for manual setup.

## Overview

The `build_custom_snpeff_dbs.sh` script is designed to automate the creation of multiple SnpEff databases from all genome folders located in this `reference_genomes` directory. It centralizes the SnpEff configuration and data into a single `reference_genomes/snpEff` directory, ensuring consistency and ease of use.

The script intelligently handles different annotation formats (GenBank and GFF) and applies specific configurations for special cases like the `C_auris_panel`.

You should run this script once during the initial project setup and again anytime you add a new reference genome.

## How to Build the Databases

The script requires you to be in the `reference_genomes` directory. This is the most reliable way to ensure SnpEff can always find its files.

**Step 1: Get the Absolute Path**

First, navigate into this directory.

```bash
cd /path/to/your/project/reference_genomes

# Make sure your conda environment is active
conda activate npm-candida

# Make build_custom_snpeff_dbs.sh executable
chmod +x build_custom_snpeff_dbs.sh

# Run the script from this directory
./build_custom_snpeff_dbs.sh
```





