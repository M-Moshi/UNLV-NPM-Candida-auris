# Reference Genomes & Database Setup

This directory contains the reference genomes, annotation files, and custom panels required for the various analysis pipelines in the UNLV NPM *Candida auris* project. Each subdirectory corresponds to a specific reference genome or amplicon panel used for a distinct analytical purpose.

---

## Available Reference Genomes

### WGS References: `candida_auris_b11205` & `candida_auris_b11221`
- **Use Case:** These directories contain the primary reference genomes for the **Whole Genome Sequencing (WGS)** phylogenetic and variant analysis.
- **Key Files:**
  - `*.fna` : The raw nucleotide genome sequence.
  - `*.gbk` : A rich annotation file in GenBank format, utilized for accurate downstream variant calling and automated annotation.

### RNA-Seq Reference: `candida_auris_b8441`
- **Use Case:** This reference genome is specifically optimized for **RNA-Seq and Differential Abundance Analysis**.
- **Key Files:**
  - `*.fasta` : The nucleotide genome sequence.
  - `*.gtf` : Gene models provided in GTF format, which is strictly required for accurately quantifying transcript abundance and mapping read features.

### Amplicon Panel: `C_auris_panel`
- **Use Case:** A custom-built reference dedicated to targeted **amplicon panel analysis**. Instead of a whole genome approach, this relies on specific gene and marker sequences of interest (e.g., antifungal resistance markers).
- **Key Folders & Files:**
  - `reference_genes/` : A subdirectory containing the individual FASTA files for all the specific genes that were combined to construct the single panel reference.
  - `C_auris_panel.fasta` : The combined nucleotide sequences of the targeted panel.
  - `C_auris.csv` : The output primer generation file. This contains the raw primer sequences, amplicon coordinates, pool assignments, and lengths used to design the panel.
  - `C_auris_panel.primer.tab` : The primer coordinate file specifically formatted for primer sequence soft-clipping/trimming using **`fgbio TrimPrimers`**.
  - `C_auris_panel.gff` : Genomic feature annotations defining the precise regions of interest within the extracted panel.
  - `C_auris_panel.tsv` : A tab-separated bed-style file providing specific genomic coordinates, utilized for calculating depth and gene coverage.
---

## Automated SnpEff Database Management

This repository utilizes an automated bash script (`build_custom_snpeff_dbs.sh`) to prepare and configure all necessary SnpEff databases. This eliminates the need for manual configuration and ensures high reproducibility across environments.

### Overview
The `build_custom_snpeff_dbs.sh` script automates the creation of multiple SnpEff databases from the genome folders located in this directory. 
* It centralizes the SnpEff configuration and resulting data into a single `reference_genomes/snpEff` directory.
* It intelligently parses different annotation formats (handling both GenBank and GFF setups).
* It applies specialized custom configurations for targeted references like the `C_auris_panel`.

**When to run:** You should execute this script once during the initial project setup, or whenever a new reference genome is added to the repository.

### Build Instructions

**Pre-requisite:** Ensure your conda environment containing `snpEff` is activated before running the script. The script relies on relative paths, so you **must** execute it directly from within the `reference_genomes` directory.

```bash
# 1. Navigate to the reference genomes directory
cd /path/to/your/project/reference_genomes

# 2. Activate the project conda environment
conda activate npm-candida

# 3. Ensure the script has execution permissions
chmod +x build_custom_snpeff_dbs.sh

# 4. Execute the database build script
./build_custom_snpeff_dbs.sh
