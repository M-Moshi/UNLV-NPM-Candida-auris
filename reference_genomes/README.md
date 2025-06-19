# Reference Genomes and SnpEff Database Preparation

The genome and annotation files in this directory are used by the analysis workflows. Before running the WGS pipeline (`mycosnp-nf`), you must build a custom SnpEff database for variant annotation.

## Build the Custom SnpEff Database

**Step 1: Create the SnpEff Directory Structure**

The configuration expects a specific directory layout. From the project's root directory, run the following command to create it:

```bash
mkdir -p 02_whole_genome_sequencing/assets/snpeff/data/Candida_auris_B11205
```

**Step 2: Create a Custom `snpEff.config` File**

Create a file named `snpEff.config` inside `02_whole_genome_sequencing/assets/snpeff/`. This file tells SnpEff about our new genome.

**File: `02_whole_genome_sequencing/assets/snpeff/snpEff.config`**
```
# --- snpEff configuration file ---
# The path to the SnpEff data directory within this structure
data.dir = ./data/

# Define our custom genome
# The key 'Candida_auris_B11205' must match the directory name under 'data/'
Candida_auris_B11205.genome : Candida auris B11205
```

**Step 3: Link Reference Files into the SnpEff Directory**

Instead of copying, we can create symbolic links (shortcuts) to the reference files. This avoids data duplication. The files must be renamed to what SnpEff expects (`genes.gtf` and `genome.fa`).

```bash
# Link the GTF file
ln -s ../../../../reference_genomes/candida_auris_b11205/GCA_016772135.1_genomic.gtf \
      02_whole_genome_sequencing/assets/snpeff/data/Candida_auris_B11205/genes.gtf

# Link the FASTA file
ln -s ../../../../reference_genomes/candida_auris_b11205/GCA_016772135.1_genomic.fna \
      02_whole_genome_sequencing/assets/snpeff/data/Candida_auris_B11205/genome.fa
```

**Step 4: Build the Database**

Finally, activate the Conda environment and run the SnpEff `build` command from the correct directory.

```bash
# Activate the environment first
conda activate cauris-publication

# Navigate into the snpEff directory to run the build command
cd 02_whole_genome_sequencing/assets/snpeff/

# Build the database using our custom config file
snpEff build -c ./snpEff.config -gtf22 -v Candida_auris_B11205

# Navigate back to the project root directory
cd ../../../..
```

After these steps, a `snpEffectDb.bin` file will be created in the `.../snpeff/data/Candida_auris_B11205/` directory. The setup is now complete, and you can run the WGS pipeline as described in the main `README.md`.
