#!/bin/bash
#
# This script builds all custom SnpEff databases into a single, centralized
# 'reference_genomes/snpEff' directory.
#
# It performs two critical actions:
#   1. Defines a custom codon table named 'Alternative_Yeast_Nuclear'.
#   2. Builds a database for each genome found within the 'reference_genomes'
#      directory using the GenBank method and this custom codon table.
#
# USAGE:
#   Run this script from the project's root directory after adding a new
#   genome directory or if a rebuild is needed.
#   Example: bash reference_genomes/build_custom_snpeff_dbs.sh
#

set -e # Exit immediately if a command exits with a non-zero status.

# --- Sanity Check: Ensure the script is run from the project root ---
if [ ! -d "reference_genomes" ]; then
    echo "Error: 'reference_genomes' directory not found."
    echo "Please run this script from the project's root directory (the parent of 'reference_genomes')."
    exit 1
fi

# --- Define Central Paths ---
# Use absolute paths to prevent issues with symbolic links.
PROJECT_ROOT=$(pwd)
SNPEFF_DIR="${PROJECT_ROOT}/reference_genomes/snpEff"
SNPEFF_DATA_DIR="${SNPEFF_DIR}/data"
SNPEFF_GENOMES_DIR="${SNPEFF_DATA_DIR}/genomes"
CONFIG_FILE="${SNPEFF_DIR}/snpEff.config"

# --- Step 1: Initialize the Central SnpEff Directory and Config File ---
echo "Initializing central SnpEff directory at: ${SNPEFF_DIR}"
mkdir -p "${SNPEFF_GENOMES_DIR}"

# Create a new config file from scratch to ensure it's clean.
# The first line MUST be the custom codon table definition.
# We use a quoted 'EOF' to prevent any shell interpretation of the codon string.
echo "Creating new master config file with custom codon table: ${CONFIG_FILE}"
cat <<'EOF' > "${CONFIG_FILE}"
# --- Custom SnpEff configuration for analyzing project genomes ---

# 1. Define the custom codon table for Candida. This MUST come before it is used.
codon.Alternative_Yeast_Nuclear             : TTT/F, TTC/F, TTA/L, TTG/L, TCT/S, TCC/S, TCA/S, TCG/S, TAT/Y, TAC/Y, TAA/*, TAG/*, TGT/C, TGC/C, TGA/*, TGG/W, CTT/L, CTC/L, CTA/L, CTG/S+, CCT/P, CCC/P, CCA/P, CCG/P, CAT/H, CAC/H, CAA/Q, CAG/Q, CGT/R, CGC/R, CGA/R, CGG/R, ATT/I, ATC/I, ATA/I, ATG/M+, ACT/T, ACC/T, ACA/T, ACG/T, AAT/N, AAC/N, AAA/K, AAG/K, AGT/S, AGC/S, AGA/R, AGG/R, GTT/V, GTC/V, GTA/V, GTG/V, GCT/A, GCC/A, GCA/A, GCG/A, GAT/D, GAC/D, GAA/E, GAG/E, GGT/G, GGC/G, GGA/G, GGG/G

# 2. Define the data directory path for SnpEff.
#    All genome definitions will be appended below by the build script.
EOF

# Append the data directory path using its absolute path for robustness.
echo "data.dir = ${SNPEFF_DATA_DIR}" >> "${CONFIG_FILE}"


# --- Step 2: Loop through each genome directory and build its database ---
# The find command now robustly excludes the 'panel' and 'snpEff' directories.
for GENOME_DIR in $(find reference_genomes -mindepth 1 -maxdepth 1 -type d ! -name "panel" ! -name "snpEff"); do

    GENOME_ID=$(basename "${GENOME_DIR}")
    echo ""
    echo "=========================================================="
    echo "===> Processing genome: ${GENOME_ID}"
    echo "=========================================================="

    # Find the required .fna and .gbk files.
    GENOME_FASTA=$(find "${GENOME_DIR}" -name "*.fna" -o -name "*.fasta" | head -n 1)
    GENOME_GBK=$(find "${GENOME_DIR}" -name "*.gbk" -o -name "*.gbff" | head -n 1)

    if [[ -z "$GENOME_FASTA" || -z "$GENOME_GBK" ]]; then
        echo "Warning: Skipping '${GENOME_ID}'. Could not find both a FASTA (.fna) and a GenBank (.gbk) file."
        continue
    fi

    echo "Found FASTA: ${GENOME_FASTA}"
    echo "Found GenBank: ${GENOME_GBK}"

    # --- Step 2a: Append this genome's definition to the master config file ---
    echo "Appending definition for '${GENOME_ID}' to master snpEff.config..."
    CHROMOSOMES=$(grep '^>' "${GENOME_FASTA}" | sed 's/>//' | awk '{print $1}')
    CHROMOSOMES_COMMA_SEPARATED=$(echo "${CHROMOSOMES}" | tr '\n' ',' | sed 's/,$//')

    {
        echo ""
        echo "# --- Definition for ${GENOME_ID} ---"
        echo "${GENOME_ID}.genome : ${GENOME_ID}"
        echo "${GENOME_ID}.chromosomes : ${CHROMOSOMES_COMMA_SEPARATED}"
        # Set the custom codon table for each chromosome.
        for CHR in ${CHROMOSOMES}; do
            echo "${GENOME_ID}.${CHR}.codonTable : Alternative_Yeast_Nuclear"
        done
    } >> "${CONFIG_FILE}"

    # --- Step 3: Link reference files into the central SnpEff directories ---
    echo "Linking reference files into the snpEff data structure..."
    # Create the directory for this specific genome's data files (e.g., genes.gbk)
    mkdir -p "${SNPEFF_DATA_DIR}/${GENOME_ID}"
    # Link the GenBank file.
    ln -sf "${PROJECT_ROOT}/${GENOME_GBK}" "${SNPEFF_DATA_DIR}/${GENOME_ID}/genes.gbk"
    # Link the FASTA file into the 'genomes' directory.
    ln -sf "${PROJECT_ROOT}/${GENOME_FASTA}" "${SNPEFF_GENOMES_DIR}/${GENOME_ID}.fa"

    # --- Step 4: Build this specific database ---
    echo "Building database for '${GENOME_ID}' using the GenBank method..."
    snpEff build -c "${CONFIG_FILE}" -genbank -v "${GENOME_ID}"

    echo "Build complete for ${GENOME_ID}."
done

echo ""
echo "=========================================================="
echo "      All SnpEff database builds finished. "
echo "=========================================================="
echo "Master config file is located at: ${CONFIG_FILE}"
echo "Run snpEff with '-c ${CONFIG_FILE}' to use these databases."
echo "=========================================================="
