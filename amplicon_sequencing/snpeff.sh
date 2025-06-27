#!/bin/bash
#
# Script to process sentinel output to annotated uwing snpeff
# Handles ivar filtered tsv to vcf conversion, vcf annotation and extraction
#
# Usage ./snpeff.sh sentinel_output/ <snpEff_genome> <path_to_snpEff.config>

# --- Argument Checks ---
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <path_to_sentinel_output_dir> <snpEff_genome_name> <path_to_snpEff.config>"
    exit 1
fi

SENTINEL_DIR=$1
GENOME_NAME=$2
CONFIG_FILE=$3

# --- Making directories ---
echo "Creating output directories..."
mkdir -p "${SENTINEL_DIR}/tsv" "${SENTINEL_DIR}/annotation" "${SENTINEL_DIR}/vcf"

# --- Convert ivar filtered tsv to vcf ---
echo "Converting TSV files to VCF..."
for i in "${SENTINEL_DIR}/filtered_tsv/"*_filtered.tsv; do
  # Check if files exist to avoid errors on empty directories
  if [ ! -f "$i" ]; then
      echo "Warning: No filtered tsv files found in ${SENTINEL_DIR}/filtered_tsv/"
      continue
  fi
  prefix=$(basename "${i/_filtered.tsv}")
  python3 ivar_variants_to_vcf.py "$i" "${SENTINEL_DIR}/vcf/${prefix}.vcf"
done

# --- Annotate vcf ---
echo "Annotating VCF files..."
for vcf in "${SENTINEL_DIR}/vcf/"*.vcf
do
    # Check if files exist
    if [ ! -f "$vcf" ]; then
        echo "Warning: No VCF files found in ${SENTINEL_DIR}/vcf/ to annotate."
        continue
    fi
    prefix=$(basename "$vcf" .vcf)
    echo "Annotating ${vcf}"

    # CORRECTED COMMAND: Use 'eff' subcommand and '-v' flag for robustness
    snpEff eff -v "${GENOME_NAME}" -c "${CONFIG_FILE}" -no-downstream -no-upstream "${vcf}" > "${SENTINEL_DIR}/annotation/${prefix}_annotated.vcf"

    # Move SnpEff output files, checking if they exist first
    if [ -f "snpEff_genes.txt" ]; then
        mv snpEff_genes.txt "${SENTINEL_DIR}/annotation/${prefix}_snpEff_genes.txt"
    fi
    if [ -f "snpEff_summary.html" ]; then
        mv snpEff_summary.html "${SENTINEL_DIR}/annotation/${prefix}_snpEff_summary.html"
    fi
done

# --- Extract VCF using SnpSift and Python ---
echo "Extracting fields from annotated VCFs..."
for vcf in "${SENTINEL_DIR}/annotation/"*_annotated.vcf; do
    # Check if files exist
    if [ ! -f "$vcf" ]; then
        echo "Warning: No annotated VCF files found in ${SENTINEL_DIR}/annotation/ to extract."
        continue
    fi
    prefix=$(basename "$vcf" _annotated.vcf)
    echo "Extracting from ${vcf}"

    # Use SnpSift to extract fields, including the header
    SnpSift extractFields "$vcf" CHROM POS ID REF ALT QUAL FILTER \
    "GEN[*].DP" "ANN[0].ALLELE" "ANN[0].EFFECT" "ANN[0].IMPACT" "ANN[0].GENE" \
    "ANN[0].GENEID" "ANN[0].FEATUREID" "ANN[0].HGVS_P" "ANN[0].HGVS_C" > "${SENTINEL_DIR}/tsv/${prefix}_snpSift.tsv"

    # Call the Python script to process the VCF and the SnpSift-generated TSV file
    echo "Processing with process_vcf.py for ${vcf}"
    python3 process_vcf.py "$vcf" "${SENTINEL_DIR}/tsv/${prefix}_snpSift.tsv" "${SENTINEL_DIR}/tsv/${prefix}_combined.tsv"

done

echo "Processing complete."
