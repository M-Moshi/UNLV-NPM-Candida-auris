#!/bin/bash
#
# This script builds all custom SnpEff databases. It is designed to be
# executed from WITHIN the 'reference_genomes' directory.
#
# It creates a centralized 'snpEff' directory here, attempts to build all
# found genomes, and provides a final report with detailed logs for any failures.
#
# USAGE:
#   cd /path/to/your/reference_genomes/
#   bash build_custom_snpeff_dbs.sh
#


# --- Define Paths ---
CWD=$(pwd)
SNPEFF_DIR="${CWD}/snpEff"
SNPEFF_DATA_DIR="${SNPEFF_DIR}/data"
SNPEFF_GENOMES_DIR="${SNPEFF_DATA_DIR}/genomes"
CONFIG_FILE="${SNPEFF_DIR}/snpEff.config"

# --- Step 1: Initialize the Central SnpEff Directory and Config File ---
echo "Initializing SnpEff directory at: ${SNPEFF_DIR}"
mkdir -p "${SNPEFF_GENOMES_DIR}"
echo "Creating new master config file: ${CONFIG_FILE}"
cat <<'EOF' > "${CONFIG_FILE}"
# --- Custom SnpEff configuration for analyzing project genomes ---
# 1. Define the custom codon table for Candida. This MUST come before it is used.
codon.Alternative_Yeast_Nuclear             : TTT/F, TTC/F, TTA/L, TTG/L, TCT/S, TCC/S, TCA/S, TCG/S, TAT/Y, TAC/Y, TAA/*, TAG/*, TGT/C, TGC/C, TGA/*, TGG/W, CTT/L, CTC/L, CTA/L, CTG/S+, CCT/P, CCC/P, CCA/P, CCG/P, CAT/H, CAC/H, CAA/Q, CAG/Q, CGT/R, CGC/R, CGA/R, CGG/R, ATT/I, ATC/I, ATA/I, ATG/M+, ACT/T, ACC/T, ACA/T, ACG/T, AAT/N, AAC/N, AAA/K, AAG/K, AGT/S, AGC/S, AGA/R, AGG/R, GTT/V, GTC/V, GTA/V, GTG/V, GCT/A, GCC/A, GCA/A, GCG/A, GAT/D, GAC/D, GAA/E, GAG/E, GGT/G, GGC/G, GGA/G, GGG/G
# 2. Define the data directory path for SnpEff.
EOF
echo "data.dir = ${SNPEFF_DATA_DIR}" >> "${CONFIG_FILE}"

# --- Create arrays to track successes and failures ---
declare -a successes=()
declare -a failures=()

# --- Step 2: Loop through each genome directory ---
for GENOME_DIR in $(find . -mindepth 1 -maxdepth 1 -type d ! -name "panel" ! -name "snpEff"); do

    GENOME_ID=$(basename "${GENOME_DIR}")
    BUILD_LOG_FILE="${SNPEFF_DIR}/${GENOME_ID}.build.log"
    echo ""
    echo "=========================================================="
    echo "===> Processing genome: ${GENOME_ID}"
    echo "=========================================================="

    # --- Find the FASTA file, with special handling for the panel ---
    GENOME_FASTA_RELATIVE=$(find "${GENOME_DIR}" -name "*.fna" -o -name "*.fasta" | head -n 1)

    if [[ -z "$GENOME_FASTA_RELATIVE" ]]; then
        echo "WARNING: Skipping '${GENOME_ID}'. Could not find a FASTA file."
        failures+=("${GENOME_ID} (FASTA file not found)")
        continue
    fi
    GENOME_FASTA_ABSOLUTE="${CWD}/${GENOME_FASTA_RELATIVE#./}"
    echo "Found FASTA file: ${GENOME_FASTA_RELATIVE}"


    # --- Append genome info to the config file (with special logic for the panel) ---
    if [[ "${GENOME_ID}" == "C_auris_panel" ]]; then
        echo "Appending simplified definition for '${GENOME_ID}' to config..."
        {
            echo ""
            echo "# --- Definition for ${GENOME_ID} (Panel) ---"
            echo "${GENOME_ID}.genome : ${GENOME_ID}"
            echo "${GENOME_ID}.codonTable : Alternative_Yeast_Nuclear"
        } >> "${CONFIG_FILE}"
    else
        echo "Appending per-chromosome definition for '${GENOME_ID}' to config..."
        CHROMOSOMES=$(grep '^>' "${GENOME_FASTA_ABSOLUTE}" | sed 's/>//' | awk '{print $1}')
        CHROMOSOMES_COMMA_SEPARATED=$(echo "${CHROMOSOMES}" | tr '\n' ',' | sed 's/,$//')
        {
            echo ""
            echo "# --- Definition for ${GENOME_ID} ---"
            echo "${GENOME_ID}.genome : ${GENOME_ID}"
            echo "${GENOME_ID}.chromosomes : ${CHROMOSOMES_COMMA_SEPARATED}"
            for CHR in ${CHROMOSOMES}; do
                echo "${GENOME_ID}.${CHR}.codonTable : Alternative_Yeast_Nuclear"
            done
        } >> "${CONFIG_FILE}"
    fi

    ln -sf "${GENOME_FASTA_ABSOLUTE}" "${SNPEFF_GENOMES_DIR}/${GENOME_ID}.fa"

    # --- Determine build method and run SnpEff ---
    GENOME_GBK_RELATIVE=$(find "${GENOME_DIR}" -name "*.gbk" -o -name "*.gbff" | head -n 1)
    if [[ -n "$GENOME_GBK_RELATIVE" ]]; then
        GENOME_GBK_ABSOLUTE="${CWD}/${GENOME_GBK_RELATIVE#./}"
        echo "Found GenBank file: ${GENOME_GBK_RELATIVE}. Using '-genbank' build method."
        mkdir -p "${SNPEFF_DATA_DIR}/${GENOME_ID}"
        ln -sf "${GENOME_GBK_ABSOLUTE}" "${SNPEFF_DATA_DIR}/${GENOME_ID}/genes.gbk"
        
        # CORRECTED FLAG: Use -noCheckProtein to skip protein validation
        if snpEff build -c "${CONFIG_FILE}" -genbank -v -noCheckProtein "${GENOME_ID}" > "${BUILD_LOG_FILE}" 2>&1; then
            successes+=("${GENOME_ID}")
            rm -f "${BUILD_LOG_FILE}" # Clean up log on success
        else
            failures+=("${GENOME_ID} (GenBank build error - see log)")
        fi
    else
        GENOME_GFF_RELATIVE=$(find "${GENOME_DIR}" -name "*.gff" -o -name "*.gtf" | head -n 1)
        if [[ -n "$GENOME_GFF_RELATIVE" ]]; then
            GENOME_GFF_ABSOLUTE="${CWD}/${GENOME_GFF_RELATIVE#./}"
            echo "Found GFF/GTF file: ${GENOME_GFF_RELATIVE}. Using '-gff3' build method."
            mkdir -p "${SNPEFF_DATA_DIR}/${GENOME_ID}"
            ln -sf "${GENOME_GFF_ABSOLUTE}" "${SNPEFF_DATA_DIR}/${GENOME_ID}/genes.gff"

            # CORRECTED FLAGS: Use valid checks to bypass protein/CDS validation, which often fails with public data.
            if snpEff build -c "${CONFIG_FILE}" -gff3 -v -noCheckProtein -noCheckCds "${GENOME_ID}" > "${BUILD_LOG_FILE}" 2>&1; then
                successes+=("${GENOME_ID}")
                rm -f "${BUILD_LOG_FILE}" # Clean up log on success
            else
                failures+=("${GENOME_ID} (GFF build error - see log)")
            fi
        else
            echo "WARNING: Skipping '${GENOME_ID}'. Could not find a .gbk, .gff, or .gtf annotation file."
            failures+=("${GENOME_ID} (No annotation file)")
            continue
        fi
    fi
done

# --- Final Report ---
echo ""
echo "=========================================================="
echo "                    FINAL BUILD REPORT"
echo "=========================================================="
echo ""
echo "Successful builds (${#successes[@]}):"
for g in "${successes[@]}"; do
    echo "  - ${g}"
done
echo ""
echo "Failed builds (${#failures[@]}):"
for g in "${failures[@]}"; do
    GENOME_ID=$(echo "$g" | awk '{print $1}')
    LOG_FILE="${SNPEFF_DIR}/${GENOME_ID}.build.log"
    echo "  - ${g}"
    if [ -f "$LOG_FILE" ]; then
        echo "    -> Check log for details: ${LOG_FILE}"
    fi
done
echo ""
echo "=========================================================="
if [ ${#failures[@]} -gt 0 ]; then
    echo "Please review the SnpEff error logs for the failed builds above."
    exit 1
fi
