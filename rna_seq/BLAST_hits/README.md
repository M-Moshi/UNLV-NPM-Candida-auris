# Candida auris Pathway Annotation Mapping

## Overview
This repository contains the BLAST hit results generated to establish robust gene mapping protocols for *Candida auris* (*C. auris*). Because neither KEGG nor Gene Ontology (GO) databases directly supported Candida Genome Database (CGD) identifiers natively, we developed a comprehensive mapping framework to translate these identifiers for downstream pathway enrichment analysis.

## Files in this Directory
* `GO-caur-BLAST.csv`: BLAST mapping results linking *C. auris* genes to GO counterparts.
* `KEGG-CGD-BLAST-res.csv`: BLAST mapping results linking *C. auris* genes to KEGG counterparts.

## Methodology
To ensure high-confidence annotations for examining biological pathways under different stress conditions, the mappings were generated using the following parameters and sources:
* **Alignment Tool:** BLAST (version 2.5.0)
* **Similarity Threshold:** A stringent threshold of **95%** was implemented for all alignments.
* **KEGG Mapping:** Nucleotide sequences were obtained directly from the KEGG database (release 111.0).
* **GO Mapping:** Gene sequences were sourced from the NCBI *C. auris* assembly ASM301371v2 (GO release 2024-09-08).

These mappings serve as the foundational inputs for Pathway enrichment analysis, which is subsequently performed using the Consensus Pathway Analysis (CPA) framework.

## Data Dictionary
The CSV files in this directory follow the standard BLAST tabular output format (format 6). The columns are defined as follows:

| Column | Name | Description |
| :--- | :--- | :--- |
| 1 | `qseqid` | Query or source sequence ID (e.g., the *C. auris* gene / CGD identifier). |
| 2 | `sseqid` | Subject or target sequence ID (the reference genome sequence in KEGG/GO). |
| 3 | `pident` | Percentage of identical positions between the query and subject. |
| 4 | `length` | Alignment length (total sequence overlap). |
| 5 | `mismatch` | Number of mismatches in the alignment. |
| 6 | `gapopen` | Number of gap openings. |
| 7 | `qstart` | Start of the alignment in the query sequence. |
| 8 | `qend` | End of the alignment in the query sequence. |
| 9 | `sstart` | Start of the alignment in the subject sequence. |
| 10 | `send` | End of the alignment in the subject sequence. |
| 11 | `evalue` | Expect value (statistical significance of the alignment). |
| 12 | `bitscore` | Bit score (normalized alignment score). |

---
*Note: These files only retain alignments that met the >=95% similarity threshold mentioned in the methodology.*
