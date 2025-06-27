# RNA-Seq and Differential Abundance Analysis Workflow

This document details the multi-stage workflow used to process RNA-Seq data, from raw reads to differential gene expression analysis. The process involves two main `nf-core` pipelines and uses specially prepared reference annotation files.

## Prerequisites

* **Environment Setup:** Ensure you have completed the main environment setup described in the root `README.md`, including activating the `npm-candida` conda environment.

## A Note on Reference Annotation Files

The final, corrected reference annotation files used in this workflow are provided in the `reference_genomes/candida_auris_b8441/` directory. The publicly available annotations required modification prior to use, as documented below for reproducibility. **End-users do not need to repeat these steps.**

* `C_auris_B8441_current_exons.gtf`: This file was generated from the original public GFF file. The key modification was the inference and addition of `exon` and `mRNA` features, which were absent in the source and are required for the `nf-core/rnaseq` pipeline to correctly quantify gene expression.

* `C_auris_B8441_current_transcripts_annotated.gtf`: This file is a further modification used for the `nf-core/differentialabundance` pipeline. The feature type `mRNA` was converted to `transcript`, as the downstream statistical pipeline specifically requires the `transcript` feature for its analysis. Gene names were also added to this file from a separate annotation source to aid in interpretation.

## Workflow Stage 1: Run `nf-core/rnaseq`

This pipeline is used for read QC, alignment, and gene quantification.
See [nf-core/rnaseq](https://github.com/nf-core/rnaseq) for more information.

**Input Samplesheet:**
Create a samplesheet CSV file (`sample,fastq_1,fastq_2,strandedness`). A template is provided in `rna_seq/samplesheet_rnaseq_template.csv`.

**Command:**
```bash
nextflow run nf-core/rnaseq \
    -profile docker \
    --input /path/to/your/rnaseq_samplesheet.csv \
    --fasta ../reference_genomes/candida_auris_b8441/C_auris_B8441_current_chromosomes.fasta \
    --gtf ../reference_genomes/candida_auris_b8441/C_auris_B8441_current_exons.gtf \
    --outdir ./rnaseq_results
```
*We use the `_exons.gtf` file here as it contains the necessary features for alignment and quantification.*


## Workflow Stage 2: Run `nf-core/differentialabundance`

This pipeline performs the statistical comparison using the count data generated in Stage 1.
See [nf-core/differentialabundance](https://github.com/nf-core/differentialabundance) for more information.


**Input Files:**
1.  **Counts Matrix:** The primary input, located in the `rnaseq` output directory at `rnaseq_results/star_salmon/salmon.merged.gene_counts.tsv`.
2.  **Design Sheet:** A CSV file linking each sample to its experimental group. See `rna_seq/design_template.csv`.
3.  **Contrasts Sheet:** A CSV file defining the comparisons to be made. See `rna_seq/contrasts_template.csv`.

**Command:**
```bash
nextflow run nf-core/differentialabundance -r 2.0.0 \
    -profile docker \
    --input /path/to/your/design.csv \
    --contrasts /path/to/your/contrasts.csv \
    --matrix ./rnaseq_results/star_salmon/salmon.merged.gene_counts.tsv \
    --transcript_length_metric ./rnaseq_results/star_salmon/salmon.merged.gene_lengths.tsv \
    --gtf ../reference_genomes/candida_auris_b8441/C_auris_B8441_current_transcripts_annotated.gtf \
    --outdir ./differential_results
```
*We use the final `_transcripts_annotated.gtf` here because it contains the `transcript` features required by this specific pipeline.*


---

## Workflow Stage 3: Consensus Pathway Analysis

To understand the biological meaning behind the differential expression results, we performed a comprehensive pathway analysis using the **R package for Consensus Pathway Analysis (RCPA)**.

This framework allowed us to identify which biological pathways were significantly impacted by the different stress conditions tested in our experiments. Because direct annotations for *Candida auris* were not available in the KEGG and Gene Ontology (GO) databases, we first mapped the genes using BLAST with a stringent similarity threshold.

We then used the `runGeneSetAnalysis()` function from the RCPA package to perform Fast Gene Set Enrichment Analysis (FGSEA) on our DESeq2 results. This revealed the distinct pathway signatures associated with each experimental condition.

-   **RCPA GitHub Repository**: [tinnlab/RCPA](https://github.com/tinnlab/RCPA](https://github.com/tinnlab/RCPA)
