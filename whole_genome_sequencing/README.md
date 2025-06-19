# WGS Analysis using `CDCgov/mycosnp-nf`

This directory contains instructions and assets for performing whole genome SNP analysis using the `CDCgov/mycosnp-nf` Nextflow pipeline.

## Prerequisites

1.  **Environment Setup:** Ensure you have completed the main environment setup described in the root `README.md`, including activating the `cauris-publication` conda environment.
2.  **SnpEff Database:** Ensure you have successfully built the custom SnpEff database by following the instructions in `reference_genomes/README.md`.

## Input File: Samplesheet

The pipeline requires a CSV samplesheet specifying the sample ID and the location of the FASTQ files.

**Format:** The file must be a CSV with the header `sample,fastq_1,fastq_2`.

* `sample`: A unique identifier for the sample.
* `fastq_1`: Full path to the forward read file (R1).
* `fastq_2`: Full path to the reverse read file (R2).

**Template:** A template file named `assets/samplesheet_wgs_template.csv` is provided.

**File: `assets/samplesheet_template.csv`**
```csv
sample,fastq_1,fastq_2
CA_SAMPLE_01,/path/to/your/fastqs/sample1_R1.fastq.gz,/path/to/your/fastqs/sample1_R2.fastq.gz
CA_SAMPLE_02,/path/to/your/fastqs/sample2_R1.fastq.gz,/path/to/your/fastqs/sample2_R2.fastq.gz
CA_SAMPLE_03,/path/to/your/fastqs/sample3_R1.fastq.gz,/path/to/your/fastqs/sample3_R2.fastq.gz
```

## Running the Pipeline

The command below will execute the MycoSNP pipeline using the specified version and parameters.

**Step 1 (Optional): Pull the pipeline and container**
```bash
nextflow pull CDCgov/mycosnp-nf -r 1.5
```

**Step 2: Run the analysis**
```bash
nextflow run CDCgov/mycosnp-nf \
    -profile docker \
    --input ./path/to/your_wgs_samplesheet.csv \
    --fasta ../reference_genomes/candida_auris_b11205/GCA_016772135.1_genomic.fna \
    --outdir ./wgs_results \
    --snpeffconfig ./assets/snpeff/snpEff.config
```

### Command Parameters Explained
* `-profile docker`: Tells Nextflow to use Docker as the execution engine. Docker must be running.
* `--input`: Path to your completed samplesheet.
* `--fasta`: Path to the reference genome.
* `--outdir`: The directory where all results will be saved.
* `--snpeffconfig`: Path to our custom SnpEff configuration file, which points to our custom *C. auris* database.

## Expected Output

The pipeline will generate a feature-rich output directory, including:
* BAM alignment files
* VCF files with raw and filtered SNP calls
* Variant annotation tables from SnpEff
* A final SNP matrix/alignment for phylogenetic analysis
* Various quality control reports.
