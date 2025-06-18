#!/bin/bash
#
# UNLV Lab of Neuogenetics and Precision Medicine Amplicon Pipeline
# Author: Michael Moshi

## Alternative manual input for CPU
CPU=8

{
   echo "UNLV Lab of Neuogenetics and Precision Medicine Amplicon Pipeline"
   echo ""
   echo "Usage: ./Amplicon_pipeline.sh <input_fastq_folder>  <output_folder>  <reference.fasta>  <primers.tab>"
   echo ""
   exit 1 # Exit script after printing help
}

# Check for dependencies
for cmd in "bwa" "fastp" "fgbio" "samtools" "ivar" "pigz" "unpigz" "parallel" "kraken2" "qualimap" "python3"; do
  if ! command -v $cmd &> /dev/null; then
    echo "Error: $cmd is not installed."
    echo "Please install $cmd before running the pipeline."
    exit 1
  fi
done


# Print helpFunction in case parameters are empty
if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi

### Check if there are files with R1_001.fastq.gz suffix in the input folder
if ! ls $1/*R1_001.fastq.gz &> /dev/null; then
  echo "Error: No files with R1_001.fastq.gz suffix found in the input folder."
  echo "Please check the input folder for fastq.gz files and try again."
  exit 1
fi



## Create run log file with version and date
# Initialize log file
out_dir=$2
log_dir=$out_dir/logs
mkdir -p $log_dir
log_file="$out_dir/pipeline.log"
fastq_files=$(ls $1/*_R1_001.fastq.gz)
number_of_files=$(ls $fastq_dir/*_R1_001.fastq.gz | wc -l)
echo "Pipeline started at $(date)" >> $log_file
echo "Log file: $log_file" >> $log_file
echo "FASTQ Directory: $1" >> $log_file
echo "Output Directory: $out_dir" >> $log_file
echo "Reference Genome: $3" >> $log_file
echo "Kraken Database: $kraken_db" >> $log_file
echo "" >> $log_file
echo "Hardware Information:" >> $log_file
echo "CPU type: $(lscpu | grep "Model name" | awk '{$1=$2=""; print $0}')" >> $log_file
echo "Number of threads available: $(nproc)" >> $log_file
echo "Number of threads used: $max_jobs" >> $log_file
echo "Memory available: $(free -h | grep Mem | awk '{print $7}')" >> $log_file
echo "Disk space available in output folder drive: $(df -h $out_dir | tail -n 1 | awk '{print $4}')" >> $log_file
echo "" >> $log_file
echo "Software Used:" >> $log_file
echo "fastqc: $(fastqc --version)" >> $log_file
fastp -v 2>&1 | head | tee -a $log_file
echo "bwa: $(bwa 2>&1 | head -n 3 | tail -n 1)" >> $log_file
echo "samtools: $(samtools --version | head -n 1)" >> $log_file
echo "ivar: $(ivar version | head -n 1)" >> $log_file
echo "fgbio: $(fgbio --version | head -n 1)" >> $log_file
echo "bcftools: $(bcftools --version | head -n 1)" >> $log_file
echo "multiqc: $(multiqc --version)" >> $log_file
echo "qualimap: $(qualimap --version | head -n 4 | tail -n 1 )" >> $log_file
echo "parallel: $(parallel --version | head -n 1)" >> $log_file
echo "kraken2: $(kraken2 --version | head -n 1)" >> $log_file
echo "bracken: $(bracken -v | head -n 1)" >> $log_file
echo "Number of paired ends reads found: $number_of_files" >> $log_file
echo "Found the following fastq files:" >> $log_file
echo "$fastq_files" >> $log_file
echo "" >> $log_file
echo "Running pipeline..."

# outdir=$(dirname ${2})

mkdir -p $2
mkdir -p $2/trimmed_fastq
mkdir -p $2/sam
mkdir -p $2/bam
mkdir -p $2/raw_tsv
mkdir -p $2/unfiltered_tsv
mkdir -p $2/filtered_bam
mkdir -p $2/cov
mkdir -p $2/histogram
mkdir -p $2/depth
mkdir -p $2/filtered_tsv
mkdir -p $2/stats
mkdir -p $2/qualimap
mkdir -p $2/kraken2
mkdir -p $2/logs
initdir=$PWD

genomedir=$(basename $(dirname ${3}))

# Trimming before alignment with 16 threads
for file in $1/*R1_001.fastq.gz
do
  read1=$file
  read2=${read1/_R1_001.fastq.gz/_R2_001.fastq.gz}
  prefix=$(basename ${read1/_R1_001.fastq.gz})
  start_time=$(date +%s)
  cut1=$2/trimmed_fastq/${prefix}_cutadapt_R1.fastq.gz
  cut2=$2/trimmed_fastq/${prefix}_cutadapt_R2.fastq.gz
  cutlog=$2/trimmed_fastq/${prefix}_cutadapt.log
  fastplog=$2/trimmed_fastq/${prefix}_fastp.html
  fastpjson=$2/trimmed_fastq/${prefix}_fastp.json
  echo $prefix ":: fastp trimming"
  fastp -i $read1 -I $read2 -o $cut1 -O $cut2 -j $fastpjson -w 16 -g
  end_time=$(date +%s)
  echo "Fastp for $prefix took $((end_time - start_time)) seconds" >> $2/logs/${prefix}.log
done
fastq_dir=$1
out_dir=$2
ref_genome=$3
export out_dir
export ref_genome

process_read() {
  read1=$1
  output_dir=$2
  reference=$3
  primers=$4

  read2=${read1/_R1_001.fastq.gz/_R2_001.fastq.gz}
  prefix=$(basename ${read1/_R1_001.fastq.gz})
  ivprefix=$output_dir/raw_tsv/${prefix}_iv
  fgprefix=$output_dir/unfiltered_tsv/${prefix}_fg
  cut1=$output_dir/trimmed_fastq/${prefix}_cutadapt_R1.fastq.gz
  cut2=$output_dir/trimmed_fastq/${prefix}_cutadapt_R2.fastq.gz
  align=$output_dir/sam/${prefix}.sam
  view=$output_dir/sam/${prefix}.sam
  sorted=$output_dir/bam/${prefix}.bam
  cutlog=$output_dir/trimmed_fastq/${prefix}_cutadapt.log
  fastplog=$output_dir/trimmed_fastq/${prefix}_fastp.html
  fastpjson=$output_dir/trimmed_fastq/${prefix}_fastp.json
  fgbam=$output_dir/bam/${prefix}_fgtrim.bam
  filteredbam=$output_dir/filtered_bam/${prefix}_filtered.bam
  filteredprefix=$output_dir/filtered_tsv/${prefix}_filtered
  variant=$output_dir/filtered_bam/bam/${prefix}.variant
  freydepth=$output_dir/filtered_bam/bam/${prefix}.depth
  demix=$output_dir/filtered_bam/bam/${prefix}.demix
  echo $prefix ":: bwa alignment"

  bwa index $3
  start_time_bwa=$(date +%s)
  bwa mem \
    -t 2 \
    -o $align \
    $3 \
    $cut1 \
    $cut2
  end_time_bwa=$(date +%s)
  echo "BWA MEM for $prefix took $((end_time_bwa - start_time_bwa)) seconds" >> $2/logs/${prefix}.log


  pigz -p 1 $align
  unpigz -c ${align}.gz |  samtools view -h -b | samtools sort -@ 1 -o $sorted
  start_time_samtools=$(date +%s)
  samtools index $sorted
  end_time_samtools=$(date +%s)
  echo "Samtools sort and index for $prefix took $((end_time_samtools - start_time_samtools)) seconds" >> $2/logs/${prefix}.log

sleep 5

  start_time_fgbio=$(date +%s)
  fgbio --tmp-dir=$2/tmp TrimPrimers \
  -i $sorted -o $fgbam -p $4 -H true
  end_time_fgbio=$(date +%s)
  echo "TrimPrimers for $prefix took $((end_time_fgbio - start_time_fgbio)) seconds" >> $2/logs/${prefix}.log
  #samtools index $fgbam

# filteredbam filter the bam for length <40 i guess.
  start_time_filter=$(date +%s)
  fgbio FilterBam \
  -i $fgbam -o $filteredbam -m 40
  end_time_filter=$(date +%s)
  echo "FilterBam for $prefix took $((end_time_filter - start_time_filter)) seconds" >> $2/logs/${prefix}.log
  #samtools index $filteredbam


############

    start_time_ivar=$(date +%s)
    samtools mpileup -aa -A -d 0 -B -Q 0 --reference $3 $sorted | ivar variants -p $ivprefix -r $3 &


    samtools mpileup -aa -A -d 0 -B -Q 0 --reference $3 $fgbam | ivar variants -p $fgprefix -r $3 &


    samtools mpileup -aa -A -d 0 -B -Q 0 --reference $3 $filteredbam | ivar variants -p $filteredprefix -r $3 

    end_time_ivar=$(date +%s)
    echo "iVar for $prefix took $((end_time_ivar - start_time_ivar)) seconds" >> $2/logs/${prefix}.log

  
    start_time_coverage=$(date +%s)
    samtools coverage $filteredbam > $2/cov/${prefix}.cov
    samtools coverage -m $filteredbam > $2/histogram/${prefix}.hist
    samtools depth -@ 1 -a -d 0 $filteredbam >$2/depth/${prefix}.udepth
    samtools stats $filteredbam > $2/stats/${prefix}.stats
    end_time_coverage=$(date +%s)
    echo "Coverage for $prefix took $((end_time_coverage - start_time_coverage)) seconds" >> $2/logs/${prefix}.log
}

export -f process_read


echo ""
echo "Now running the pipeline in parallel."
echo "This can take a while, be patient"
echo ""

find $1 -type f -name "*R1_001.fastq.gz" | parallel -j $CPU --bar process_read {} $2 $3 $4 $6 > /dev/null
wait
echo "Finished processing reads in parallel at $(date)" >> $2/pipeline.log


### Run Qaulimap
mkdir -p $2/qualimap
for file in $2/filtered_bam/*_filtered.bam
do
    prefix=$(basename ${file/_filtered.bam})
    start_time_qualimap=$(date +%s)
    echo "running qualimap on $prefix"
    qualimap bamqc -bam $file -outdir $2/qualimap/${prefix} -nt 8
    end_time_qualimap=$(date +%s)
    echo "Qualimap for $prefix took $((end_time_qualimap - start_time_qualimap)) seconds" >> $2/logs/${prefix}.log
done
echo "Finished running Qualimap at $(date)" >> $2/logs/pipeline.log
cd $initdir
cd $2/filtered_tsv

## Run MultiQC


# IF function for $6 is provided, then use grep and se commands
if [ -z "$6" ]
then
  echo "No GFF file provided. Skipping amino acid translation."
    grep -HP '^[A-Z][A-Z][_0123456789]' *.tsv | cut -sf 1-5,8,11-14 | perl -pe 's/_S[0123456789]+_filtered.tsv:NC_045512.2//g' > all_filtered.tsv
    sed -i "1i\sample	POS	REF	ALT	REF_DP	ALT_DP	ALT_FREQ	TOTAL_DP	PVAL	PASS" all_filtered.tsv
else
  echo "GFF file provided. Amino acid translation will be performed."
    grep -HP '^[A-Z][A-Z][_0123456789]' *.tsv | cut -sf 1-5,8,11-20 | perl -pe 's/_S[0123456789]+_filtered.tsv:NC_045512.2//g' > all_filtered.tsv
    sed -i "1i\sample	POS	REF	ALT	REF_DP	ALT_DP	ALT_FREQ	TOTAL_DP	PVAL	PASS	GFF_FEATURE	REF_CODON	REF_AA	ALT_CODON	ALT_AA	POS_AA" all_filtered.tsv
fi

# Merge all filtered tsv files into one
grep -HP '^[A-Z][A-Z][_0123456789]' *.tsv | cut -sf 1-5,8,11-14 | perl -pe 's/_S[0123456789]+_filtered.tsv:NC_045512.2//g' > all_filtered.tsv
sed -i "1i\sample	POS	REF	ALT	REF_DP	ALT_DP	ALT_FREQ	TOTAL_DP	PVAL	PASS" all_filtered.tsv
mv $2/filtered_tsv/all_filtered.tsv $2

echo "Finished running pipeline at $(date)" >> $2/logs/pipeline.log
cd $initdir
echo "Running integrated Python statistics script..."
python3 - "$2" << 'EOF'
import os
import glob
import pandas as pd
import sys

# Use command-line arguments for output_dir
# sys.argv[0] is the script name ('-'), sys.argv[1] is the first argument passed
if len(sys.argv) < 2:
    print("Error: Output directory not provided to Python script.")
    sys.exit(1)
output_dir = sys.argv[1]

depth_dir = os.path.join(output_dir, 'depth')
cov_dir = os.path.join(output_dir, 'cov')
stat_dir = os.path.join(output_dir, 'stats') # Corrected from 'stat' to match your bash script

# Ensure the statistics directory exists
os.makedirs(stat_dir, exist_ok=True)

# Print all files in the depth directory (for verification)
# print(f"Files in {depth_dir}:")
# for file in os.listdir(depth_dir):
#     print(file)

# List to store all output_data DataFrames
all_data = []

# Process each file in the depth directory
for depth_file in glob.glob(os.path.join(depth_dir, '*.udepth')):
    # Improved sample name extraction to be more robust
    sample_name = os.path.basename(depth_file).replace('.udepth', '')
    
    # Initialize placeholders for data
    genome = 'N/A'
    length = 0
    mean_depth = 0.0
    median_depth = 0.0
    numreads = 0
    covered1 = 0.0
    covered10 = 0.0
    coveredby30 = 0.0
    covered50 = 0.0
    covered100 = 0.0
    
    # Process depth files
    try:
        depth_data = pd.read_csv(depth_file, sep='\t', names=['genome', 'pos', 'depth'])
        if not depth_data.empty:
            genome = depth_data['genome'].iloc[0]
            length = len(depth_data)
            mean_depth = depth_data['depth'].mean()
            median_depth = depth_data['depth'].median()
            # Calculate coverage percentages
            if length > 0:
                covered1 = (depth_data['depth'] >= 1).sum() / length * 100
                covered10 = (depth_data['depth'] >= 10).sum() / length * 100
                coveredby30 = (depth_data['depth'] >= 30).sum() / length * 100
                covered50 = (depth_data['depth'] >= 50).sum() / length * 100
                covered100 = (depth_data['depth'] >= 100).sum() / length * 100
    except pd.errors.EmptyDataError:
        print(f"Warning: Depth file {depth_file} is empty.")

    # Match coverage files for each sample
    matching_cov_files = glob.glob(os.path.join(cov_dir, f'{sample_name}.cov'))
    
    # Process coverage files (if any)
    if matching_cov_files:
        cov_file = matching_cov_files[0]
        try:
            # samtools coverage output has a header line, so we skip it.
            cov_data = pd.read_csv(cov_file, sep='\t', header=0)
            if not cov_data.empty and 'numreads' in cov_data.columns:
                numreads = cov_data['numreads'].iloc[0]
        except (pd.errors.EmptyDataError, IndexError):
            print(f"Warning: Could not process coverage file {cov_file}.")


    # Prepare the output DataFrame for this sample
    output_data = pd.DataFrame({
        'Sample Name': [sample_name],
        'Genome': [genome],
        'Length': [length],
        'Total Reads': [numreads],
        'Mean Depth': [f"{mean_depth:.2f}"],
        'Median Depth': [median_depth],
        'Coveredby1x': [f"{covered1:.2f}%"],
        'Coveredby10x': [f"{covered10:.2f}%"],
        'Coveredby30x': [f"{coveredby30:.2f}%"],
        'Coveredby50x': [f"{covered50:.2f}%"],
        'Coveredby100x': [f"{covered100:.2f}%"],
    })

    # Save individual sample statistics
    output_file = os.path.join(stat_dir, f'{sample_name}.stats.tsv')
    output_data.to_csv(output_file, sep='\t', index=False)
    # print(f"Saved output for sample {sample_name} to {output_file}")

    # Append to the list for combined output later
    all_data.append(output_data)

# Combine all individual sample DataFrames into one and save
if all_data:
    combined_data = pd.concat(all_data, ignore_index=True)
    combined_file = os.path.join(output_dir, 'combined_stats.tsv')
    combined_data.to_csv(combined_file, sep='\t', index=False)
    print(f"Successfully generated combined statistics report at: {combined_file}")
else:
    print("Warning: No data was processed to generate a combined statistics report.")

EOF
