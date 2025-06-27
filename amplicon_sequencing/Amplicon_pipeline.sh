#!/bin/bash
#
# UNLV Lab of Neuogenetics and Precision Medicine Amplicon Pipeline
# Author: Michael Moshi

## Alternative manual input for CPU
CPU=8

{
   echo "UNLV Lab of Neuogenetics and Precision Medicine Amplicon Pipeline"
   echo ""
   echo "Usage: ./Amplicon_pipeline.sh <input_fastq_folder>  <output_folder>  <reference.fasta>  <primers.tab>  <genes.tsv>"
   echo ""
   exit 1 # Exit script after printing help
}

# Check for dependencies
for cmd in "bwa" "fastp" "fgbio" "samtools" "ivar" "pigz" "unpigz" "parallel" "qualimap" "python3"; do
  if ! command -v $cmd &> /dev/null; then
    echo "Error: $cmd is not installed."
    echo "Please install $cmd before running the pipeline."
    exit 1
  fi
done


# Print helpFunction in case parameters are empty
if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi

### Check if there are files with R1_001.fastq.gz suffix in the input folder
if ! ls $1/*.fastq.gz &> /dev/null; then
  echo "Error: No files with .fastq.gz suffix found in the input folder."
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
fastq_files=$(ls $1/*_1.fastq.gz)
number_of_files=$(ls $fastq_dir/*_R1_001.fastq.gz | wc -l)
number_of_files=$(ls $fastq_dir/*_1.fastq.gz | wc -l)
echo "Pipeline started at $(date)" >> $log_file
echo "Log file: $log_file" >> $log_file
echo "FASTQ Directory: $1" >> $log_file
echo "Output Directory: $out_dir" >> $log_file
echo "Reference Genome: $3" >> $log_file
echo "" >> $log_file
echo "Hardware Information:" >> $log_file
echo "CPU type: $(lscpu | grep "Model name" | awk '{$1=$2=""; print $0}')" >> $log_file
echo "Number of threads available: $(nproc)" >> $log_file
echo "Number of threads used: $max_jobs" >> $log_file
echo "Memory available: $(free -h | grep Mem | awk '{print $7}')" >> $log_file
echo "Disk space available in output folder drive: $(df -h $out_dir | tail -n 1 | awk '{print $4}')" >> $log_file
echo "" >> $log_file
echo "Software Used:" >> $log_file
fastp -v 2>&1 | head | tee -a $log_file
echo "bwa: $(bwa 2>&1 | head -n 3 | tail -n 1)" >> $log_file
echo "samtools: $(samtools --version | head -n 1)" >> $log_file
echo "ivar: $(ivar version | head -n 1)" >> $log_file
echo "fgbio: $(fgbio --version | head -n 1)" >> $log_file
echo "bcftools: $(bcftools --version | head -n 1)" >> $log_file
echo "qualimap: $(qualimap --version | head -n 4 | tail -n 1 )" >> $log_file
echo "parallel: $(parallel --version | head -n 1)" >> $log_file
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
mkdir -p $2/logs
initdir=$PWD

genomedir=$(basename $(dirname ${3}))
bwa index $3

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
  wait
  unpigz -c ${align}.gz |  samtools view -h -b | samtools sort -@ 1 -o $sorted
  wait
  start_time_samtools=$(date +%s)
  samtools index $sorted
  wait
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
python3 gene_coverage.py $2 $5
