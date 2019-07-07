#!/bin/bash
# Initialize our own variables:
output=""
in1=""
in2=""
run_type=""
run_name=""
parent_dir=""
genome=""
bbtools_dir="$BSCRATCH/bin/bbmap"
dropseq_dir="$BSCRATCH/bin/Drop-seq_tools-2.3.0"
# Fetch parameters

for i in "$@"
do
case $i in
    -o=*|--output_file=*)
    output="${i#*=}"
    shift # past argument=value
    ;;
    -a=*|--in1=*)
    in1="${i#*=}"
    shift # past argument=value
    ;;
    -b=*|--in2=*)
    in2="${i#*=}"
    shift # past argument=value
    ;;
    -r=*|--run_type=*)
    run_type="${i#*=}"
    shift # past argument=value
    ;;
    -n=*|--run_name=*)
    run_name="${i#*=}"
    shift # past argument=value
    ;;
    -g=*|--genome=*)
    genome="${i#*=}"
    shift # past argument=value
    ;;
    -p=*|--parent_dir=*)
    parent_dir="${i#*=}"
    shift # past argument=value
    ;;
    *)
esac
done

genomedir=$BSCRATCH/at.sc.db/scratch/at10
genomefa=${genomedir}/${genome}.fa

echo "Setting up directory ${parent_dir}/${run_name}"

mkdir $parent_dir
mkdir ${parent_dir}\/${run_name}

file_prefix=${parent_dir}/${run_name}/${run_name}

python3 $BSCRATCH/at.sc.db/src/map/tag_fastq_file.py --in1=${in1} --in2=${in2} --out=${file_prefix}_step01.sam --run=${run_name}

export JAVA_HOME=$BSCRATCH/bin/jdk-12.0.1
export PATH=$JAVA_HOME/bin:$PATH
PICARD='$BSCRATCH/bin/picard/picard.jar'
alias picard="java -jar $PICARD"

picard SortSam \
  INPUT=${file_prefix}_step01.sam \
  OUTPUT=${file_prefix}_step02.bam \
  SORT_ORDER=queryname \
  TMP_DIR=${file_prefix}_tmp

STAR \
  --genomeDir ${genomedir} \
  --readFilesIn ${in2} \
  --outFileNamePrefix ${file_prefix}_step03_ \
  --runThreadN 4

picard SortSam \
  INPUT=${file_prefix}_step03_Aligned.out.sam \
  OUTPUT=${file_prefix}_step04.bam \
  SORT_ORDER=queryname \
  TMP_DIR=${file_prefix}_tmp

picard MergeBamAlignment \
 REFERENCE_SEQUENCE=${genomedir}/${genome}.fa \
  UNMAPPED_BAM=${file_prefix}_step02.bam \
  ALIGNED_BAM=${file_prefix}_step04.bam \
  INCLUDE_SECONDARY_ALIGNMENTS=false \
  PAIRED_RUN=false \
  CLIP_ADAPTERS=false \
  TMP_DIR=${file_prefix}_tmp \
  OUTPUT=${file_prefix}_step06.bam

${dropseq_dir}/TagReadWithInterval \
  I=${file_prefix}_step06.bam \
  O=${file_prefix}_step07.bam \
  TMP_DIR=${file_prefix}_tmp \
  INTERVALS=${genomedir}/${genome}.genes.intervals \
  TAG=XG

${dropseq_dir}/TagReadWithGeneFunction \
  ANNOTATIONS_FILE=${genomedir}/${genome}.refFlat \
  INPUT=${file_prefix}_step07.bam \
  O=${file_prefix}_step08.bam

${dropseq_dir}/DetectBeadSubstitutionErrors \
  INPUT=${file_prefix}_step08.bam \
  OUTPUT=${file_prefix}_step09.bam \
  TMP_DIR=${file_prefix}_tmp \
  MIN_UMIS_PER_CELL=20 \
  OUTPUT_REPORT=${file_prefix}_substitution_error_report.txt

${dropseq_dir}/DetectBeadSynthesisErrors \
  INPUT=${file_prefix}_step09.bam \
  MIN_UMIS_PER_CELL=20 \
  OUTPUT_STATS=${file_prefix}_synthesis_error_stats.txt \
  SUMMARY=${file_prefix}_synthesis_error_summary.txt \
  REPORT=${file_prefix}_synthesis_error_report.txt \
  CREATE_INDEX=true \
  TMP_DIR=${file_prefix}_tmp \
  OUTPUT=${file_prefix}_final.bam


