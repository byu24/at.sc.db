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
dropseq_dir="$BSCRATCH/Drop-seq_tools-2.3.0"
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

genomedir=$BSCRATCH/at.sc.db/scratch/genomes/${genome}
genomefa=${genomedir}/${genome}.fa

echo "Setting up directory ${parent_dir}/${run_name}"

mkdir $parent_dir
mkdir ${parent_dir}\/${run_name}

file_prefix=${parent_dir}/${run_name}/${run_name}

#python3 src/map/tag_fastq_file.py --in1=${in1} --in2=${in2} --out=${file_prefix}_step01.sam --run=${run_name}

picard FastqToSam \
  FASTQ=${in1} \
  FASTQ2=${in2} \
  OUTPUT=${file_prefix}_step01.bam \
  SAMPLE_NAME=${run_name}

${dropseq_dir}/TagBamWithReadSequenceExtended \
  SUMMARY=${file_prefix}_unaligned_tagged_Cellular.bam_summary.txt \
  BASE_RANGE=1-12 \
  BASE_QUALITY=10 \
  BARCODED_READ=1 \
  DISCARD_READ=false \
  TAG_NAME=XC \
  NUM_BASES_BELOW_QUALITY=1 \
  INPUT=${file_prefix}_step01.bam \
  OUTPUT=${file_prefix}_step02.bam

${dropseq_dir}/TagBamWithReadSequenceExtended \
  SUMMARY=${file_prefix}_unaligned_tagged_Molecular.bam_summary.txt.txt \
  BASE_RANGE=13-20 \
  BASE_QUALITY=10 \
  BARCODED_READ=1 \
  DISCARD_READ=true \
  TAG_NAME=XM \
  NUM_BASES_BELOW_QUALITY=1 \
  INPUT=${file_prefix}_step02.bam \
  OUTPUT=${file_prefix}_step03.bam

${dropseq_dir}/FilterBam \
  TAG_REJECT=XQ \
  INPUT=${file_prefix}_step03.bam \
  OUTPUT=${file_prefix}_step04.bam

${dropseq_dir}/TrimStartingSequence \
  OUTPUT_SUMMARY=${file_prefix}_adapter_trimming_report.txt \
  SEQUENCE=AAGCAGTGGTATCAACGCAGAGTGAATGGG \
  MISMATCHES=0 \
  NUM_BASES=5 \
  INPUT=${file_prefix}_step04.bam \
  OUTPUT=${file_prefix}_step05.bam

${dropseq_dir}/PolyATrimmer \
  OUTPUT=${file_prefix}_step06.bam \
  OUTPUT_SUMMARY=${file_prefix}_polyA_trimming_report.txt \
  MISMATCHES=0 \
  NUM_BASES=6 \
  NEW=true \
  INPUT=${file_prefix}_step05.bam

picard SamToFastq \
  INPUT=${file_prefix}_step06.bam \
  FASTQ=${file_prefix}_step07.fastq

STAR \
  --genomeDir ${genomedir} \
  --readFilesIn ${file_prefix}_step07.fastq \
  --outFileNamePrefix ${file_prefix}_step08_ \
  --runThreadN 8

picard SortSam \
  INPUT=${file_prefix}_step08_Aligned.out.sam \
  OUTPUT=${file_prefix}_step09.bam \
  SORT_ORDER=queryname \
  TMP_DIR=${file_prefix}_tmp

picard MergeBamAlignment \
 REFERENCE_SEQUENCE=${genomedir}/${genome}.fa \
  UNMAPPED_BAM=${file_prefix}_step06.bam \
  ALIGNED_BAM=${file_prefix}_step09.bam \
  INCLUDE_SECONDARY_ALIGNMENTS=false \
  PAIRED_RUN=false \
  CLIP_ADAPTERS=false \
  TMP_DIR=${file_prefix}_tmp \
  OUTPUT=${file_prefix}_step10.bam

${dropseq_dir}/TagReadWithInterval \
  I=${file_prefix}_step10.bam \
  O=${file_prefix}_step11.bam \
  TMP_DIR=${file_prefix}_tmp \
  INTERVALS=${genomedir}/${genome}.genes.intervals \
  TAG=XG

${dropseq_dir}/TagReadWithGeneFunction \
  ANNOTATIONS_FILE=${genomedir}/${genome}.refFlat \
  INPUT=${file_prefix}_step11.bam \
  O=${file_prefix}_step12.bam

${dropseq_dir}/DetectBeadSubstitutionErrors \
  INPUT=${file_prefix}_step12.bam \
  OUTPUT=${file_prefix}_step13.bam \
  TMP_DIR=${file_prefix}_tmp \
  MIN_UMIS_PER_CELL=20 \
  OUTPUT_REPORT=${file_prefix}_substitution_error_report.txt

${dropseq_dir}/DetectBeadSynthesisErrors \
  INPUT=${file_prefix}_step13.bam \
  MIN_UMIS_PER_CELL=20 \
  OUTPUT_STATS=${file_prefix}_synthesis_error_stats.txt \
  SUMMARY=${file_prefix}_synthesis_error_summary.txt \
  REPORT=${file_prefix}_synthesis_error_report.txt \
  CREATE_INDEX=true \
  TMP_DIR=${file_prefix}_tmp \
  OUTPUT=${file_prefix}_final.bam
