#!/bin/bash
prefix=$1
technology=$2
bbmap_path=$BSCRATCH/bin/bbmap

echo $prefix 
echo $technology

if [[ $technology == "DropSeq" ]]; then
  cb_len=12
  umi_len=8
elif [[ $technology == "10x_V2" ]]; then
  cb_len=16
  umi_len=10
elif [[$technology == "10x_V3" ]]; then
  cb_len=16
  umi_len=12
fi

file_prefix=$BSCRATCH/at.sc.db/scratch/${prefix}/${prefix}

#${bbmap_path}/bbduk.sh \
#  in1=${file_prefix}_1.fastq \
#  in2=${file_prefix}_2.fastq \
#  out1=${file_prefix}_rRNA_filtered_1.fastq \
#  out2=${file_prefix}_rRNA_filtered_2.fastq \
#  k=21 edist=2 ref=$BSCRATCH/at.sc.db/scratch/at10/rrna_seqs.fa

${bbmap_path}/bbduk.sh \
  in1=${file_prefix}_1.fastq \
  in2=${file_prefix}_2.fastq \
  out1=${file_prefix}_R1.fastq \
  out2=${file_prefix}_R2.fastq \
  ref=${bbmap_path}/resources/adapters.fa \
  k=21 mink=11 minlen=$((cb_len + umi_len)) edist=2 ktrim=r qtrim=10
