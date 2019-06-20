#!/bin/bash
prefix=$1
technology=$2

if [[ $technology eq "ds" ]] then
  cb_len=12
  umi_len=8
elif [[ $technology eq "10x.V2" ]] then
  cb_len=16
  umi_len=10
elif [[$technology eq "10x.V3" ]] then
  cb_len=16
  umi_len=12
fi

file_prefix=$BSCRATCH/at.sc.db/scratch/${prefix}/${file_prefix}

$BSCRATCH/bin/bbmap/bbduk.sh \
  in1=${file_prefix}_1.fastq \
  in2=${file_prefix}_2.fastq \
  out1=${file_prefix}_rRNA_filtered_1.fastq \
  out2=${file_prefix}_rRNA_filtered_2.fastq \
  k=21 edist=2 ref=$BSCRATCH/at.sc.db/scratch/at10/rrna_seqs.fa

$BSCRATCH/bin/bbmap/bbduk.sh \
  in1=${file_prefix}_rRNA_filtered_1.fastq \
  in2=${file_prefix}_rRNA_filtered_2.fastq \
  out1=${file_prefix}_rRNA_adapter_filtered_1.fastq \
  out2=${file_prefix}_R2.fastq \
  ref=$BSCRATCH/bin/bbmap/resources/adapters.fa \
  k=21 mink=11 minlen=$((cb_len + umi_len - 1)) edist=2 ktrim=r qtrim=30
