#!/bin/bash
prefix=$1 #sample number according to your naming convention from the metadata CSV
technology=$2 #scRNAseq technique (DropSeq or 10xGenomics) from the metadata CSV
bbmap_path=$BSCRATCH/bin/bbmap #change to path for BBtools which includes bbmap

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

${bbmap_path}/reformat.sh in1=${file_prefix}_1.fastq in2=${file_prefix}_2.fastq out=stdout.fq addslash=t | \
  ${bbmap_path}/shuffle2.sh in=stdin.fq out=stdout.fq -Xmx4G -da int=t | \
  ${bbmap_path}/bbduk.sh in=stdin.fq out=stdout.fq ref=${bbmap_path}/resources/adapters.fa k=23 mink=11 hdist=1 ktrim=r qtrim=10 -Xmx1g t=16 int=t | \
  ${bbmap_path}/bbduk.sh in=stdin.fq out=stdout.fq k=23 mink=11 literal=AAGCAGTGGTATCAACGCAGAGTGAATGGG hdist=1 ktrim=l -Xmx1g t=16 int=t | \
  ${bbmap_path}/reformat.sh in=stdin.fq out1=${file_prefix}_tmp_1.fastq out2=${file_prefix}_R2.fastq int=t minlen=$((cb_len + umi_len))

${bbmap_path}/bbduk.sh in=${file_prefix}_tmp_1.fastq \
  out=${file_prefix}_R1.fastq \
  ftr=$((cb_len+umi_len-1)) \
  -Xmx4g t=1
