#!/bin/bash
#SBATCH --mail-user=brendayu@lbl.gov
#SBATCH --mail-type=ALL
#SBATCH -A gtrnd
#SBATCH -q genepool_shared
#SBATCH -J filter_dc_019
#SBATCH -t 12:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --output=filter_dc_019.out

prefix=$1
genome=$2
version=$3
genomedir=$BSCRATCH/at.sc.db/scratch/${prefix}/${genome}


## Trim and filter fastq files

$BSCRATCH/bin/bbmap/bbduk.sh \
  in1=$BSCRATCH/at.sc.db/scratch/dc_019_1.fastq \
  in2=$BSCRATCH/at.sc.db/scratch/dc_019_2.fastq \
  out1=$BSCRATCH/at.sc.db/scratch/filtered/dc_019_1_filtered.fastq \
  out2=$BSCRATCH/at.sc.db/scratch/filtered/dc_019_2_filtered.fastq \
  k=21 edist=2 ref=$BSCRATCH/at.sc.db/scratch/at10/rrna_seqs.fa \
  
$BSCRATCH/bin/bbmap/bbduk.sh \
  in1=$BSCRATCH/at.sc.db/scratch/filtered/dc_019_1_filtered.fastq \
  in2=$BSCRATCH/at.sc.db/scratch/filtered/dc_019_2_filtered.fastq \
  out1=$BSCRATCH/at.sc.db/scratch/filtered/dc_019_1_R1_polyA.fastq \
  out2=$BSCRATCH/at.sc.db/scratch/filtered/dc_019_2_R1_polyA.fastq \
  k=21 mink=11 minlen=30 \
  ref=$BSCRATCH/bin/bbmap/resources/adapters.fa \
  edist=2 ktrim=r qtrim=30 \

$BSCRATCH/bin/bbmap/bbduk.sh \
  in1=$BSCRATCH/at.sc.db/scratch/filtered/dc_019_1_filtered.fastq \
  in2=$BSCRATCH/at.sc.db/scratch/filtered/dc_019_2_filtered.fastq \
  out1=$BSCRATCH/at.sc.db/scratch/filtered/dc_019_1_R2.fastq \
  out2=$BSCRATCH/at.sc.db/scratch/filtered/dc_019_2_R2.fastq \
  k=21 mink=11 minlen=30 \
  ref=$BSCRATCH/bin/bbmap/resources/adapters.fa \
  edist=2 ktrim=r qtrim=30 \

$BSCRATCH/bin/bbmap/bbduk.sh \
  in1=$BSCRATCH/at.sc.db/scratch/filtered/dc_019_1_R1_polyA.fastq \
  in2=$BSCRATCH/at.sc.db/scratch/filtered/dc_019_2_R1_polyA.fastq \
  out1=$BSCRATCH/at.sc.db/scratch/filtered/dc_019_1_R1.fastq \
  out2=$BSCRATCH/at.sc.db/scratch/filtered/dc_019_2_R1.fastq \
  ftr=$((cb_len + umi_len - 1)) \


