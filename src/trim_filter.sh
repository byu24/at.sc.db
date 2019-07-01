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

#!/bin/bash
prefix=$1
technology=$2
bbmap_path=$BSCRATCH/bin/bbmap

echo $prefix 
echo $technology

if [[ $technology -eq "DropSeq" ]]; then
  cb_len=12
  umi_len=8
elif [[ $technology -eq "TenxV2" ]]; then
  cb_len=16
  umi_len=10
elif [[$technology -eq "TenxV3" ]]; then
  cb_len=16
  umi_len=12
fi

file_prefix=$BSCRATCH/at.sc.db/scratch/${prefix}

${bbmap_path}/bbduk.sh \
  in1=${file_prefix}_1.fastq \
  in2=${file_prefix}_2.fastq \
  out1=${file_prefix}_rRNA_filtered_1.fastq \
  out2=${file_prefix}_rRNA_filtered_2.fastq \
  k=21 edist=2 ref=$BSCRATCH/at.sc.db/scratch/at10/rrna_seqs.fa

${bbmap_path}/bbduk.sh \
  in1=${file_prefix}_rRNA_filtered_1.fastq \
  in2=${file_prefix}_rRNA_filtered_2.fastq \
  out1=${file_prefix}_rRNA_adapter_filtered_1.fastq \
  out2=${file_prefix}_R2.fastq \
  ref=${bbmap_path}/resources/adapters.fa \
  k=21 mink=11 minlen=$((cb_len + umi_len - 1)) edist=2 ktrim=r qtrim=30
