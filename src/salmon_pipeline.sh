#!/bin/bash
prefix=$1
technology=$2

echo $prefix 
echo $technology

salmon_path=$BSCRATCH/bin/salmon/bin
index_path=$BSCRATCH/at.sc.db/scratch/at10
file_prefix=$BSCRATCH/at.sc.db/scratch/${prefix}
data_path=$BSCRATCH/at.sc.db/data

export PATH=$PATH:$BSCRATCH/bin/salmon/bin

if [[ $technology -eq "DropSeq" ]]; then
  ${salmon_path}/salmon alevin -l ISR -1 ${file_prefix}_rRNA_adapter_filtered_1.fastq -2 ${file_prefix}_R2.fastq --dropseq -i ${index_path}/at10_index -p 10 -o ${data_path}/${prefix} --tgMap ${index_path}/at10_tgMap.tsv --dumpMtx ${data_path}/${prefix}
elif [[ $technology -eq "10xV2" ]]; then
  ${salmon_path}/salmon alevin -l ISR -1 ${file_prefix}_rRNA_adapter_filtered_1.fastq -2 ${file_prefix}_R2.fastq --chromium -i ${index_path}/at10_index -p 10 -o ${data_path}/${prefix} --tgMap ${index_path}/at10_tgMap.tsv --dumpMtx ${data_path}/${prefix}
elif [[$technology -eq "10xV3" ]]; then
  ${salmon_path}/salmon alevin -l ISR -1 ${file_prefix}_rRNA_adapter_filtered_1.fastq -2 ${file_prefix}_R2.fastq --chromium -i ${index_path}/at10_index -p 10 -o ${data_path}/${prefix} --tgMap ${index_path}/at10_tgMap.tsv --dumpMtx ${data_path}/${prefix}
fi


