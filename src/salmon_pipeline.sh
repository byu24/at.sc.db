#!/bin/bash
prefix=$1
technology=$2

echo $prefix 
echo $technology

salmon_path=$BSCRATCH/bin/salmon/bin
index_path=$BSCRATCH/bin/at.sc.db/scratch/at10
file_prefix=$BSCRATCH/at.sc.db/scratch/${prefix}
data_path=$BSCRATCH/at.sc.db/data

if [[ $technology -eq "DropSeq" ]]; then
  protocol=dropseq
elif [[ $technology -eq "10x_V2" ]]; then
  protocol=chromium
elif [[$technology -eq "10x_V3" ]]; then
  protocol=chromium
fi

${salmon_path}/salmon alevin -l ISR -1 ${file_prefix}_rRNA_adapter_filtered_1.fastq -2 ${file_prefix}_R2.fastq --$(protocol) -i ${index_path}/at10_index -p 10 -o ${data_path}/${file_prefix} --tgMap ${index_path}/at10_tgMap.tsv --dumpMtx ${data_path}/${file_prefix}
