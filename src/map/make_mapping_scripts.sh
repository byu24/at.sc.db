cd $BSCRATCH/at.sc.db/

lib_names=($(awk -F, '{if(NR>1) print $1}' data/sample_metadata.csv))
lib_types=($(awk -F, '{if(NR>1) print $8}' data/sample_metadata.csv))

len=${#lib_names[@]}
i=0
echo "" > src/map/launch_mapping_scripts.sh

while [ "$i" -lt "$len" ]; do
  lib_name=${lib_names[$i]}
  lib_type=${lib_types[$i]}

  echo "#!/bin/bash" > src/map/${lib_name}_run.bs
  echo "#SBATCH -A gtrnd" >> src/map/${lib_name}_run.bs
  echo "#SBATCH -q genepool_shared" >> src/map/${lib_name}_run.bs
  echo "#SBATCH -J ${lib_name}_run" >> src/map/${lib_name}_run.bs
  echo "#SBATCH -t 6:00:00" >> src/map/${lib_name}_run.bs
  echo "#SBATCH --mem-per-cpu=8000" >> src/map/${lib_name}_run.bs
  echo "#SBATCH --ntasks=1" >> src/map/${lib_name}_run.bs
  echo "#SBATCH --cpus-per-task=4" >> src/map/${lib_name}_run.bs
  echo "#SBATCH --output=$BSCRATCH/at.sc.db/log/${lib_name}_run.out" >> src/map/${lib_name}_run.bs
  echo "" >> src/map/${lib_name}_run.bs
 
  echo "cd $BSCRATCH/at.sc.db" >> src/map/${lib_name}_run.bs
  echo "mkdir $BSCRATCH/at.sc.db/scratch/${lib_name}" >> src/map/${lib_name}_run.bs
  echo "" >> src/map/${lib_name}_run.bs

  if [[ $lib_type == "DropSeq" ]]; then
    echo "module load python3" >> src/map/${lib_name}_run.bs
	echo "source activate /global/projectb/scratch/bjcole/env_STARsolo" >> src/map/${lib_name}_run.bs
	echo "sh src/map/map_dropseq.sh \\" >> src/map/${lib_name}_run.bs
    echo "  --in1=scratch/${lib_name}_R1.fastq \\" >> src/map/${lib_name}_run.bs
    echo "  --in2=scratch/${lib_name}_R2.fastq \\" >> src/map/${lib_name}_run.bs
    echo "  --run_type="DropSeq" \\" >> src/map/${lib_name}_run.bs
    echo "  --run_name="${lib_name}" \\" >> src/map/${lib_name}_run.bs
    echo "  --genome=at10 \\" >> src/map/${lib_name}_run.bs
    echo "  --parent_dir=scratch/" >> src/map/${lib_name}_run.bs
    echo "" >> src/map/${lib_name}_run.bs
    echo "awk '{if(NR > 1) {print \$1}}' scratch/${lib_name}/${lib_name}_synthesis_error_stats.txt > \\" >> src/map/${lib_name}_run.bs
    echo "  scratch/${lib_name}/${lib_name}_whitelist.csv" >> src/map/${lib_name}_run.bs
    echo "" >> src/map/${lib_name}_run.bs
  fi

  if [[ $lib_type == "TenxV2" ]]; then
    echo "cat $BSCRATCH/at.sc.db/data/whitelist.txt > \\" >> src/map/${lib_name}_run.bs
    echo "  $BSCRATCH/at.sc.db/scratch/${lib_name}/${lib_name}_whitelist.csv" >> src/map/${lib_name}_run.bs
	echo "module load python3" >> src/map/${lib_name}_run.bs	
	echo "source activate $BSCRATCH/bin/env_STAR" >> src/map/${lib_name}_run.bs
  fi

  echo "export PATH=$PATH:$BSCRATCH/bin/salmon/bin" >> src/map/${lib_name}_run.bs
  echo "salmon alevin -l ISR \\" >> src/map/${lib_name}_run.bs
  echo "  -1 scratch/${lib_name}_R1.fastq \\" >> src/map/${lib_name}_run.bs
  echo "  -2 scratch/${lib_name}_R2.fastq \\" >> src/map/${lib_name}_run.bs
  
  if [[ $lib_type == "DropSeq" ]]; then
    echo "  --dropseq  \\" >> src/map/${lib_name}_run.bs
  fi
  if [[ $lib_type == "TenxV2" ]]; then
    echo "  --chromium \\" >> src/map/${lib_name}_run.bs
  fi

  echo "  -i $BSCRATCH/at.sc.db/scratch/at10 \\" >> src/map/${lib_name}_run.bs 
  echo "  -p 2 -o scratch/${lib_name}/ \\" >> src/map/${lib_name}_run.bs
  echo "  --dumpMtx \\" >> src/map/${lib_name}_run.bs
  echo "  --tgMap $BSCRATCH/at.sc.db/scratch/at10/at10_tgMap.tsv \\" >> src/map/${lib_name}_run.bs
  echo "  --mrna $BSCRATCH/at.sc.db/scratch/at10/at10_ptGenes.tsv \\" >> src/map/${lib_name}_run.bs
  echo "  --rrna $BSCRATCH/at.sc.db/scratch/at10/at10_riboGenes.tsv \\" >> src/map/${lib_name}_run.bs
  echo "  --whitelist scratch/${lib_name}/${lib_name}_whitelist.csv" >> src/map/${lib_name}_run.bs

  echo "sbatch $BSCRATCH/at.sc.db/src/map/${lib_name}_run.bs &" >> $BSCRATCH/at.sc.db/src/map/launch_mapping_scripts.sh

  i=$(($i + 1))
 
done
