cd $BSCRATCH/at.sc.db/

lib_names=($(awk -F, '{if(NR>1) print $1}' data/sample_metadata.csv))
lib_types=($(awk -F, '{if(NR>1) print $8}' data/sample_metadata.csv))
references=($(awk -F, '{if(NR>1) print $11}' data/sample_metadata.csv))

len=${#lib_names[@]}
i=0
cb_len=0
umi_len=0
echo "" > src/map/launch_mapping_scripts.sh

while [ "$i" -lt "$len" ]; do
  lib_name=${lib_names[$i]}
  lib_type=${lib_types[$i]}
  reference=${references[$i]}

  echo "#!/bin/bash" > src/map/${lib_name}_run.bs
  echo "#SBATCH -A gtrnd" >> src/map/${lib_name}_run.bs
  echo "#SBATCH -q genepool_shared" >> src/map/${lib_name}_run.bs
  echo "#SBATCH -J ${lib_name}_run" >> src/map/${lib_name}_run.bs
  echo "#SBATCH -t 24:00:00" >> src/map/${lib_name}_run.bs
  echo "#SBATCH --mem-per-cpu=4000" >> src/map/${lib_name}_run.bs
  echo "#SBATCH --ntasks=1" >> src/map/${lib_name}_run.bs
  echo "#SBATCH --cpus-per-task=8" >> src/map/${lib_name}_run.bs
  echo "#SBATCH --output=$BSCRATCH/at.sc.db/log/${lib_name}_run.out" >> src/map/${lib_name}_run.bs
  echo "" >> src/map/${lib_name}_run.bs
  echo "module load python3" >> src/map/${lib_name}_run.bs
  echo "source activate $BSCRATCH/bin/env_STARsolo" >> src/map/${lib_name}_run.bs

  echo "cd $BSCRATCH/at.sc.db/" >> src/map/${lib_name}_run.bs
  echo "mkdir scratch/${lib_name}" >> src/map/${lib_name}_run.bs
  echo "" >> src/map/${lib_name}_run.bs

  if [[ $lib_type == "DropSeq" ]]; then
    echo "" >> src/map/${lib_name}_run.bs
    echo "sh $BSCRATCH/at.sc.db/src/map/map_dropseq.sh \\" >> src/map/${lib_name}_run.bs
    echo "  --in1=scratch/${lib_name}/${lib_name}_R1.fastq \\" >> src/map/${lib_name}_run.bs
    echo "  --in2=scratch/${lib_name}/${lib_name}_R2.fastq \\" >> src/map/${lib_name}_run.bs
    echo "  --run_type="DropSeq" \\" >> src/map/${lib_name}_run.bs
    echo "  --run_name="${lib_name}" \\" >> src/map/${lib_name}_run.bs
    echo "  --genome=${reference} \\" >> src/map/${lib_name}_run.bs
    echo "  --parent_dir=scratch/" >> src/map/${lib_name}_run.bs
    echo "" >> src/map/${lib_name}_run.bs
    echo "awk '{if(NR > 1) {print \$1}}' scratch/${lib_name}/${lib_name}_synthesis_error_stats.txt > \\" >> src/map/${lib_name}_run.bs
    echo "  scratch/${lib_name}/${lib_name}_whitelist.csv" >> src/map/${lib_name}_run.bs
    echo "" >> src/map/${lib_name}_run.bs
    cb_len=12
    umi_len=8
  fi

  if [[ $lib_type == "10x_V2" ]]; then
    echo "cat $BSCRATCH/at.sc.db/data/whitelist.txt > \\" >> src/map/${lib_name}_run.bs
    echo "  $BSCRATCH/at.sc.db/scratch/${lib_name}/${lib_name}_whitelist.csv" >> src/map/${lib_name}_run.bs
    cb_len=16
    umi_len=10
  fi

  echo "STAR \\" >> src/map/${lib_name}_run.bs
  echo "  --genomeDir scratch/genomes/${reference} \\" >> src/map/${lib_name}_run.bs
  echo "  --readFilesIn scratch/${lib_name}/${lib_name}_R2.fastq scratch/${lib_name}/${lib_name}_R1.fastq \\" >> src/map/${lib_name}_run.bs
  echo "  --soloCBwhitelist scratch/${lib_name}/${lib_name}_whitelist.csv \\" >> src/map/${lib_name}_run.bs
  echo "  --outFileNamePrefix scratch/${lib_name}/${lib_name}_star. \\" >> src/map/${lib_name}_run.bs
  echo "  --soloType Droplet \\" >> src/map/${lib_name}_run.bs
  echo "  --soloCBlen $cb_len \\" >> src/map/${lib_name}_run.bs
  echo "  --soloUMIlen $umi_len \\" >> src/map/${lib_name}_run.bs
  echo "  --soloUMIstart $((cb_len + 1)) \\" >> src/map/${lib_name}_run.bs
  echo "  --runThreadN 8" >> src/map/${lib_name}_run.bs
 
  echo "sbatch $BSCRATCH/at.sc.db/src/map/${lib_name}_run.bs &" >> src/map/launch_mapping_scripts.sh

  i=$(($i + 1))
 
done
