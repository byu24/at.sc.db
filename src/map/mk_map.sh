cd $BSCRATCH/at.sc.db/

lib_names=($(awk -F, '{if(NR>1) print $1}' data/sample_metadata.csv))
lib_types=($(awk -F, '{if(NR>1) print $8}' data/sample_metadata.csv))
references=($(awk -F, '{if(NR>1) print $11}' data/sample_metadata.csv))

len=${#lib_names[@]}
i=0
cb_len=0
umi_len=0
echo "" > src/map/launch_map.sh

while [ "$i" -lt "$len" ]; do
  lib_name=${lib_names[$i]}
  lib_type=${lib_types[$i]}
  reference=${references[$i]}

  echo "#!/bin/bash" > src/map/map_${lib_name}.bs
  echo "#SBATCH -A gtrnd" >> src/map/map_${lib_name}.bs
  echo "#SBATCH -q genepool_shared" >> src/map/map_${lib_name}.bs
  echo "#SBATCH -J map_${lib_name}" >> src/map/map_${lib_name}.bs
  echo "#SBATCH -t 24:00:00" >> src/map/map_${lib_name}.bs
  echo "#SBATCH --mem-per-cpu=4000" >> src/map/map_${lib_name}.bs
  echo "#SBATCH --ntasks=1" >> src/map/map_${lib_name}.bs
  echo "#SBATCH --cpus-per-task=8" >> src/map/map_${lib_name}.bs
  echo "#SBATCH --output=$BSCRATCH/at.sc.db/log/map_${lib_name}.out" >> src/map/map_${lib_name}.bs
  echo "" >> src/map/map_${lib_name}.bs
  echo "module load python/3.7-anaconda-2019.07" >> src/map/map_${lib_name}.bs
  echo "source activate /global/projectb/scratch/bjcole/env_STARsolo" >> src/map/map_${lib_name}.bs

  echo "cd $BSCRATCH/at.sc.db/" >> src/map/map_${lib_name}.bs
  echo "mkdir scratch/${lib_name}" >> src/map/map_${lib_name}.bs
  echo "" >> src/map/map_${lib_name}.bs

  if [[ $lib_type == "DropSeq" ]]; then
    echo "" >> src/map/map_${lib_name}.bs
    echo "sh $BSCRATCH/at.sc.db/src/map/map_dropseq.sh \\" >> src/map/map_${lib_name}.bs
    echo "  --in1=scratch/${lib_name}/${lib_name}_R1.fastq \\" >> src/map/map_${lib_name}.bs
    echo "  --in2=scratch/${lib_name}/${lib_name}_R2.fastq \\" >> src/map/map_${lib_name}.bs
    echo "  --run_type="DropSeq" \\" >> src/map/map_${lib_name}.bs
    echo "  --run_name="${lib_name}" \\" >> src/map/map_${lib_name}.bs
    echo "  --genome=${reference} \\" >> src/map/map_${lib_name}.bs
    echo "  --parent_dir=scratch/" >> src/map/map_${lib_name}.bs
    echo "" >> src/map/map_${lib_name}.bs
    echo "awk '{if(NR > 1) {print \$1}}' scratch/${lib_name}/${lib_name}_synthesis_error_stats.txt > \\" >> src/map/map_${lib_name}.bs
    echo "  scratch/${lib_name}/${lib_name}_whitelist.csv" >> src/map/map_${lib_name}.bs
    echo "" >> src/map/map_${lib_name}.bs
    cb_len=12
    umi_len=8
  fi

  if [[ $lib_type == "10x_V2" ]]; then
    echo "cat $BSCRATCH/at.sc.db/data/whitelist.txt > \\" >> src/map/map_${lib_name}.bs
    echo "  $BSCRATCH/at.sc.db/scratch/${lib_name}/${lib_name}_whitelist.csv" >> src/map/map_${lib_name}.bs
    cb_len=16
    umi_len=10
  fi

  echo "STAR \\" >> src/map/map_${lib_name}.bs
  echo "  --genomeDir $BSCRATCH/at.sc.db/scratch/genomes/${reference}/STAR \\" >> src/map/map_${lib_name}.bs
  echo "  --readFilesIn $BSCRATCH/at.sc.db/scratch/${lib_name}/${lib_name}_R2.fastq scratch/${lib_name}/${lib_name}_R1.fastq \\" >> src/map/map_${lib_name}.bs
  echo "  --soloCBwhitelist $BSCRATCH/at.sc.db/scratch/${lib_name}/${lib_name}_whitelist.csv \\" >> src/map/map_${lib_name}.bs
  echo "  --outFileNamePrefix $BSCRATCH/at.sc.db/scratch/${lib_name}/${lib_name}_star. \\" >> src/map/map_${lib_name}.bs
  echo "  --soloType Droplet \\" >> src/map/map_${lib_name}.bs
  echo "  --soloCBlen $cb_len \\" >> src/map/map_${lib_name}.bs
  echo "  --soloUMIlen $umi_len \\" >> src/map/map_${lib_name}.bs
  echo "  --soloUMIstart $((cb_len + 1)) \\" >> src/map/map_${lib_name}.bs
  echo "  --runThreadN 8" >> src/map/map_${lib_name}.bs
 
  echo "sbatch $BSCRATCH/at.sc.db/src/map/map_${lib_name}.bs &" >> src/map/launch_map.sh

  i=$(($i + 1))
 
done
