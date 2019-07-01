cd $BSCRATCH/at.sc.db/

srr_names=($(awk -F, '{if(NR>1) print $9}' data/sample_metadata.csv))
set_names=($(awk -F, '{if(NR>1) print $1}' data/sample_metadata.csv))
run_names=($(awk -F, '{if(NR>1) print $8}' data/sample_metadata.csv))

len=${#srr_names[@]}
i=0

while [ "$i" -lt "$len" ]; do
  set_name=${set_names[$i]}
  srr_name=${srr_names[$i]}
  run_name=${run_names[$i]}

  echo "#!/bin/bash" > src/alevin/alevin_${set_name}.bs
  echo "#SBATCH --mail-user=brendayu@lbl.gov" >> src/alevin/alevin_${set_name}.bs
  echo "#SBATCH --mail-type=ALL" >> src/alevin/alevin_${set_name}.bs
  echo "#SBATCH -A gtrnd" >> src/alevin/alevin_${set_name}.bs
  echo "#SBATCH -q genepool_shared" >> src/alevin/alevin_${set_name}.bs
  echo "#SBATCH -J alevin_SRA_${set_name}" >> src/alevin/alevin_${set_name}.bs
  echo "#SBATCH -t 12:00:00" >> src/alevin/alevin_${set_name}.bs
  echo "#SBATCH --mem-per-cpu=4000" >> src/alevin/alevin_${set_name}.bs
  echo "#SBATCH --ntasks=1" >> src/alevin/alevin_${set_name}.bs
  echo "#SBATCH --cpus-per-task=4" >> src/alevin/alevin_${set_name}.bs
  echo "#SBATCH --output=$BSCRATCH/at.sc.db/log/alevin_SRA_${set_name}.out" >> src/alevin/alevin_${set_name}.bs
  echo "" >> src/alevin/alevin_${set_name}.bs
  echo "cd $BSCRATCH/at.sc.db/" >> src/alevin/alevin_${set_name}.bs

  echo "" >> src/alevin/alevin_${set_name}.bs
  echo "sh src/salmon_pipeline.sh ${set_name} ${run_name}" >> src/alevin/alevin_${set_name}.bs

  i=$(($i + 1))
done

awk -F, '{if(NR > 1) print "sbatch $BSCRATCH/at.sc.db/src/alevin/alevin_"$1".bs &"}' data/sample_metadata.csv > src/alevin/launch_alevin.sh
