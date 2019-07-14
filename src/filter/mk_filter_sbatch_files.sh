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

  echo "#!/bin/bash" > src/filter/filter_${set_name}.bs
  echo "#SBATCH -A gtrnd" >> src/filter/filter_${set_name}.bs
  echo "#SBATCH -q genepool" >> src/filter/filter_${set_name}.bs
  echo "#SBATCH -J filter_SRA_${set_name}" >> src/filter/filter_${set_name}.bs
  echo "#SBATCH -t 12:00:00" >> src/filter/filter_${set_name}.bs
  echo "#SBATCH --mem-per-cpu=2000" >> src/filter/filter_${set_name}.bs
  echo "#SBATCH --ntasks=1" >> src/filter/filter_${set_name}.bs
  echo "#SBATCH --cpus-per-task=16" >> src/filter/filter_${set_name}.bs
  echo "#SBATCH --output=$BSCRATCH/at.sc.db/log/filter_SRA_${set_name}.out" >> src/filter/filter_${set_name}.bs
  echo "" >> src/filter/filter_${set_name}.bs
  echo "module load python3" >> src/filter/filter_${set_name}.bs
  echo "source activate \$BSCRATCH/bin/env_STARsolo" >> src/filter/filter_${set_name}.bs

  echo "cd \$BSCRATCH/at.sc.db/" >> src/filter/filter_${set_name}.bs

  echo "" >> src/filter/filter_${set_name}.bs
  echo "sh $BSCRATCH/at.sc.db/src/trim_filter.sh ${set_name} ${run_name}" >> src/filter/filter_${set_name}.bs

  i=$(($i + 1))
done

awk -F, '{if(NR > 1) print "sbatch $BSCRATCH/at.sc.db/src/filter/filter_"$1".bs &"}' data/sample_metadata.csv > $BSCRATCH/at.sc.db/src/filter/launch_filter.sh
