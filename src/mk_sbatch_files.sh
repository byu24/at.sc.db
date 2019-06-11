srr_names=($(awk -F, '{if(NR>1) print $9}' $BSCRATCH/at.sc.db/data/sample_metadata.csv))
set_names=($(awk -F, '{if(NR>1) print $1}' $BSCRATCH/at.sc.db/data/sample_metadata.csv))

#for i in {1 .. ${#srr_names[@]}}; do
#  echo ${#srr_names[@]}
#  #echo "$BSCRATCH/sratoolkit.2.9.6-1-ubuntu64/bin/prefetch.2.9.3 "$i -o tmp -t http -p 1 -v 
#done

len=${#srr_names[@]}
i=0
while [ "$i" -lt "$len" ]; do
  set_name=${set_names[$i]}
  srr_name=${srr_names[$i]}

  echo "#!/bin/bash" > download_${set_name}.bs
  echo "#SBATCH --mail-user=bjcole@lbl.gov" >> download_${set_name}.bs
  echo "#SBATCH --mail-type=ALL" >> download_${set_name}.bs
  echo "#SBATCH -A gtrnd" >> download_${set_name}.bs
  echo "#SBATCH -q genepool_shared" >> download_${set_name}.bs
  echo "#SBATCH -J download_SRA_${set_name}" >> download_${set_name}.bs
  echo "#SBATCH -t 12:00:00" >> download_${set_name}.bs
  echo "#SBATCH --mem-per-cpu=4000" >> download_${set_name}.bs
  echo "#SBATCH --ntasks=1" >> download_${set_name}.bs
  echo "#SBATCH --cpus-per-task=4" >> download_${set_name}.bs
  echo "#SBATCH --output=download_SRA_${set_name}.out" >> download_${set_name}.bs
  echo "" >> download_${set_name}.bs

  echo "$BSCRATCH/sratoolkit.2.9.6-1-ubuntu64/bin/prefetch.2.9.3 ${srr_name} -o $BSCRATCH/at.sc.db/scratch/${set_name} -t http" >> download_${set_name}.bs
  echo "$BSCRATCH/sratoolkit.2.9.6-1-ubuntu64/bin/fasterq-dump $BSCRATCH/at.sc.db/scratch/${set_name} \\" >> download_${set_name}.bs
  echo "  -o ${set_name}.fastq \\" >> download_${set_name}.bs
  echo "  -O $BSCRATCH/at.sc.db/scratch \\" >> download_${set_name}.bs
  echo "  -t $BSCRATCH/at.sc.db/scratch \\" >> download_${set_name}.bs
  echo "  -sS --include-technical" >> download_${set_name}.bs

  i=$(($i + 1))
done
