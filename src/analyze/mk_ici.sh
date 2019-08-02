cd $BSCRATCH/at.sc.db/

lib_names=($(awk -F, '{if(NR>1) print $1}' data/sample_metadata.csv))

len=${#lib_names[@]}
i=0
cb_len=0
umi_len=0
echo "" > src/analyze/launch_ici.sh

while [ "$i" -lt "$len" ]; do
  lib_name=${lib_names[$i]}
  
  echo "source('"/global/projectb/scratch/byu24/at.sc.db/src/analyze/ici.R"')" > src/analyze/${lib_name}_ici.r
  
  echo "setwd('"/global/projectb/scratch/byu24/at.sc.db/scratch"')" >> src/analyze/${lib_name}_ici.r
  echo "library('"Seurat"')" >> src/analyze/${lib_name}_ici.r
  echo "library(ggplot2)" >> src/analyze/${lib_name}_ici.r
  echo "library(tidyverse)" >> src/analyze/${lib_name}_ici.r
  echo "library(furrr)" >> src/analyze/${lib_name}_ici.r
  echo "library(future)" >> src/analyze/${lib_name}_ici.r
  echo "options(future.globals.maxSize = 100000 * 1024^2)" >> src/analyze/${lib_name}_ici.r
  echo "" >> src/analyze/${lib_name}_ici.r

  echo "${lib_name}<-readRDS(file = '"/global/projectb/scratch/byu24/at.sc.db/scratch/robjects/${lib_name}.rds"')" >> src/analyze/${lib_name}_ici.r
  echo "ici_${lib_name}<- get_ICI(${lib_name})" >> src/analyze/${lib_name}_ici.r
  echo "saveRDS(ici_${lib_name}, file = '"/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/robjects/ici_${lib_name}.rds"')" >> src/analyze/${lib_name}_ici.r
  echo "" >> src/analyze/${lib_name}_ici.r

  echo "#!/bin/bash" > src/analyze/${lib_name}_ici.sh
  echo "#SBATCH -C skylake" >> src/analyze/${lib_name}_ici.sh
  echo "#SBATCH -A gtrnd" >> src/analyze/${lib_name}_ici.sh
  echo "#SBATCH -q jgi_exvivo" >> src/analyze/${lib_name}_ici.sh
  echo "#SBATCH -J ici_${lib_name}" >> src/analyze/${lib_name}_ici.sh
  echo "#SBATCH -t 12:00:00" >> src/analyze/${lib_name}_ici.sh
  echo "#SBATCH --mem-per-cpu=2000" >> src/analyze/${lib_name}_ici.sh
  echo "#SBATCH --ntasks=1" >> src/analyze/${lib_name}_ici.sh
  echo "#SBATCH --cpus-per-task=8" >> src/analyze/${lib_name}_ici.sh
  echo "#SBATCH --output=$BSCRATCH/at.sc.db/log/{lib_name}_ici.out" >> src/analyze/${lib_name}_ici.sh
  echo "" >> src/analyze/${lib_name}_ici.sh
  
  echo "module load python/3.7-anaconda-2019.07" >> src/analyze/${lib_name}_ici.sh
  echo "source activate /global/projectb/scratch/bjcole/env/scRNAseq2" >> src/analyze/${lib_name}_ici.sh
  echo "cd $BSCRATCH/at.sc.db/" >> src/analyze/${lib_name}_ici.sh
  echo "" >> src/analyze/${lib_name}_ici.sh
  
  echo "Rscript --verbose $BSCRATCH/at.sc.db/src/analyze/${lib_name}_ici.r >> $BSCRATCH/at.sc.db/log/{lib_name}_ici.Rout" >> src/analyze/${lib_name}_ici.sh

  echo "sbatch $BSCRATCH/at.sc.db/src/analyze/${lib_name}_ici.sh" >> src/analyze/launch_ici.sh

  i=$(($i + 1))
 
done
