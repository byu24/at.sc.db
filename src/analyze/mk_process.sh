cd $BSCRATCH/at.sc.db/

lib_names=($(awk -F, '{if(NR>1) print $1}' data/sample_metadata.csv))
cell_num=($(awk -F, '{if(NR>1) print $13}' data/sample_metadata.csv))

len=${#lib_names[@]}
i=0
cb_len=0
umi_len=0
echo "" > src/analyze/launch_process.sh

while [ "$i" -lt "$len" ]; do
  lib_name=${lib_names[$i]}
  cell_num=${cell_num[$i]}
  
  echo "if(!is.null(dev.list())) dev.off()" > src/analyze/process_${lib_name}.r
  echo "cat('"\014"')" > src/analyze/process_${lib_name}.r
  echo "rm(list=ls())" > src/analyze/process_${lib_name}.r
  echo "source('"/global/projectb/scratch/byu24/at.sc.db/src/analyze/process_cluster.r"')" > src/analyze/process_${lib_name}.r
  
  echo "setwd('"/global/projectb/scratch/byu24/at.sc.db/scratch"')" >> src/analyze/process_${lib_name}.r
  echo "library('"Seurat"')" >> src/analyze/process_${lib_name}.r
  echo "library(ggplot2)" >> src/analyze/process_${lib_name}.r
  echo "library(tidyverse)" >> src/analyze/process_${lib_name}.r
  echo "library(furrr)" >> src/analyze/process_${lib_name}.r
  echo "library(future)" >> src/analyze/process_${lib_name}.r
  echo "options(future.globals.maxSize = 100000 * 1024^2) #Extends memory size allowed in R" >> src/analyze/process_${lib_name}.r
  echo "" >> src/analyze/process_${lib_name}.r

  echo "${lib_name}.read = get_dge('"${lib_name}"')" >> src/analyze/process_${lib_name}.r
  echo "${lib_name}.data = get_dge_stats(${lib_name}.read)" >> src/analyze/process_${lib_name}.r
  echo "${lib_name}.filtered = filter_dge(${lib_name}.read,${lib_name}.data, expected_cells = ${cell_num})" >> src/analyze/process_${lib_name}.r
  echo "${lib_name} = get_sobj(${lib_name}.filtered, ${lib_name}.data, group = '"${lib_name}"')" >> src/analyze/process_${lib_name}.r
  echo "${lib_name} = trans_dge(${lib_name})" >> src/analyze/process_${lib_name}.r
  echo "saveRDS(${lib_name}, file = "/global/projectb/scratch/byu24/at.sc.db/scratch/robjects/${lib_name}.rds")" >> src/analyze/process_${lib_name}.r
  echo "${lib_name}" >> src/analyze/process_${lib_name}.r

  echo "DimPlot(${lib_name}, reduction = "umap")" >> src/analyze/process_${lib_name}.r
  echo "ggsave(file="/global/projectb/scratch/byu24/at.sc.db/reports/samples/umap_${lib_name}.png", width = 20, height = 20, units = "cm")" >> src/analyze/process_${lib_name}.r
  echo "ElbowPlot(${lib_name})" >> src/analyze/process_${lib_name}.r
  echo "ggsave(file="/global/projectb/scratch/byu24/at.sc.db/reports/samples/elbow_${lib_name}.png", width = 20, height = 20, units = "cm")" >> src/analyze/process_${lib_name}.r
  
  echo "source('"/global/projectb/scratch/byu24/at.sc.db/src/analyze/process_ici.R"')" > src/analyze/process_${lib_name}.r
  echo "ici_${lib_name}<- get_ICI(${lib_name})" >> src/analyze/process_${lib_name}.r
  echo "ici_${lib_name}<- summarize_ici(${lib_name})" >> src/analyze/process_${lib_name}.r
  echo "saveRDS(ici_${lib_name}, file = '"/global/projectb/scratch/byu24/at.sc.db/scratch/robjects/ici_${lib_name}.rds"')" >> src/analyze/process_${lib_name}.r
  echo "" >> src/analyze/process_${lib_name}.r

  echo "#!/bin/bash" > src/analyze/process_${lib_name}.sh
  echo "#SBATCH -C skylake" >> src/analyze/process_${lib_name}.sh
  echo "#SBATCH -A gtrnd" >> src/analyze/process_${lib_name}.sh
  echo "#SBATCH -q jgi_exvivo" >> src/analyze/process_${lib_name}.sh
  echo "#SBATCH -J ${lib_name}_process" >> src/analyze/process_${lib_name}.sh
  echo "#SBATCH -t 12:00:00" >> src/analyze/process_${lib_name}.sh
  echo "#SBATCH --mem-per-cpu=4000" >> src/analyze/process_${lib_name}.sh
  echo "#SBATCH --ntasks=1" >> src/analyze/process_${lib_name}.sh
  echo "#SBATCH --cpus-per-task=8" >> src/analyze/process_${lib_name}.sh
  echo "#SBATCH --output=$BSCRATCH/at.sc.db/log/process_${lib_name}.out" >> src/analyze/process_${lib_name}.sh
  echo "" >> src/analyze/process_${lib_name}.sh
  
  echo "module load python/3.7-anaconda-2019.07" >> src/analyze/process_${lib_name}.sh
  echo "source activate /global/projectb/scratch/bjcole/env/scRNAseq2" >> src/analyze/process_${lib_name}.sh
  echo "cd $BSCRATCH/at.sc.db/" >> src/analyze/process_${lib_name}.sh
  echo "" >> src/analyze/process_${lib_name}.sh
  
  echo "Rscript --verbose $BSCRATCH/at.sc.db/src/analyze/process_${lib_name}.r >> $BSCRATCH/at.sc.db/log/${lib_name}_ici.Rout" >> src/analyze/process_${lib_name}.sh

  echo "cd $BSCRATCH/at.sc.db/" >> src/analyze/launch_process.sh
  echo "mkdir scratch/robjects" >> src/analyze/launch_process.sh
  echo "mkdir reports/samples" >> src/analyze/launch_process.sh
  echo "sbatch $BSCRATCH/at.sc.db/src/analyze/process_${lib_name}.sh" >> src/analyze/launch_process.sh

  i=$(($i + 1))
 
done
