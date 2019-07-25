cd $BSCRATCH/at.sc.db/

lib_names=($(awk -F, '{if(NR>1) print $1}' data/sample_metadata.csv))

len=${#lib_names[@]}
i=0
cb_len=0
umi_len=0
echo "" > src/analyze/launch_analyze_scripts.sh
echo "module load python3" >> src/analyze/launch_analyze_scripts.sh
echo "source activate /global/projectb/scratch/bjcole/env/scRNAseq2" >> src/analyze/launch_analyze_scripts.sh

while [ "$i" -lt "$len" ]; do
  lib_name=${lib_names[$i]}
  
  echo "#!/usr/bin/env" > src/analyze/${lib_name}_analyze.r
  echo "" >> src/analyze/${lib_name}_analyze.r
  
  echo "setwd("/global/projectb/scratch/byu24/at.sc.db/scratch")" >> src/analyze/${lib_name}_analyze.r
  echo "library("Seurat")" >> src/analyze/${lib_name}_analyze.r
  echo "library(ggplot2)" >> src/analyze/${lib_name}_analyze.r
  echo "" >> src/analyze/${lib_name}_analyze.r

  echo "${lib_name}.data<-Read10X(data.dir = "/global/projectb/scratch/byu24/at.sc.db/scratch/${lib_name}/${lib_name}_star.Solo.out")" >> src/analyze/${lib_name}_analyze.r
  echo "${lib_name}<- CreateSeuratObject(counts = ${lib_name}.data, project = "${lib_name}")" >> src/analyze/${lib_name}_analyze.r
  echo "${lib_name} <- NormalizeData(${lib_name}, normalization.method = "LogNormalize", scale.factor = 10000)" >> src/analyze/${lib_name}_analyze.r
  echo "" >> src/analyze/${lib_name}_analyze.r

  echo "${lib_name} <- FindVariableFeatures(${lib_name}, selection.method = "vst", nfeatures = 2000)" >> src/analyze/${lib_name}_analyze.r
  echo "${lib_name} <- ScaleData(${lib_name})" >> src/analyze/${lib_name}_analyze.r
  echo "" >> src/analyze/${lib_name}_analyze.r

  echo "${lib_name} <- RunPCA(${lib_name}, features = VariableFeatures(object = ${lib_name}))" >> src/analyze/${lib_name}_analyze.r
  echo "print(lib_name[["pca"]], dims = 1:5, nfeatures = 5)" >> src/analyze/${lib_name}_analyze.r
  echo "" >> src/analyze/${lib_name}_analyze.r

  echo "png(file="/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/${lib_name}_vizdim.png")" >> src/analyze/${lib_name}_analyze.r
  echo "VizDimLoadings(${lib_name}, dims = 1:2, reduction = "pca")" >> src/analyze/${lib_name}_analyze.r
  echo "dev.off()" >> src/analyze/${lib_name}_analyze.r
  echo "" >> src/analyze/${lib_name}_analyze.r

  echo "png(file="/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/${lib_name}_elbow.png")" >> src/analyze/${lib_name}_analyze.r
  echo "ElbowPlot(${lib_name})" >> src/analyze/${lib_name}_analyze.r
  echo "dev.off()" >> src/analyze/${lib_name}_analyze.r
  echo "" >> src/analyze/${lib_name}_analyze.r

  echo "${lib_name} <- RunUMAP(${lib_name}, dims = 1:10)" >> src/analyze/${lib_name}_analyze.r
  echo "png(file="/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/${lib_name}_dimplot.png")" >> src/analyze/${lib_name}_analyze.r
  echo "DimPlot({lib_name}, reduction = "umap")" >> src/analyze/${lib_name}_analyze.r
  echo "dev.off()" >> src/analyze/${lib_name}_analyze.r
  echo "" >> src/analyze/${lib_name}_analyze.r

  echo "saveRDS(${lib_name}, file = "/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/${lib_name}.rds")" >> src/analyze/${lib_name}_analyze.r

  echo "nohup R CMD BATCH $BSCRATCH/at.sc.db/src/analyze/${lib_name}_analyze.r $BSCRATCH/at.sc.db/log/analyze_${lib_name}.Rout &" >> src/analyze/launch_analyze_scripts.sh

  i=$(($i + 1))
 
done
