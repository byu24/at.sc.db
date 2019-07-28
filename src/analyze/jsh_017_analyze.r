#!/usr/bin/env

setwd('/global/projectb/scratch/byu24/at.sc.db/scratch')
library('Seurat')
library(ggplot2)

jsh_017.data<-Read10X(data.dir = '/global/projectb/scratch/byu24/at.sc.db/scratch/jsh_017/jsh_017_star.Solo.out')
jsh_017<- CreateSeuratObject(counts = jsh_017.data, project = 'jsh_017')
jsh_017 <- NormalizeData(jsh_017, normalization.method = 'LogNormalize', scale.factor = 10000)

jsh_017 <- FindVariableFeatures(jsh_017, selection.method = 'vst', nfeatures = 2000)
jsh_017 <- ScaleData(jsh_017)

jsh_017 <- RunPCA(jsh_017, features = VariableFeatures(object = jsh_017))
print(jsh_017[['pca']], dims = 1:5, nfeatures = 5)

png(file='/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/jsh_017_vizdim.png')
VizDimLoadings(jsh_017, dims = 1:2, reduction = 'pca')
dev.off()

png(file='/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/jsh_017_elbow.png')
ElbowPlot(jsh_017)
dev.off()

jsh_017 <- RunUMAP(jsh_017, dims = 1:10)
png(file='/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/jsh_017_dimplot.png')
DimPlot(jsh_017, reduction = 'umap')
dev.off()

saveRDS(jsh_017, file = '/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/jsh_017.rds')
