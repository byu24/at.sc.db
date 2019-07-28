#!/usr/bin/env

setwd('/global/projectb/scratch/byu24/at.sc.db/scratch')
library('Seurat')
library(ggplot2)

ss_003.data<-Read10X(data.dir = '/global/projectb/scratch/byu24/at.sc.db/scratch/ss_003/ss_003_star.Solo.out')
ss_003<- CreateSeuratObject(counts = ss_003.data, project = 'ss_003')
ss_003 <- NormalizeData(ss_003, normalization.method = 'LogNormalize', scale.factor = 10000)

ss_003 <- FindVariableFeatures(ss_003, selection.method = 'vst', nfeatures = 2000)
ss_003 <- ScaleData(ss_003)

ss_003 <- RunPCA(ss_003, features = VariableFeatures(object = ss_003))
print(ss_003[['pca']], dims = 1:5, nfeatures = 5)

png(file='/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/ss_003_vizdim.png')
VizDimLoadings(ss_003, dims = 1:2, reduction = 'pca')
dev.off()

png(file='/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/ss_003_elbow.png')
ElbowPlot(ss_003)
dev.off()

ss_003 <- RunUMAP(ss_003, dims = 1:10)
png(file='/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/ss_003_dimplot.png')
DimPlot(ss_003, reduction = 'umap')
dev.off()

saveRDS(ss_003, file = '/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/ss_003.rds')
