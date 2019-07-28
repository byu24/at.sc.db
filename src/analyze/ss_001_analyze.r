#!/usr/bin/env

setwd('/global/projectb/scratch/byu24/at.sc.db/scratch')
library('Seurat')
library(ggplot2)

ss_001.data<-Read10X(data.dir = '/global/projectb/scratch/byu24/at.sc.db/scratch/ss_001/ss_001_star.Solo.out')
ss_001<- CreateSeuratObject(counts = ss_001.data, project = 'ss_001')
ss_001 <- NormalizeData(ss_001, normalization.method = 'LogNormalize', scale.factor = 10000)

ss_001 <- FindVariableFeatures(ss_001, selection.method = 'vst', nfeatures = 2000)
ss_001 <- ScaleData(ss_001)

ss_001 <- RunPCA(ss_001, features = VariableFeatures(object = ss_001))
print(ss_001[['pca']], dims = 1:5, nfeatures = 5)

png(file='/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/ss_001_vizdim.png')
VizDimLoadings(ss_001, dims = 1:2, reduction = 'pca')
dev.off()

png(file='/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/ss_001_elbow.png')
ElbowPlot(ss_001)
dev.off()

ss_001 <- RunUMAP(ss_001, dims = 1:10)
png(file='/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/ss_001_dimplot.png')
DimPlot(ss_001, reduction = 'umap')
dev.off()

saveRDS(ss_001, file = '/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/ss_001.rds')
