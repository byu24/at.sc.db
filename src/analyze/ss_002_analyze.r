#!/usr/bin/env

setwd('/global/projectb/scratch/byu24/at.sc.db/scratch')
library('Seurat')
library(ggplot2)

ss_002.data<-Read10X(data.dir = '/global/projectb/scratch/byu24/at.sc.db/scratch/ss_002/ss_002_star.Solo.out')
ss_002<- CreateSeuratObject(counts = ss_002.data, project = 'ss_002')
ss_002 <- NormalizeData(ss_002, normalization.method = 'LogNormalize', scale.factor = 10000)

ss_002 <- FindVariableFeatures(ss_002, selection.method = 'vst', nfeatures = 2000)
ss_002 <- ScaleData(ss_002)

ss_002 <- RunPCA(ss_002, features = VariableFeatures(object = ss_002))
print(ss_002[['pca']], dims = 1:5, nfeatures = 5)

png(file='/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/ss_002_vizdim.png')
VizDimLoadings(ss_002, dims = 1:2, reduction = 'pca')
dev.off()

png(file='/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/ss_002_elbow.png')
ElbowPlot(ss_002)
dev.off()

ss_002 <- RunUMAP(ss_002, dims = 1:10)
png(file='/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/ss_002_dimplot.png')
DimPlot(ss_002, reduction = 'umap')
dev.off()

saveRDS(ss_002, file = '/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/ss_002.rds')
