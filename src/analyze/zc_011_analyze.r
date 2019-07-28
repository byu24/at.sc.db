#!/usr/bin/env

setwd('/global/projectb/scratch/byu24/at.sc.db/scratch')
library('Seurat')
library(ggplot2)

zc_011.data<-Read10X(data.dir = '/global/projectb/scratch/byu24/at.sc.db/scratch/zc_011/zc_011_star.Solo.out')
zc_011<- CreateSeuratObject(counts = zc_011.data, project = 'zc_011')
zc_011 <- NormalizeData(zc_011, normalization.method = 'LogNormalize', scale.factor = 10000)

zc_011 <- FindVariableFeatures(zc_011, selection.method = 'vst', nfeatures = 2000)
zc_011 <- ScaleData(zc_011)

zc_011 <- RunPCA(zc_011, features = VariableFeatures(object = zc_011))
print(zc_011[['pca']], dims = 1:5, nfeatures = 5)

png(file='/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/zc_011_vizdim.png')
VizDimLoadings(zc_011, dims = 1:2, reduction = 'pca')
dev.off()

png(file='/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/zc_011_elbow.png')
ElbowPlot(zc_011)
dev.off()

zc_011 <- RunUMAP(zc_011, dims = 1:10)
png(file='/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/zc_011_dimplot.png')
DimPlot(zc_011, reduction = 'umap')
dev.off()

saveRDS(zc_011, file = '/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/zc_011.rds')
