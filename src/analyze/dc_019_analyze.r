#!/usr/bin/env

setwd('/global/projectb/scratch/byu24/at.sc.db/scratch')
library('Seurat')
library(ggplot2)

dc_019.data<-Read10X(data.dir = '/global/projectb/scratch/byu24/at.sc.db/scratch/dc_019/dc_019_star.Solo.out')
dc_019<- CreateSeuratObject(counts = dc_019.data, project = 'dc_019')
dc_019 <- NormalizeData(dc_019, normalization.method = 'LogNormalize', scale.factor = 10000)

dc_019 <- FindVariableFeatures(dc_019, selection.method = 'vst', nfeatures = 2000)
dc_019 <- ScaleData(dc_019)

dc_019 <- RunPCA(dc_019, features = VariableFeatures(object = dc_019))
print(dc_019[['pca']], dims = 1:5, nfeatures = 5)

png(file='/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/dc_019_vizdim.png')
VizDimLoadings(dc_019, dims = 1:2, reduction = 'pca')
dev.off()

png(file='/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/dc_019_elbow.png')
ElbowPlot(dc_019)
dev.off()

dc_019 <- RunUMAP(dc_019, dims = 1:10)
png(file='/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/dc_019_dimplot.png')
DimPlot(dc_019, reduction = 'umap')
dev.off()

saveRDS(dc_019, file = '/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/dc_019.rds')
