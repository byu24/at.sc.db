#!/usr/bin/env

setwd('/global/projectb/scratch/byu24/at.sc.db/scratch')
library('Seurat')
library(ggplot2)

dc_018.data<-Read10X(data.dir = '/global/projectb/scratch/byu24/at.sc.db/scratch/dc_018/dc_018_star.Solo.out')
dc_018<- CreateSeuratObject(counts = dc_018.data, project = 'dc_018')
dc_018 <- NormalizeData(dc_018, normalization.method = 'LogNormalize', scale.factor = 10000)

dc_018 <- FindVariableFeatures(dc_018, selection.method = 'vst', nfeatures = 2000)
dc_018 <- ScaleData(dc_018)

dc_018 <- RunPCA(dc_018, features = VariableFeatures(object = dc_018))
print(dc_018[['pca']], dims = 1:5, nfeatures = 5)

png(file='/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/dc_018_vizdim.png')
VizDimLoadings(dc_018, dims = 1:2, reduction = 'pca')
dev.off()

png(file='/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/dc_018_elbow.png')
ElbowPlot(dc_018)
dev.off()

dc_018 <- RunUMAP(dc_018, dims = 1:10)
png(file='/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/dc_018_dimplot.png')
DimPlot(dc_018, reduction = 'umap')
dev.off()

saveRDS(dc_018, file = '/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/dc_018.rds')
