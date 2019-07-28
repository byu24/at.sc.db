#!/usr/bin/env

setwd('/global/projectb/scratch/byu24/at.sc.db/scratch')
library('Seurat')
library(ggplot2)

sc_009.data<-Read10X(data.dir = '/global/projectb/scratch/byu24/at.sc.db/scratch/sc_009/sc_009_star.Solo.out')
sc_009<- CreateSeuratObject(counts = sc_009.data, project = 'sc_009')
sc_009 <- NormalizeData(sc_009, normalization.method = 'LogNormalize', scale.factor = 10000)

sc_009 <- FindVariableFeatures(sc_009, selection.method = 'vst', nfeatures = 2000)
sc_009 <- ScaleData(sc_009)

sc_009 <- RunPCA(sc_009, features = VariableFeatures(object = sc_009))
print(sc_009[['pca']], dims = 1:5, nfeatures = 5)

png(file='/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/sc_009_vizdim.png')
VizDimLoadings(sc_009, dims = 1:2, reduction = 'pca')
dev.off()

png(file='/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/sc_009_elbow.png')
ElbowPlot(sc_009)
dev.off()

sc_009 <- RunUMAP(sc_009, dims = 1:10)
png(file='/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/sc_009_dimplot.png')
DimPlot(sc_009, reduction = 'umap')
dev.off()

saveRDS(sc_009, file = '/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/sc_009.rds')
