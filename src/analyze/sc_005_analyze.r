#!/usr/bin/env

setwd('/global/projectb/scratch/byu24/at.sc.db/scratch')
library('Seurat')
library(ggplot2)

sc_005.data<-Read10X(data.dir = '/global/projectb/scratch/byu24/at.sc.db/scratch/sc_005/sc_005_star.Solo.out')
sc_005<- CreateSeuratObject(counts = sc_005.data, project = 'sc_005')
sc_005 <- NormalizeData(sc_005, normalization.method = 'LogNormalize', scale.factor = 10000)

sc_005 <- FindVariableFeatures(sc_005, selection.method = 'vst', nfeatures = 2000)
sc_005 <- ScaleData(sc_005)

sc_005 <- RunPCA(sc_005, features = VariableFeatures(object = sc_005))
print(sc_005[['pca']], dims = 1:5, nfeatures = 5)

png(file='/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/sc_005_vizdim.png')
VizDimLoadings(sc_005, dims = 1:2, reduction = 'pca')
dev.off()

png(file='/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/sc_005_elbow.png')
ElbowPlot(sc_005)
dev.off()

sc_005 <- RunUMAP(sc_005, dims = 1:10)
png(file='/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/sc_005_dimplot.png')
DimPlot(sc_005, reduction = 'umap')
dev.off()

saveRDS(sc_005, file = '/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/sc_005.rds')
