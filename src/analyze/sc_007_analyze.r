#!/usr/bin/env

setwd('/global/projectb/scratch/byu24/at.sc.db/scratch')
library('Seurat')
library(ggplot2)

sc_007.data<-Read10X(data.dir = '/global/projectb/scratch/byu24/at.sc.db/scratch/sc_007/sc_007_star.Solo.out')
sc_007<- CreateSeuratObject(counts = sc_007.data, project = 'sc_007')
sc_007 <- NormalizeData(sc_007, normalization.method = 'LogNormalize', scale.factor = 10000)

sc_007 <- FindVariableFeatures(sc_007, selection.method = 'vst', nfeatures = 2000)
sc_007 <- ScaleData(sc_007)

sc_007 <- RunPCA(sc_007, features = VariableFeatures(object = sc_007))
print(sc_007[['pca']], dims = 1:5, nfeatures = 5)

png(file='/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/sc_007_vizdim.png')
VizDimLoadings(sc_007, dims = 1:2, reduction = 'pca')
dev.off()

png(file='/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/sc_007_elbow.png')
ElbowPlot(sc_007)
dev.off()

sc_007 <- RunUMAP(sc_007, dims = 1:10)
png(file='/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/sc_007_dimplot.png')
DimPlot(sc_007, reduction = 'umap')
dev.off()

saveRDS(sc_007, file = '/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/sc_007.rds')
