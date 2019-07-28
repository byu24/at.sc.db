#!/usr/bin/env

setwd('/global/projectb/scratch/byu24/at.sc.db/scratch')
library('Seurat')
library(ggplot2)

sc_008.data<-Read10X(data.dir = '/global/projectb/scratch/byu24/at.sc.db/scratch/sc_008/sc_008_star.Solo.out')
sc_008<- CreateSeuratObject(counts = sc_008.data, project = 'sc_008')
sc_008 <- NormalizeData(sc_008, normalization.method = 'LogNormalize', scale.factor = 10000)

sc_008 <- FindVariableFeatures(sc_008, selection.method = 'vst', nfeatures = 2000)
sc_008 <- ScaleData(sc_008)

sc_008 <- RunPCA(sc_008, features = VariableFeatures(object = sc_008))
print(sc_008[['pca']], dims = 1:5, nfeatures = 5)

png(file='/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/sc_008_vizdim.png')
VizDimLoadings(sc_008, dims = 1:2, reduction = 'pca')
dev.off()

png(file='/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/sc_008_elbow.png')
ElbowPlot(sc_008)
dev.off()

sc_008 <- RunUMAP(sc_008, dims = 1:10)
png(file='/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/sc_008_dimplot.png')
DimPlot(sc_008, reduction = 'umap')
dev.off()

saveRDS(sc_008, file = '/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/sc_008.rds')
