#!/usr/bin/env

setwd('/global/projectb/scratch/byu24/at.sc.db/scratch')
library('Seurat')
library(ggplot2)

sc_010.data<-Read10X(data.dir = '/global/projectb/scratch/byu24/at.sc.db/scratch/sc_010/sc_010_star.Solo.out')
sc_010<- CreateSeuratObject(counts = sc_010.data, project = 'sc_010')
sc_010 <- NormalizeData(sc_010, normalization.method = 'LogNormalize', scale.factor = 10000)

sc_010 <- FindVariableFeatures(sc_010, selection.method = 'vst', nfeatures = 2000)
sc_010 <- ScaleData(sc_010)

sc_010 <- RunPCA(sc_010, features = VariableFeatures(object = sc_010))
print(sc_010[['pca']], dims = 1:5, nfeatures = 5)

png(file='/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/sc_010_vizdim.png')
VizDimLoadings(sc_010, dims = 1:2, reduction = 'pca')
dev.off()

png(file='/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/sc_010_elbow.png')
ElbowPlot(sc_010)
dev.off()

sc_010 <- RunUMAP(sc_010, dims = 1:10)
png(file='/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/sc_010_dimplot.png')
DimPlot(sc_010, reduction = 'umap')
dev.off()

saveRDS(sc_010, file = '/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/sc_010.rds')
