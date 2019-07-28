#!/usr/bin/env

setwd('/global/projectb/scratch/byu24/at.sc.db/scratch')
library('Seurat')
library(ggplot2)

sc_004.data<-Read10X(data.dir = '/global/projectb/scratch/byu24/at.sc.db/scratch/sc_004/sc_004_star.Solo.out')
sc_004<- CreateSeuratObject(counts = sc_004.data, project = 'sc_004')
sc_004 <- NormalizeData(sc_004, normalization.method = 'LogNormalize', scale.factor = 10000)

sc_004 <- FindVariableFeatures(sc_004, selection.method = 'vst', nfeatures = 2000)
sc_004 <- ScaleData(sc_004)

sc_004 <- RunPCA(sc_004, features = VariableFeatures(object = sc_004))
print(sc_004[['pca']], dims = 1:5, nfeatures = 5)

png(file='/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/sc_004_vizdim.png')
VizDimLoadings(sc_004, dims = 1:2, reduction = 'pca')
dev.off()

png(file='/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/sc_004_elbow.png')
ElbowPlot(sc_004)
dev.off()

sc_004 <- RunUMAP(sc_004, dims = 1:10)
png(file='/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/sc_004_dimplot.png')
DimPlot(sc_004, reduction = 'umap')
dev.off()

saveRDS(sc_004, file = '/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/sc_004.rds')
