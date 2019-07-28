#!/usr/bin/env

setwd('/global/projectb/scratch/byu24/at.sc.db/scratch')
library('Seurat')
library(ggplot2)

sc_006.data<-Read10X(data.dir = '/global/projectb/scratch/byu24/at.sc.db/scratch/sc_006/sc_006_star.Solo.out')
sc_006<- CreateSeuratObject(counts = sc_006.data, project = 'sc_006')
sc_006 <- NormalizeData(sc_006, normalization.method = 'LogNormalize', scale.factor = 10000)

sc_006 <- FindVariableFeatures(sc_006, selection.method = 'vst', nfeatures = 2000)
sc_006 <- ScaleData(sc_006)

sc_006 <- RunPCA(sc_006, features = VariableFeatures(object = sc_006))
print(sc_006[['pca']], dims = 1:5, nfeatures = 5)

png(file='/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/sc_006_vizdim.png')
VizDimLoadings(sc_006, dims = 1:2, reduction = 'pca')
dev.off()

png(file='/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/sc_006_elbow.png')
ElbowPlot(sc_006)
dev.off()

sc_006 <- RunUMAP(sc_006, dims = 1:10)
png(file='/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/sc_006_dimplot.png')
DimPlot(sc_006, reduction = 'umap')
dev.off()

saveRDS(sc_006, file = '/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/sc_006.rds')
