#!/usr/bin/env

setwd('/global/projectb/scratch/byu24/at.sc.db/scratch')
library('Seurat')
library(ggplot2)

rc_013.data<-Read10X(data.dir = '/global/projectb/scratch/byu24/at.sc.db/scratch/rc_013/rc_013_star.Solo.out')
rc_013<- CreateSeuratObject(counts = rc_013.data, project = 'rc_013')
rc_013 <- NormalizeData(rc_013, normalization.method = 'LogNormalize', scale.factor = 10000)

rc_013 <- FindVariableFeatures(rc_013, selection.method = 'vst', nfeatures = 2000)
rc_013 <- ScaleData(rc_013)

rc_013 <- RunPCA(rc_013, features = VariableFeatures(object = rc_013))
print(rc_013[['pca']], dims = 1:5, nfeatures = 5)

png(file='/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/rc_013_vizdim.png')
VizDimLoadings(rc_013, dims = 1:2, reduction = 'pca')
dev.off()

png(file='/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/rc_013_elbow.png')
ElbowPlot(rc_013)
dev.off()

rc_013 <- RunUMAP(rc_013, dims = 1:10)
png(file='/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/rc_013_dimplot.png')
DimPlot(rc_013, reduction = 'umap')
dev.off()

saveRDS(rc_013, file = '/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/rc_013.rds')
