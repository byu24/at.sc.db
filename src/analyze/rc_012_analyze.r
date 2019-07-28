#!/usr/bin/env

setwd('/global/projectb/scratch/byu24/at.sc.db/scratch')
library('Seurat')
library(ggplot2)

rc_012.data<-Read10X(data.dir = '/global/projectb/scratch/byu24/at.sc.db/scratch/rc_012/rc_012_star.Solo.out')
rc_012<- CreateSeuratObject(counts = rc_012.data, project = 'rc_012')
rc_012 <- NormalizeData(rc_012, normalization.method = 'LogNormalize', scale.factor = 10000)

rc_012 <- FindVariableFeatures(rc_012, selection.method = 'vst', nfeatures = 2000)
rc_012 <- ScaleData(rc_012)

rc_012 <- RunPCA(rc_012, features = VariableFeatures(object = rc_012))
print(rc_012[['pca']], dims = 1:5, nfeatures = 5)

png(file='/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/rc_012_vizdim.png')
VizDimLoadings(rc_012, dims = 1:2, reduction = 'pca')
dev.off()

png(file='/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/rc_012_elbow.png')
ElbowPlot(rc_012)
dev.off()

rc_012 <- RunUMAP(rc_012, dims = 1:10)
png(file='/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/rc_012_dimplot.png')
DimPlot(rc_012, reduction = 'umap')
dev.off()

saveRDS(rc_012, file = '/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/rc_012.rds')
