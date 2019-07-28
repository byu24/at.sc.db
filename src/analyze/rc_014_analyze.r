#!/usr/bin/env

setwd('/global/projectb/scratch/byu24/at.sc.db/scratch')
library('Seurat')
library(ggplot2)

rc_014.data<-Read10X(data.dir = '/global/projectb/scratch/byu24/at.sc.db/scratch/rc_014/rc_014_star.Solo.out')
rc_014<- CreateSeuratObject(counts = rc_014.data, project = 'rc_014')
rc_014 <- NormalizeData(rc_014, normalization.method = 'LogNormalize', scale.factor = 10000)

rc_014 <- FindVariableFeatures(rc_014, selection.method = 'vst', nfeatures = 2000)
rc_014 <- ScaleData(rc_014)

rc_014 <- RunPCA(rc_014, features = VariableFeatures(object = rc_014))
print(rc_014[['pca']], dims = 1:5, nfeatures = 5)

png(file='/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/rc_014_vizdim.png')
VizDimLoadings(rc_014, dims = 1:2, reduction = 'pca')
dev.off()

png(file='/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/rc_014_elbow.png')
ElbowPlot(rc_014)
dev.off()

rc_014 <- RunUMAP(rc_014, dims = 1:10)
png(file='/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/rc_014_dimplot.png')
DimPlot(rc_014, reduction = 'umap')
dev.off()

saveRDS(rc_014, file = '/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/rc_014.rds')
