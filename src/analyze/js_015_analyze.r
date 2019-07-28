#!/usr/bin/env

setwd('/global/projectb/scratch/byu24/at.sc.db/scratch')
library('Seurat')
library(ggplot2)

js_015.data<-Read10X(data.dir = '/global/projectb/scratch/byu24/at.sc.db/scratch/js_015/js_015_star.Solo.out')
js_015<- CreateSeuratObject(counts = js_015.data, project = 'js_015')
js_015 <- NormalizeData(js_015, normalization.method = 'LogNormalize', scale.factor = 10000)

js_015 <- FindVariableFeatures(js_015, selection.method = 'vst', nfeatures = 2000)
js_015 <- ScaleData(js_015)

js_015 <- RunPCA(js_015, features = VariableFeatures(object = js_015))
print(js_015[['pca']], dims = 1:5, nfeatures = 5)

png(file='/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/js_015_vizdim.png')
VizDimLoadings(js_015, dims = 1:2, reduction = 'pca')
dev.off()

png(file='/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/js_015_elbow.png')
ElbowPlot(js_015)
dev.off()

js_015 <- RunUMAP(js_015, dims = 1:10)
png(file='/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/js_015_dimplot.png')
DimPlot(js_015, reduction = 'umap')
dev.off()

saveRDS(js_015, file = '/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/js_015.rds')
