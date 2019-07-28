#!/usr/bin/env

setwd('/global/projectb/scratch/byu24/at.sc.db/scratch')
library('Seurat')
library(ggplot2)

js_016.data<-Read10X(data.dir = '/global/projectb/scratch/byu24/at.sc.db/scratch/js_016/js_016_star.Solo.out')
js_016<- CreateSeuratObject(counts = js_016.data, project = 'js_016')
js_016 <- NormalizeData(js_016, normalization.method = 'LogNormalize', scale.factor = 10000)

js_016 <- FindVariableFeatures(js_016, selection.method = 'vst', nfeatures = 2000)
js_016 <- ScaleData(js_016)

js_016 <- RunPCA(js_016, features = VariableFeatures(object = js_016))
print(js_016[['pca']], dims = 1:5, nfeatures = 5)

png(file='/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/js_016_vizdim.png')
VizDimLoadings(js_016, dims = 1:2, reduction = 'pca')
dev.off()

png(file='/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/js_016_elbow.png')
ElbowPlot(js_016)
dev.off()

js_016 <- RunUMAP(js_016, dims = 1:10)
png(file='/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/js_016_dimplot.png')
DimPlot(js_016, reduction = 'umap')
dev.off()

saveRDS(js_016, file = '/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/js_016.rds')
