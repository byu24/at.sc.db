#!/usr/bin/env Rscript

#Load necessary packages

install.packages('Seurat')
library("Seurat")
library(tidyverse)
library(ggplot2)
library(cowplot)

setwd("/global/projectb/scratch/byu24/at.sc.db/scratch")

dc_018.data<-Read10X(data.dir = paste0("/global/projectb/scratch/byu24/at.sc.db/scratch", dataset,"/",dataset,"_star.Solo.out"))

dc_018 <- CreateSeuratObject(counts = dc_018.data, project = "dc_018")
dc_018 <- NormalizeData(dc_018, normalization.method = "LogNormalize", scale.factor = 10000)

dc_018 <- FindVariableFeatures(dc_018, selection.method = "vst", nfeatures = 2000)
dc_018 <- ScaleData(dc_018)

dc_018 <- RunPCA(dc_018, features = VariableFeatures(object = dc_018))
print(dc_018[["pca"]], dims = 1:5, nfeatures = 5)

png(file="../analysis/dc_018_vizdim.png")
VizDimLoadings(dc_018, dims = 1:2, reduction = "pca")
dev.off()

png(file="/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/dc_018_elbow.png")
ElbowPlot(dc_018)
dev.off()

#Long step
dc_018 <- RunUMAP(dc_018, dims = 1:10)
png(file="/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/dc_018_dimplot.png")
DimPlot(dc_018, reduction = "umap")

saveRDS(dc_018, file = "/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/dc_018.rds")

#at.big<-merge(dc_018, y=c(dc_019,sc_004),add.cell.ids=c("dc018", "dc019","sc004"),project = "at3")
#unique(sapply(X = strsplit(colnames(at.big), split = "_"), FUN = "[", 1))
