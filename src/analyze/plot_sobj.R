
library("Seurat")
library(ggplot2)
library("tidyverse")
library(future)
options(future.globals.maxSize = 100000 * 1024^2)

#Load all R objects
at.integrated = readRDS(file = "/global/projectb/scratch/byu24/at.sc.db/scratch/robjects/at_integrated.rds")

at.integrated
head(at.integrated@metadata)

png(file="/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/at_integrated.png")
DimPlot(at.integrated, reduction = "umap")
dev.off()