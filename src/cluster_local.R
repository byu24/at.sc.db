#!/usr/bin/env Rscript

#Load necessary packages
install.packages('reticulate')

library('reticulate')
library("Seurat")
library(ggplot2)
library(tidyverse)
library(umap)
use_virtualenv("scRNAseq")

setwd("C:/Users/BYU24/Desktop/JGI 2019/at.sc.db/scratch/")

sample_metadata = read.csv("C:/Users/BYU24/Desktop/JGI 2019/at.sc.db/data/sample_metadata.csv")
dataset = sample_metadata$Name
protoplast_loci = read_csv("C:/Users/BYU24/Desktop/JGI 2019/at.sc.db/data/protoplast_loci.csv")%>%
  pull(Locus)

rm(list = ls(pattern = "sc_005"))
# Functions ------------------------------------------------------------------
#Load dataset
get_dge <- function(dataset) {
  dge <- Seurat::Read10X(data.dir = paste0("C:/Users/BYU24/Desktop/JGI 2019/at.sc.db/scratch/",dataset,"/",dataset,"_star.Solo.out")) %>%
    `colnames<-`(c(paste0(colnames(.), "_", dataset)))
}

#Pull cell statistics
get_dge_stats <- function(dge) {
  nAt <- Matrix::colSums(dge[rownames(dge)[str_detect(rownames(dge), "AT.G.....")],]) %>%
    enframe("Cell","nAt")
  nGene <- Matrix::colSums(dge[rownames(dge)[str_detect(rownames(dge), "AT.G.....")],] > 0) %>% 
    enframe("Cell", "nGene")
  nMa <- Matrix::colSums(dge[rownames(dge)[str_detect(rownames(dge), "AT.G.....")],]) %>%
    enframe("Cell", "nMa")
  nCM <- Matrix::colSums(dge[rownames(dge)[str_detect(rownames(dge), "AT[MC]G.....")],]) %>%
    enframe("Cell", "nCM")
  full_join(nAt, nMa) %>%
    full_join(nCM) %>%
    full_join(nGene) %>% 
    mutate(
      pAt = nAt/(nMa + nAt),
      pCM = nCM/nAt) %>%
    column_to_rownames("Cell")
}

#Filter dataset
filter_dge <- function(dge, dge_stats, thresh_Gene = 200) {
  dge_stats <- dge_stats %>% 
    as_tibble(rownames = "Cell") %>% 
    filter(nGene >= thresh_Gene)
  
  nn_percentile <- dge_stats %>%
    summarize(mincell = round(n()*0.01)) %>%
    pull(mincell)
  
  threshold <- dge_stats %>% 
    top_n(nn_percentile, nAt) %>%
    summarize(threshold = min(nAt)*0.05) %>%
    pull(threshold)
  
  keep_cells <- dge_stats %>% 
    filter(nAt >= threshold) %>%
    pull(Cell)
  
  dge <- dge[,keep_cells]
  
  protoplast_loci = read_csv("C:/Users/BYU24/Desktop/JGI 2019/at.sc.db/data/protoplast_loci.csv")%>%
    pull(Locus)
  
  keep_genes <- setdiff(rownames(dge)[str_detect(rownames(dge), "AT.G.....")], protoplast_loci)
  dge <- dge[keep_genes,]
}

#Seurat Object
get_sobj <- function(dge_filtered, dge_stats, group) {
  sobj <- Seurat::CreateSeuratObject(
    counts = dge_filtered,
    project = group,
    meta.data = dge_stats) 
}

#Transform data
trans_dge <- function(sobj) {
  sobj <- sobj %>%
    SCTransform(vars.to.regress = "nCM", variable.features.n = 2000) %>%
    RunPCA(npcs=50) %>%
    RunUMAP(dims=1:30)
  }


# SS_001 ------------------------------------------------------------------
ss_001.read = get_dge("ss_001")
ss_001.data = get_dge_stats(ss_001.read)
ss_001.filtered = filter_dge(ss_001.read,ss_001.data)
ss_001 = get_sobj(ss_001.filtered, ss_001.data, group = "ss_001")
ss_001 = trans_dge(ss_001)

png(file="C:/Users/BYU24/Desktop/JGI 2019/at.sc.db/scratch/analysis/ss_001_umap.png")
DimPlot(ss_001, reduction = "umap")
dev.off()

png(file="C:/Users/BYU24/Desktop/JGI 2019/at.sc.db/scratch/analysis/ss_001_elbow.png")
ElbowPlot(ss_001)
dev.off()

saveRDS(ss_001, file = "C:/Users/BYU24/Desktop/JGI 2019/at.sc.db/scratch/robjects/ss_001.rds")

FindNeighbors(ss_001, dims = 1:30)
FindClusters(ss_001)
FindAllMarkers(ss_001)

#Find cluster biomarkers
#cluster1.markers <- FindMarkers(dc_018, ident.1 = 1, min.pct = 0.25)
#head(cluster1.markers, n = 5)

# SS_002 ------------------------------------------------------------------
ss_002.read = get_dge("ss_002")
ss_002.data = get_dge_stats(ss_002.read)
ss_002.filtered = filter_dge(ss_002.read,ss_002.data)
ss_002 = get_sobj(ss_002.filtered, ss_002.data, group = "ss_002")
ss_002 = trans_dge(ss_002)

png(file="C:/Users/BYU24/Desktop/JGI 2019/at.sc.db/scratch/analysis/ss_002_umap.png")
DimPlot(ss_002, reduction = "umap")
dev.off()

png(file="C:/Users/BYU24/Desktop/JGI 2019/at.sc.db/scratch/analysis/ss_002_elbow.png")
ElbowPlot(ss_002)
dev.off()

saveRDS(ss_002, file = "C:/Users/BYU24/Desktop/JGI 2019/at.sc.db/scratch/robjects/ss_002.rds")

# SS_003 ------------------------------------------------------------------
ss_003.read = get_dge("ss_003")
ss_003.data = get_dge_stats(ss_003.read)
ss_003.filtered = filter_dge(ss_003.read,ss_003.data)
ss_003 = get_sobj(ss_003.filtered, ss_003.data, group = "ss_003")
ss_003 = trans_dge(ss_003)

png(file="C:/Users/BYU24/Desktop/JGI 2019/at.sc.db/scratch/analysis/ss_003_umap.png")
DimPlot(ss_003, reduction = "umap")
dev.off()

png(file="C:/Users/BYU24/Desktop/JGI 2019/at.sc.db/scratch/analysis/ss_003_elbow.png")
ElbowPlot(ss_003)
dev.off()

saveRDS(ss_003, file = "C:/Users/BYU24/Desktop/JGI 2019/at.sc.db/scratch/robjects/ss_003.rds")

# SC_004 ------------------------------------------------------------------
sc_004.read = get_dge("sc_004")
sc_004.data = get_dge_stats(sc_004.read)
sc_004.filtered = filter_dge(sc_004.read,sc_004.data)
sc_004 = get_sobj(sc_004.filtered, sc_004.data, group = "sc_004")
sc_004 = trans_dge(sc_004)

png(file="C:/Users/BYU24/Desktop/JGI 2019/at.sc.db/scratch/analysis/sc_004_umap.png")
DimPlot(sc_004, reduction = "umap")
dev.off()

png(file="C:/Users/BYU24/Desktop/JGI 2019/at.sc.db/scratch/analysis/sc_004_elbow.png")
ElbowPlot(sc_004)
dev.off()

saveRDS(sc_004, file = "C:/Users/BYU24/Desktop/JGI 2019/at.sc.db/scratch/robjects/sc_004.rds")

# sc_005 ------------------------------------------------------------------
sc_005.read = get_dge("sc_005")
sc_005.data = get_dge_stats(sc_005.read)
sc_005.filtered = filter_dge(sc_005.read,sc_005.data)
sc_005 = get_sobj(sc_005.filtered, sc_005.data, group = "sc_005")
sc_005 = trans_dge(sc_005)

png(file="C:/Users/BYU24/Desktop/JGI 2019/at.sc.db/scratch/analysis/sc_005_umap.png")
DimPlot(sc_005, reduction = "umap")
dev.off()

png(file="C:/Users/BYU24/Desktop/JGI 2019/at.sc.db/scratch/analysis/sc_005_elbow.png")
ElbowPlot(sc_005)
dev.off()

saveRDS(sc_005, file = "C:/Users/BYU24/Desktop/JGI 2019/at.sc.db/scratch/robjects/sc_005.rds")

# SC_006 ------------------------------------------------------------------
sc_006.read = get_dge("sc_006")
sc_006.data = get_dge_stats(sc_006.read)
sc_006.filtered = filter_dge(sc_006.read,sc_006.data)
sc_006 = get_sobj(sc_006.filtered, sc_006.data, group = "sc_006")
sc_006 = trans_dge(sc_006)

png(file="C:/Users/BYU24/Desktop/JGI 2019/at.sc.db/scratch/analysis/sc_006_umap.png")
DimPlot(sc_006, reduction = "umap")
dev.off()

png(file="C:/Users/BYU24/Desktop/JGI 2019/at.sc.db/scratch/analysis/sc_006_elbow.png")
ElbowPlot(sc_006)
dev.off()

saveRDS(sc_006, file = "C:/Users/BYU24/Desktop/JGI 2019/at.sc.db/scratch/robjects/sc_006.rds")

# SC_007 ------------------------------------------------------------------
sc_007.read = get_dge("sc_007")
sc_007.data = get_dge_stats(sc_007.read)
sc_007.filtered = filter_dge(sc_007.read,sc_007.data)
sc_007 = get_sobj(sc_007.filtered, sc_007.data, group = "sc_007")
sc_007 = trans_dge(sc_007)

png(file="C:/Users/BYU24/Desktop/JGI 2019/at.sc.db/scratch/analysis/sc_007_umap.png")
DimPlot(sc_007, reduction = "umap")
dev.off()

png(file="C:/Users/BYU24/Desktop/JGI 2019/at.sc.db/scratch/analysis/sc_007_elbow.png")
ElbowPlot(sc_007)
dev.off()

saveRDS(sc_007, file = "C:/Users/BYU24/Desktop/JGI 2019/at.sc.db/scratch/robjects/sc_007.rds")

# SC_008 ------------------------------------------------------------------
sc_008.read = get_dge("sc_008")
sc_008.data = get_dge_stats(sc_008.read)
sc_008.filtered = filter_dge(sc_008.read,sc_008.data)
sc_008 = get_sobj(sc_008.filtered, sc_008.data, group = "sc_008")
sc_008 = trans_dge(sc_008)

png(file="C:/Users/BYU24/Desktop/JGI 2019/at.sc.db/scratch/analysis/sc_008_umap.png")
DimPlot(sc_008, reduction = "umap")
dev.off()

png(file="C:/Users/BYU24/Desktop/JGI 2019/at.sc.db/scratch/analysis/sc_008_elbow.png")
ElbowPlot(sc_008)
dev.off()

saveRDS(sc_008, file = "C:/Users/BYU24/Desktop/JGI 2019/at.sc.db/scratch/robjects/sc_008.rds")

# SC_009 ------------------------------------------------------------------
sc_009.read = get_dge("sc_009")
sc_009.data = get_dge_stats(sc_009.read)
sc_009.filtered = filter_dge(sc_009.read,sc_009.data)
sc_009 = get_sobj(sc_009.filtered, sc_009.data, group = "sc_009")
sc_009 = trans_dge(sc_009)

png(file="C:/Users/BYU24/Desktop/JGI 2019/at.sc.db/scratch/analysis/sc_009_umap.png")
DimPlot(sc_009, reduction = "umap")
dev.off()

png(file="C:/Users/BYU24/Desktop/JGI 2019/at.sc.db/scratch/analysis/sc_009_elbow.png")
ElbowPlot(sc_009)
dev.off()

saveRDS(sc_009, file = "C:/Users/BYU24/Desktop/JGI 2019/at.sc.db/scratch/robjects/sc_009.rds")

# SC_010 ------------------------------------------------------------------
sc_010.read = get_dge("sc_010")
sc_010.data = get_dge_stats(sc_010.read)
sc_010.filtered = filter_dge(sc_010.read,sc_010.data)
sc_010 = get_sobj(sc_010.filtered, sc_010.data, group = "sc_010")
sc_010 = trans_dge(sc_010)

png(file="C:/Users/BYU24/Desktop/JGI 2019/at.sc.db/scratch/analysis/sc_010_umap.png")
DimPlot(sc_010, reduction = "umap")
dev.off()

png(file="C:/Users/BYU24/Desktop/JGI 2019/at.sc.db/scratch/analysis/sc_010_elbow.png")
ElbowPlot(sc_010)
dev.off()

saveRDS(sc_010, file = "C:/Users/BYU24/Desktop/JGI 2019/at.sc.db/scratch/robjects/sc_010.rds")

# ZC_011 ------------------------------------------------------------------
zc_011.read = get_dge("zc_011")
zc_011.data = get_dge_stats(zc_011.read)
zc_011.filtered = filter_dge(zc_011.read,zc_011.data)
zc_011 = get_sobj(zc_011.filtered, zc_011.data, group = "zc_011")
zc_011 = trans_dge(zc_011)

png(file="C:/Users/BYU24/Desktop/JGI 2019/at.sc.db/scratch/analysis/zc_011_umap.png")
DimPlot(zc_011, reduction = "umap")
dev.off()

png(file="C:/Users/BYU24/Desktop/JGI 2019/at.sc.db/scratch/analysis/zc_011_elbow.png")
ElbowPlot(zc_011)
dev.off()

saveRDS(zc_011, file = "C:/Users/BYU24/Desktop/JGI 2019/at.sc.db/scratch/robjects/zc_011.rds")

# RC_012 ------------------------------------------------------------------
rc_012.read = get_dge("rc_012")
rc_012.data = get_dge_stats(rc_012.read)
rc_012.filtered = filter_dge(rc_012.read,rc_012.data)
rc_012 = get_sobj(rc_012.filtered, rc_012.data, group = "rc_012")
rc_012 = trans_dge(rc_012)

png(file="C:/Users/BYU24/Desktop/JGI 2019/at.sc.db/scratch/analysis/rc_012_umap.png")
DimPlot(rc_012, reduction = "umap")
dev.off()

png(file="C:/Users/BYU24/Desktop/JGI 2019/at.sc.db/scratch/analysis/rc_012_elbow.png")
ElbowPlot(rc_012)
dev.off()

saveRDS(rc_012, file = "C:/Users/BYU24/Desktop/JGI 2019/at.sc.db/scratch/robjects/rc_012.rds")

# RC_013 ------------------------------------------------------------------
rc_013.read = get_dge("rc_013")
rc_013.data = get_dge_stats(rc_013.read)
rc_013.filtered = filter_dge(rc_013.read,rc_013.data)
rc_013 = get_sobj(rc_013.filtered, rc_013.data, group = "rc_013")
rc_013 = trans_dge(rc_013)

png(file="C:/Users/BYU24/Desktop/JGI 2019/at.sc.db/scratch/analysis/rc_013_umap.png")
DimPlot(rc_013, reduction = "umap")
dev.off()

png(file="C:/Users/BYU24/Desktop/JGI 2019/at.sc.db/scratch/analysis/rc_013_elbow.png")
ElbowPlot(rc_013)
dev.off()

saveRDS(rc_013, file = "C:/Users/BYU24/Desktop/JGI 2019/at.sc.db/scratch/robjects/rc_013.rds")

# RC_014 ------------------------------------------------------------------
rc_014.read = get_dge("rc_014")
rc_014.data = get_dge_stats(rc_014.read)
rc_014.filtered = filter_dge(rc_014.read,rc_014.data)
rc_014 = get_sobj(rc_014.filtered, rc_014.data, group = "rc_014")
rc_014 = trans_dge(rc_014)

png(file="C:/Users/BYU24/Desktop/JGI 2019/at.sc.db/scratch/analysis/rc_014_umap.png")
DimPlot(rc_014, reduction = "umap")
dev.off()

png(file="C:/Users/BYU24/Desktop/JGI 2019/at.sc.db/scratch/analysis/rc_014_elbow.png")
ElbowPlot(rc_014)
dev.off()

saveRDS(rc_014, file = "C:/Users/BYU24/Desktop/JGI 2019/at.sc.db/scratch/robjects/rc_014.rds")

# JS_015 ------------------------------------------------------------------
js_015.read = get_dge("js_015")
js_015.data = get_dge_stats(js_015.read)
js_015.filtered = filter_dge(js_015.read,js_015.data)
js_015 = get_sobj(js_015.filtered, js_015.data, group = "js_015")
js_015 = trans_dge(js_015)

png(file="C:/Users/BYU24/Desktop/JGI 2019/at.sc.db/scratch/analysis/js_015_umap.png")
DimPlot(js_015, reduction = "umap")
dev.off()

png(file="C:/Users/BYU24/Desktop/JGI 2019/at.sc.db/scratch/analysis/js_015_elbow.png")
ElbowPlot(js_015)
dev.off()

saveRDS(js_015, file = "C:/Users/BYU24/Desktop/JGI 2019/at.sc.db/scratch/robjects/js_015.rds")

# JS_016 ------------------------------------------------------------------
js_016.read = get_dge("js_016")
js_016.data = get_dge_stats(js_016.read)
js_016.filtered = filter_dge(js_016.read,js_016.data)
js_016 = get_sobj(js_016.filtered, js_016.data, group = "js_016")
js_016 = trans_dge(js_016)

png(file="C:/Users/BYU24/Desktop/JGI 2019/at.sc.db/scratch/analysis/js_016_umap.png")
DimPlot(js_016, reduction = "umap")
dev.off()

png(file="C:/Users/BYU24/Desktop/JGI 2019/at.sc.db/scratch/analysis/js_016_elbow.png")
ElbowPlot(js_016)
dev.off()

saveRDS(js_016, file = "C:/Users/BYU24/Desktop/JGI 2019/at.sc.db/scratch/robjects/js_016.rds")

# JSH_017 ------------------------------------------------------------------
jsh_017.read = get_dge("jsh_017")
jsh_017.data = get_dge_stats(jsh_017.read)
jsh_017.filtered = filter_dge(jsh_017.read,jsh_017.data)
jsh_017 = get_sobj(jsh_017.filtered, jsh_017.data, group = "jsh_017")
jsh_017 = trans_dge(jsh_017)

png(file="C:/Users/BYU24/Desktop/JGI 2019/at.sc.db/scratch/analysis/jsh_017_umap.png")
DimPlot(jsh_017, reduction = "umap")
dev.off()

png(file="C:/Users/BYU24/Desktop/JGI 2019/at.sc.db/scratch/analysis/jsh_017_elbow.png")
ElbowPlot(jsh_017)
dev.off()

saveRDS(jsh_017, file = "C:/Users/BYU24/Desktop/JGI 2019/at.sc.db/scratch/robjects/jsh_017.rds")

# DC_018 ------------------------------------------------------------------
dc_018.read = get_dge("dc_018")
dc_018.data = get_dge_stats(dc_018.read)
dc_018.filtered = filter_dge(dc_018.read,dc_018.data)
dc_018 = get_sobj(dc_018.filtered, dc_018.data, group = "dc_018")
dc_018 = trans_dge(dc_018)

png(file="C:/Users/BYU24/Desktop/JGI 2019/at.sc.db/scratch/analysis/dc_018_umap.png")
DimPlot(dc_018, reduction = "umap")
dev.off()

png(file="C:/Users/BYU24/Desktop/JGI 2019/at.sc.db/scratch/analysis/dc_018_elbow.png")
ElbowPlot(dc_018)
dev.off()

saveRDS(dc_018, file = "C:/Users/BYU24/Desktop/JGI 2019/at.sc.db/scratch/robjects/dc_018.rds")

# DC_019 ------------------------------------------------------------------
dc_019.read = get_dge("dc_019")
dc_019.data = get_dge_stats(dc_019.read)
dc_019.filtered = filter_dge(dc_019.read,dc_019.data)
dc_019 = get_sobj(dc_019.filtered, dc_019.data, group = "dc_019")
dc_019 = trans_dge(dc_019)

png(file="C:/Users/BYU24/Desktop/JGI 2019/at.sc.db/scratch/analysis/dc_019_umap.png")
DimPlot(dc_019, reduction = "umap")
dev.off()

png(file="C:/Users/BYU24/Desktop/JGI 2019/at.sc.db/scratch/analysis/dc_019_elbow.png")
ElbowPlot(dc_019)
dev.off()

saveRDS(dc_019, file = "C:/Users/BYU24/Desktop/JGI 2019/at.sc.db/scratch/robjects/dc_019.rds")
