#!/usr/bin/env Rscript

#Load necessary packages

install.packages('reticulate', repos="https://cran.cnr.berkeley.edu/")
install.packages('future', repos="https://cran.cnr.berkeley.edu/")
install.packages('furrr', repos="https://cran.cnr.berkeley.edu/")
install.packages("Seurat",repos="https://cran.cnr.berkeley.edu/")
install.packages("ggplot2",repos="https://cran.cnr.berkeley.edu/")

library("Seurat")
library(ggplot2)
library(tidyverse)

#set working directory
setwd("/global/projectb/scratch/byu24/at.sc.db/scratch")

#replace sample_metadata filepath with actual filepath
sample_metadata = read.csv("/global/projectb/scratch/byu24/at.sc.db/data/sample_metadata.csv")
#replace protoplast_loci filepath with actual filepath
protoplast_loci = read_csv("/global/projectb/scratch/byu24/at.sc.db/data/protoplast_loci.csv")%>%
  pull(Locus)

# Functions ------------------------------------------------------------------
#Load dataset into seurat matrix
get_dge <- function(dataset) {
  dge <- Seurat::Read10X(data.dir = paste0("/global/projectb/scratch/byu24/at.sc.db/scratch/",dataset,"/",dataset,"_star.Solo.out")) %>%
    `colnames<-`(c(paste0(colnames(.), "_", dataset)))
}

#Pull cell statistics from seurat matrix
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

#Filter dataset with threshold gene statistics
filter_dge <- function(dge, dge_stats, thresh_Gene = 200, expected_cells) {
  dge_stats <- dge_stats %>% 
    as_tibble(rownames = "Cell") %>% 
    filter(nGene >= thresh_Gene)
  
  nn_percentile <- expected_cells*0.01
  
  threshold <- dge_stats %>% 
    top_n(nn_percentile, nAt) %>%
    summarize(threshold = min(nAt)*0.05) %>%
    pull(threshold)
  
  keep_cells <- dge_stats %>% 
    filter(nAt >= threshold) %>%
    pull(Cell)
  
  dge <- dge[,keep_cells]
  
  #replace protoplast_loci filepath with actual filepath
  protoplast_loci = read_csv("/global/projectb/scratch/byu24/at.sc.db/data/protoplast_loci.csv")%>% 
    pull(Locus)
  
  keep_genes <- setdiff(rownames(dge)[str_detect(rownames(dge), "AT.G.....")], protoplast_loci)
  dge <- dge[keep_genes,]
}

#Creates Seurat Object after filtering
get_sobj <- function(dge_filtered, dge_stats, group) {
  sobj <- Seurat::CreateSeuratObject(
    counts = dge_filtered,
    project = group,
    meta.data = dge_stats) 
}

#Transforms Seurat object using PCA and UMAP
trans_dge <- function(sobj) {
  sobj <- sobj %>%
    SCTransform(vars.to.regress = "nCM", variable.features.n = 2000) %>%
    RunPCA(npcs=50) %>%
    RunUMAP(dims=1:30) %>%
    FindNeighbors(reduction = "pca", dims = 1:30) %>%
    FindClusters(resolution = 0.8)
}