library(tidyverse)
get_dge <- function(dataset) {
  dge <- Seurat::Read10X(data.dir = paste0("scratch/",dataset,"/",dataset,"_star.Solo.out")) %>%
    `colnames<-`(c(paste0(colnames(.), "_", dataset)))
}

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

get_dge_group <- function(group, sample_metadata) {
  datasets <- dplyr::filter(sample_metadata, Group %in% group) %>%
    dplyr::pull(Name)
  dge <- map(datasets, get_dge) %>% reduce(cbind)
}

get_sobj <- function(dge, dge_stats, group) {
  dge_stats <- get_dge_stats(dge)
  sobj <- Seurat::CreateSeuratObject(
    counts = dge,
    project = group,
    meta.data = dge_stats) %>%
    SCTransform(vars.to.regress = "nCM", variable.features.n = 2000) %>%
    RunPCA(npcs=50) %>%
    RunUMAP(dims=1:30)
}

filter_dge <- function(dge, thresh_Gene = 200, max_umi = 50000, thresh_pAt = 0.98, thresh_pCM = 0.1, protoplast_loci) {
  expressed_genes <- Matrix::rowSums(dge > 0) %>%
    enframe("Locus", "nCells") %>%
    filter(nCells >= 1) %>%
    pull(Locus)
  dge <- dge[expressed_genes,]
  dge_stats <- get_dge_stats(dge) %>%
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
    filter(pAt >= 0.98 & nAt > max(1000, threshold) & nAt <= max_umi) %>%
    pull(Cell)
  sprintf("Threshold UMI set at %g; keeping %g cells", max(1000, threshold), length(keep_cells))
  keep_genes <- setdiff(rownames(dge)[str_detect(rownames(dge), "AT.G.....")], protoplast_loci)
  dge <- dge[keep_genes,]
  dge <- dge[,keep_cells]
}
