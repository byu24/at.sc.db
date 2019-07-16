library(tidyverse)
get_dge <- function(dataset) {
  dge <- Seurat::Read10X(data.dir = paste0("scratch/",dataset,"/",dataset,"_star.Solo.out")) %>%
    `colnames<-`(c(paste0(colnames(.), "_", dataset)))
}

get_dge_stats <- function(dge) {
  nAt <- Matrix::colSums(dge[rownames(dge)[str_detect(rownames(dge), "AT.G.....")],]) %>%
    as_tibble(rownames = "Cell") %>%
    dplyr::rename(nAt = value)
  nGene <- Matrix::colSums(dge[rownames(dge)[str_detect(rownames(dge), "AT.G.....")],] > 0) %>% 
    as_tibble(rownames = "Cell") %>% 
    dplyr::rename(nGene = value)
  nMa <- Matrix::colSums(dge[rownames(dge)[!str_detect(rownames(dge), "AT.G.....")],]) %>%
    as_tibble(rownames = "Cell") %>%
    dplyr::rename(nMa = value)
  nCM <- Matrix::colSums(dge[rownames(dge)[str_detect(rownames(dge), "AT[MC]G.....")],]) %>%
    as_tibble(rownames = "Cell") %>%
    dplyr::rename(nCM = value)
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
    Seurat::SCTransform(variable.features.n = 3000)
}

filter_dge <- function(dge, dge_stats, thresh_Gene = 800, thresh_pAt = 0.98, thresh_pCM = 0.1) {
  keep_cells <- dge_stats %>% 
    as_tibble(rownames = "Cell") %>% 
    filter(pAt > 0.98 & nGene > 800 & pCM < 0.1, .preserve = T) %>% 
    pull(Cell)
  dge[,keep_cells]
}
