library("Seurat")
library(ggplot2)
library(tidyverse)
library(furrr)
library(future)

get_ICI <- function(df, ci = NULL, sig = TRUE, ma.thresh = 0, do.par = F) {
  ## If ci is not set, load from Birnbaum data
  if(is.null(ci)) ci <- read.csv("/global/projectb/scratch/byu24/at.sc.db/data/root_ci.csv")
  
  ## Move locus identifier to row.names (makes subsetting easier)
  if("Locus" %in% colnames(ci)) ci <- column_to_rownames(ci, "Locus")
  universe <- sort(row.names(ci))
  df = df %>%
    GetAssayData(assay = "RNA", slot = "counts") %>% 
    as.data.frame() %>% 
    rownames_to_column("Locus")
  
  df <- right_join(df, tibble(Locus = universe)) %>%
    as.data.frame() %>%
    column_to_rownames("Locus")
  
  ## Filter markers to only those with > 0.15 scores. Set NA values in df ot 0.
  ci[ci < ma.thresh] = 0
  df[is.na(df)] <-  0
  
  ## Re-sort all matrices so that the locus order is the same
  df <- df[universe,]
  ci <- ci[universe,]
  
  clean_output <- function(x, cn) {
    do.call(what = rbind, x) %>%
      as.data.frame() %>%
      rownames_to_column("Cell") %>%
      tidyr::gather("cell_type", !!as.name(cn), -Cell)
  }
  future::plan(future.callr::callr, workers = 71)
  pvals <- furrr::future_map(df, function(x) compute_all_scores(x, ci, 'pval')) %>%
    clean_output(cn = 'p.val')
  
  scores <- furrr::future_map(df, function(x) compute_all_scores(x, ci, 'score')) %>%
    clean_output(cn = 'ici.score')
  
  ici_scores <- dplyr::full_join(pvals, scores) %>%
    mutate(p.adj = p.adjust(p.val, method = "BH")) %>%
    group_by(Cell) %>%
    mutate(ici.score.norm = ici.score/sum(ici.score)) %>%
    ungroup()
  
  return(ici_scores)
}
