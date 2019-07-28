#' @title Compute ICI scores from Arabidopsis data sets
#'
#' @description Takes in scRNAseq data, and computes ICI (index of cell
#' identity) scores,
#'
#' @usage get_ICI(df, ci = spec, ma = markers, sig = TRUE, ma.thresh = 0.15)
#'
#' @param df a data.frame containing a column, "Locus" with Arabidopsis locus
#' identifiers (e.g. AT1G30420). All other columns are STAMP barcodes containing
#' expression data for each cell for all loci.
#'
#' @param ci a data.frame containing spec scores for all loci in all cell types
#' considered. If this is not set (or set to NULL), the function will use spec
#' scores and cell types described in Elfroni, et al., 2016.
#'
#' @param ma a data.frame containing marker genes and their importance in
#' specifying each cell type considered. If this is not set (or set to NULL),
#' the function will use spec scores and cell types described in Elfroni, et
#' al., 2016.
#'
#' @param sig logical, whether or not to compute p-values for cell type calls.
#'
#' @param ma.thresh numeric, the marker importance threshold for ICI computation
#'
#' @details ICI scores are computed as previously escribed in
#' \href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4354993/}{Elfroni, et al.
#' (2016)}. The cell-type specific data was derived from previous FACS RNA-seq
#' or microarray data, and is contained in two data files: one file containing
#' the specificity score ('spec') and another containing marker genes and their
#' importance in each cell type. You can supply your own spec (ci) and marker
#' (ma) data to this function, however, if new data becomes available. Note,
#' setting sig = TRUE will randomly permute which genes are considered markers
#' to derive p-values (otherwise not computed) and will cause the function to
#' take a long time to execute.
#'
#' @import dplyr
#' @import pbapply
#'
#' @export

get_ICI <- function(df, ci = NULL, sig = TRUE, ma.thresh = 0.15, do.par = F) {
  ## If ci is not set, load from Birnbaum data
  if(is.null(ci)) ci <- root_ci

  ## Move locus identifier to row.names (makes subsetting easier)
  if("Locus" %in% colnames(ci)) ci <- column_to_rownames(ci, "Locus")
  universe <- sort(row.names(ci))

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

#' @title Summarize ICI score matrix
#'
#' @description Given an ICI score matrix (ICI scores for all cells for all
#' cell types), this function finds the most likely cell type, and returns the
#' cell type, p.value, and adjusted p-value associated with that call. The
#' function also collapses vascular-related cell types into a single category
#' ("Stele"), storing this identity in the cell_type_simple column.
#'
#' @param ici_scores a data.frame of ICI scores for all cell types for all
#' cells, as computed by get_ICI().
#'
#' @param p_thresh a p-value threshold
#'
#' @import tidyverse
#'
#' @export
summarize_ici <- function(ici_scores, p_thresh = 0.1) {
  stele_types <- c("Meri_Xylem", "Late_PPP", "Protoxylem", "Protophloem",
                   "Phloem_CC", "Pericycle", "Phloem", "Late_XPP")

  if(sum(!is.na(ici_scores$p.val)) > 0) {
    ici_summary <-  ici_scores %>%
      mutate(cell_type_simple = if_else(
        (cell_type %in% stele_types),
        "Stele",
        as.character(cell_type))) %>%
      group_by(Cell) %>%
      arrange(p.val) %>%
      summarize(
        cell_type_mix = if_else(
          sum(p.adj < p_thresh) > 0,
          paste(sort(cell_type[p.adj < p_thresh]), collapse = "/"),
          "Not Significant"),
        cell_type_simple_mix = if_else(
          sum(p.adj < p_thresh) > 0,
          paste(sort(unique(cell_type_simple[p.adj < p_thresh])), collapse = "/"),
          "Not Significant"),
        cell_type = cell_type[1],
        cell_type_simple = cell_type_simple[1],
        p.val = p.val[1],
        p.adj = p.adj[1],
        ici.score = ici.score[1],
        ici.score.norm = ici.score.norm[1]
      )
  } else {
    ici_summary <- ici_scores %>%
      group_by(Cell) %>%
      summarize(cell_type = cell_type[which.max(ici.score)],
                p.adj = NA,
                ici.score = max(ici.score),
                ici.score.norm = ici.score.norm[which.min(p.val)],
                p.val = NA) %>%
      mutate(cell_type_simple = case_when(
        cell_type %in% stele_types,
        "Stele",
        as.character(cell_type))) %>%
      mutate(cell_type_sig = NA,
             cell_type_simple_sig = NA)
  }
  return(ici_summary)
}

#' Compute individual ICI score and associated p-value. (Internal)
#'
#' @description Internal function fo computing ICI score. Computes ICI score
#' using method described by Birnbaum (Efroni et al., 2016)
#'
#' @param x numeric vector, single-cell transcript values
#' @param ci_ct numeric vector, spec scores for given cell type
#'
#' @return a vector of ICI scores
compute_ici_score <- function(x, ci_ct) {
  nt <- sum(ci_ct > 0)
  score.1 <- sum(x*ci_ct)/nt
  score.2 <- sum(x[ci_ct > 0] > 0)/nt
  return(score.1*score.2)
}

#' Compute random permutations of ICI scores for a given cell type. (Internal)
#'
#' @param x numeric vector, a single-cell transcriptome
#' @param ci_ct numeric vector, spec scores for a single cell type
#'
#' @return numeric vector of random ICI scores
compute_ici_pval <- function(x, ci_ct) {
  score <- compute_ici_score(x, ci_ct)
  p <- lapply(1:1000, function(y) {
    compute_ici_score(x = x, ci_ct = sample(ci_ct))
  })
  p <- sum(as.numeric(p) > score)/1000
  return(p)
}

#' Compute p-values for a single cell over all cell types (Internal)
#'
#' @param x numeric vector, a single cell transcriptome
#'
#' @return a vector of p-values over all cell types for a given cell.
compute_all_scores <- function(x, ci, type = 'pval') {
  apply(ci, 2, function(ci_ct) {
    if(type == 'pval') return(compute_ici_pval(x, ci_ct))
    if(type == 'score') return(compute_ici_score(x, ci_ct))
  })
}

