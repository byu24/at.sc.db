library(tidyverse)
library(drake)
library(Seurat)
library(future.batchtools)

source("src/analyze/process_dge_functions.R")
source("src/analyze/ici.R")

md <- readr::read_csv("data/sample_metadata.csv")

runs <- md %>%
  pull(Name) %>%
  unique()

my_plan <- drake_plan(
  protoplast_loci = target(
    readr::read_csv(file_in("data/protoplast_loci.csv")) %>%
      pull(Locus),
    hpc = FALSE),

  sample_metadata = target(readr::read_csv(file_in("data/sample_metadata.csv")), hpc = FALSE),

  root_ci = target(readr::read_csv(file_in("data/root_ci.csv")), hpc = FALSE),

  dge = target(
    get_dge(dataset = x),
    transform = map(x = (!!runs)),
    hpc = FALSE),

  stats = target(
    get_dge_stats(filtered),
    transform = map(filtered),
    hpc = FALSE),
  
  filtered = target(
    filter_dge(dge, protoplast_loci = protoplast_loci), 
    transform = map(dge, .id = dge),
    hpc = FALSE),

  sobj = target(
    get_sobj(filtered, dge_stats = stats, group = group),
    transform = map(
      filtered, 
      stats = !!syms(paste0("stats_dge_", runs)), 
      group = !!runs, .id = group),
    resources = list(ncpus = 1, memory = 4000, walltime = 7200)),

  norm = target(
    GetAssayData(sobj, assay = "SCT", slot = "data") %>% as_tibble(rownames = "Locus"),
    transform = map(sobj),
    hpc = FALSE),

  ICI = target(
    get_ICI(norm, ci = root_ci, sig = TRUE, do.par = F),
    transform = map(norm)
  ),

  sobj_list_raw = target(
    list(sobj),
    transform = combine(sobj),
    hpc = FALSE),

  sobj_features = SelectIntegrationFeatures(object.list = sobj_list_raw, nfeatures = 2000),

  sobj_list = PrepSCTIntegration(object.list = sobj_list_raw, anchor.features = sobj_features),

  sobj_anchors = FindIntegrationAnchors(
    object.list = sobj_list,
    normalization.method = "SCT",
    dims=1:50,
    anchor.features = sobj_features),

  sobj_integrated = IntegrateData(anchorset = sobj_anchors, dims = 1:50, normalization.method = "SCT") %>%
    RunPCA(npcs=100) %>%
    RunUMAP(dims=1:50)
)
options(future.globals.maxSize= 32000 * 1024^2)
make(my_plan, memory_strategy = "memory", target = "sobj_integrated")
make(my_plan, memory_strategy = "memory")
