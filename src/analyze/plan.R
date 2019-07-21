library(drake)
library(dplyr)
library(ggplot2)

source("$BSCRATCH/at.sc.db/src/analyze/process_dge_functions.R")

groups <- readr::read_csv("$BSCRATCH/at.sc.db/data/sample_metadata.csv") %>% 
  pull(Group) %>% 
  unique()

my_plan <- drake_plan(
  sample_metadata = readr::read_csv(file_in("$BSCRATCH/at.sc.db/data/sample_metadata.csv")),
  dge = target(
    get_dge_group(group = x, sample_metadata = sample_metadata),
    transform = map(x = (!!groups))
  ),
  stats = target(
    get_dge_stats(dge),
    transform = map(dge)),
  
  filtered = target(
    filter_dge(dge, dge_stats = stats), 
    transform = map(dge, stats = !!syms(paste0("stats_dge_", groups)), .id = dge)
  ),

  sobj = target(
    get_sobj(filtered, dge_stats = stats, group = group),
    transform = map(
      filtered, 
      stats = !!syms(paste0("stats_dge_", groups)), 
      group = !!groups, .id = group
    )
  ),
  
  sobj_list = target(
    list(sobj),
    transform = combine(sobj)
  )
)

make(my_plan)
