# Clear plots
if(!is.null(dev.list())) dev.off()
# Clear console
cat("\014") 
# Clean workspace
rm(list=ls())

setwd('/global/projectb/scratch/byu24/at.sc.db/scratch')
library('Seurat')
library(ggplot2)
library(tidyverse)
library(furrr)
library(future)
options(future.globals.maxSize = 100000 * 1024^2)

at.final<-readRDS(file = "/global/projectb/scratch/byu24/at.sc.db/data/at.final.rds") #change the path to correct directory
# Plots cluster map of final database ------------------------------------------------------------------
DimPlot(at.final, group.by = "cell_type_simple")
ggsave(file="/global/projectb/scratch/byu24/at.sc.db/reports/final_simplecell.png", width = 30, height = 20, units = "cm")
DimPlot(at.final, reduction = "pca")
ggsave(file="/global/projectb/scratch/byu24/at.sc.db/reports/final_pca.png", width = 30, height = 20, units = "cm")

# Isolate specific cell types ------------------------------------------------------------------
#Pulls QC only data out
umap_coords <- at.final %>% Embeddings(reduction = "umap") %>% 
	as_tibble(rownames="Cell")
cell_data <- at.final@meta.data %>% 
	as_tibble(rownames="Cell") %>%
	left_join(umap_coords)

#Plots only QC cells
qc = cell_data %>%
	filter(cell_type == "QC")
ggplot(qc, aes(x=UMAP_1, y=UMAP_2, color=orig.ident)) + geom_point()
ggsave(file="/global/projectb/scratch/byu24/at.sc.db/reports/QC.png", width = 30, height = 20, units = "cm")

# Visualize cell identity composition per cluster ------------------------------------------------------------------
#Creates barplots of cell identity composition per cluster
ggplot(data=at.final@meta.data, aes(x=seurat_clusters, y=seurat_clusters, fill=cell_type_simple)) + geom_col()
  geom_bar(stat="Identity")
ggsave(file="/global/projectb/scratch/byu24/at.sc.db/reports/idents.png", width = 30, height = 20, units = "cm")

#Reclassify cluster IDs to group similar cell types together. Can be calculated manually (shown below) or with R.
new.cluster.ids <- c("10", "19", "5", "0", "13", "1", "24", "23", "22", "21", "20", "12", "3", "6", "25", "14", "7", "2", "8", "15", "26", "4", "16", "17", "18", "9", "11") #Numerical order for sorting.
names(new.cluster.ids) <- levels(at.final) #Attaches names to the levels of Seurat object
at.reorder <- RenameIdents(at.final, new.cluster.ids) #assigns new names to clusters
at.reorder@active.ident <- forcats::fct_relevel(at.reorder@active.ident, "0") #Reorders in numerical order. Save to new seurat object if desired.
ggplot(data=at.reorder@meta.data, aes(x=at.reorder@active.ident, y=seurat_clusters, fill=cell_type_simple)) +
  geom_bar(stat="Identity") #Creates new barplot in reordered form
ggsave(file="/global/projectb/scratch/byu24/at.sc.db/reports/idents_reorder.png", width = 30, height = 20, units = "cm")

#Plots new cluster map of final database with reclassified cluster IDs
DimPlot(at.reorder, group.by = "cell_type_simple") #can be plotted using different metadata if desired
ggsave(file="/global/projectb/scratch/byu24/at.sc.db/reports/final_simplecell_reorder.png", width = 30, height = 20, units = "cm")
DimPlot(at.reorder, reduction = "pca") #can be plotted using different metadata if desired
ggsave(file="/global/projectb/scratch/byu24/at.sc.db/reports/final_pca_reorder.png", width = 30, height = 20, units = "cm")

#If desired, save 'at.reorder' to a new R object. Not recommended overwriting 'at.final' if new cluster ordering is desired.
