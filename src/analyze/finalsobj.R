# Clear plots
if(!is.null(dev.list())) dev.off()
# Clear console
cat("\014") 
# Clean workspace
rm(list=ls())

source("/global/projectb/scratch/byu24/at.sc.db/src/analyze/ici.R")
setwd('/global/projectb/scratch/byu24/at.sc.db/scratch')
library('Seurat')
library(ggplot2)
library(tidyverse)
library(furrr)
library(future)
options(future.globals.maxSize = 100000 * 1024^2)

sum_df <- function(ici_df) {
	ici_df <-readRDS(file = paste0("/global/projectb/scratch/byu24/at.sc.db/scratch/robjects/",ici_df,".rds"))
	
	{
		sum_df = summarize_ici(ici_df)
	}
	
	sum_df<- sum_df %>% as.data.frame() %>% column_to_rownames("Cell")
	}

ici_ss_001<-sum_df("ici_ss_001")
ici_ss_002<-sum_df("ici_ss_002")
ici_ss_003<-sum_df("ici_ss_003")
ici_sc_004<-sum_df("ici_sc_004")
ici_sc_005<-sum_df("ici_sc_005")
ici_sc_006<-sum_df("ici_sc_006")
ici_sc_007<-sum_df("ici_sc_007")
ici_sc_008<-sum_df("ici_sc_008")
ici_sc_009<-sum_df("ici_sc_009")
ici_sc_010<-sum_df("ici_sc_010")
ici_zc_011<-sum_df("ici_zc_011")
ici_rc_012<-sum_df("ici_rc_012")
ici_rc_013<-sum_df("ici_rc_013")
ici_rc_014<-sum_df("ici_rc_014")
ici_js_015<-sum_df("ici_js_015")
ici_js_016<-sum_df("ici_js_016")
ici_jsh_017<-sum_df("ici_jsh_017")
ici_dc_018<-sum_df("ici_dc_018")
ici_dc_019<-sum_df("ici_dc_019")

ls_ici <- list(ici_ss_001, ici_ss_002, ici_ss_003, ici_sc_004, ici_sc_005, ici_sc_006, ici_sc_007, ici_sc_008, ici_sc_009, ici_sc_010,
                    ici_zc_011, ici_rc_012, ici_rc_013, ici_rc_014, ici_js_015, ici_js_016, ici_jsh_017, ici_dc_018, ici_dc_019)

all_ici<-do.call(rbind, ls_ici)

at.integrated<-readRDS(file = "/global/projectb/scratch/byu24/at.sc.db/scratch/robjects/at_integrated.rds")	

at.integrated <- at.integrated %>% 
	FindNeighbors(reduction = "pca", dims = 1:30) %>%
    FindClusters(resolution = 0.60)
	
at.final<-AddMetaData(object=at.integrated, metadata=all_ici)
saveRDS(at.final, file = "/global/projectb/scratch/byu24/at.sc.db/scratch/robjects/at.final.rds")

DimPlot(at.final, group.by = "cell_type_simple")
ggsave(file="/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/at_final_cts.png", width = 30, height = 20, units = "cm")
DimPlot(at.final, reduction = "pca")
ggsave(file="/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/at_final_pca.png", width = 30, height = 20, units = "cm")
at.final@meta.data

at.final<-readRDS(file = "/global/projectb/scratch/byu24/at.sc.db/scratch/robjects/at.final.rds")

#Pulls QC only data out
umap_coords <- at.final %>% Embeddings(reduction = "umap") %>% 
	as_tibble(rownames="Cell")
cell_data <- at.final@meta.data %>% 
	as_tibble(rownames="Cell") %>%
	left_join(umap_coords)

qconly = cell_data %>%
	filter(cell_type == "QC")
ggplot(qconly, aes(x=UMAP_1, y=UMAP_2, color=orig.ident)) + geom_point()
ggsave(file="/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/QC.png", width = 30, height = 20, units = "cm")


ggplot(data=at.final@meta.data, aes(x=Idents(at.final), y=seurat_clusters, fill=cell_type_simple)) + geom_col()
  geom_bar(stat="Identity")
ggsave(file="/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/idents.png", width = 30, height = 20, units = "cm")


new.cluster.ids <- c("Endodermis_1", "Trichoblast_1", "Cortex_1", "Atrichoblast_1", "Stele_1", "Atrichoblast_2", "Trichoblast_2", "Trichoblast_3", "Trichoblast_4", "Trichoblast_5", "Trichoblast_6", "Endodermis_2", "Columella_1", "Cortex_2", "Trichoblast_7", "Stele_2", "Cortex_3", "Atrichoblast_3", "Cortex_4", "Stele_3", "Trichoblast_8", "Columella_2", "Stele_4", "Stele_5", "Stele_6", "Cortex_5", "Endodermis_3")

new.cluster.ids <- c("10", "19", "5", "0", "13", "1", "24", "23", "22", "21", "20", "12", "3", "6", "25", "14", "7", "2", "8", "15", "26", "4", "16", "17", "18", "9", "11")
names(new.cluster.ids) <- levels(at.final)
test <- RenameIdents(at.final, new.cluster.ids)
test@active.ident <- reorder(test@active.ident)
ggplot(data=test@meta.data, aes(x=test@active.ident, y=seurat_clusters, fill=cell_type_simple)) +
  geom_bar(stat="Identity")
ggsave(file="/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/reorder.png", width = 30, height = 20, units = "cm")


ggplot(at.final@meta.data, aes(x = reorder(cell_type_simple, -perc), y = perc)) + geom_bar(stat = "identity")

test@active.ident <- forcats::fct_relevel(test@active.ident, "0")

