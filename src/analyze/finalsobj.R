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
	FindNeighbors(dims = 1:30) %>%
    FindClusters(resolution = 0.8)
	
at.final<-AddMetaData(object=at.integrated, metadata=all_ici)
saveRDS(at.final, file = "/global/projectb/scratch/byu24/at.sc.db/scratch/robjects/at.final.rds")
at.final@metadata

png(file="/global/projectb/scratch/byu24/at.sc.db/scratch/analysis/at_final.png")
DimPlot(at_final, group.by = "cell_type_simple")
dev.off()

