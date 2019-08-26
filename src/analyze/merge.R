# Clear plots
if(!is.null(dev.list())) dev.off()
# Clear console
cat("\014") 
# Clean workspace
rm(list=ls())

source("/global/projectb/scratch/byu24/at.sc.db/src/analyze/process_ici.R")
setwd('/global/projectb/scratch/byu24/at.sc.db/scratch')
library('Seurat')
library(ggplot2)
library(tidyverse)
library(furrr)
library(future)
options(future.globals.maxSize = 100000 * 1024^2)

# Merge samples to Seurat object ------------------------------------------------------------------
#Load all R objects
ss_001 = readRDS(file = "/global/projectb/scratch/byu24/at.sc.db/scratch/robjects/ss_001.rds")
ss_002 = readRDS(file = "/global/projectb/scratch/byu24/at.sc.db/scratch/robjects/ss_002.rds")
ss_003 = readRDS(file = "/global/projectb/scratch/byu24/at.sc.db/scratch/robjects/ss_003.rds")
sc_004 = readRDS(file = "/global/projectb/scratch/byu24/at.sc.db/scratch/robjects/sc_004.rds")
sc_005 = readRDS(file = "/global/projectb/scratch/byu24/at.sc.db/scratch/robjects/sc_005.rds")
sc_006 = readRDS(file = "/global/projectb/scratch/byu24/at.sc.db/scratch/robjects/sc_006.rds")
sc_007 = readRDS(file = "/global/projectb/scratch/byu24/at.sc.db/scratch/robjects/sc_007.rds")
sc_008 = readRDS(file = "/global/projectb/scratch/byu24/at.sc.db/scratch/robjects/sc_008.rds")
sc_009 = readRDS(file = "/global/projectb/scratch/byu24/at.sc.db/scratch/robjects/sc_009.rds")
sc_010 = readRDS(file = "/global/projectb/scratch/byu24/at.sc.db/scratch/robjects/sc_010.rds")
zc_011 = readRDS(file = "/global/projectb/scratch/byu24/at.sc.db/scratch/robjects/zc_011.rds")
rc_012 = readRDS(file = "/global/projectb/scratch/byu24/at.sc.db/scratch/robjects/rc_012.rds")
rc_013 = readRDS(file = "/global/projectb/scratch/byu24/at.sc.db/scratch/robjects/rc_013.rds")
rc_014 = readRDS(file = "/global/projectb/scratch/byu24/at.sc.db/scratch/robjects/rc_014.rds")
js_015 = readRDS(file = "/global/projectb/scratch/byu24/at.sc.db/scratch/robjects/js_015.rds")
js_016 = readRDS(file = "/global/projectb/scratch/byu24/at.sc.db/scratch/robjects/js_016.rds")
jsh_017 = readRDS(file = "/global/projectb/scratch/byu24/at.sc.db/scratch/robjects/jsh_017.rds")
dc_018 = readRDS(file = "/global/projectb/scratch/byu24/at.sc.db/scratch/robjects/dc_018.rds")
dc_019 = readRDS(file = "/global/projectb/scratch/byu24/at.sc.db/scratch/robjects/dc_019.rds")

#Merge all Seurat objects into a list
at.list <- list(ss_001, ss_002, ss_003, sc_004, sc_005, sc_006, sc_007, sc_008, sc_009, sc_010,
                    zc_011, rc_012, rc_013, rc_014, js_015, js_016, jsh_017, dc_018, dc_019)

at.features <- SelectIntegrationFeatures(object.list = at.list, nfeatures = 2000)
at.list <- PrepSCTIntegration(object.list = at.list, anchor.features = at.features, 
                                    verbose = FALSE)

at.anchors <- FindIntegrationAnchors(object.list = at.list, normalization.method = "SCT", 
                                           anchor.features = at.features, verbose = FALSE)
at.integrated <- IntegrateData(anchorset = at.anchors, normalization.method = "SCT", 
                                     verbose = FALSE)

at.integrated <- at.integrated %>% 
	RunPCA(verbose = FALSE) %>%
	RunUMAP(dims = 1:30) %>%
	FindNeighbors(reduction = "pca", dims = 1:30) %>%
    FindClusters(resolution = 0.6) #resolution adjusted based on desired clarity

at.integrated
at.integrated@meta.data

# Merge ICI scores to Seurat ------------------------------------------------------------------
# Convert summary ICI scores to dataframe
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

#List ICI scores into one object
ls_ici <- list(ici_ss_001, ici_ss_002, ici_ss_003, ici_sc_004, ici_sc_005, ici_sc_006, ici_sc_007, ici_sc_008, ici_sc_009, ici_sc_010,
                    ici_zc_011, ici_rc_012, ici_rc_013, ici_rc_014, ici_js_015, ici_js_016, ici_jsh_017, ici_dc_018, ici_dc_019)

all_ici<-do.call(rbind, ls_ici)

#Add ICI scores to Seurat object
at.final<-AddMetaData(object=at.integrated, metadata=all_ici)

saveRDS(at.final, file = "/global/projectb/scratch/byu24/at.sc.db/data/at.final.rds")
at.final
at.final@meta.data
