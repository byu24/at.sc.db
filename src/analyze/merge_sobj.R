
library("Seurat")
library(ggplot2)
library("tidyverse")
library(future)
options(future.globals.maxSize = 100000 * 1024^2)

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

at.integrated <- RunPCA(at.integrated, verbose = FALSE)
at.integrated <- RunUMAP(at.integrated, dims = 1:30)
saveRDS(at.integrated, file = "/global/projectb/scratch/byu24/at.sc.db/scratch/robjects/at_integrated.rds")



