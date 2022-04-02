library("Seurat")
library("ggplot2")
library("patchwork")
library("dplyr")
library("hdf5r")
source("/share/pub/dengcy/Singlecell/CRC_space/src/Functions.r")
setwd("/share/pub/dengcy/Singlecell/CRC_space/5.spatial_progress")
source("/share/pub/dengcy/Singlecell/CRC_space/src/library.r")
Spatial_RST2bei<-readRDS("Spatial_RST2bei.rds")
Spatial_RNST2<-readRDS("Spatial_RNST2.rds")
Spatial_CST1<-readRDS("Spatial_CST1.rds")
Spatial_CNST1<-readRDS("Spatial_CNST1.rds")

Spatial_RST2bei<-cluster_pca_umap(Spatial_RST2bei, assay = "SCT",reduction="pca",cluster_res = 0.3)
Spatial_RNST2<-cluster_pca_umap(Spatial_RNST2, assay = "SCT",reduction="pca",cluster_res = 0.3)
Spatial_CST1<-cluster_pca_umap(Spatial_CST1, assay = "SCT",reduction="pca",cluster_res = 0.3)
Spatial_CNST1<-cluster_pca_umap(Spatial_CNST1, assay = "SCT",reduction="pca",cluster_res = 0.3)

saveRDS(Spatial_RST2bei,file="Spatial_RST2bei.rds")
saveRDS(Spatial_RNST2,file="Spatial_RNST2.rds")
saveRDS(Spatial_CST1,file="Spatial_CST1.rds")
saveRDS(Spatial_CNST1,file="Spatial_CNST1.rds")

plan("multiprocess", workers = N_WORKERS)

  Spatial_RST2bei_Allmarkers <- parallelFindAllMarkers(Spatial_RST2bei)
  save(Spatial_RST2bei_Allmarkers,file="Spatial_RST2bei_Allmarkers.RData")

  Spatial_RNST2_Allmarkers <- parallelFindAllMarkers(Spatial_RNST2)
  save(Spatial_RNST2_Allmarkers,file="Spatial_RNST2_Allmarkers.RData")

  Spatial_CST1_Allmarkers <- parallelFindAllMarkers(Spatial_CST1)
  save(Spatial_CST1_Allmarkers,file="Spatial_CST1_Allmarkers.RData")

  Spatial_CNST1_Allmarkers <- parallelFindAllMarkers(Spatial_CNST1)
  save(Spatial_CNST1_Allmarkers,file="Spatial_CNST1_Allmarkers.RData")
library("celldex")
library("SingleR")
library("pracma")
ref <- HumanPrimaryCellAtlasData() 

#######################################################
pred.brca_RST2bei <- singler_func(Spatial_RST2bei,ref)
pred.brca_RNST2 <- singler_func(Spatial_RNST2,ref)
pred.brca_CST1 <- singler_func(Spatial_CST1,ref)
pred.brca_CNST1 <- singler_func(Spatial_CNST1,ref)

pred_final_RST2bei<-pred_final(pred.brca=pred.brca_RST2bei,CRC_scrna=Spatial_RST2bei)
pred_final_RNST2<-pred_final(pred.brca=pred.brca_RNST2,CRC_scrna=Spatial_RNST2)
pred_final_CST1<-pred_final(pred.brca=pred.brca_CST1,CRC_scrna=Spatial_CST1)
pred_final_CNST1<-pred_final(pred.brca=pred.brca_CNST1,CRC_scrna=Spatial_CNST1)

Spatial_RST2bei$singler <- pred.brca_RST2bei$labels
Spatial_RNST2$singler <- pred.brca_RNST2$labels
Spatial_CST1$singler <- pred.brca_CST1$labels
Spatial_CNST1$singler <- pred.brca_CNST1$labels

saveRDS(Spatial_RST2bei,file="Spatial_RST2bei.rds")
saveRDS(Spatial_RNST2,file="Spatial_RNST2.rds")
saveRDS(Spatial_CST1,file="Spatial_CST1.rds")
saveRDS(Spatial_CNST1,file="Spatial_CNST1.rds")