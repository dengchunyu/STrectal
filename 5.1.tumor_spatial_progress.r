library(ggpubr)
library(Seurat)
library(SeuratData)
library(dplyr)
library(ggtext)

source("/share/pub/dengcy/Singlecell/CC_space/src/Functions.r")
setwd("/share/pub/dengcy/Singlecell/CC_space/5.1.tumor_spatial_progress")
#source("/share/pub/dengcy/Singlecell/CC_space/src/library.r")
#############
#Epithelial_cells stem cell
############
#load("/share/pub/dengcy/Singlecell/CC_space/1.scRNA_seq_prepare/CC_scrna.RData")
load("/share/pub/dengcy/Singlecell/CC_space/3.Epi_cells_stage/cancer_scrna_stage1.RData")
load("/share/pub/dengcy/Singlecell/CC_space/3.Epi_cells_stage/cancer_scrna_stage2.RData")
###################
#将肿瘤细胞的各个cluster映射到不同的空间位置
################
Spatial_RST2bei<-readRDS("/share/pub/dengcy/Singlecell/CC_space/5.spatial_progress/Spatial_RST2bei.rds")
#Spatial_RNST2<-readRDS("Spatial_RNST2.rds")
Spatial_CST1<-readRDS("/share/pub/dengcy/Singlecell/CC_space/5.spatial_progress/Spatial_CST1.rds")
#Spatial_CNST1<-readRDS("Spatial_CNST1.rds")

####################
#空间图像的聚类映射位置
###################
##导入函数
files<-list.files("/share/pub/dengcy/Singlecell/CC_space/spotlight/R")
for(i in 1:length(files)){
	source(paste0("/share/pub/dengcy/Singlecell/CC_space/spotlight/R/",files[i]))
}
##########################################
############CC1
####1.1state
stage1_markers <- Seurat::FindAllMarkers(object = cancer_scrna_stage1,
                                              assay = "RNA",
                                              slot = "data",
                                              verbose = TRUE,
                                              only.pos = TRUE)
saveRDS(object = stage1_markers, file = here::here("stage1_markers.rds"))
###
spotlight_cancer_cc1 <- spotlight_deconvolution(
  se_sc = cancer_scrna_stage1,
  counts_spatial = Spatial_CST1@assays$Spatial@counts,
  clust_vr = "State", # Variable in sc_seu containing the cell-type annotation
  cluster_markers = stage1_markers, # Dataframe with the marker genes
  cl_n = 100, # number of cells per cell type to use
  hvg = 2000, # Number of HVG to use
  ntop = NULL, # How many of the marker genes to use (by default all)
  transf = "uv", # Perform unit-variance scaling per cell and spot prior to factorzation and NLS
  method = "nsNMF", # Factorization method
  min_cont = 0 # Remove those cells contributing to a spot below a certain threshold
  )

saveRDS(object = spotlight_cancer_cc1, file = here::here("spotlight_cancer_cc1.RDS"))

nmf_mod1 <- spotlight_cancer_cc1[[1]]
decon_mtrx1 <- spotlight_cancer_cc1[[2]]
#Assess deconvolution
h <- NMF::coef(nmf_mod1[[1]])
rownames(h) <- paste("Topic", 1:nrow(h), sep = "_")
topic_profile_plts <- dot_plot_profiles_fun(
  h = h,
  train_cell_clust = nmf_mod1[[2]])

topic_profile_plts[[2]] + ggplot2::theme(
  axis.text.x = ggplot2::element_text(angle = 90),
  axis.text = ggplot2::element_text(size = 12))

#Visualization
# This is the equivalent to setting min_cont to 0.04
decon_mtrx_sub <- decon_mtrx1[, colnames(decon_mtrx1) != "res_ss"]
decon_mtrx_sub[decon_mtrx_sub < 0.0001] <- 0
decon_mtrx1 <- cbind(decon_mtrx_sub, "res_ss" = decon_mtrx1[, "res_ss"])
rownames(decon_mtrx1) <- colnames(Spatial_CST1)

decon_df <- decon_mtrx1 %>%
  data.frame() %>%
  tibble::rownames_to_column("barcodes")

Spatial_CST1@meta.data <- Spatial_CST1@meta.data %>%
  tibble::rownames_to_column("barcodes") %>%
  dplyr::left_join(decon_df, by = "barcodes") %>%
  tibble::column_to_rownames("barcodes")

#Specific cell-types
#we can use the standard Seurat::SpatialFeaturePlot to view predicted celltype proportions one at a time.
colnames(Spatial_CST1@meta.data)[15:19]<-paste0("State",1:5)
pdf("spotlight_cancer_stage1_FeaturePlot.pdf",width=7,height=7)
Seurat::SpatialFeaturePlot(
  object = Spatial_CST1,
  features = c("State1", "State2", "State3", "State4", "State5"),
  image.alpha=0.8,
   stroke=0.3,
  alpha = c(0.2, 1))+
 scale_colour_manual(name = "stageI", values =c("#f38181","#fce38a","#eaffd0","#95e1d3","#64958f"))
dev.off()

Spatial_CST1@meta.data$State1<-Spatial_CST1@meta.data$X1
Spatial_CST1@meta.data$State2<-Spatial_CST1@meta.data$X2
Spatial_CST1@meta.data$State3<-Spatial_CST1@meta.data$X3
Spatial_CST1@meta.data$State4<-Spatial_CST1@meta.data$X4
Spatial_CST1@meta.data$State5<-Spatial_CST1@meta.data$X5

Spatial_CST1@meta.data$X1<-NULL

#Spatial scatterpiesPlot spot composition of all the spots.
#cell_types_all <- colnames(decon_mtrx1)[which(colnames(decon_mtrx1) != "res_ss")]
pdf("SPOTlight_spatial_scatterpie_cancer_CC1.pdf")
spatial_scatterpie(se_obj = Spatial_CST1,
                              cell_types_all = c("State1", "State2", "State3", "State4", "State5"),
                              img_path = "/share/pub/dengcy/Singlecell/CC_space/gcy_analysis/Analysis_200860-01/CST1/2.Basic_analysis/2.3.h5_files/spatial/tissue_lowres_image.png",
                              pie_scale = 0.4)
dev.off()



####1.2cluster
Idents(cancer_scrna_stage1)<-cancer_scrna_stage1$Cluster
stage1.2_markers <- Seurat::FindAllMarkers(object = cancer_scrna_stage1,
                                              assay = "RNA",
                                              slot = "data",
                                              verbose = TRUE,
                                              only.pos = TRUE)
saveRDS(object = stage1.2_markers, file = here::here("stage1.2_markers.rds"))
stage1.2_markers<-readRDS("stage1.2_markers.rds")
###
spotlight_cancer_cc1_cluster <- spotlight_deconvolution(
  se_sc = cancer_scrna_stage1,
  counts_spatial = Spatial_CST1@assays$Spatial@counts,
  clust_vr = "Cluster", # Variable in sc_seu containing the cell-type annotation
  cluster_markers = stage1.2_markers, # Dataframe with the marker genes
  cl_n = 500, # number of cells per cell type to use
  hvg = 2000, # Number of HVG to use
  ntop = NULL, # How many of the marker genes to use (by default all)
  transf = "uv", # Perform unit-variance scaling per cell and spot prior to factorzation and NLS
  method = "nsNMF", # Factorization method
  min_cont = 0 # Remove those cells contributing to a spot below a certain threshold
  )

saveRDS(object = spotlight_cancer_cc1_cluster, file = here::here("spotlight_cancer_cc1_cluster.RDS"))
spotlight_cancer_cc1<-readRDS("spotlight_cancer_cc1_cluster.RDS")
nmf_mod1 <- spotlight_cancer_cc1[[1]]
decon_mtrx1 <- spotlight_cancer_cc1[[2]]
#Assess deconvolution
h <- NMF::coef(nmf_mod1[[1]])
rownames(h) <- paste("Topic", 1:nrow(h), sep = "_")
topic_profile_plts <- dot_plot_profiles_fun(
  h = h,
  train_cell_clust = nmf_mod1[[2]])

topic_profile_plts[[2]] + ggplot2::theme(
  axis.text.x = ggplot2::element_text(angle = 90),
  axis.text = ggplot2::element_text(size = 12))

#Visualization
# This is the equivalent to setting min_cont to 0.04
decon_mtrx_sub <- decon_mtrx1[, colnames(decon_mtrx1) != "res_ss"]
decon_mtrx_sub[decon_mtrx_sub < 0.0001] <- 0
decon_mtrx1 <- cbind(decon_mtrx_sub, "res_ss" = decon_mtrx1[, "res_ss"])
rownames(decon_mtrx1) <- colnames(Spatial_CST1)

decon_df <- decon_mtrx1 %>%
  data.frame() %>%
  tibble::rownames_to_column("barcodes")

Spatial_CST1@meta.data <- Spatial_CST1@meta.data %>%
  tibble::rownames_to_column("barcodes") %>%
  dplyr::left_join(decon_df, by = "barcodes") %>%
  tibble::column_to_rownames("barcodes")

#Specific cell-types
#we can use the standard Seurat::SpatialFeaturePlot to view predicted celltype proportions one at a time.
#colnames(Spatial_CST1@meta.data)[15:19]<-paste0("State",1:5)
color_CRC<-c("#7FBCD2","#A5F1E9","#FFEEAF","#7FB77E","#F7F6DC","#B1D7B4","#FFC090",
  "#F3E0B5","#FFAE6D","#FECD70","#E3C770","#94B49F","#CEE5D0","#ECB390","#F2D7D9",
  "#D3CEDF","#B2C8DF","#C4D7E0","#F8F9D7","#76BA99","#ADCF9F",
  "#CED89E","#FFDCAE","#FFE3A9","#FFC3C3","#FF8C8C","#97C4B8","#F1F0C0","#B1BCE6")
pdf("spotlight_cancer_stage1_cluster_FeaturePlot.pdf",width=7,height=7)
Seurat::SpatialFeaturePlot(
  object = Spatial_CST1,
  features = c('X1','X2','X3','X4','X5','X6','X7'),
  image.alpha=0.8,
   stroke=0.3,
  alpha = c(0.2, 1))+
 scale_colour_manual(name = "stageI", values =color_CRC[1:7])
dev.off()

Spatial_CST1@meta.data$cluster1<-Spatial_CST1@meta.data$X1
Spatial_CST1@meta.data$cluster2<-Spatial_CST1@meta.data$X2
Spatial_CST1@meta.data$cluster3<-Spatial_CST1@meta.data$X3
Spatial_CST1@meta.data$cluster4<-Spatial_CST1@meta.data$X4
Spatial_CST1@meta.data$cluster5<-Spatial_CST1@meta.data$X5
Spatial_CST1@meta.data$cluster6<-Spatial_CST1@meta.data$X6
Spatial_CST1@meta.data$cluster7<-Spatial_CST1@meta.data$X7


packageurl='https://cran.r-project.org/src/contrib/Archive/randomForest/randomForest_4.6-14.tar.gz'
install.packages(packageurl, repos=NULL, type="source")
#Spatial scatterpiesPlot spot composition of all the spots.
#cell_types_all <- colnames(decon_mtrx1)[which(colnames(decon_mtrx1) != "res_ss")]
pdf("SPOTlight_spatial_scatterpie_cancer_CC1.pdf")
spatial_scatterpie(se_obj = Spatial_CST1,
                              cell_types_all = c('X1','X2','X3','X4','X5','X6','X7'),
                              img_path = "/share/pub/dengcy/Singlecell/CC_space/gcy_analysis/Analysis_200860-01/CST1/2.Basic_analysis/2.3.h5_files/spatial/tissue_lowres_image.png",
                              pie_scale = 0.4)
dev.off()

saveRDS(object = Spatial_CST1, file = "/share/pub/dengcy/Singlecell/CC_space/5.spatial_progress/Spatial_CST1.rds")
#######################CC2
####2.1stage

stage2_markers <- Seurat::FindAllMarkers(object = cancer_scrna_stage2,
                                              assay = "RNA",
                                              slot = "data",
                                              verbose = TRUE,
                                              only.pos = TRUE)
saveRDS(object = stage2_markers, file = here::here("stage2_markers.RDS"))
stage2_markers <-readRDS("/share/pub/dengcy/Singlecell/CC_space/5.1.tumor_spatial_progress/stage2_markers.RDS")

cancer_scrna_stage2_13<-subset(cancer_scrna_stage2,State %in% c(1,3))
stage2_markers_13 <- Seurat::FindAllMarkers(object = cancer_scrna_stage2_13,
                                              assay = "RNA",
                                              slot = "data",
                                              verbose = TRUE,
                                              only.pos = TRUE)

cancer_scrna_stage2_12<-subset(cancer_scrna_stage2,State %in% c(1,2))
stage2_markers_13 <- Seurat::FindAllMarkers(object = cancer_scrna_stage2_12,
                                              assay = "RNA",
                                              slot = "data",
                                              verbose = TRUE,
                                              only.pos = TRUE)


###这一步有点慢,半个多小时。
spotlight_cancer_cc2_13 <- spotlight_deconvolution(
  se_sc = cancer_scrna_stage2_13,
  counts_spatial = Spatial_RST2bei@assays$Spatial@counts,
  clust_vr = "State", # Variable in sc_seu containing the cell-type annotation
  cluster_markers = stage2_markers_13, # Dataframe with the marker genes
  cl_n = 3000, # number of cells per cell type to use
  hvg = 5000, # Number of HVG to use
  ntop = NULL, # How many of the marker genes to use (by default all)
  transf = "uv", # Perform unit-variance scaling per cell and spot prior to factorzation and NLS
  method = "nsNMF", # Factorization method
  min_cont = 0 # Remove those cells contributing to a spot below a certain threshold
  )

saveRDS(object = spotlight_cancer_cc2_13, file = "/share/pub/dengcy/Singlecell/CC_space/5.1.tumor_spatial_progress/spotlight_cancer_cc2_13.RDS")
nmf_mod2 <- spotlight_cancer_cc2_13[[1]]
decon_mtrx2 <- spotlight_cancer_cc2_13[[2]]
#Assess deconvolution
h <- NMF::coef(nmf_mod2[[1]])
rownames(h) <- paste("Topic", 1:nrow(h), sep = "_")
topic_profile_plts <- dot_plot_profiles_fun(
  h = h,
  train_cell_clust = nmf_mod2[[2]])

topic_profile_plts[[2]] + ggplot2::theme(
  axis.text.x = ggplot2::element_text(angle = 90),
  axis.text = ggplot2::element_text(size = 12))

#Visualization
# This is the equivalent to setting min_cont to 0.04
decon_mtrx_sub <- decon_mtrx2[, colnames(decon_mtrx2) != "res_ss"]
decon_mtrx_sub[decon_mtrx_sub < 0.0001] <- 0
decon_mtrx2 <- cbind(decon_mtrx_sub, "res_ss" = decon_mtrx2[, "res_ss"])
rownames(decon_mtrx2) <- colnames(Spatial_RST2bei)

decon_df <- decon_mtrx2 %>%
  data.frame() %>%
  tibble::rownames_to_column("barcodes")

Spatial_RST2bei@meta.data <- Spatial_RST2bei@meta.data %>%
  tibble::rownames_to_column("barcodes") %>%
  dplyr::left_join(decon_df, by = "barcodes") %>%
  tibble::column_to_rownames("barcodes")

#Specific cell-types
#we can use the standard Seurat::SpatialFeaturePlot to view predicted celltype proportions one at a time.
#colnames(Spatial_RST2bei@meta.data)[15:17]<-paste0("State",1:3)

pdf("spotlight_cancer_CC2_13_FeaturePlot.pdf",width=7,height=7)
Seurat::SpatialFeaturePlot(
  object = Spatial_RST2bei,
  features = c("X1.y" ,"X3.y"),
  image.alpha=0.8,
   stroke=0.3,
  alpha = c(0.2, 1))
dev.off()

spotlight_cancer_cc2_12 <- spotlight_deconvolution(
  se_sc = cancer_scrna_stage2_12,
  counts_spatial = Spatial_RST2bei@assays$Spatial@counts,
  clust_vr = "State", # Variable in sc_seu containing the cell-type annotation
  cluster_markers = stage2_markers_13, # Dataframe with the marker genes
  cl_n = 3000, # number of cells per cell type to use
  hvg = 5000, # Number of HVG to use
  ntop = NULL, # How many of the marker genes to use (by default all)
  transf = "uv", # Perform unit-variance scaling per cell and spot prior to factorzation and NLS
  method = "nsNMF", # Factorization method
  min_cont = 0 # Remove those cells contributing to a spot below a certain threshold
  )

saveRDS(object = spotlight_cancer_cc2_12, file = "/share/pub/dengcy/Singlecell/CC_space/5.1.tumor_spatial_progress/spotlight_cancer_cc2_12.RDS")

nmf_mod2 <- spotlight_cancer_cc2_13[[1]]
decon_mtrx2 <- spotlight_cancer_cc2_13[[2]]
#Assess deconvolution
h <- NMF::coef(nmf_mod2[[1]])
rownames(h) <- paste("Topic", 1:nrow(h), sep = "_")
topic_profile_plts <- dot_plot_profiles_fun(
  h = h,
  train_cell_clust = nmf_mod2[[2]])

topic_profile_plts[[2]] + ggplot2::theme(
  axis.text.x = ggplot2::element_text(angle = 90),
  axis.text = ggplot2::element_text(size = 12))

#Visualization
# This is the equivalent to setting min_cont to 0.04
decon_mtrx_sub <- decon_mtrx2[, colnames(decon_mtrx2) != "res_ss"]
decon_mtrx_sub[decon_mtrx_sub < 0.0001] <- 0
decon_mtrx2 <- cbind(decon_mtrx_sub, "res_ss" = decon_mtrx2[, "res_ss"])
rownames(decon_mtrx2) <- colnames(Spatial_RST2bei)

decon_df <- decon_mtrx2 %>%
  data.frame() %>%
  tibble::rownames_to_column("barcodes")

Spatial_RST2bei@meta.data <- Spatial_RST2bei@meta.data %>%
  tibble::rownames_to_column("barcodes") %>%
  dplyr::left_join(decon_df, by = "barcodes") %>%
  tibble::column_to_rownames("barcodes")

#Specific cell-types
#we can use the standard Seurat::SpatialFeaturePlot to view predicted celltype proportions one at a time.
#colnames(Spatial_RST2bei@meta.data)[15:17]<-paste0("State",1:3)

pdf("spotlight_cancer_CC2_12_FeaturePlot.pdf",width=7,height=7)
Seurat::SpatialFeaturePlot(
  object = Spatial_RST2bei,
  features = c("X1.z" ,"X2"),
  image.alpha=0.8,
   stroke=0.3,
  alpha = c(0.2, 1))
dev.off()

####2.2cluster
####2.1stage
cancer_scrna_stage2<-subset(cancer_scrna_stage2,)
stage2_markers <- Seurat::FindAllMarkers(object = cancer_scrna_stage2,
                                              assay = "RNA",
                                              slot = "data",
                                              verbose = TRUE,
                                              only.pos = TRUE)
saveRDS(object = stage2_markers, file = here::here("stage2_markers.RDS"))
stage2_markers <-readRDS("/share/pub/dengcy/Singlecell/CC_space/5.1.tumor_spatial_progress/stage2_markers.RDS")

cancer_scrna_stage2_13<-subset(cancer_scrna_stage2,State %in% c(1,3))
stage2_markers_13 <- Seurat::FindAllMarkers(object = cancer_scrna_stage2_13,
                                              assay = "RNA",
                                              slot = "data",
                                              verbose = TRUE,
                                              only.pos = TRUE)

cancer_scrna_stage2_12<-subset(cancer_scrna_stage2,State %in% c(1,2))
stage2_markers_13 <- Seurat::FindAllMarkers(object = cancer_scrna_stage2_12,
                                              assay = "RNA",
                                              slot = "data",
                                              verbose = TRUE,
                                              only.pos = TRUE)


###这一步有点慢,半个多小时。
spotlight_cancer_cc2_13 <- spotlight_deconvolution(
  se_sc = cancer_scrna_stage2_13,
  counts_spatial = Spatial_RST2bei@assays$Spatial@counts,
  clust_vr = "State", # Variable in sc_seu containing the cell-type annotation
  cluster_markers = stage2_markers_13, # Dataframe with the marker genes
  cl_n = 3000, # number of cells per cell type to use
  hvg = 5000, # Number of HVG to use
  ntop = NULL, # How many of the marker genes to use (by default all)
  transf = "uv", # Perform unit-variance scaling per cell and spot prior to factorzation and NLS
  method = "nsNMF", # Factorization method
  min_cont = 0 # Remove those cells contributing to a spot below a certain threshold
  )

saveRDS(object = spotlight_cancer_cc2_13, file = "/share/pub/dengcy/Singlecell/CC_space/5.1.tumor_spatial_progress/spotlight_cancer_cc2_13.RDS")
nmf_mod2 <- spotlight_cancer_cc2_13[[1]]
decon_mtrx2 <- spotlight_cancer_cc2_13[[2]]
#Assess deconvolution
h <- NMF::coef(nmf_mod2[[1]])
rownames(h) <- paste("Topic", 1:nrow(h), sep = "_")
topic_profile_plts <- dot_plot_profiles_fun(
  h = h,
  train_cell_clust = nmf_mod2[[2]])

topic_profile_plts[[2]] + ggplot2::theme(
  axis.text.x = ggplot2::element_text(angle = 90),
  axis.text = ggplot2::element_text(size = 12))

#Visualization
# This is the equivalent to setting min_cont to 0.04
decon_mtrx_sub <- decon_mtrx2[, colnames(decon_mtrx2) != "res_ss"]
decon_mtrx_sub[decon_mtrx_sub < 0.0001] <- 0
decon_mtrx2 <- cbind(decon_mtrx_sub, "res_ss" = decon_mtrx2[, "res_ss"])
rownames(decon_mtrx2) <- colnames(Spatial_RST2bei)

decon_df <- decon_mtrx2 %>%
  data.frame() %>%
  tibble::rownames_to_column("barcodes")

Spatial_RST2bei@meta.data <- Spatial_RST2bei@meta.data %>%
  tibble::rownames_to_column("barcodes") %>%
  dplyr::left_join(decon_df, by = "barcodes") %>%
  tibble::column_to_rownames("barcodes")

#Specific cell-types
#we can use the standard Seurat::SpatialFeaturePlot to view predicted celltype proportions one at a time.
#colnames(Spatial_RST2bei@meta.data)[15:17]<-paste0("State",1:3)

pdf("spotlight_cancer_CC2_13_FeaturePlot.pdf",width=7,height=7)
Seurat::SpatialFeaturePlot(
  object = Spatial_RST2bei,
  features = c("X1.y" ,"X3.y"),
  image.alpha=0.8,
   stroke=0.3,
  alpha = c(0.2, 1))
dev.off()

spotlight_cancer_cc2_12 <- spotlight_deconvolution(
  se_sc = cancer_scrna_stage2_12,
  counts_spatial = Spatial_RST2bei@assays$Spatial@counts,
  clust_vr = "State", # Variable in sc_seu containing the cell-type annotation
  cluster_markers = stage2_markers_13, # Dataframe with the marker genes
  cl_n = 3000, # number of cells per cell type to use
  hvg = 5000, # Number of HVG to use
  ntop = NULL, # How many of the marker genes to use (by default all)
  transf = "uv", # Perform unit-variance scaling per cell and spot prior to factorzation and NLS
  method = "nsNMF", # Factorization method
  min_cont = 0 # Remove those cells contributing to a spot below a certain threshold
  )

saveRDS(object = spotlight_cancer_cc2_12, file = "/share/pub/dengcy/Singlecell/CC_space/5.1.tumor_spatial_progress/spotlight_cancer_cc2_12.RDS")

nmf_mod2 <- spotlight_cancer_cc2_13[[1]]
decon_mtrx2 <- spotlight_cancer_cc2_13[[2]]
#Assess deconvolution
h <- NMF::coef(nmf_mod2[[1]])
rownames(h) <- paste("Topic", 1:nrow(h), sep = "_")
topic_profile_plts <- dot_plot_profiles_fun(
  h = h,
  train_cell_clust = nmf_mod2[[2]])

topic_profile_plts[[2]] + ggplot2::theme(
  axis.text.x = ggplot2::element_text(angle = 90),
  axis.text = ggplot2::element_text(size = 12))

#Visualization
# This is the equivalent to setting min_cont to 0.04
decon_mtrx_sub <- decon_mtrx2[, colnames(decon_mtrx2) != "res_ss"]
decon_mtrx_sub[decon_mtrx_sub < 0.0001] <- 0
decon_mtrx2 <- cbind(decon_mtrx_sub, "res_ss" = decon_mtrx2[, "res_ss"])
rownames(decon_mtrx2) <- colnames(Spatial_RST2bei)

decon_df <- decon_mtrx2 %>%
  data.frame() %>%
  tibble::rownames_to_column("barcodes")

Spatial_RST2bei@meta.data <- Spatial_RST2bei@meta.data %>%
  tibble::rownames_to_column("barcodes") %>%
  dplyr::left_join(decon_df, by = "barcodes") %>%
  tibble::column_to_rownames("barcodes")

#Specific cell-types
#we can use the standard Seurat::SpatialFeaturePlot to view predicted celltype proportions one at a time.
#colnames(Spatial_RST2bei@meta.data)[15:17]<-paste0("State",1:3)

pdf("spotlight_cancer_CC2_12_FeaturePlot.pdf",width=7,height=7)
Seurat::SpatialFeaturePlot(
  object = Spatial_RST2bei,
  features = c("X1.z" ,"X2"),
  image.alpha=0.8,
   stroke=0.3,
  alpha = c(0.2, 1))
dev.off()




#Spatial scatterpiesPlot spot composition of all the spots.
##cell_types_all <- colnames(decon_mtrx2)[which(colnames(decon_mtrx2) != "res_ss")]
pdf("SPOTlight_spatial_scatterpie_cancer_CC2.pdf")
spatial_scatterpie(se_obj = Spatial_RST2bei,
                              cell_types_all = c("State1" ,"State2", "State3"),
                              img_path = "/share/pub/dengcy/Singlecell/CC_space/gcy_analysis/Analysis_200860-01/RST2bei/2.Basic_analysis/2.3.h5_files/spatial/tissue_lowres_image.png",
                              pie_scale = 0.5)
dev.off()




