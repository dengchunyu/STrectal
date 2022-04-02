source("/share/pub/dengcy/Singlecell/CC_space/src/Functions.r")
setwd("/share/pub/dengcy/Singlecell/CC_space/5.1.tumor_spatial_progress")
source("/share/pub/dengcy/Singlecell/CC_space/src/library.r")
#############
#提取干细胞结果
############
load("/share/pub/dengcy/Singlecell/CC_space/1.scRNA_seq_prepare/CC_scrna.RData")
load("/share/pub/dengcy/Singlecell/CC_space/3.Epi_cells_stage/cancer_scrna_cc1.RData")
load("/share/pub/dengcy/Singlecell/CC_space/3.Epi_cells_stage/cancer_scrna_cc2.RData")

stem_scrna<-subset(CC_scrna,subset= annotation %in% c("Cancer_stem_cells","Epithelial_cells"))
stem_scrna_cc1<-subset(stem_scrna,subset= samples %in% "CC1")
stem_scrna_cc2<-subset(stem_scrna,subset= samples %in% "CC2")

 stem_scrna_cc1<-cluster_pca_umap(stem_scrna_cc1, assay = "SCT",reduction="pca",cluster_res = 0.3)
 stem_scrna_cc2<-cluster_pca_umap(stem_scrna_cc2, assay = "SCT",reduction="pca",cluster_res = 0.3)

Idents(stem_scrna_cc1)<-stem_scrna_cc1$annotation
Idents(stem_scrna_cc2)<-stem_scrna_cc2$annotation

save(stem_scrna_cc1,file="stem_scrna_cc1.RData")
save(stem_scrna_cc2,file="stem_scrna_cc2.RData")

stemcc1_markers <- Seurat::FindAllMarkers(object = stem_scrna_cc1,
                                          #features ="seurat_clusters",
                                              assay = "SCT",
                                              slot = "data",
                                              verbose = TRUE,
                                              only.pos = TRUE)
stemcc2_markers <- Seurat::FindAllMarkers(object = stem_scrna_cc2,
                                              assay = "SCT",
                                              slot = "data",
                                              verbose = TRUE,
                                              only.pos = TRUE)

####################
#1.分别展示不同样本的cluster
#######################
pdf(file="Malignant_cc1_umap_0.3.pdf",height=6,width=6)
DimPlot(cancer_scrna_cc1,pt.size=0.9,label = TRUE, repel=TRUE,label.size = 6)+ 
umap_theme()+ labs(x="UMAP",y="")+
ggtitle("CC1 Malignant+stem cells")+
        scale_colour_manual(name = "Cluster", values = color_scanpy_viridis28[c(1,4,19)]) +
        theme(aspect.ratio=1)+
        geom_path(data=embedded1.1, aes(x=UMAP_1, y=UMAP_2), size=0.7,color="black",alpha=1)+
        geom_path(data=embedded1.2, aes(x=UMAP_1, y=UMAP_2), size=0.7,color="black",alpha=1)

dev.off()

#Idents(cancer_scrna_cc2) <- factor(cancer_scrna_cc2$annotation,levels=c("CT3","CT4","CT5","CT6"))
pdf(file="Malignant_cc2_umap_0.3.pdf",height=6,width=6)
DimPlot(cancer_scrna_cc2,pt.size=0.9,label = TRUE, repel=TRUE,label.size = 6)+ 
umap_theme()+ labs(x="UMAP",y="")+
ggtitle("CC2 Malignant cells")+
        scale_colour_manual(name = "Cluster", values = color_scanpy_viridis28[c(20,19,7,11)]) +
        theme(aspect.ratio=1)+
      geom_path(data=embedded2.1, aes(x=UMAP_1, y=UMAP_2), size=0.7,color="black",alpha=1)+
       geom_path(data=embedded2.2, aes(x=UMAP_1, y=UMAP_2), size=0.7,color="black",alpha=1)

dev.off()
####################stem
pdf(file="stem_cc1_umap_0.3.pdf",height=6,width=6)
DimPlot(stem_scrna_cc1,pt.size=0.9,label = TRUE, repel=TRUE,label.size = 6)+ 
umap_theme()+ labs(x="UMAP",y="")+
ggtitle("CC1 stem cells")+
        scale_colour_manual(name = "Cluster", values = color_scanpy_viridis28[c(18,19)]) +
        theme(aspect.ratio=1)
dev.off()

#Idents(cancer_scrna_cc2) <- factor(cancer_scrna_cc2$annotation,levels=c("CT3","CT4","CT5","CT6"))
pdf(file="stem_cc2_umap_0.3.pdf",height=6,width=6)
DimPlot(stem_scrna_cc2,pt.size=0.9,label = TRUE, repel=TRUE,label.size = 6)+ 
umap_theme()+ labs(x="UMAP",y="")+
ggtitle("CC2 stem cells")+
        scale_colour_manual(name = "Cluster", values = color_scanpy_viridis28[c(18,19)]) +
        theme(aspect.ratio=1)

dev.off()
###################
#将肿瘤细胞的各个cluster映射到不同的空间位置
################
Spatial_RST2bei<-readRDS("/share/pub/dengcy/Singlecell/CC_space/5.spatial_progress/Spatial_RST2bei.rds")
#Spatial_RNST2<-readRDS("Spatial_RNST2.rds")
Spatial_CST1<-readRDS("/share/pub/dengcy/Singlecell/CC_space/5.spatial_progress/Spatial_CST1.rds")
#Spatial_CNST1<-readRDS("Spatial_CNST1.rds")

#load("/share/pub/dengcy/Singlecell/CC_space/3.Epi_cells_stage/Slingshot_scrna1.RData")
#load("/share/pub/dengcy/Singlecell/CC_space/3.Epi_cells_stage/Slingshot_scrna2.RData")
#CC1_fortify_can <- fortify.Seurat(cancer_stem_scrna_cc1)
#CC2_fortify_can <- fortify.Seurat(cancer_stem_scrna_cc1)
####################
#空间图像的聚类映射位置
###################
##导入函数
files<-list.files("/share/pub/dengcy/Singlecell/CC_space/spotlight/R")
for(i in 1:length(files)){
	source(paste0("/share/pub/dengcy/Singlecell/CC_space/spotlight/R/",files[i]))
}
############CC1
cc1_markers <- Seurat::FindAllMarkers(object = cancer_scrna_cc1,
                                              assay = "SCT",
                                              slot = "data",
                                              verbose = TRUE,
                                              only.pos = TRUE)
saveRDS(object = cc1_markers, file = here::here("cc1_markers.RDS"))
###这一步有点慢,半个多小时。
spotlight_cancer_cc1 <- spotlight_deconvolution(
  se_sc = cancer_scrna_cc1,
  counts_spatial = Spatial_CST1@assays$Spatial@counts,
  clust_vr = "annotation", # Variable in sc_seu containing the cell-type annotation
  cluster_markers = cc1_markers, # Dataframe with the marker genes
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
pdf("spotlight_cancer_CC1_FeaturePlot.pdf",width=7,height=7)
Seurat::SpatialFeaturePlot(
  object = Spatial_CST1,
  features = c("CT1","CT2"),
  image.alpha=0.8,
   stroke=0.3,
  alpha = c(0.2, 1))+
 scale_colour_manual(name = "CC1", values =c("#F4E185","#F4E185","#CD1818"))
dev.off()

#Spatial scatterpiesPlot spot composition of all the spots.
cell_types_all <- colnames(decon_mtrx1)[which(colnames(decon_mtrx1) != "res_ss")]
pdf("SPOTlight_spatial_scatterpie_cancer_CC1.pdf")
spatial_scatterpie(se_obj = Spatial_CST1,
                              cell_types_all = cell_types_all,
                              img_path = "/share/pub/dengcy/Singlecell/CC_space/gcy_analysis/Analysis_200860-01/CST1/2.Basic_analysis/2.3.h5_files/spatial/tissue_lowres_image.png",
                              pie_scale = 0.4)
dev.off()
########################################
##############################stem

#saveRDS(object = stemcc1_markers, file = here::here("cc1_markers.RDS"))
###这一步有点慢,半个多小时。
spotlight_stem_cc1 <- spotlight_deconvolution(
  se_sc = stem_scrna_cc1,
  counts_spatial = Spatial_CST1@assays$Spatial@counts,
  clust_vr = "annotation", # Variable in sc_seu containing the cell-type annotation
  cluster_markers = stemcc1_markers, # Dataframe with the marker genes
  cl_n = 100, # number of cells per cell type to use
  hvg = 1500, # Number of HVG to use
  ntop = NULL, # How many of the marker genes to use (by default all)
  transf = "uv", # Perform unit-variance scaling per cell and spot prior to factorzation and NLS
  method = "nsNMF", # Factorization method
  min_cont = 0 # Remove those cells contributing to a spot below a certain threshold
  )

saveRDS(object = spotlight_stem_cc1, file = here::here("spotlight_stem_cc1.RDS"))

nmf_mod1 <- spotlight_stem_cc1[[1]]
decon_mtrx1 <- spotlight_stem_cc1[[2]]
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
decon_mtrx_sub[decon_mtrx_sub < 0.000001] <- 0
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
pdf("spotlight_stem_CC1_FeaturePlot.pdf",width=7,height=7)
Seurat::SpatialFeaturePlot(
  object = Spatial_CST1,
  features = c("Epithelial_cells"),#"Cancer_stem_cells",
  image.alpha=0.8,
   stroke=0.3,
  alpha = c(0.2, 1))+
 scale_colour_manual(name = "CC1", values =c("#F4E185","#F4E185"))#,"#CD1818"
dev.off()

#######################CC2
cc2_markers <- Seurat::FindAllMarkers(object = cancer_scrna_cc2,
                                              assay = "SCT",
                                              slot = "data",
                                              verbose = TRUE,
                                              only.pos = TRUE)
saveRDS(object = cc2_markers, file = here::here("cc2_markers.RDS"))
###这一步有点慢,半个多小时。
spotlight_cancer_cc2 <- spotlight_deconvolution(
  se_sc = cancer_scrna_cc2,
  counts_spatial = Spatial_RST2bei@assays$Spatial@counts,
  clust_vr = "annotation", # Variable in sc_seu containing the cell-type annotation
  cluster_markers = cc2_markers, # Dataframe with the marker genes
  cl_n = 100, # number of cells per cell type to use
  hvg = 2000, # Number of HVG to use
  ntop = NULL, # How many of the marker genes to use (by default all)
  transf = "uv", # Perform unit-variance scaling per cell and spot prior to factorzation and NLS
  method = "nsNMF", # Factorization method
  min_cont = 0 # Remove those cells contributing to a spot below a certain threshold
  )

saveRDS(object = spotlight_cancer_cc2, file = here::here("spotlight_cancer_cc2.RDS"))

nmf_mod2 <- spotlight_cancer_cc2[[1]]
decon_mtrx2 <- spotlight_cancer_cc2[[2]]
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
pdf("spotlight_cancer_CC2_FeaturePlot.pdf",width=7,height=7)
Seurat::SpatialFeaturePlot(
  object = Spatial_RST2bei,
  features = c("CT4","CT5","CT6"),
  image.alpha=0.8,
   stroke=0.3,
  alpha = c(0.2, 1))
dev.off()

#Spatial scatterpiesPlot spot composition of all the spots.
cell_types_all <- colnames(decon_mtrx2)[which(colnames(decon_mtrx2) != "res_ss")]
pdf("SPOTlight_spatial_scatterpie_cancer_CC2.pdf")
spatial_scatterpie(se_obj = Spatial_RST2bei,
                              cell_types_all = c("CT4","CT5","CT6"),
                              img_path = "/share/pub/dengcy/Singlecell/CC_space/gcy_analysis/Analysis_200860-01/RST2bei/2.Basic_analysis/2.3.h5_files/spatial/tissue_lowres_image.png",
                              pie_scale = 0.5)
dev.off()




