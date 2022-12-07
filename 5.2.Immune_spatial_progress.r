##################################
#免疫细胞的空间展示
############################

source("/share/pub/dengcy/Singlecell/CC_space/src/Functions.r")
setwd("/share/pub/dengcy/Singlecell/CC_space/5.2.Immune_spatial_progress")
source("/share/pub/dengcy/Singlecell/CC_space/src/library.r")

Spatial_RST2bei<-readRDS("/share/pub/dengcy/Singlecell/CC_space/5.spatial_progress/Spatial_RST2bei.rds")
#Spatial_RNST2<-readRDS("Spatial_RNST2.rds")
Spatial_CST1<-readRDS("/share/pub/dengcy/Singlecell/CC_space/5.spatial_progress/Spatial_CST1.rds")
#Spatial_CNST1<-readRDS("Spatial_CNST1.rds")

#load("/share/pub/dengcy/Singlecell/CC_space/3.2.Immune_cells_stage/Immune_scrna.RData")
load("/share/pub/dengcy/Singlecell/CC_space/3.2.Immune_cells_stage/Immune_markers2.0.RData")

##导入函数
files<-list.files("/share/pub/dengcy/Singlecell/CC_space/spotlight/R")
for(i in 1:length(files)){
	source(paste0("/share/pub/dengcy/Singlecell/CC_space/spotlight/R/",files[i]))
}
###########################
############CC1
###
load("/share/pub/dengcy/Singlecell/CC_space/6.cellcom/scrna_stage1.RData")
#Immune_scrna_cc1<- subset(Immune_scrna,samples == "CC1" )
#Immune_scrna_cc1<-subset(scrna_cc1,annotation != "NK.cell")

Immune_scrna1_T <- subset(scrna_stage1,annotation %in% c("Teff","Tex","Th","Treg"))
Immune_scrna1_M <- subset(scrna_stage1,annotation %in% c("B","Plasma","Neutrophils","TAM"))
rm(scrna_stage1)
################################################旧的代码
spotlight_Immune_magliant<-function(Spatial_data=Spatial_CST1,
  Immune_scrna=Immune_scrna_cc1,
  Immunecells="Cytotoxic.CD8.T",
  tumorcells=c("CT1","CT2")
  ){

  Immune_scrna1_cell <- subset(Immune_scrna,annotation %in% c(tumorcells,Immunecells))
cc1markers_cell <- Seurat::FindAllMarkers(object = Immune_scrna1_cell,
                                              assay = "RNA",
                                              slot = "data",
                                              verbose = TRUE,
                                              only.pos = TRUE)
#saveRDS(object = cc1markers_t, file = here::here("cc1markers_t.RDS"))
spotlight_cc1 <- spotlight_deconvolution(
  se_sc = Immune_scrna1_cell,
  counts_spatial = Spatial_data@assays$Spatial@counts,
  clust_vr = "annotation", # Variable in sc_seu containing the cell-type annotation
  cluster_markers = cc1markers_cell, # Dataframe with the marker genes
  cl_n = 100, # number of cells per cell type to use
  hvg = 6000, # Number of HVG to use
  ntop = NULL, # How many of the marker genes to use (by default all)
  transf = "uv", # Perform unit-variance scaling per cell and spot prior to factorzation and NLS
  method = "nsNMF", # Factorization method
  min_cont = 0 # Remove those cells contributing to a spot below a certain threshold
  )

#saveRDS(object = spotlight_immunect1_cc1, file = here::here("spotlight_immunect1_cc1.RDS"))

Spatial_data<-spotlight_fuc(Spatial_data=Spatial_data,spotlight_immune_cc1=spotlight_cc1)
return(Spatial_data)
}
##################

Spatial_data_cd8t<-spotlight_Immune_magliant(Spatial_data=Spatial_CST1,
  Immune_scrna=Immune_scrna_cc1,
  Immunecells=c("Cytotoxic.CD8.T"),
  tumorcells=c("LatePseudotime_TumerCells"))
save(Spatial_data_cd8t,file="Spatial_data_cd8t.RData")

pdf("scatterpie_cd8_CC1.pdf",width=8)
scatterpie_plot(se_obj= Spatial_data_cd8t,
                            cell_types_all=c("CT1","CT2","Cytotoxic.CD8.T"),
                            slice = NULL,
                            scatterpie_alpha = 1,
                            cell_types_interest = NULL,
                            pie_scale = 0.4,colors=NA,fillcolor=c("#D3ECA7","#A1B57D","#FC4F4F"))
dev.off()

Spatial_data_tumor<-spotlight_Immune_magliant(Spatial_data=Spatial_CST1,
  Immune_scrna=Immune_scrna_cc1,
  Immunecells=NULL,
  #fillcolor=c("#D3ECA7","#A1B57D","#FC4F4F","#2EB086","#11468F"),
  tumorcells=c("CT1","CT2"))
#save(Spatial_data_cd8t,file="Spatial_data_cd8t.RData")

pdf("scatterpie_tumor_CC1.pdf",width=8)
scatterpie_plot(se_obj= Spatial_data_tumor,
                            cell_types_all=c("CT1","CT2"),
                            slice = NULL,
                            scatterpie_alpha = 1,
                            cell_types_interest = NULL,
                            pie_scale = 0.4,colors=NA,fillcolor=color_scanpy_viridis28[c(1,4)])
dev.off()



Spatial_data_et<-spotlight_Immune_magliant(Spatial_data=Spatial_CST1,
  Immune_scrna=Immune_scrna_cc1,
  Immunecells=c("Terminally.Exhausted CD8T"),
  #fillcolor=c("#D3ECA7","#A1B57D","#FC4F4F","#2EB086","#11468F"),
  tumorcells=c("CT1","CT2"))
save(Spatial_data_et,file="Spatial_data_et.RData")

pdf("scatterpie_et_CC1.pdf",width=8)
scatterpie_plot(se_obj= Spatial_data_et,
                            cell_types_all=c("CT1","CT2","Terminally.Exhausted.CD8T"),
                            slice = NULL,
                            scatterpie_alpha = 1,
                            cell_types_interest = NULL,
                            pie_scale = 0.4,colors=NA,fillcolor=c("#D3ECA7","#A1B57D","#2EB086"))
dev.off()

Spatial_data_tam<-spotlight_Immune_magliant(Spatial_data=Spatial_CST1,
  Immune_scrna=Immune_scrna_cc1,
  Immunecells=c("TAM"),
  #fillcolor=c("#D3ECA7","#A1B57D","#FC4F4F","#2EB086","#11468F"),
  tumorcells=c("CT1","CT2"))
save(Spatial_data_tam,file="Spatial_data_tam.RData")


  #filenames="spotlight_immune_Ct1_FeaturePlot.pdf",
  filenames_pie="scatterpie_immune_ct1_CC1.pdf"

pdf("scatterpie_immune_CC1.pdf",width=8)
scatterpie_plot(se_obj= Spatial_data,
                            cell_types_all=c("CT1","CT2","Cytotoxic.CD8.T","Terminally.Exhausted.CD8T","TAM"),
                            slice = NULL,
                            scatterpie_alpha = 1,
                            cell_types_interest = NULL,
                            pie_scale = 0.4,colors=NA,fillcolor=c("#D3ECA7","#A1B57D","#FC4F4F","#2EB086","#11468F"))
dev.off()

######################################################################################2.20以前
cc1markers_t <- Seurat::FindAllMarkers(object = Immune_scrna1_T,
                                              assay = "RNA",
                                              slot = "data",
                                              verbose = TRUE,
                                              only.pos = TRUE)
saveRDS(object = cc1markers_t, file = here::here("cc1markers_t.RDS"))

cc1markers_2 <- Seurat::FindAllMarkers(object = Immune_scrna1_M,
                                              assay = "RNA",
                                              slot = "data",
                                              verbose = TRUE,
                                              only.pos = TRUE)
saveRDS(object = cc1markers_2, file = here::here("cc1markers_M.RDS"))

spotlight_immuneT_cc1 <- spotlight_deconvolution(
  se_sc = Immune_scrna1_T,
  counts_spatial = Spatial_CST1@assays$Spatial@counts,
  clust_vr = "annotation", # Variable in sc_seu containing the cell-type annotation
  cluster_markers = cc1markers_t, # Dataframe with the marker genes
  cl_n = 500, # number of cells per cell type to use
  hvg = 6000, # Number of HVG to use
  ntop = NULL, # How many of the marker genes to use (by default all)
  transf = "uv", # Perform unit-variance scaling per cell and spot prior to factorzation and NLS
  method = "nsNMF", # Factorization method
  min_cont = 0 # Remove those cells contributing to a spot below a certain threshold
  )

saveRDS(object = spotlight_immuneT_cc1, file = here::here("spotlight_immuneT_cc1.RDS"))

Spatial_CST1<-spotlight_fuc(Spatial_data=Spatial_CST1,spotlight_immune_cc1=spotlight_immuneT_cc1)

df<-Spatial_CST1@meta.data
p1 <- ggplot(df,aes(x =State4 ,
                y =Treg))+
          geom_point(size=0.5,color="#8785a2")+
          geom_smooth(se = F,
                      color = "#8E0505")+
          labs(x = "Pseudotime-State4",y="Treg",title="StageI/II  Treg")+
          theme_classic()+
          stat_cor(data =df,
           method = "pearson")

p2 <- ggplot(df,aes(x =State4 ,
                y =Tex))+
          geom_point(size=0.5)+
          geom_smooth(se = F,
                      color = "#8E0505")+
          labs(x = "Pseudotime-State4",y="Tex",title="StageI/II  Tex")+
          theme_classic()+
          stat_cor(data =df,
           method = "pearson")
p3 <- ggplot(df,aes(x =State4 ,
                y =TAM))+
          geom_point(size=0.5)+
          geom_smooth(se = F,
                      color = "#8E0505")+
          labs(x = "Pseudotime-State4",y="TAM",title="StageI/II  TAM")+
          theme_classic()+
          stat_cor(data =df,
           method = "pearson")
          p3.1 <- ggplot(df,aes(x =State5 ,
                y =TAM))+
          geom_point(size=0.5)+
          geom_smooth(se = F,
                      color = "#8E0505")+
          labs(x = "Pseudotime-State5",y="TAM",title="StageI/II  TAM")+
          theme_classic()+
          stat_cor(data =df,
           method = "pearson")
p4 <- ggplot(df,aes(x =State4 ,
                y =Neutrophils))+
          geom_point(size=0.5)+
          geom_smooth(se = F,
                      color = "#8E0505")+
          labs(x = "Pseudotime-State5",y="Neutrophils",title="StageI/II  Neutrophils")+
          theme_classic()+
          stat_cor(data =df,
           method = "pearson")
pdf("Immune_emt_score_stage1.pdf",width=10,height=9)
ggarrange(p1,p2,p3,p3.1,p4, ncol = 2)
dev.off()


pdf("spotlight_Tcell_CC1_FeaturePlot.pdf",width=7,height=7)
Seurat::SpatialFeaturePlot(
  object = Spatial_CST1,
  features = c("Teff","Tex","Th","Treg"),
  image.alpha=0.8,
   stroke=0.3,
  alpha = c(0.2, 1))
dev.off()
#Spatial scatterpiesPlot spot composition of all the spots.

###########B细胞和中性粒细胞

spotlight_immuneB_cc1 <- spotlight_deconvolution(
  se_sc = Immune_scrna1_M,
  counts_spatial = Spatial_CST1@assays$Spatial@counts,
  clust_vr = "annotation", # Variable in sc_seu containing the cell-type annotation
  cluster_markers = cc1markers_2, # Dataframe with the marker genes
  cl_n = 500, # number of cells per cell type to use
  hvg = 6000, # Number of HVG to use
  ntop = NULL, # How many of the marker genes to use (by default all)
  transf = "uv", # Perform unit-variance scaling per cell and spot prior to factorzation and NLS
  method = "nsNMF", # Factorization method
  min_cont = 0 # Remove those cells contributing to a spot below a certain threshold
  )

saveRDS(object = spotlight_immuneB_cc1, file = here::here("spotlight_immuneB_cc1.RDS"))

Spatial_CST1<-spotlight_fuc(Spatial_data=Spatial_CST1,spotlight_immune_cc1=spotlight_immuneB_cc1)

pdf("spotlight_Bcell_CC1_FeaturePlot.pdf",width=9,height=7)
Seurat::SpatialFeaturePlot(
  object = Spatial_CST1,
  features = unique(as.vector(Immune_scrna1_M$annotation)),
  image.alpha=0.8,
   stroke=0.3,
  alpha = c(0.2, 1))
dev.off()

###############################
#CC2
########################
load("/share/pub/dengcy/Singlecell/CC_space/6.cellcom/scrna_stage2.RData")
#Immune_scrna_cc2<- subset(Immune_scrna,samples == "CC2" )
#Immune_scrna_cc1<-subset(Immune_scrna_cc1,annotation != "NK.cell")
Immune_scrna2_T <- subset(scrna_stage2,annotation %in% c("Teff","Tex","Th","Treg"))
Immune_scrna2_2 <- subset(scrna_stage2,annotation %in% c("B","Plasma","Neutrophils","TAM"))


cc2markers_t <- Seurat::FindAllMarkers(object = Immune_scrna2_T,
                                              assay = "RNA",
                                              slot = "data",
                                              verbose = TRUE,
                                              only.pos = TRUE)
saveRDS(object = cc2markers_t, file = here::here("cc2markers_t.RDS"))

cc2markers_B <- Seurat::FindAllMarkers(object = Immune_scrna2_2,
                                              assay = "RNA",
                                              slot = "data",
                                              verbose = TRUE,
                                              only.pos = TRUE)
saveRDS(object = cc2markers_B, file = here::here("cc2markers_B.RDS"))

spotlight_immuneT_cc2 <- spotlight_deconvolution(
  se_sc = Immune_scrna2_T,
  counts_spatial = Spatial_RST2bei@assays$Spatial@counts,
  clust_vr = "annotation", # Variable in sc_seu containing the cell-type annotation
  cluster_markers = cc2markers_t, # Dataframe with the marker genes
  cl_n = 500, # number of cells per cell type to use
  hvg = 6000, # Number of HVG to use
  ntop = NULL, # How many of the marker genes to use (by default all)
  transf = "uv", # Perform unit-variance scaling per cell and spot prior to factorzation and NLS
  method = "nsNMF", # Factorization method
  min_cont = 0 # Remove those cells contributing to a spot below a certain threshold
  )

saveRDS(object = spotlight_immuneT_cc2, file = here::here("spotlight_immuneT_cc2.RDS"))

Spatial_RST2bei<-spotlight_fuc(Spatial_data=Spatial_RST2bei,spotlight_immune_cc1=spotlight_immuneT_cc2)

df<-Spatial_CST1@meta.data
p1 <- ggplot(df,aes(x =State4 ,
                y =Treg))+
          geom_point(size=0.5)+
          geom_smooth(se = F,
                      color = "#8E0505")+
          labs(x = "Pseudotime-State4",y="Treg",title="StageI/II  Treg")+
          theme_classic()+
          stat_cor(data =df,
           method = "pearson")
p2 <- ggplot(df,aes(x =State4 ,
                y =Tex))+
          geom_point(size=0.5)+
          geom_smooth(se = F,
                      color = "#8E0505")+
          labs(x = "Pseudotime-State4",y="Tex",title="StageI/II  Tex")+
          theme_classic()+
          stat_cor(data =df,
           method = "pearson")
p3 <- ggplot(df,aes(x =State5 ,
                y =TAM))+
          geom_point(size=0.5)+
          geom_smooth(se = F,
                      color = "#8E0505")+
          labs(x = "Pseudotime-State5",y="TAM",title="StageI/II  TAM")+
          theme_classic()+
          stat_cor(data =df,
           method = "pearson")
p4 <- ggplot(df,aes(x =State5 ,
                y =Neutrophils))+
          geom_point(size=0.5)+
          geom_smooth(se = F,
                      color = "#8E0505")+
          labs(x = "Pseudotime-State5",y="Neutrophils",title="StageI/II  Neutrophils")+
          theme_classic()+
          stat_cor(data =df,
           method = "pearson")
pdf("Immune_emt_score.pdf",width=10,height=10)
ggarrange(p1,p2,p3,p4, ncol = 2)
dev.off()


pdf("spotlight_Tcell_CC2_FeaturePlot.pdf",width=7,height=7)
Seurat::SpatialFeaturePlot(
  object = Spatial_RST2bei,
  features = c("Teff","Tex","Th","Treg"),
  image.alpha=0.8,
   stroke=0.3,
  alpha = c(0.2, 1))
dev.off()
#Spatial scatterpiesPlot spot composition of all the spots.

###########B细胞和中性粒细胞
spotlight_immuneB_cc2 <- spotlight_deconvolution(
  se_sc = Immune_scrna2_2,
  counts_spatial = Spatial_RST2bei@assays$Spatial@counts,
  clust_vr = "annotation", # Variable in sc_seu containing the cell-type annotation
  cluster_markers = cc2markers_B, # Dataframe with the marker genes
  cl_n = 500, # number of cells per cell type to use
  hvg = 6000, # Number of HVG to use
  ntop = NULL, # How many of the marker genes to use (by default all)
  transf = "uv", # Perform unit-variance scaling per cell and spot prior to factorzation and NLS
  method = "nsNMF", # Factorization method
  min_cont = 0 # Remove those cells contributing to a spot below a certain threshold
  )

saveRDS(object = spotlight_immuneB_cc2, file = here::here("spotlight_immuneB_cc2.RDS"))

Spatial_RST2bei<-spotlight_fuc(Spatial_data=Spatial_RST2bei,spotlight_immune_cc1=spotlight_immuneB_cc2)

pdf("spotlight_Bcell_CC2_FeaturePlot.pdf",width=7,height=7)
Seurat::SpatialFeaturePlot(
  object = Spatial_RST2bei,
  features = unique(as.vector(Immune_scrna2_2$annotation)),
  image.alpha=0.8,
   stroke=0.3,
  alpha = c(0.2, 1))
dev.off()

#######
spotlight_fuc<-function(Spatial_data=Spatial_RST2bei,spotlight_immune_cc1=spotlight_immuneT_cc1){
	nmf_mod1 <- spotlight_immune_cc1[[1]]
decon_mtrx1 <- spotlight_immune_cc1[[2]]
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
rownames(decon_mtrx1) <- colnames(Spatial_data)

decon_df <- decon_mtrx1 %>%
  data.frame() %>%
  tibble::rownames_to_column("barcodes")
#Spatial_data@meta.data$barcodes<-rownames(Spatial_data@meta.data)
#Spatial_data@meta.data<-Spatial_data@meta.data[,1:10]

Spatial_data@meta.data <- Spatial_data@meta.data %>%
  tibble::rownames_to_column("barcodes") %>%
  dplyr::left_join(decon_df, by = "barcodes") %>%
  tibble::column_to_rownames("barcodes")

  cell_types_all <- colnames(decon_mtrx1)[which(colnames(decon_mtrx1) != "res_ss")]

return(Spatial_data)
}

#Specific cell-types
#we can use the standard Seurat::SpatialFeaturePlot to view predicted celltype proportions one at a time.

