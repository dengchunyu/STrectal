
library("Seurat")
library("ggplot2")
library("patchwork")
library("dplyr")
library("hdf5r")
source("/share/pub/dengcy/Singlecell/CC_space/src/Functions.r")
setwd("/share/pub/dengcy/Singlecell/CC_space/5.spatial_progress")
#source("/share/pub/dengcy/Singlecell/CC_space/src/library.r")
#tar czvf CRC_RT1.tar.gz /share/pub/dengcy/Singlecell/CC_space/gcy_analysis/Analysis_200860-01/CST1/2.Basic_analysis/2.2.filtered_feature_bc_matrix

#tar czvf CRC_RT2.tar.gz /share/pub/dengcy/Singlecell/CC_space/gcy_analysis/Analysis_200860-01/RST2bei/2.Basic_analysis/2.2.filtered_feature_bc_matrix

######RST2bei
data_dir <- "/share/pub/dengcy/Singlecell/CC_space/gcy_analysis/Analysis_200860-01/RST2bei/2.Basic_analysis/2.3.h5_files"
file_name <- "raw_feature_bc_matrix.h5"
Spatial_RST2bei <- Load10X_Spatial(data.dir=data_dir,filename=file_name,slice="RST2bei")
Spatial_RST2bei <- normalize_spacial(Spatial_data=Spatial_RST2bei,lowqspot=0.02,mitper=25,geneExprMin=10,spot_meta="RST2bei")
saveRDS(Spatial_RST2bei,file="Spatial_RST2bei.rds")
#get the VariableFeatures 20
#top20<- head(VariableFeatures(Spatial_RST2bei),20)
######RNST2
data_dir <- "/share/pub/dengcy/Singlecell/CC_space/gcy_analysis/Analysis_200860-01/RNST2/2.Basic_analysis/2.3.h5_files"
file_name <- "raw_feature_bc_matrix.h5"
Spatial_RNST2<-Load10X_Spatial(data.dir=data_dir,filename=file_name,slice="RNST2")
Spatial_RNST2 <- normalize_spacial(Spatial_data=Spatial_RNST2,lowqspot=0.02,mitper=25,geneExprMin=10,spot_meta="RNST2")
saveRDS(Spatial_RNST2,file="Spatial_RNST2.rds")
######CST1
data_dir <- "/share/pub/dengcy/Singlecell/CC_space/gcy_analysis/Analysis_200860-01/CST1/2.Basic_analysis/2.3.h5_files"
file_name <- "raw_feature_bc_matrix.h5"
Spatial_CST1<-Load10X_Spatial(data.dir=data_dir,filename=file_name,slice="CST1")
Spatial_CST1 <- normalize_spacial(Spatial_data=Spatial_CST1,lowqspot=0.02,mitper=25,geneExprMin=10,spot_meta="CST1")
saveRDS(Spatial_CST1,file="Spatial_CST1.rds")
######CNST1
data_dir <- "/share/pub/dengcy/Singlecell/CC_space/gcy_analysis/Analysis_200860-01/CNST1/2.Basic_analysis/2.3.h5_files"
file_name <- "raw_feature_bc_matrix.h5"
Spatial_CNST1<-Load10X_Spatial(data.dir=data_dir,filename=file_name,slice="CNST1")
Spatial_CNST1 <- normalize_spacial(Spatial_data=Spatial_CNST1,lowqspot=0.02,mitper=25,geneExprMin=10,spot_meta="CNST1")
saveRDS(Spatial_CNST1,file="Spatial_CNST1.rds")
#################
#2.cluster
############################################
Spatial_RST2bei<-readRDS("Spatial_RST2bei.rds")
#Spatial_RNST2<-readRDS("Spatial_RNST2.rds")
Spatial_CST1<-readRDS("Spatial_CST1.rds")
#Spatial_CNST1<-readRDS("Spatial_CNST1.rds")

Spatial_RST2bei<-cluster_pca_umap(Spatial_RST2bei, assay = "SCT",reduction="pca",cluster_res = 0.2)
#plot UMAP
#p1 <- DimPlot(Spatial_RST2bei, reduction = "umap")
#ggsave("Spatial_RST2bei-firstCluster-r0.2-UMAP.pdf", plot=p1, width = 6, height =5)

#Spatial_RNST2<-cluster_pca_umap(Spatial_RNST2, assay = "SCT",reduction="pca",cluster_res = 0.3)
Spatial_CST1<-cluster_pca_umap(Spatial_CST1, assay = "SCT",reduction="pca",cluster_res = 0.2)
#Spatial_CNST1<-cluster_pca_umap(Spatial_CNST1, assay = "SCT",reduction="pca",cluster_res = 0.3)

Spatial_RST2bei@meta.data$samples="CC2"
write.csv(Spatial_RST2bei@meta.data,file="SupplementaryTable13.csv")
Spatial_CST1@meta.data$samples="CC1"
write.csv(Spatial_CST1@meta.data,file="SupplementaryTable14.csv")

saveRDS(Spatial_RST2bei,file="Spatial_RST2bei.rds")
#saveRDS(Spatial_RNST2,file="Spatial_RNST2.rds")
saveRDS(Spatial_CST1,file="Spatial_CST1.rds")
#saveRDS(Spatial_CNST1,file="Spatial_CNST1.rds")
#######
Spatial_RST2bei<-readRDS("Spatial_RST2bei.rds")
#Spatial_RNST2<-readRDS("Spatial_RNST2.rds")
Spatial_CST1<-readRDS("Spatial_CST1.rds")
#Spatial_CNST1<-readRDS("Spatial_CNST1.rds")

color_scanpy_patient <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22","#D9DD6B","#A45D5D","#D54C4C","#8D2828")

plot1 <- VlnPlot(Spatial_RST2bei, features = "nCount_Spatial", pt.size = 0.3) + NoLegend()
plot2 <- SpatialFeaturePlot(Spatial_RST2bei, features = "nCount_Spatial",pt.size.factor = 2.5) + theme(legend.position = "right")
pdf("Spatial_RST2bei_count.pdf",width=12)
wrap_plots(plot1, plot2)
dev.off()

plot1 <- VlnPlot(Spatial_CST1, features = "nCount_Spatial", pt.size = 0.3) + NoLegend()
plot2 <- SpatialFeaturePlot(Spatial_CST1, features = "nCount_Spatial",pt.size.factor = 2.5) + theme(legend.position = "right")
pdf("Spatial_CST1_count.pdf",width=12)
wrap_plots(plot1, plot2)
dev.off()

#
plot1 <- VlnPlot(Spatial_CNST1, features = "nCount_Spatial", pt.size = 0.3) + NoLegend()
plot2 <- SpatialFeaturePlot(Spatial_CNST1, features = "nCount_Spatial",pt.size.factor = 2.5) + theme(legend.position = "right")
pdf("Spatial_CNST1_count.pdf",width=12)
wrap_plots(plot1, plot2)
dev.off()
#
p1<-DimPlot(Spatial_RST2bei, reduction = "umap", label = TRUE)+ scale_colour_manual(name = "Cluster", values = color_scanpy_patient)+
SpatialDimPlot(Spatial_RST2bei, label = TRUE, label.size = 5,cols =color_scanpy_patient,pt.size.factor = 2)

#p2 <- DimPlot(Spatial_RNST2, reduction = "umap", label = TRUE)+scale_colour_manual(name = "Cluster", values = color_scanpy_patient)+
#SpatialDimPlot(Spatial_RNST2, label = TRUE, label.size = 5,cols =color_scanpy_patient,pt.size.factor = 2)

p3 <- DimPlot(Spatial_CST1, reduction = "umap", label = TRUE)+scale_colour_manual(name = "Cluster", values = color_scanpy_patient)+
SpatialDimPlot(Spatial_CST1, label = TRUE, label.size = 5,cols =color_scanpy_patient,pt.size.factor = 2)

#p4 <- DimPlot(Spatial_CNST1, reduction = "umap", label = TRUE)+scale_colour_manual(name = "Cluster", values = color_scanpy_patient)+
#SpatialDimPlot(Spatial_CNST1, label = TRUE, label.size = 5,cols =color_scanpy_patient,pt.size.factor = 2)

pdf("dimplot_RST2bei.pdf",width=12)
p1
dev.off()

pdf("dimplot_CST1.pdf",width=12)
p3
dev.off()


pdf("Spatial_RST2bei.pdf")
SpatialDimPlot(Spatial_RST2bei, label = F,pt.size.factor = 0)
dev.off()

pdf("Spatial_CST.pdf")
SpatialDimPlot(Spatial_CST1, label = F,pt.size.factor = 0)
dev.off()

######################
##get marker genes
plan("multiprocess", workers = N_WORKERS)

  Spatial_RST2bei_Allmarkers <- parallelFindAllMarkers(Spatial_RST2bei)
  save(Spatial_RST2bei_Allmarkers,file="Spatial_RST2bei_Allmarkers.RData")

  Spatial_RNST2_Allmarkers <- parallelFindAllMarkers(Spatial_RNST2)
  save(Spatial_RNST2_Allmarkers,file="Spatial_RNST2_Allmarkers.RData")

  Spatial_CST1_Allmarkers <- parallelFindAllMarkers(Spatial_CST1)
  save(Spatial_CST1_Allmarkers,file="Spatial_CST1_Allmarkers.RData")

  Spatial_CNST1_Allmarkers <- parallelFindAllMarkers(Spatial_CNST1)
  save(Spatial_CNST1_Allmarkers,file="Spatial_CNST1_Allmarkers.RData")

#get high VariableFeatures
RST2beitop20<- head(VariableFeatures(Spatial_RST2bei),20)
RNST2top20<- head(VariableFeatures(Spatial_RNST2),20)
CST1top20<- head(VariableFeatures(Spatial_CST1),20)
CNST1top20<- head(VariableFeatures(Spatial_CNST1),20)

Spatial_RST2bei_avrage <- AverageExpression(Spatial_RST2bei, return.seurat = TRUE, group.by= "seurat_clusters")
CellScatter(cluster.averages, cell1 = "CD8_T_rep1", cell2 = "CD8_T_rep2")
Idents(Spatial_RST2bei_avrage)<-c(0,1,2,3,4,5)
pdf("test2.pdf")
CellScatter(Spatial_RST2bei_avrage, cell1 = "1", cell2 = "0")
dev.off()
pdf("Spatial_RST2bei_avrage_top20.pdf")
DoHeatmap(Spatial_RST2bei_avrage , features = RST2beitop20)
dev.off()

###########################
#emt invation score
#######################
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
install.packages("hrbrthemes", repo=site)
library("ggsignif")
library(hrbrthemes)

load("/share/pub/dengcy/Singlecell/CC_space/4.stage_visualize/emt_gene.RData")

Spatial_RST2bei<-AddModuleScore(Spatial_RST2bei,list(stage_upgenes),name=c("Invasion_score"))
Spatial_CST1<-AddModuleScore(Spatial_CST1,list(stage_upgenes),name=c("Invasion_score"))

saveRDS(Spatial_RST2bei,file="Spatial_RST2bei.rds")
saveRDS(Spatial_CST1,file="Spatial_CST1.rds")

write.csv(Spatial_RST2bei@meta.data,file="SupplementaryData2.csv")
write.csv(Spatial_CST1@meta.data,file="SupplementaryData3.csv")

 #VlnPlot(Spatial_RST2bei,features ="Invasion_up_score1",sort="increasing")+ggtitle("Invasion_up_score")+

ggRST2bei<-Spatial_RST2bei@meta.data
ggRST2bei<-ggRST2bei[which(ggRST2bei$seurat_clusters %in% c("0","1","2","3","4")),]
ggRST2bei<-na.omit(ggRST2bei)
ggRST2bei$seurat_clusters<-as.vector(ggRST2bei$seurat_clusters)
ggRST2bei$seurat_clusters<-factor(ggRST2bei$seurat_clusters, levels =c("2","3","0","1","4"))

#sort(tapply(ggRST2bei$Invasion_up_score1,ggRST2bei$seurat_clusters,mean))

  plots_RST2bei <- ggplot(ggRST2bei,aes(x=seurat_clusters,y=Invasion_score1,fill=seurat_clusters,color=seurat_clusters))+
      geom_boxplot(outlier.shape = NA,alpha=0.6)+
      stat_compare_means(method ="wilcox.test",
    	comparisons=list(c("3","1"),c("3","4"),c("2","4"),c('2','1')),
                       hide.ns=T,label="p.signif",
                       color="black",size=4,
                       symnum.args= list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                         symbols = c("****", "***", "**", "*", "ns")))+
    
    #theme_classic()+
      #  geom_jitter(color="black" ,size=0.3,alpha=0.5)+theme(legend.position="none")+ #不需要图例
    theme_ipsum_rc(base_family = "" )+
    scale_fill_manual(name = "cluster", values =color_scanpy_patient[c(3,4,1,2,5)])+
    scale_color_manual(name = "", values =color_scanpy_patient[c(3,4,1,2,5)])

    #theme(legend.position="none",aspect.ratio=1)
pdf("Spatial_CC2_Invasion_score_vinplot2.pdf",width=7,height=5)
print(plots_RST2bei)
dev.off()

#####
ggCST1<-Spatial_CST1@meta.data
ggCST1<-ggCST1[which(ggCST1$seurat_clusters %in% c("0","1","2")),]
ggCST1<-na.omit(ggCST1)
ggCST1$seurat_clusters<-as.vector(ggCST1$seurat_clusters)
ggCST1$seurat_clusters<-factor(ggCST1$seurat_clusters, levels =c("2","0","1"))

  plots_ggCST1 <- ggplot(ggCST1,aes(x=seurat_clusters,y=Invasion_score1,fill=seurat_clusters,color=seurat_clusters))+
      geom_boxplot(outlier.shape = NA,alpha=0.6)+
      stat_compare_means(method ="wilcox.test",
        comparisons=list(c('2','0'),c('1','2'),c('1','0')),
                       hide.ns=T,label="p.signif",
                       color="black",size=4,
                       symnum.args= list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                         symbols = c("****", "***", "**", "*", "ns")))+
    theme_ipsum_rc(base_family = "" )+
    scale_fill_manual(name = "cluster", values =color_scanpy_patient[c(3,1,2)])+
    scale_color_manual(name = "", values =color_scanpy_patient[c(3,1,2)])

pdf("Spatial_CC1_Invasion_score_vinplot2.pdf",width=5,height=5)
print(plots_ggCST1)
dev.off()

Spatial_CST1$Invasion_score<-Spatial_CST1$Invasion_score1
Spatial_CST1$Invasion_score[Spatial_CST1$Invasion_score1>1.5]<-1.5
p_1.1 <-SpatialFeaturePlot(Spatial_CST1, features="Invasion_score")
# p_1.2 <-SpatialFeaturePlot(Spatial_CST1, features="Invasion_down_score2")
pdf("Spatial_CC1_Invasion_score_dimplot.pdf",width=6)
p_1.1
dev.off()


p_2.1 <-SpatialFeaturePlot(Spatial_RST2bei, features="Invasion_score1")
# p_2.2 <-SpatialFeaturePlot(Spatial_RST2bei, features="Invasion_down_score2")
pdf("Spatial_CC2_Invasion_score_dimplot.pdf",width=12)
p_2.1
dev.off()

##########EMT
load("/share/pub/dengcy/Singlecell/CC_space/4.stage_visualize/genes.by.pathway.h.RData")

Spatial_RST2bei<-AddModuleScore(Spatial_RST2bei,list(c(genes.by.pathway.h[["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]])),name="EMT")
Spatial_CST1<-AddModuleScore(Spatial_CST1,list(c(genes.by.pathway.h[["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]])),name="EMT")

  plots <- VlnPlot(Spatial_RST2bei,features ="EMT1",sort="increasing")+ggtitle("EMT")+
    stat_compare_means(method ="wilcox.test",
                       hide.ns=T,label="p.signif",
                       color="red",size=4,
                       symnum.args= list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                         symbols = c("****", "***", "**", "*", "ns")))+
    scale_fill_manual(name = "T", values =color_scanpy_patient) +
    theme(legend.position="none",aspect.ratio=1)
pdf("Spatial_RST2bei_emt_vinplot.pdf")
print(plots)
dev.off()

  plots <- VlnPlot(Spatial_CST1,features ="EMT1",sort="increasing")+ggtitle("EMT")+
    stat_compare_means(method ="wilcox.test",
                       hide.ns=T,label="p.signif",
                       color="red",size=4,
                       symnum.args= list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                         symbols = c("****", "***", "**", "*", "ns")))+
    scale_fill_manual(name = "T", values =color_scanpy_patient) +
    theme(legend.position="none",aspect.ratio=1)
pdf("Spatial_CST1_emt_vinplot.pdf")
print(plots)
dev.off()


######################################
#######################################################
pred.brca_RST2bei <- singler_func(Spatial_RST2bei,ref)
pred.brca_RNST2 <- singler_func(Spatial_RNST2,ref)
pred.brca_CST1 <- singler_func(Spatial_CST1,ref)
pred.brca_CNST1 <- singler_func(Spatial_CNST1,ref)

pred_final_RST2bei<-pred_final(pred.brca=pred.brca_RST2bei,CC_scrna=Spatial_RST2bei)
pred_final_RNST2<-pred_final(pred.brca=pred.brca_RNST2,CC_scrna=Spatial_RNST2)
pred_final_CST1<-pred_final(pred.brca=pred.brca_CST1,CC_scrna=Spatial_CST1)
pred_final_CNST1<-pred_final(pred.brca=pred.brca_CNST1,CC_scrna=Spatial_CNST1)

Spatial_RST2bei$singler <- pred.brca_RST2bei$labels
Spatial_RNST2$singler <- pred.brca_RNST2$labels
Spatial_CST1$singler <- pred.brca_CST1$labels
Spatial_CNST1$singler <- pred.brca_CNST1$labels
#######################################################

saveRDS(Spatial_RST2bei,file="Spatial_RST2bei.rds")
#saveRDS(Spatial_RNST2,file="Spatial_RNST2.rds")
saveRDS(Spatial_CST1,file="Spatial_CST1.rds")
#saveRDS(Spatial_CNST1,file="Spatial_CNST1.rds")
#
#################################################
##stage
########################################

#load("/share/pub/dengcy/Singlecell/CC_space/4.stage_visualize/cor_p.RData")
#Spatial_merged <- merge(Spatial_RST2bei, y = Spatial_CST1, project = "merged", merge.data = TRUE)

Spatial_RST2bei$Invasion_score<-Spatial_RST2bei$Invasion_score1
plots1 <- SpatialFeaturePlot(Spatial_RST2bei, features = "Invasion_score",pt.size.factor = 2)+
scale_fill_gradient2(limits=range(Spatial_RST2bei$Invasion_score),low="#2C2891",mid="white",high="#E02401",midpoint= median(Spatial_RST2bei$Invasion_score))+ggtitle("StageIII")

plots3 <- SpatialFeaturePlot(Spatial_CST1, features = "Invasion_score",pt.size.factor = 2)+
scale_fill_gradient2(limits=range(Spatial_CST1$Invasion_score),low="#2C2891",mid="white",high="#E02401",midpoint= median(Spatial_CST1$Invasion_score))+ggtitle("StageI/II")

pdf("Spatial_Invasion_score.pdf",width=15,height=15)
wrap_plots(plots1,plots3)
dev.off()
