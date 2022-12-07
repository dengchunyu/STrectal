###################################################
#step3:Plot heatmap of the PG clusters' mean/median marker expressions
#################################################
library(testSctpa)
library(Seurat)
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(rtracklayer)
library(sctransform)
library(clustree)
library(harmony)
library(cowplot)
library(stringr)
options(future.globals.maxSize = 10000 * 1024^2)

source("/share/pub/dengcy/Singlecell/CC_space/src/Functions.r")
source("/share/pub/dengcy/Singlecell/CC_space/src/library.r")
setwd("/share/pub/dengcy/Singlecell/CC_space/2.MarkerHeatmap")
load("/share/pub/dengcy/Singlecell/CC_space/1.scRNA_seq_prepare/CRC_merged_all.RData")
load("/share/pub/dengcy/Singlecell/CC_space/1.scRNA_seq_prepare/CRC.markers.RData")
#load("/share/pub/dengcy/Singlecell/CC_space/1.scRNA_seq_prepare/Allmarkers.RData")

Idents(CRC_merged_all)<-CRC_merged_all$Initial_annotation
a1<-unlist(lapply(levels(Idents(CRC_merged_all))[c(1:5,8)],function(x){
    a<-colnames(CRC_merged_all)[Idents(CRC_merged_all)==x]
    a<-sample(a,3000)
    return(a)
  }))
a2<-unlist(lapply(levels(Idents(CRC_merged_all))[c(6,7,9)],function(x){
    a<-colnames(CRC_merged_all)[Idents(CRC_merged_all)==x]
    a<-sample(a,300)
    return(a)
  }))
CRC_subset<-CRC_merged_all[,c(a1,a2)]
################################################
#计算不同细胞类型的差异marker 
markers_annotation = FindAllMarkers(CRC_merged_all,logfc.threshold = 0)
save(markers_annotation,file="markers_annotation.RData")
load("markers_annotation.RData")

color_CRC<-c("#7FBCD2","#A5F1E9","#FFEEAF","#7FB77E","#F7F6DC","#B1D7B4","#FFC090",
  "#F3E0B5","#FFAE6D","#FECD70","#E3C770","#94B49F","#CEE5D0","#ECB390","#F2D7D9",
  "#D3CEDF","#B2C8DF","#C4D7E0","#F8F9D7","#76BA99","#ADCF9F",
  "#CED89E","#FFDCAE","#FFE3A9","#FFC3C3","#FF8C8C","#97C4B8","#F1F0C0","#B1BCE6")

markers = markers_annotation %>% group_by(cluster) %>% top_n(n=5,wt=avg_log2FC)  #For Seurat>4.0, set "wt=avg_logFC"
Idents(CRC_merged_all)<-CRC_merged_all$Initial_annotation
pdf("step3_DoHeatmap_marker_for_annotation.pdf",width=7,height=8)
DoHeatmap(CRC_merged_all,features=markers$gene,group.colors=color_CRC,size = 2)
dev.off()

#############
##########Sctpa绘制通路热图
library(dplyr)
#counts = CC_scrna
### you could load your own data:
#counts = read.table("folder/expression_counts.txt.gz",
#                     header = TRUE,
#                     sep = '\t',
#                     row.names = 1)
#se_oj = CreateSeuratObject(GetAssayData(object=CRC_merged_all,assay="RNA",slot="data"),
#  meta.data=CRC_merged_all@meta.data)
#CRC_subset<-CRC_merged_all[,sample(1:ncol(CRC_merged_all),50000)]
#rm(CRC_merged_all)

CC_scrna_pas = cal_PAS(seurat_object = CRC_subset,
              tool = 'AUCell',   ## GSVA, ssGSEA, plage, zscore or Vision
              normalize = "log",
              species = 'human', 
              pathway='kegg')
Idents(CC_scrna_pas)<-CC_scrna_pas$Initial_annotation

CC_scrna_pas = FindVariableFeatures(CC_scrna_pas, verbose = FALSE)
CC_scrna_pas = ScaleData(CC_scrna_pas)
CC_scrna_pas = RunPCA(CC_scrna_pas)
CC_scrna_pas = RunUMAP(CC_scrna_pas,dims = 1:8)

markers_pathway = FindAllMarkers(CC_scrna_pas,logfc.threshold = 0)
save(markers_pathway,file="annatation_markers_pathway.RData")
color_CRC<-c("#7FBCD2","#A5F1E9","#FFEEAF","#7FB77E","#F7F6DC","#B1D7B4","#FFC090",
  "#F3E0B5","#FFAE6D","#FECD70","#E3C770","#94B49F","#CEE5D0","#ECB390","#F2D7D9",
  "#D3CEDF","#B2C8DF","#C4D7E0","#F8F9D7","#76BA99","#ADCF9F",
  "#CED89E","#FFDCAE","#FFE3A9","#FFC3C3","#FF8C8C","#97C4B8","#F1F0C0","#B1BCE6")
pathways = markers_pathway %>% group_by(Initial_annotation) %>% top_n(n=5,wt=avg_log2FC)  #For Seurat>4.0, set "wt=avg_logFC"

pdf("step3_DoHeatmap_pathway_for_annotation.pdf",width=7,height=8)
DoHeatmap(CC_scrna_pas,features=pathways$gene,group.colors=color_CRC,size = 2)
dev.off()


############

top_markers1 = markers %>% group_by(cluster) %>% top_n(n=1,wt=avg_log2FC)  #For Seurat>4.0, set "wt=avg_logFC"

p_list<-lapply(1:nrow(top_markers1),function(i){
#pdf(file=paste0("marker_",top_markers1$gene[i],"_UMAP.pdf"), width=5, height=5)
p<-FeaturePlot(CC_scrna, features=top_markers1$gene[i],reduction ="tsne", pt.size = 0.6) +
        featureplot_theme() +
        theme(aspect.ratio=1) +
        scale_color_viridis(name = "Expression", direction = -1, trans = "sqrt") +
        theme(legend.text=element_text(size=3),
              legend.title=element_text(size=3),
              legend.key.size = unit(0.8, "cm"),
        plot.title = element_text(face = "italic"))
#dev.off()
return(p)
  })
pdf("step3_marker_top_markers_UMAP.pdf", width=15, height=15)
ggarrange(p_list[[1]],p_list[[2]],p_list[[3]],p_list[[4]],p_list[[5]],p_list[[6]],p_list[[7]],p_list[[8]],ncol =3)
dev.off()
#########################
#两个样本之间得所有细胞得比例图：
#
pdf("cellproportion_in_samples.pdf",width=3,height=5)
stack_plot(CC_scrna,options=c("annotation","samples"),colors= color_scanpy_8)
dev.off()

  ggplot_df<-cancer_scrna@meta.data[,c("seurat_clusters","samples")]
  
  a<-tapply(ggplot_df$seurat_clusters,factor(ggplot_df$samples),function(index2){table(index2) })
  freq_s<-data.frame(C1=a$C1,R2=a$R2)
  freq_s<-freq_s[,c(1,2,4)]
  colnames(freq_s)[1]<-"Cluster"
  write.csv(freq_s,file = "freq_sample_epicells.csv")


