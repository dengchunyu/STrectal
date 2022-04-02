###################################################
#Plot heatmap of the PG clusters' mean/median marker expressions
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
load("/share/pub/dengcy/Singlecell/CC_space/1.scRNA_seq_prepare/CC_scrna.RData")
load("/share/pub/dengcy/Singlecell/CC_space/1.scRNA_seq_prepare/All_cluster_markers.RData")
#load("/share/pub/dengcy/Singlecell/CC_space/1.scRNA_seq_prepare/Allmarkers.RData")
#plan("multiprocess", workers = 10)
CC_scrna$annotation<-factor(CC_scrna$annotation)
Idents(CC_scrna)<-CC_scrna$annotation

################################
#marker hearmap
markers= FindAllMarkers(CC_scrna,logfc.threshold = 0)
save(markers,file="annotation_markers.RData")
load("annotation_markers.RData")
top_markers = markers %>% group_by(cluster) %>% top_n(n=5,wt=avg_log2FC)  #For Seurat>4.0, set "wt=avg_logFC"
pdf("DoHeatmap_markers.pdf",width=6)
DoHeatmap(CC_scrna,features=top_markers$gene,group.colors=color_scanpy_8,size = 3)
dev.off()

######################################################
##########Sctpa绘制通路热图
library(dplyr)
#counts = CC_scrna
### you could load your own data:
#counts = read.table("folder/expression_counts.txt.gz",
#                     header = TRUE,
#                     sep = '\t',
#                     row.names = 1)
se_oj = CreateSeuratObject(GetAssayData(object=CC_scrna,assay="RNA",slot="data"),
  meta.data=CC_scrna@meta.data[,c("annotation","samples")])

CC_scrna_pas = cal_PAS(seurat_object = se_oj,
              tool = 'AUCell',   ## GSVA, ssGSEA, plage, zscore or Vision
              normalize = "log",
              species = 'human', 
              pathway='kegg')
Idents(CC_scrna_pas)<-CC_scrna_pas$annotation

CC_scrna_pas = FindVariableFeatures(CC_scrna_pas, verbose = FALSE)
CC_scrna_pas = ScaleData(CC_scrna_pas)
CC_scrna_pas = RunPCA(CC_scrna_pas)
CC_scrna_pas = RunUMAP(CC_scrna_pas,dims = 1:8)

markers_pathway = FindAllMarkers(CC_scrna_pas,logfc.threshold = 0)
save(markers_pathway,file="annatation_markers_pathway.RData")

pathways = markers_pathway %>% group_by(cluster) %>% top_n(n=5,wt=avg_log2FC)  #For Seurat>4.0, set "wt=avg_logFC"
pdf("DoHeatmap_pathway.pdf",width=7,height=8)
DoHeatmap(CC_scrna_pas,features=pathways$gene,group.colors=color_scanpy_8,size = 2)
dev.off()

################################################
#

celltype_marker=c("PCGF1","LGR5","POU5F1","PROM1")
pdf(file="ClusterMarker_vlnplot.pdf")
VlnPlot(CC_scrna,features = celltype_marker,pt.size = 0,ncol = 2,group.by = "annotation")#这里定义10为干细胞群
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
pdf("marker_top_markers_UMAP.pdf", width=15, height=15)
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


