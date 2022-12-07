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
#source("/share/pub/dengcy/Singlecell/CC_space/src/library.r")
#"E:/OneDrive/colorectal_sc/CC_space/src/Functions.r"
#"E:/OneDrive/colorectal_sc/CC_space/3.Epi_cells_stage"
setwd("/share/pub/dengcy/Singlecell/CC_space/3.2.Immune_cells_stage")
load("/share/pub/dengcy/Singlecell/CC_space/1.scRNA_seq_prepare/CRC_merged_all.RData")

Idents(CRC_merged_all)<-CRC_merged_all$Initial_annotation

Immune_scrna<-subset(CRC_merged_all,Initial_annotation != "Epithelial_cells")
rm(CRC_merged_all)
gc()
Immune_scrna<-cluster_pca_umap(Immune_scrna, assay = "RNA",reduction="pca",cluster_res = 0.2)
Immune_scrna$stage<-"StageIII"
Immune_scrna$stage[Immune_scrna$samples=="CC1"]<-"StageI"
Immune_scrna$stage[Immune_scrna$orig.ident %in% c("KUL01-B","KUL01-T","KUL28-B","KUL28-T","KUL30-B","KUL30-T",
  "SMC01-T","SMC05-T","SMC09-T","SMC10-T","SMC15-T","SMC18-T","scrEXT001","scrEXT002","scrEXT021",
  "scrEXT022","scrEXT027","scrEXT028")]<-"StageII"
Immune_scrna$stage[Immune_scrna$orig.ident %in% c("KUL31-B","KUL31-T","SMC07-T","SMC24-T","scrEXT018","scrEXT019")]<-"StageI"

###过滤低细胞群体
cluster_a1<- table(Immune_scrna$seurat_clusters)
Immune_scrna<- subset(Immune_scrna, subset = seurat_clusters %in% names(cluster_a1)[which(cluster_a1>50)])
Immune_scrna$seurat_clusters <- as.numeric(Immune_scrna$seurat_clusters)
Immune_scrna[,sample(1:ncol(Immune_scrna),20000)]
markers= FindAllMarkers(Immune_scrna[,sample(1:ncol(Immune_scrna),10000)],logfc.threshold = 0)
markers<-markers[!grepl("MT-",markers$gene,ignore.case=F),]

save(Immune_scrna,file="Immune_scrna.RData")
save(markers,file="Immune_scrna_seucluster_sample.RData")

cluster_de<- lapply(0:10,function(x) markers[markers$avg_log2FC>0 & markers$cluster==x ,])


library("Seurat")
library("dplyr")
library("celldex")
library("SingleR")
library("pracma")
library("SingleCellExperiment")
library(scuttle)
setwd("/share/pub/dengcy/Singlecell/CC_space/3.2.Immune_cells_stage")
load("Immune_scrna.RData")
ref <- HumanPrimaryCellAtlasData() 
#############
common <- intersect(rownames(ref), rownames(Immune_scrna))
Immune_scrna.singler <- Immune_scrna[common,]
ref <- ref[common,]
Immune_scrna.singler.sce <-SingleCellExperiment(assays = List(counts = GetAssayData(object =Immune_scrna.singler, slot = "counts")))
Immune_scrna.singler.sce <- logNormCounts(Immune_scrna.singler.sce)

pred.brca_immune <- SingleR(test = Immune_scrna.singler.sce, 
  ref = ref, assay.type.test=1,labels = ref$label.main)
save(pred.brca_immune ,file="singler_pred.brca_immune.RData")

singler.results <- merge(data.frame(cell = rownames(pred.brca_immune), singler = pred.brca_immune$labels), 
                         data.frame(cell = rownames(Immune_scrna@meta.data), 
                                    cluster = Immune_scrna@meta.data$seurat_clusters), 
                         by = "cell", 
                         all.y = FALSE)
singler.results$cell <- NULL
singler.results$count <- 1
singler.results <- aggregate(count ~ ., singler.results, FUN = sum)
singler.final <- singler.results %>% group_by(cluster) %>% top_n(n = 1, wt = count)
write.csv(singler.final,file="singler.final.immune.csv")

#####绘制注释后的细胞图
#
seurat_clusters<-as.vector(Immune_scrna$seurat_clusters)
annotation<-as.vector(seurat_clusters)
#"Terminally.Exhausted.T"
#"B.cell.activate.memory"

annotation[seurat_clusters == 1]<-"T"
annotation[seurat_clusters == 2]<-"T"
annotation[seurat_clusters == 3]<-"Macrophage"
annotation[seurat_clusters == 4]<-"B"
annotation[seurat_clusters == 5]<-"Fibroblasts"
#"NK.cell"
annotation[seurat_clusters == 6]<-"B"
#annotation[seurat_clusters == 1]<-"Neutrophils"
annotation[seurat_clusters == 7] <-"Edothelial_cells"
annotation[seurat_clusters == 8] <-"T"
annotation[seurat_clusters == 9]<-"Epi"
annotation[seurat_clusters == 10] <-"Tissue_stem_cells"
annotation[seurat_clusters == 11]<-"T"
annotation[seurat_clusters == 12]<-"B"
annotation[seurat_clusters == 13]<-"Neutrophils"
annotation[seurat_clusters == 14]<-"CMP"
annotation[seurat_clusters == 15]<-"Monocyte"


Immune_scrna$annotation<-annotation
#删除无用的细胞
Immune_scrna<-subset(Immune_scrna,annotation != "Epi")
Immune_scrna<-subset(Immune_scrna,annotation != "Edothelial_cells")
Immune_scrna<-subset(Immune_scrna,annotation != "Monocyte")
Immune_scrna<-subset(Immune_scrna,annotation != "Fibroblasts")
Immune_scrna<-subset(Immune_scrna,annotation != "Tissue_stem_cells")
Immune_scrna<-subset(Immune_scrna,annotation != "CMP")
Immune_scrna<-cluster_pca_umap(Immune_scrna, assay = "RNA",reduction="pca",cluster_res = 0.2)
save(Immune_scrna,file="Immune_scrna.RData")

pdf(file="Immune_cells_TSNE_0.2.pdf",height=8,width=8)
DimPlot(Immune_scrna,group.by = "annotation",pt.size=0.8,label = TRUE,reduction="tsne", repel=TRUE)+ 
umap_theme()+ ggtitle("Immune cells")+labs(x="TSNE",y="")+
        scale_colour_manual(name = "Cluster", values = color_CRC[14:29]) +
        theme(aspect.ratio=1)
dev.off()

pdf(file="Immune_seurat_clusters_TSNE_0.2.pdf",height=8,width=8)
DimPlot(Immune_scrna,group.by = "seurat_clusters",pt.size=0.8,label = TRUE,reduction="tsne", repel=TRUE)+ 
umap_theme()+ ggtitle("seurat_clusters")+labs(x="TSNE",y="")+
        scale_colour_manual(name = "seurat_clusters", values = color_CRC[14:29]) +
        theme(aspect.ratio=1)
dev.off()
############################
#####2.聚类
Immune_scrna$annotation<- factor(Immune_scrna$annotation,levels= c("CD4.TCM","Cytotoxic.CD8.T","Terminally.Exhausted CD8T","B.cell","B.cell.activite.memory","TAM","NK.cell","Neutrophils"))
Idents(Immune_scrna) <- Immune_scrna$annotation
color_CRC<-c("#7FBCD2","#A5F1E9","#FFEEAF","#7FB77E","#F7F6DC","#B1D7B4","#FFC090",
  "#F3E0B5","#FFAE6D","#FECD70","#E3C770","#94B49F","#CEE5D0","#ECB390","#F2D7D9",
  "#D3CEDF","#B2C8DF","#C4D7E0","#F8F9D7","#76BA99","#ADCF9F",
  "#CED89E","#FFDCAE","#FFE3A9","#FFC3C3","#FF8C8C","#97C4B8","#F1F0C0","#B1BCE6")
pdf(file="Immune_cells_TSNE_0.2.pdf",height=8,width=8)
DimPlot(Immune_scrna,group.by = "annotation",pt.size=0.8,label = TRUE,reduction="tsne", repel=TRUE)+ 
umap_theme()+ ggtitle("Immune cells")+labs(x="TSNE",y="")+
        scale_colour_manual(name = "Cluster", values = color_CRC[14:29]) +
        theme(aspect.ratio=1)
dev.off()



Lineage<-c("CD3E","CD4","CD8A")
 Naive_memory<- c("SELL","CCR7","IL7R","FAS","CD27","CD28","ITGAE","KLRB1")
 Inhibitory<-c("PDCD1","TIGIT","HAVCR2","LAG3","CTLA4","TNFRSF14","BTLA","CD244","CD160")
Ectonucleases<-c("CD38","ENTPD1","NT5E")
 Activation<-c("CD69","IL2RA","ICOS","TNFRSF4","TNFRSF9","HLA-DRA","FASLG","CD40LG")
Cytolytic<-c("GZMA","GZMB","GZMH","GZMK","GZMM","PRF1","NKG7","GNLY","LAMP1")
 Cytokines<-c("IL2","IFNG","TNF")


####################
#根据文献注释
cluster_de<- lapply(0:10,function(x) markers[markers$avg_log2FC>0 & markers$p_val<0.01 & markers$cluster==x ,])
top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
write.csv(top10,file="/share/pub/dengcy/Singlecell/CC_space/1.scRNA_seq_prepare/immune_cluster_markers.csv")

#49074 
Immune_scrna<-subset(Immune_scrna,seurat_clusters !=  7)
Immune_scrna<-subset(Immune_scrna,seurat_clusters !=  8)
Immune_scrna<-subset(Immune_scrna,seurat_clusters !=  10)

seurat_clusters<-as.vector(Immune_scrna$seurat_clusters)
annotation<-as.vector(seurat_clusters)
#"Terminally.Exhausted.T"
#"B.cell.activate.memory"
annotation[seurat_clusters == 0]<-"Th" 
annotation[seurat_clusters == 1]<-"TAM"
annotation[seurat_clusters == 2]<-"Teff"
annotation[seurat_clusters == 3]<-"B"
annotation[seurat_clusters == 4]<-"Tex"
annotation[seurat_clusters == 5]<-"Plasma"
#"NK.cell"
annotation[seurat_clusters == 6]<-"Treg"
annotation[seurat_clusters == 9]<-"MDSC"
#T_Helper_Cells
#Effector_T_Cells
#Exhausted_T_Cells
#Regulatory_T_cells
Immune_scrna$annotation<-annotation
Idents(Immune_scrna)<-annotation

#####################
#可视化一些功能marker

#T_scrna<-subset(Immune_scrna,annotation %in% c("Th","Teff","Tex","Treg"))
Anti_inflammatory<-c("TIGIT","IDO1","LGALS3","PDCD1","FOXP3","ENTPD1","CD274","CSF2","CTLA4","CXCL12","CXCL5","IL8","MIF","PTGS2","VEGFA")
Pro_inflammatory<-c("IL1A","IL1B","TNF","IFNG","TBX21","CCL3","CCL4","PRF1","GZMA","GZMB","GZMK","GZMH","CD8A","FASLG","CCL2","CCL20","IL2","IL6","IL12a","IL17a","IL23a","PTGS2","TLR4","TNF")

Anti_inflammatory<-intersect(Anti_inflammatory,rownames(Immune_scrna))
Pro_inflammatory<-intersect(Pro_inflammatory,rownames(Immune_scrna))

Immune_scrna<-AddModuleScore(Immune_scrna,list(Anti_inflammatory,Pro_inflammatory,Naive_memory,Inhibitory,Activation,Cytolytic),name=c("Anti_inflammatory","Pro_inflammatory","Naive_memory","Inhibitory","Activation","Cytolytic"))

########
#可视化新注释的细胞类型

pdf(file="Immune_cells_TSNE_anotation.pdf",height=8,width=8)
DimPlot(Immune_scrna,group.by = "annotation",pt.size=0.5,label = TRUE,reduction="tsne", repel=TRUE)+ 
umap_theme()+ ggtitle("Immune cells")+labs(x="TSNE",y="")+
        scale_colour_manual(name = "Cluster", values = color_CRC[14:29]) +
        theme(aspect.ratio=1)
dev.off()

##############
#通路观察
##############
load("Immune_scrna.RData")

library(doRNG)
library(testSctpa)

se_oj = CreateSeuratObject(GetAssayData(object=Immune_scrna,assay="RNA",slot="data"))
Immune_scrna_pas = cal_PAS(seurat_object = se_oj,
              tool = 'AUCell',   ## GSVA, ssGSEA, plage, zscore or Vision
              normalize = "log",
              species = 'human', 
              pathway='Immu')

Idents(Immune_scrna_pas)<-Immune_scrna$annotation
Immune_scrna_pas = FindVariableFeatures(Immune_scrna_pas, verbose = FALSE)
Immune_scrna_pas = ScaleData(Immune_scrna_pas)
save(Immune_scrna_pas,file="Immune_scrna_pas.RData")

markers_Immu = FindAllMarkers(Immune_scrna_pas,logfc.threshold = 0)
save(markers_Immu,file="markers_Immu.RData")

pathways = markers_Immu %>% group_by(cluster) %>% top_n(n=3,wt=avg_log2FC)  #For Seurat>4.0, set "wt=avg_logFC"

Idents(Immune_scrna_pas) <- factor(as.vector(Idents(Immune_scrna_pas)),levels=sort(unique(as.vector(Idents(Immune_scrna_pas)))))

pdf("DoHeatmap_immunecells_immu.pdf",width=12,height=6)
DoHeatmap(Immune_scrna_pas,features=pathways$gene,,
  group.colors=color_CRC,size = 2)+ggtitle("ImmunePathways")
dev.off()

####reactome

Immune_scrna_reactome = cal_PAS(seurat_object = se_oj,
              tool = 'AUCell',   ## GSVA, ssGSEA, plage, zscore or Vision
              normalize = "log",
              species = 'human', 
              pathway='reactome')

Idents(Immune_scrna_reactome)<-Immune_scrna$annotation
Immune_scrna_reactome = FindVariableFeatures(Immune_scrna_reactome, verbose = FALSE)
Immune_scrna_reactome = ScaleData(Immune_scrna_reactome)
save(Immune_scrna_reactome,file="Immune_scrna_reactome.RData")

markers_reactome = FindAllMarkers(Immune_scrna_reactome,logfc.threshold = 0)
save(markers_reactome,file="markers_reactome.RData")

pathways = markers_reactome %>% group_by(cluster) %>% top_n(n=3,wt=avg_log2FC)  #For Seurat>4.0, set "wt=avg_logFC"

Idents(Immune_scrna_reactome) <- factor(as.vector(Idents(Immune_scrna_reactome)),
	levels=sort(unique(as.vector(Idents(Immune_scrna_reactome)))))

pdf("DoHeatmap_immunecells_reactome.pdf",width=12,height=6)
DoHeatmap(Immune_scrna_reactome,features=pathways$gene,,
  group.colors=color_CRC,size = 2)+ggtitle("ReactomePathways")
dev.off()

#####################绘制免疫细胞比例图
load("Immune_scrna.RData")
load("Immune_markers.RData")

library(plyr)
library(ggplot2)


pdf("Immunecellproportion_in_stage.pdf",width=8,height=3)
stack_plot(obj=Immune_scrna,options=c("annotation","stage"),colors= color_CRC[14:29],angle=FALSE)
dev.off()

save(Immune_scrna,file="Immune_scrna.RData")
##########################促炎基因和抗炎基因集合在不同样本，不同免疫细胞中的表达情况

#################两个得分的umap图

T_scrna<-subset(Immune_scrna,annotation %in% c("Teff","Tex","Th","Treg"))
Anti_inflammatory<-c("TIGIT","IDO1","LGALS3","PDCD1","FOXP3","ENTPD1","CD274","CSF2","CTLA4","CXCL12","CXCL5","IL8","MIF","PTGS2","VEGFA")
Pro_inflammatory<-c("IL1A","IL1B","TNF","IFNG","TBX21","CCL3","CCL4","PRF1","GZMA","GZMB","GZMK","GZMH","CD8A","FASLG","CCL2","CCL20","IL2","IL6","IL12a","IL17a","IL23a","PTGS2","TLR4","TNF")

Anti_inflammatory<-intersect(Anti_inflammatory,rownames(Immune_scrna))
Pro_inflammatory<-intersect(Pro_inflammatory,rownames(Immune_scrna))

T_scrna<-AddModuleScore(T_scrna,list(Anti_inflammatory,Pro_inflammatory),name=c("Anti_inflammatory","Pro_inflammatory"))

T_fortify_can<-fortify.Seurat(T_scrna)
T_fortify_can$Anti_inflammatory1[T_fortify_can$Anti_inflammatory1>0.1]<-0.1
T_fortify_can$Anti_inflammatory1[T_fortify_can$Anti_inflammatory1< -0.6]<- -0.6

T_fortify_can$Pro_inflammatory2[T_fortify_can$Pro_inflammatory2>1]<-1
T_fortify_can$Pro_inflammatory2[T_fortify_can$Pro_inflammatory2< -1]<- -1
library(ggtext)
library(ggpubr)
Immune_anti_plot <- ggplot() +
        geom_point(data = T_fortify_can,
                   aes(x = TSNE_1, y = TSNE_2,color =Anti_inflammatory1), size = 0.2, alpha = 1) +
        umap_theme() +
        scale_colour_gradient(low="#FFFCDC",high="#516BEB")+
        theme(aspect.ratio=1) +
        theme(legend.text=element_markdown(size=14),
              legend.title=element_text(size=14)) +
        guides(colour = guide_legend(override.aes = list(size=3)))  +
        #ggtitle("G2M score")

        pdf(file="Immune_Anti_inflammatory_plot.pdf",width = 5, height = 5)
        print(Immune_anti_plot)
        dev.off()

Immune_pro_plot <- ggplot() +
        geom_point(data = T_fortify_can,
                   aes(x = TSNE_1, y = TSNE_2,color =Pro_inflammatory2), size = 0.2, alpha = 1) +
        umap_theme() +
        scale_colour_gradient(low="#FFFCDC",high="#FF1700")+
        theme(aspect.ratio=1) +
        theme(legend.text=element_markdown(size=14),
              legend.title=element_text(size=14)) +
        guides(colour = guide_legend(override.aes = list(size=3)))  +
        #ggtitle("G2M score")

        pdf(file="Immune_Pro_inflammatory_plot.pdf",width = 5, height = 5)
        print(Immune_pro_plot)
        dev.off()

####
M_scrna<-subset(Immune_scrna,annotation %in% c("TAM","Neutrophils"))

pdf("M_CXCR1_mdsc.pdf",width=14)
FeaturePlot(M_scrna, features = c("CXCR1","CXCR2"),reduction = "tsne")
dev.off()

M_scrna<-AddModuleScore(M_scrna,list(Anti_inflammatory,Pro_inflammatory),name=c("Anti_inflammatory","Pro_inflammatory"))

M_fortify_can<-fortify.Seurat(M_scrna)
M_fortify_can$Anti_inflammatory1[M_fortify_can$Anti_inflammatory1>0.1]<-0.1
M_fortify_can$Anti_inflammatory1[M_fortify_can$Anti_inflammatory1< -4]<- -4

M_fortify_can$Pro_inflammatory2[M_fortify_can$Pro_inflammatory2>50]<-50
M_fortify_can$Pro_inflammatory2[M_fortify_can$Pro_inflammatory2< -1]<- -1
library(ggtext)
library(ggpubr)
Immune_anti_plot <- ggplot() +
        geom_point(data = M_fortify_can,
                   aes(x = TSNE_1, y = TSNE_2,color =Anti_inflammatory1), size = 0.2, alpha = 1) +
        umap_theme() +
        scale_colour_gradient(low="#FFFCDC",high="#516BEB")+
        theme(aspect.ratio=1) +
        theme(legend.text=element_markdown(size=14),
              legend.title=element_text(size=14)) +
        guides(colour = guide_legend(override.aes = list(size=3)))  +
        #ggtitle("G2M score")

        pdf(file="M_Anti_inflammatory_plot.pdf",width = 5, height = 5)
        print(Immune_anti_plot)
        dev.off()

Immune_pro_plot <- ggplot() +
        geom_point(data = M_fortify_can,
                   aes(x = TSNE_1, y = TSNE_2,color =Pro_inflammatory2), size = 0.2, alpha = 1) +
        umap_theme() +
        scale_colour_gradient(low="#FFFCDC",high="#FF1700")+
        theme(aspect.ratio=1) +
        theme(legend.text=element_markdown(size=14),
              legend.title=element_text(size=14)) +
        guides(colour = guide_legend(override.aes = list(size=3)))  +
        #ggtitle("G2M score")

        pdf(file="M_Pro_inflammatory_plot.pdf",width = 5, height = 5)
        print(Immune_pro_plot)
        dev.off()



#############
#差异基因热图
###################
Idents(Immune_scrna)<-Immune_scrna$annotation
markers= FindAllMarkers(Immune_scrna,logfc.threshold = 0)
markers<-markers[!grepl("MT-",markers$gene,ignore.case=F),]

save(markers,file="Immune_markers2.0.RData")
write.csv(markers,file="Immune_markers.csv")
#

top_markers = markers %>% group_by(cluster) %>% top_n(n=5,wt=avg_log2FC)  #For Seurat>4.0, set "wt=avg_logFC"
setwd("/share/pub/dengcy/Singlecell/CC_space/3.2.Immune_cells_stage")

pdf("DotPlot_immuneMarkers.pdf",width=15,height=4)
DotPlot(Immune_scrna, features =top_markers$gene, cols = c("lightgrey", "#CD1818"))+ RotatedAxis()
dev.off()
