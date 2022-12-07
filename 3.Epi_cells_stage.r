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
#library(slingshot)

options(future.globals.maxSize = 10000 * 1024^2)
source("/share/pub/dengcy/Singlecell/CC_space/src/Functions.r")
source("/share/pub/dengcy/Singlecell/CC_space/src/library.r")
#"E:/OneDrive/colorectal_sc/CC_space/src/Functions.r"
#"E:/OneDrive/colorectal_sc/CC_space/3.Epi_cells_stage"
setwd("/share/pub/dengcy/Singlecell/CC_space/3.Epi_cells_stage")

load("/share/pub/dengcy/Singlecell/CC_space/1.scRNA_seq_prepare/CRC_merged_all.RData")
Idents(CRC_merged_all)<-CRC_merged_all$Initial_annotation


############################
####1.get the epi cells
cancer_scrna<-subset(CRC_merged_all,Initial_annotation %in% c("Epithelial_cells"))

rm(CRC_merged_all)

cancer_scrna$stage<-"StageIII"
cancer_scrna$stage[cancer_scrna$samples=="CC1"]<-"StageI"
cancer_scrna$stage[cancer_scrna$orig.ident %in% c("KUL01-B","KUL01-T","KUL28-B","KUL28-T","KUL30-B","KUL30-T",
  "SMC01-T","SMC05-T","SMC09-T","SMC10-T","SMC15-T","SMC18-T","scrEXT001","scrEXT002","scrEXT021",
  "scrEXT022","scrEXT027","scrEXT028")]<-"StageII"
cancer_scrna$stage[cancer_scrna$orig.ident %in% c("KUL31-B","KUL31-T","SMC07-T","SMC24-T","scrEXT018","scrEXT019")]<-"StageI"

cancer_scrna_stage1<-subset(cancer_scrna,stage %in% c("StageI","StageII"))
cancer_scrna_stage2<-subset(cancer_scrna,stage== "StageIII")

cancer_scrna_stage1<- cluster_pca_umap(cancer_scrna_stage1, assay = "RNA",reduction="harmony",cluster_res = 0.1)
cancer_scrna_stage2<- cluster_pca_umap(cancer_scrna_stage2, assay = "RNA",reduction="harmony",cluster_res = 0.1)
save(cancer_scrna_stage1,file="cancer_scrna_stage1.RData")
save(cancer_scrna_stage2,file="cancer_scrna_stage2.RData")

color_CRC<-c("#7FBCD2","#A5F1E9","#FFEEAF","#7FB77E","#F7F6DC","#B1D7B4","#FFC090",
  "#F3E0B5","#FFAE6D","#FECD70","#E3C770","#94B49F","#CEE5D0","#ECB390","#F2D7D9",
  "#D3CEDF","#B2C8DF","#C4D7E0","#F8F9D7","#76BA99","#ADCF9F",
  "#CED89E","#FFDCAE","#FFE3A9","#FFC3C3","#FF8C8C","#97C4B8","#F1F0C0","#B1BCE6")


############################
#####2. cluster
pdf(file="Malignant_cells_umap_0.3.pdf",height=8,width=8)
DimPlot(cancer_scrna,pt.size=0.8,label = TRUE, repel=TRUE,label.size = 6)+ 
umap_theme()+ labs(x="UMAP",y="")+
ggtitle("Malignant cells")+
        scale_colour_manual(name = "Cluster", values = color_scanpy_viridis28[c(1,4,5,7,20,11)]) +
        theme(aspect.ratio=1)
dev.off()

pdf(file="Malignant_cells_tsne_0.3.pdf",height=8,width=8)
DimPlot(cancer_scrna,pt.size=0.8,label = TRUE,reduction="tsne", repel=TRUE,label.size = 6)+ 
umap_theme()+ labs(x="TSNE",y="")+
ggtitle("Malignant cells")+
        scale_colour_manual(name = "Cluster", values = color_scanpy_viridis28[c(1,4,5,7,20,11)]) +
        theme(aspect.ratio=1)
dev.off()

pdf(file="Malignant_cells_samples_0.3.pdf",height=8,width=8)
DimPlot(cancer_scrna,group.by = "samples",pt.size=0.8,label = TRUE, repel=TRUE)+ umap_theme()+ 
ggtitle("Malignant cells in samples")+labs(x="UMAP",y="")+
        scale_colour_manual(name = "Samples", values =c("#E05D5D","#00A19D") ) +
        theme(aspect.ratio=1)
dev.off()

####################
# inferCNV
#######################
setwd("/share/pub/dengcy/Singlecell/CC_space/3.Epi_cells_stage")
load("/share/pub/dengcy/Singlecell/CC_space/1.scRNA_seq_prepare/CRC_merged_all.RData")
library(AnnoProbe)
library(SeuratDisk)
library(Seurat)

Idents(CRC_merged_all)<-CRC_merged_all$Initial_annotation

a1<-unlist(lapply(c("Macrophage","Fibroblasts"),function(x){
    a<-colnames(CRC_merged_all)[Idents(CRC_merged_all)==x]
    a<-sample(a,2000)
    return(a)
  }))
CRC_immune_subset<-CRC_merged_all[,a1]
rm(CRC_merged_all)
load("cancer_scrna_stage1.RData")
CC_infercnv<-merge(cancer_scrna_stage1,CRC_immune_subset)
CC_infercnv$cnv_annotation<-c(cancer_scrna_stage1$seurat_clusters,CRC_immune_subset$Initial_annotation)

a<-GetAssayData(CC_infercnv,assay = "RNA",slot="counts")
a <- CreateSeuratObject(counts = a,
                              min.cells = 3,
                              meta.data =CC_infercnv@meta.data,
                              min.features = 200)
DefaultAssay(a) <- "RNA"
SaveH5Seurat(a, "cancer_scrna_stage1.h5seurat")
Convert("cancer_scrna_stage1.h5seurat", dest="h5ad")

###
load("cancer_scrna_stage2.RData")
CC_infercnv<-merge(cancer_scrna_stage2,CRC_immune_subset)
CC_infercnv$cnv_annotation<-c(cancer_scrna_stage2$seurat_clusters,CRC_immune_subset$Initial_annotation)

rm(cancer_scrna_stage2,CRC_immune_subset)
a<-GetAssayData(CC_infercnv,assay = "RNA",slot="counts")
a <- CreateSeuratObject(counts = a,
                              min.cells = 3,
                              meta.data =CC_infercnv@meta.data,
                              min.features = 200)
DefaultAssay(a) <- "RNA"
SaveH5Seurat(a, "cancer_scrna_stage2.h5seurat")
Convert("cancer_scrna_stage2.h5seurat", dest="h5ad")
########################
###remove mormal cells
#####################
load("cancer_scrna_stage1.RData")
##16115个细胞
stage1_t<-read.csv("/share/pub/dengcy/Singlecell/CC_space/3.Epi_cells_stage/cancer_infercnv_stage1.csv")
tn<-intersect(colnames(cancer_scrna_stage1),stage1_t$X[stage1_t$cnv_status=="tumor"])
cancer_scrna_stage1<-cancer_scrna_stage1[,tn]
save(cancer_scrna_stage1,file="cancer_scrna_stage1.RData")
##11826个细胞
load("cancer_scrna_stage2.RData")
#18267
stage2_t<-read.csv("/share/pub/dengcy/Singlecell/CC_space/3.Epi_cells_stage/cancer_infercnv_stage2.csv")
tn<-intersect(colnames(cancer_scrna_stage2),stage2_t$X[stage2_t$cnv_status=="tumor"])
cancer_scrna_stage2<-cancer_scrna_stage2[,tn]
#14026
save(cancer_scrna_stage2,file="cancer_scrna_stage2.RData")

cancer_scrna_stage1<- cluster_pca_umap(cancer_scrna_stage1, assay = "RNA",reduction="harmony",cluster_res = 0.2)
cancer_scrna_stage2<- cluster_pca_umap(cancer_scrna_stage2, assay = "RNA",reduction="harmony",cluster_res = 0.2)

pdf(file="stage1_magliant_clusters_umap_0.2.pdf",height=5,width=5)
DimPlot(cancer_scrna_stage1,group.by="seurat_clusters",pt.size=0.3,label = TRUE, repel=TRUE,label.size = 6)+ 
umap_theme()+ labs(x="UMAP",y="")+
ggtitle("StageI/II Epithelial cells")+
        scale_colour_manual(name = "seurat_clusters", values = color_CRC) +
        theme(aspect.ratio=1)
dev.off()

pdf(file="stage2_magliant_clusters_umap_0.2.pdf",height=5,width=5)
DimPlot(cancer_scrna_stage2,group.by="seurat_clusters",pt.size=0.3,label = TRUE, repel=TRUE,label.size = 6)+ 
umap_theme()+ labs(x="UMAP",y="")+
ggtitle("StageIII Epithelial cells")+
        scale_colour_manual(name = "seurat_clusters", values = color_CRC) +
        theme(aspect.ratio=1)
dev.off()
save(cancer_scrna_stage2,file="cancer_scrna_stage2.RData")
save(cancer_scrna_stage1,file="cancer_scrna_stage1.RData")


############################
#4.monocle analysis for each stage
########################
setwd("/share/pub/dengcy/Singlecell/CC_space/3.Epi_cells_stage")
#source("/share/pub/dengcy/Singlecell/CC_space/src/STrectal/Functions.r")
load("cancer_scrna_stage1.RData")
#BiocManager::install("monocle")
output="/share/pub/dengcy/Singlecell/CC_space/3.Epi_cells_stage/monocle/"
library(monocle)
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

seuTocds<-function(seu_obj=cancer_scrna_stage1,
  output="/share/pub/dengcy/Singlecell/CC_space/3.Epi_cells_stage/monocle/",
  index=""){
        expr_mat1<-as(as.matrix(seu_obj@assays$RNA@counts),'sparseMatrix')
        p_data<-seu_obj@meta.data
        p_data$celltype<-seu_obj$Initial_annotation
        f_data<-data.frame(gene_short_name=rownames(seu_obj),row.names=rownames(seu_obj))
        pd<-new('AnnotatedDataFrame',data=p_data)
        pf<-new('AnnotatedDataFrame',data=f_data)
        cds<-newCellDataSet(expr_mat1,phenoData=pd,featureData=pf,lowerDetectionLimit=0.5,expressionFamily=negbinomial.size())
        cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr = 0.1)
#聚类分群
disp_table <- dispersionTable(cds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
cds <- setOrderingFilter(cds, unsup_clustering_genes$gene_id)
cds<- reduceDimension(cds, max_components = 2, num_dim = 8, reduction_method = 'tSNE', verbose = T,check_duplicates = FALSE)
cds<- clusterCells(cds, num_clusters = 8)

expressed_genes <- row.names(subset(fData(cds),  num_cells_expressed >= 10))
clustering_DEG_genes <- differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr = '~Cluster',cores = 8)
ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:800] #排序后取前800个基因
cds <- setOrderingFilter(cds,ordering_genes = ordering_genes)
cds <- reduceDimension(cds, method = 'DDRTree')
cds <- orderCells(cds)
return(cds)
}


output="/share/pub/dengcy/Singlecell/CC_space/3.Epi_cells_stage/monocle/"
index=1
pdf(paste0(output,index,"_cell_trajectoryPseudotime.pdf"))
plot_cell_trajectory(cds1, cell_size=1,color_by = "Pseudotime")+
scale_color_gradient(low = "#f4e2d8",high = "#ba5370")
dev.off()
pdf(paste0(output,index,"_cell_trajectoryCluster.pdf"))
plot_cell_trajectory(cds1, cell_size=2,color_by = "Cluster")+
scale_colour_manual(name = "Cluster", values = color_CRC[1:7] )
dev.off()
pdf(paste0(output,index,"_cell_trajectoryState.pdf"))
plot_cell_trajectory(cds1, cell_size=2,color_by = "State")+
scale_colour_manual(name = "State", values = c("#f38181","#fce38a","#eaffd0","#95e1d3","#64958f") )
dev.off()

cds1<-seuTocds(seu_obj=cancer_scrna_stage1,index="1")
save(cds1,file=paste0(output,"cds_cancer_scrna_stage1.RData"))
load(paste0(output,"cds_cancer_scrna_stage1.RData"))
#########################
###############得到emt相关的基因做交集
###############
#2）.hallmark genes
  x <- readLines("/share/pub/dengcy/Singlecell/CC_space/3.Epi_cells_stage/h.all.v7.4.symbols.gmt")
  res <- strsplit(x, "\t")
  names(res) <- vapply(res, function(y) y[1], character(1))
  genes.by.pathway.h <- lapply(res, "[", -c(1:2))
  save(genes.by.pathway.h,file="genes.by.pathway.h.RData")
  genes.pathway.h <- unique(unlist(genes.by.pathway.h))
  load("/share/pub/dengcy/Singlecell/CC_space/4.stage_visualize/genes.by.pathway.h.RData")
#4）EMT
  EMT_genes<-read.delim("/share/pub/dengcy/Singlecell/CC_space/4.stage_visualize/EMTgene.txt")
#intersecton
  dehallgene<-intersect(unique(unlist(genes.by.pathway.h)),EMT_genes$x)

###提取轨迹上的差异基因
#expressed_genes <- row.names(subset(fData(cds1),  num_cells_expressed >= 10))
dehallgene<-intersect(dehallgene,expressed_genes)
cds1_sub <- cds1[,pData(cds1)$State %in% c('1','3','4','5')] #提取State1、2、6上的细胞
pseudotime_de1 <- differentialGeneTest(cds1_sub[dehallgene,], fullModelFormulaStr = "~sm.ns(Pseudotime)") #求拟时序上基因基因的显著性
pseudotime_de1.1 <- subset(pseudotime_de1, qval < 1e-4) #取q值小于0.0001的基因
save(pseudotime_de1,pseudotime_de1.1,file=paste0(output,"pseudotime_de1.RData"))
library(RColorBrewer)

#做热图
pdf("/share/pub/dengcy/Singlecell/CC_space/3.Epi_cells_stage/monocle/cds1_pseudotime_heatmap.pdf")
plot_pseudotime_heatmap(cds1_sub[pseudotime_de1.1$gene_short_name,],
  return_heatmap=T,
  num_clusters=4,
  hmcols = colorRampPalette(rev(brewer.pal(11, "Spectral")))(70),
  show_rownames = F)
dev.off()


load("cancer_scrna_stage2.RData")
cds2<-seuTocds(seu_obj=cancer_scrna_stage2,index="2")
output="/share/pub/dengcy/Singlecell/CC_space/3.Epi_cells_stage/monocle/"
index=2
pdf(paste0(output,index,"_cell_trajectoryPseudotime.pdf"))
plot_cell_trajectory(cds2, cell_size=2,color_by = "Pseudotime")+
scale_color_gradient(low = "#f4e2d8",high = "#ba5370")
dev.off()
pdf(paste0(output,index,"_cell_trajectoryCluster.pdf"))
plot_cell_trajectory(cds2, cell_size=2,color_by = "Cluster")+
scale_colour_manual(name = "Cluster", values = color_CRC[1:7] )
dev.off()
pdf(paste0(output,index,"_cell_trajectoryState.pdf"))
plot_cell_trajectory(cds2, cell_size=2,color_by = "State")+
scale_colour_manual(name = "State", values = c("#f38181","#fce38a","#95e1d3","#64958f") )
dev.off()

save(cds2,file=paste0(output,"cds_cancer_scrna_stage2.RData"))

load(paste0(output,"cds_cancer_scrna_stage2.RData"))
#expressed_genes <- row.names(subset(fData(cds2),  num_cells_expressed >= 10))
#cds2_sub <- cds2[,pData(cds2)$State %in% c('1','3','4','5')] #提取State1、2、6上的细胞
pseudotime_de2 <- differentialGeneTest(cds2[dehallgene,], fullModelFormulaStr = "~sm.ns(Pseudotime)") #求拟时序上基因基因的显著性
pseudotime_de2.1 <- subset(pseudotime_de2, qval < 1e-4) #取q值小于0.0001的基因
save(pseudotime_de2,pseudotime_de2.1,file=paste0(output,"pseudotime_de2.RData"))
#做热图
pdf("/share/pub/dengcy/Singlecell/CC_space/3.Epi_cells_stage/monocle/cds2_pseudotime_heatmap.pdf")
p = plot_pseudotime_heatmap(cds2[pseudotime_de2.1$gene_short_name,],
  return_heatmap=T,
  num_clusters=4,
  hmcols = colorRampPalette(rev(brewer.pal(11, "Spectral")))(70),
  show_rownames = F)
dev.off()

load(paste0(output,"pseudotime_de2.RData"))

write.csv(pseudotime_de2,file="/share/pub/dengcy/Singlecell/CC_space/SupplementFiles10.26/pseudotime_de2.csv")
load(paste0(output,"pseudotime_de1.RData"))
write.csv(pseudotime_de1,file="/share/pub/dengcy/Singlecell/CC_space/SupplementFiles10.26/pseudotime_de1.csv")
###################
#选择两组基因的交集
emt_pseudotim_gene<-intersect(pseudotime_de2.1$gene_short_name,pseudotime_de1.1$gene_short_name)
#365
save(emt_pseudotim_gene,file="emt_pseudotim_gene.RData")
###############
#Pathway heatmap for stage
################
library(Seurat)
library(SeuratObject)
library("dplyr")
library("foreach")
library("ggplot2")
library(reshape2)  # melt
library(gridExtra) # arrangeGrob
library(ggpubr)    # as_ggplot
library(cowplot)   # draw_plot_label
library(ggtext)
library(testSctpa)
library(VISION)
library(AUCell)

setwd("/share/pub/dengcy/Singlecell/CC_space/3.Epi_cells_stage")
load("cancer_scrna_stage1.RData")
counts = load_counts()
se_oj = CreateSeuratObject(counts)
se_oj = cal_PAS(seurat_object = se_oj,
              tool = 'AUCell',   ## GSVA, ssGSEA, plage, zscore or Vision
              normalize = 'log',
              species = 'mouse', 
              pathway='kegg')

DefaultAssay(object = cancer_scrna_stage1) <- "RNA"

gmt<-"/share/pub/dengcy/Singlecell/CC_space/3.Epi_cells_stage/h.all.v7.5.1.symbols.gmt"
cancer_scrna_stage1 = ScaleData(cancer_scrna_stage1)
cancer_scrna_stage1 = cal_PAS(seurat_object = cancer_scrna_stage1,
                       tool = 'AUCell',   ## GSVA, ssGSEA, plage, zscore or Vision
                       normalize = "log",
                       gmt_file=gmt)

library(dplyr)
Idents(cancer_scrna_stage1)<-cancer_scrna_stage1$State

markers = FindAllMarkers(cancer_scrna_stage1,logfc.threshold = 0)
pathways = markers %>% group_by(cluster) %>% top_n(n=5,wt=avg_log2FC)  #For Seurat>4.0, set "wt=avg_logFC"
write.csv(pathways,file="/share/pub/dengcy/Singlecell/CC_space/SupplementFiles10.26/pathways_stage1_State.csv")
pdf("state_hallmark_pathway.pdf")
DoHeatmap(cancer_scrna_stage1,features=pathways$gene,slot="data",
  group.colors =  c("#f38181","#fce38a","#eaffd0","#95e1d3","#64958f"))
dev.off()

Idents(cancer_scrna_stage1)<-cancer_scrna_stage1$Cluster
markers = FindAllMarkers(cancer_scrna_stage1,logfc.threshold = 0)
pathways = markers %>% group_by(cluster) %>% top_n(n=5,wt=avg_log2FC)  #For Seurat>4.0, set "wt=avg_logFC"
write.csv(pathways,file="/share/pub/dengcy/Singlecell/CC_space/SupplementFiles10.26/pathways_stage1_Cluster.csv")

pdf("cluster_hallmark_pathway.pdf")
DoHeatmap(cancer_scrna_stage1,features=pathways$gene,slot="data",
  group.colors = color_CRC[1:7])
dev.off()

###############
load("cancer_scrna_stage2.RData")
cancer_scrna_stage2 = ScaleData(cancer_scrna_stage2)
DefaultAssay(object = cancer_scrna_stage2) <- "RNA"
cancer_scrna_stage2 = cal_PAS(seurat_object = cancer_scrna_stage2,
                       tool = 'AUCell',   ## GSVA, ssGSEA, plage, zscore or Vision
                       normalize = "log",
                       gmt_file=gmt)

library(dplyr)
Idents(cancer_scrna_stage2)<-cancer_scrna_stage2$State

markers = FindAllMarkers(cancer_scrna_stage2,logfc.threshold = 0)
pathways = markers %>% group_by(cluster) %>% top_n(n=5,wt=avg_log2FC)  #For Seurat>4.0, set "wt=avg_logFC"
write.csv(pathways,file="/share/pub/dengcy/Singlecell/CC_space/SupplementFiles10.26/pathways_stage3_State.csv")

pdf("state_hallmark_pathway2.pdf")
DoHeatmap(cancer_scrna_stage2,features=pathways$gene,slot="data",
  group.colors =  c("#f38181","#fce38a","#eaffd0","#95e1d3","#64958f"))
dev.off()

Idents(cancer_scrna_stage2)<-cancer_scrna_stage2$Cluster
markers = FindAllMarkers(cancer_scrna_stage2,logfc.threshold = 0)
pathways = markers %>% group_by(cluster) %>% top_n(n=5,wt=avg_log2FC)  #For Seurat>4.0, set "wt=avg_logFC"
write.csv(pathways,file="/share/pub/dengcy/Singlecell/CC_space/SupplementFiles10.26/pathways_stage3_cluster.csv")

pdf("cluster_hallmark_pathway2.pdf")
DoHeatmap(cancer_scrna_stage2,features=pathways$gene,slot="data",
  group.colors = color_CRC[1:7])
dev.off()


##########
#infer cnv output
infercnv_stage1<-read.csv("/share/pub/dengcy/Singlecell/CC_space/3.Epi_cells_stage/cancer_infercnv_stage1.csv")
infercnv_stage1<-infercnv_stage1[,c("X","cnv_annotation","cnv_leiden","cnv_score","cnv_status")]

infercnv_stage2<-read.csv("/share/pub/dengcy/Singlecell/CC_space/3.Epi_cells_stage/cancer_infercnv_stage2.csv")
infercnv_stage2<-infercnv_stage2[,c("X","cnv_annotation","cnv_leiden","cnv_score","cnv_status")]

infercnv_stage<-rbind(infercnv_stage1,infercnv_stage2)
write.csv(infercnv_stage,file="/share/pub/dengcy/Singlecell/CC_space/SupplementFiles10.26/infercnv_stage_result.csv")

