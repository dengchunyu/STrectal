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
library(slingshot)
library(infercnv)
options(future.globals.maxSize = 10000 * 1024^2)
source("/share/pub/dengcy/Singlecell/CC_space/src/Functions.r")
source("/share/pub/dengcy/Singlecell/CC_space/src/library.r")
#"E:/OneDrive/colorectal_sc/CC_space/src/Functions.r"
#"E:/OneDrive/colorectal_sc/CC_space/3.Epi_cells_stage"
setwd("/share/pub/dengcy/Singlecell/CC_space/3.Epi_cells_stage")

if(F){
load("/share/pub/dengcy/Singlecell/CC_space/1.scRNA_seq_prepare/CC_scrna.RData")

#CC_scrna$annotation<-factor(CC_scrna$annotation)
#Idents(CC_scrna)<-CC_scrna$annotation
############################
####1.获得epi细胞
cancer_scrna<-subset(CC_scrna,annotation %in% c("Epithelial_cells"))

#cancerstem_scrna<-subset(CC_scrna,annotation=="Cancer_stem_cells")

#cancer_scrna <- RunPCA(object = cancer_scrna, assay = "SCT", npcs = 50)

cancer_scrna<- cluster_pca_umap(cancer_scrna, assay = "SCT",reduction="pca",cluster_res = 0.3)

pdf(file="Malignant_seurat_clusters_umap_0.3.pdf",height=8,width=8)
DimPlot(cancer_scrna,group.by="seurat_clusters",pt.size=0.8,label = TRUE, repel=TRUE,label.size = 6)+ 
umap_theme()+ labs(x="UMAP",y="")+
ggtitle("Malignant cells")+
        scale_colour_manual(name = "seurat_clusters", values = color_scanpy_viridis28[c(1,4,5,7,20,11,13,15,16,17)]) +
        theme(aspect.ratio=1)
dev.off()


cluster1<-table(cancer_scrna@meta.data$seurat_clusters,cancer_scrna@meta.data$annotation)

sc<-as.vector(cancer_scrna@meta.data$seurat_clusters)
sc[sc==3]<- "CT1"
sc[sc==6]<- "CT1"
sc[sc==2]<- "CT2"
sc[sc==5]<- "NA"
sc[sc==4]<- "CT3"
sc[sc==0]<- "CT4"
sc[sc==1]<- "CT5"
sc[sc==7]<- "CT6"

cancer_scrna@meta.data$annotation<-sc
#Idents(cancer_scrna)<-as.factor(sc)
cancer_scrna<-subset(cancer_scrna,subset = annotation %in% c("CT1","CT2","CT3","CT4","CT5","CT6"))
Idents(cancer_scrna)<- factor(cancer_scrna$annotation,levels=c("CT1","CT2","CT3","CT4","CT5","CT6"))
save(cancer_scrna,file="cancer_scrna.RData")

############################
#####2.聚类
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
########################
######3.展示G2M得分
all_fortify_can <- fortify.Seurat(cancer_scrna)

CC_g2m_plot <- ggplot() +
        geom_point(data = all_fortify_can,
                   aes(x = UMAP_1, y = UMAP_2,color =G2M.Score), size = 0.2, alpha = 1) +
        umap_theme() +
        scale_colour_gradient(low="#003638",high="#FFF338")+
        theme(aspect.ratio=1) +
        theme(legend.text=element_markdown(size=14),
              legend.title=element_text(size=14)) +
        guides(colour = guide_legend(override.aes = list(size=3)))  +
        ggtitle("G2M score")

        pdf(file="Epi_G2M_gradient_UMAP.pdf",width = 5, height = 5)
        print(CC_g2m_plot)
        dev.off()

CC_S.Score_plot <- ggplot() +
        geom_point(data = all_fortify_can,
                   aes(x = UMAP_1, y = UMAP_2,color =S.Score), size = 0.2, alpha = 1) +
        umap_theme() +
        scale_colour_gradient(low="#003638",high="#FFF338")+
        theme(aspect.ratio=1) +
        theme(legend.text=element_markdown(size=14),
              legend.title=element_text(size=14)) +
        guides(colour = guide_legend(override.aes = list(size=3)))  +
        ggtitle("S.Score")

        pdf(file="Epi_S.Score_gradient_UMAP.pdf",width = 5, height = 5)
        print(CC_S.Score_plot)
        dev.off()

CC_Phase_plot <- ggplot() +
        geom_point(data = all_fortify_can,
                   aes(x = UMAP_1, y = UMAP_2,color =Phase), size = 0.2, alpha = 1) +
        umap_theme() +
        scale_colour_manual(values=c("#01937C","#B6C867","#FFC074"))+
        theme(aspect.ratio=1) +
        theme(legend.text=element_markdown(size=14),
              legend.title=element_text(size=14)) +
        guides(colour = guide_legend(override.aes = list(size=3)))  +
        ggtitle("Phase")

        pdf(file="Epi_Phase_gradient_UMAP.pdf",width = 5, height = 5)
        print(CC_Phase_plot)
        dev.off()

####################
# inferCNV
#######################
load("/share/pub/dengcy/Singlecell/CC_space/1.scRNA_seq_prepare/CC_scrna.RData")
load("cancer_scrna.RData")
CC_infercnv <- subset(CC_scrna, subset = annotation %in% c("Neutrophils","Macrophage"))
cancer_scrna$annotation <- cancer_scrna$seurat_clusters_tumor
CC_infercnv<-merge(cancer_scrna,CC_infercnv)
#table(CC_scrna_infercnv$seurat_clusters)
#celltype <- cbind(colnames(CC_scrna_infercnv),as.character(CC_scrna_infercnv@meta.data$annotation))

#head(celltype)
mtx <-  as.data.frame(CC_infercnv@assays$RNA@counts)
groupinfo= data.frame(cellId = colnames(mtx))
groupinfo$cellType = CC_infercnv$annotation
# 第三文件
library(AnnoProbe)
geneInfor=annoGene(rownames(mtx),"SYMBOL",'human')
geneInfor=geneInfor[with(geneInfor, order(chr, start)),c(1,4:6)]
geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
##

expFile='expFile.txt'
write.table(mtx ,file = expFile,sep = '\t',quote = F)

groupFiles='groupFiles.txt'
write.table(groupinfo,file = groupFiles,sep = '\t',quote = F,col.names = F,row.names = F)

head(geneInfor)
geneFile='geneFile.txt'
write.table(geneInfor,file = geneFile,sep = '\t',quote = F,col.names = F,row.names = F)
}

library(future)
expFile='expFile.txt' 
groupFiles='groupFiles.txt'  
geneFile='geneFile.txt'
infercnv_obj = infercnv::CreateInfercnvObject(raw_counts_matrix=expFile,
                                    annotations_file=groupFiles,
                                    delim="\t",
                                    gene_order_file= geneFile,
                                    ref_group_names=c("Macrophage","Neutrophils")) # 如果有正常细胞的话，把正常细胞的分组填进去
future::plan("multiprocess",workers=20)# 多核并行处理
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir='infercnv_out/', 
                             cluster_by_groups=TRUE, 
                             denoise=TRUE,
                             HMM=F)# HMM添加后需要跑很长的时间，这里要注意一一下

if(F){
####################
#marker可视化
##################
load("cancer_scrna.RData")

markers= FindAllMarkers(cancer_scrna,logfc.threshold = 0)
#load("annotation_markers.RData")
markers<-markers[!grepl("MT-",markers$gene,ignore.case=F),]
save(markers,file="cancer_markers.RData")
write.csv(markers,file="Epicell_markers.csv")
top_markers = markers %>% group_by(cluster) %>% top_n(n=5,wt=avg_log2FC)  #For Seurat>4.0, set "wt=avg_logFC"
top_markers$cluster<- as.character(top_markers$cluster)
top_markers<-top_markers[order(top_markers$cluster),]

cancer_scrna <- FindVariableFeatures(cancer_scrna,nfeatures = 10000)

pdf("DoHeatmap_Magliant_markers.pdf",width=7)
DoHeatmap(cancer_scrna,features=top_markers$gene,
  angle =0,
  group.colors=color_scanpy_viridis28[c(1,4,5,7,20,11,13,15)],
  size = 3)
dev.off()
###############
#通路观察
################
library(doRNG)
library(testSctpa)

se_oj = CreateSeuratObject(GetAssayData(object=cancer_scrna,assay="RNA",slot="data"))
cancer_scrna_pas = cal_PAS(seurat_object = se_oj,
              tool = 'AUCell',   ## GSVA, ssGSEA, plage, zscore or Vision
              normalize = "log",
              species = 'human', 
              pathway='hallmarker')

Idents(cancer_scrna_pas)<- cancer_scrna$seurat_clusters_tumor
cancer_scrna_pas = FindVariableFeatures(cancer_scrna_pas, verbose = FALSE)
cancer_scrna_pas = ScaleData(cancer_scrna_pas)
#save(cancer_scrna_pas,file="cancer_scrna_hallmarker.RData")
markers_Hallmarker = FindAllMarkers(cancer_scrna_pas,logfc.threshold = 0)
save(markers_Hallmarker,file="cancer_markers_Hallmarker.RData")

pathways = markers_Hallmarker %>% group_by(cluster) %>% top_n(n=3,wt=avg_log2FC)  #For Seurat>4.0, set "wt=avg_logFC"
pdf("DoHeatmap_cancer_hallmarker.pdf",width=8,height=6)
DoHeatmap(cancer_scrna_pas,features=pathways$gene,,
  group.colors=color_scanpy_viridis28[c(1,4,5,7,20,11,13,15)],size = 2)+ggtitle("Hallmarker")
dev.off()

####reactome
cancer_scrna_pas = cal_PAS(seurat_object = se_oj,
              tool = 'AUCell',   ## GSVA, ssGSEA, plage, zscore or Vision
              normalize = "log",
              species = 'human', 
              pathway='reactome')

Idents(cancer_scrna_pas)<-cancer_scrna$seurat_clusters_tumor
cancer_scrna_pas = FindVariableFeatures(cancer_scrna_pas, verbose = FALSE)
cancer_scrna_pas = ScaleData(cancer_scrna_pas)
save(cancer_scrna_pas,file="cancer_scrna_reactome.RData")

markers_reactome = FindAllMarkers(cancer_scrna_pas,logfc.threshold = 0)
save(markers_reactome,file="cancer_markers_reactome.RData")

pathways = markers_reactome %>% group_by(cluster) %>% top_n(n=5,wt=avg_log2FC)  #For Seurat>4.0, set "wt=avg_logFC"
pdf("DoHeatmap_cancer_reactomer.pdf",width=10,height=8)
DoHeatmap(cancer_scrna_pas,features=pathways$gene,,
  group.colors=color_scanpy_viridis28[c(1,4,5,7,20,11,13,15)],size = 2)+ggtitle("Reactome")
dev.off()
############################
#4.slingshot为时间分析，分样本进行分析
########################
save(cancer_scrna,file="cancer_scrna.RData")
 cancer_scrna_cc1<-subset(cancer_scrna,subset= samples %in% "CC1")
 cancer_scrna_cc1<-subset(cancer_scrna_cc1,subset= annotation %in% c("CT1","CT2"))

 cancer_scrna_cc2<-subset(cancer_scrna,subset= samples %in% "CC2")
 cancer_scrna_cc2<-subset(cancer_scrna_cc2,subset= annotation %in% c("CT4","CT5","CT6"))

 cancer_scrna_cc1<-cluster_pca_umap(cancer_scrna_cc1, assay = "SCT",reduction="pca",cluster_res = 0.3)
 cancer_scrna_cc2<-cluster_pca_umap(cancer_scrna_cc2, assay = "SCT",reduction="pca",cluster_res = 0.3)

Idents(cancer_scrna_cc1) <- factor(cancer_scrna_cc1$annotation,levels=c("CT1","CT2"))
pdf(file="Malignant_cc1_umap_0.3.pdf",height=6,width=6)
DimPlot(cancer_scrna_cc1,pt.size=0.9,label = TRUE, repel=TRUE,label.size = 6)+ 
umap_theme()+ labs(x="UMAP",y="")+
ggtitle("CC1 Malignant cells")+
        scale_colour_manual(name = "Cluster", values = color_scanpy_viridis28[c(1,4)]) +
        theme(aspect.ratio=1)
dev.off()

Idents(cancer_scrna_cc2) <- factor(cancer_scrna_cc2$annotation,levels=c("CT4","CT5","CT6"))
pdf(file="Malignant_cc2_umap_0.3.pdf",height=6,width=6)
DimPlot(cancer_scrna_cc2,pt.size=0.9,label = TRUE, repel=TRUE,label.size = 6)+ 
umap_theme()+ labs(x="UMAP",y="")+
ggtitle("CC2 Malignant cells")+
        scale_colour_manual(name = "Cluster", values = color_scanpy_viridis28[c(7,20,11)]) +
        theme(aspect.ratio=1)
dev.off()

save(cancer_scrna_cc1,file="cancer_scrna_cc1.RData")
save(cancer_scrna_cc2,file="cancer_scrna_cc2.RData")

all_fortify_can1 <- fortify.Seurat(cancer_scrna_cc1)
all_fortify_can2 <- fortify.Seurat(cancer_scrna_cc2)

library(slingshot)
library(SingleCellExperiment)
cancer_count1<-SingleCellExperiment(assays = List(counts = GetAssayData(object =cancer_scrna_cc1, slot = "counts")))
cancer_count2<-SingleCellExperiment(assays = List(counts = GetAssayData(object =cancer_scrna_cc2, slot = "counts")))

###标准化处理

FQnorm <- function(counts){
  rk <- apply(counts,2,rank,ties.method='min')
  counts.sort <- apply(counts,2,sort)
  refdist <- apply(counts.sort,1,median)
  norm <- apply(rk,2,function(r){
    refdist[r]
  })
  rownames(norm) <- rownames(counts)
  return(norm)
}
assays(cancer_count1)$norm <- FQnorm(assays(cancer_count1)$counts)
assays(cancer_count2)$norm <- FQnorm(assays(cancer_count2)$counts)

rd1<-all_fortify_can1[,c("UMAP_1" , "UMAP_2")]
rd2<-all_fortify_can2[,c("UMAP_1" , "UMAP_2")]

#rd1<-all_fortify_can[,c("UMAP_1" , "UMAP_2")]
library(mclust, quietly = TRUE)
cl1 <- Mclust(rd1)$classification
colData(cancer_count1)$GMM <- cl1

cl2 <- Mclust(rd2)$classification
colData(cancer_count2)$GMM <- cl2

cancer_count1<-  SingleCellExperiment(assays = List(counts = assays(cancer_count1)$counts,
                                                   norm=assays(cancer_count1)$norm),
                                     reducedDims = SimpleList(UMAP = as.matrix(rd1)),
                                     colData = data.frame(clus = cl1))
colData(cancer_count1)$GMM <- cl1

cancer_count2<-  SingleCellExperiment(assays = List(counts = assays(cancer_count2)$counts,
                                                   norm=assays(cancer_count2 )$norm),
                                     reducedDims = SimpleList(UMAP = as.matrix(rd2)),
                                     colData = data.frame(clus = cl2))
colData(cancer_count2)$GMM <- cl2
##运行slingshot
  Slingshot_scrna1 <- slingshot(cancer_count1, 
                  clusterLabels = 'GMM', 
                  reducedDim = 'UMAP'
                  )

  save(Slingshot_scrna1,file = "Slingshot_scrna1.RData")
  Slingshot_scrna2 <- slingshot(cancer_count2, 
                  clusterLabels = 'GMM', 
                  reducedDim = 'UMAP'
                  )

  save(Slingshot_scrna2,file = "Slingshot_scrna2.RData")
###可视化
  library(grDevices)
  library(scater)
  library(scran)
  library(RColorBrewer)
  library(ggthemes)

  colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
  pdf("Slingshot_cc1.plot.pdf")
  plotcol <- colors[cut(Slingshot_scrna1$slingPseudotime_1, breaks=100)]
  plot(reducedDims(Slingshot_scrna1)$UMAP, col = plotcol, pch=16, asp = 1)
  lines(SlingshotDataSet(Slingshot_scrna1), lwd=2, col='black')
  dev.off()

  pdf("Slingshot_cc2.plot.pdf")
  plotcol <- colors[cut(Slingshot_scrna2$slingPseudotime_1, breaks=100)]
  plot(reducedDims(Slingshot_scrna2)$UMAP, col = plotcol, pch=16, asp = 1)
  lines(SlingshotDataSet(Slingshot_scrna2), lwd=2, col='black')
  dev.off()

  embedded1 <- embedCurves(Slingshot_scrna1, "UMAP")
  embedded1.1 <- slingCurves(embedded1)[[1]]
  embedded1.1 <- data.frame(embedded1.1$s[embedded1.1$ord,])
  embedded1.2 <- slingCurves(embedded1)[[2]]
  embedded1.2 <- data.frame(embedded1.2$s[embedded1.2$ord,])


  embedded2 <- embedCurves(Slingshot_scrna2, "UMAP")
  embedded2.1 <- slingCurves(embedded2)[[1]]
  embedded2.1 <- data.frame(embedded2.1$s[embedded2.1$ord,])

  embedded2.2 <- slingCurves(embedded2)[[2]]
  embedded2.2 <- data.frame(embedded2.2$s[embedded2.2$ord,])

  save(embedded1.1,embedded1.2,embedded2.1,embedded2.2,file="/share/pub/dengcy/Singlecell/CC_space/3.Epi_cells_stage/Slingshot_embedded.RData")

  gg_UMAP1 <- as.data.frame(reducedDim(Slingshot_scrna1, "UMAP"))
  colnames(gg_UMAP1)<-c("UMAP1","UMAP2")
  gg_UMAP1$Type<-all_fortify_can1$seurat_clusters_new1
  
  df<-data.frame(slingPseudotime_1=Slingshot_scrna1$slingPseudotime_1,
  slingPseudotime_2=Slingshot_scrna1$slingPseudotime_2,
  slingPseudotime_3=Slingshot_scrna1$slingPseudotime_3,
  slingPseudotime_4=Slingshot_scrna1$slingPseudotime_4)

  df1<-apply(df, 1, function(x){mean(x,na.rm =T)})
  gg_UMAP1$Pseudotime<- df1
  cancer_scrna_cc1$Pseudotime<-df1
  save(cancer_scrna_cc1,file="cancer_scrna_cc1.RData")

  library(ggplot2)
  require(ggpubr)
  
  pt_plot<-ggplot()+
    geom_point(data = gg_UMAP1, aes(x = UMAP1, y = UMAP2,color=Pseudotime),size = .5) +
    umap_theme() +
    #scale_colour_manual(name = "Timepoint", values = color_scanpy_default) +
    theme(aspect.ratio=1) +
    theme(legend.text=element_markdown(size=16),
          legend.title=element_text(size=16)) +
    guides(colour = guide_legend(override.aes = list(size=3))) +
    #new_scale_color() +
    scale_color_gradient(low = "#EBE645",high = "#000957")+
  geom_path(data=embedded1.1, aes(x=UMAP_1, y=UMAP_2), size=0.8,color="black",alpha=1)
    #geom_path(data=embedded1.2, aes(x=UMAP_1, y=UMAP_2), size=0.8,color="black",alpha=1)
      #geom_path(data=embedded3, aes(x=UMAP_1, y=UMAP_2), size=0.8,color="black",alpha=1)

  pdf("/share/pub/dengcy/Singlecell/CC_space/3.Epi_cells_stage/ptime_umap_cc1.pdf",height = 5,width = 5)
  pt_plot
  dev.off()

  gg_UMAP2 <- as.data.frame(reducedDim(Slingshot_scrna2, "UMAP"))
  colnames(gg_UMAP2)<-c("UMAP1","UMAP2")
  gg_UMAP2$Type<-all_fortify_can2$seurat_clusters_new1
  
  df<-data.frame(slingPseudotime_1=Slingshot_scrna2$slingPseudotime_1,
    slingPseudotime_2=Slingshot_scrna2$slingPseudotime_2,
    slingPseudotime_3=Slingshot_scrna2$slingPseudotime_3
    )
  df2<-apply(df, 1, function(x){mean(x,na.rm =T)})
  gg_UMAP2$Pseudotime<- df2
  cancer_scrna_cc2$Pseudotime<-df2
  save(cancer_scrna_cc2,file="cancer_scrna_cc2.RData")

  pt_plot2<-ggplot()+
    geom_point(data = gg_UMAP2, aes(x = UMAP1, y = UMAP2,color=Pseudotime),size = .5) +
    umap_theme() +
    #scale_colour_manual(name = "Timepoint", values = color_scanpy_default) +
    theme(aspect.ratio=1) +
    theme(legend.text=element_markdown(size=16),
          legend.title=element_text(size=16)) +
    guides(colour = guide_legend(override.aes = list(size=3))) +
    #new_scale_color() +
    scale_color_gradient(low = "#EBE645",high = "#000957")+
  geom_path(data=embedded2.1, aes(x=UMAP_1, y=UMAP_2), size=0.8,color="black",alpha=1)+
    geom_path(data=embedded2.2, aes(x=UMAP_1, y=UMAP_2), size=0.8,color="black",alpha=1)
      #geom_path(data=embedded3, aes(x=UMAP_1, y=UMAP_2), size=0.8,color="black",alpha=1)

  pdf("/share/pub/dengcy/Singlecell/CC_space/3.Epi_cells_stage/ptime_umap_cc2.pdf",height = 5,width =5)
  pt_plot2
  dev.off()

 #######4.展示不同样本区域
####重定义聚类名字
###########################################remove
if(F){
all_fortify_can$seurat_clusters_new1<- as.character(as.vector(all_fortify_can$seurat_clusters))
all_fortify_can$seurat_clusters_new1[all_fortify_can$seurat_clusters_new1 == "5"]<-"CT1"
all_fortify_can$seurat_clusters_new1[all_fortify_can$seurat_clusters_new1 == "8"]<-"CT2"
all_fortify_can$seurat_clusters_new1[all_fortify_can$seurat_clusters_new1 == "2"]<-"CT3"
cancer_scrna$seurat_clusters_new1<-all_fortify_can$seurat_clusters_new1
 save(cancer_scrna,file="/share/pub/dengcy/Singlecell/CC_space/3.Epi_cells_stage/cancer_scrna.RData")
 CC_Onlyc1_plot <- ggplot() +
   geom_point(data = all_fortify_can[!(all_fortify_can$seurat_clusters_new1 %in% c("CT1","CT2","CT3")),],
              aes(x = UMAP_1, y = UMAP_2), size = 0.2, alpha = 0.8, color = "gray") +
   umap_theme() +
   #scale_colour_manual(name = "Timepoint", values = color_scanpy_default) +
   theme(aspect.ratio=1) +
   theme(legend.text=element_markdown(size=14),
         legend.title=element_text(size=14)) +
   guides(colour = guide_legend(override.aes = list(size=3))) +
   new_scale_color() +
   
   geom_point(data = all_fortify_can[(all_fortify_can$seurat_clusters_new1 %in% c("CT1","CT2","CT3")),], aes(x = UMAP_1, y = UMAP_2,color=seurat_clusters_new1), size = .2) +
   scale_colour_manual(name = "Only C1", values = color_scanpy_13[1:3], guide = guide_legend(override.aes = list(size=3), order = 1)) +
   theme(aspect.ratio=1,
         legend.text=element_markdown(size=16),
         legend.title=element_text(size=16)) +
   new_scale_color()+
   ggtitle("Only C1 cluster")
 
 pdf(file="Epi_Onlyc1_highlight_UMAP.pdf",width = 5, height = 5)
 print(CC_Onlyc1_plot)
 dev.off()
 
 CC_OnlyR2_plot <- ggplot() +
   geom_point(data = all_fortify_can[!(all_fortify_can$seurat_clusters_new1 %in% c("RT1","RT2","RT3","RT4","RT5")),],
              aes(x = UMAP_1, y = UMAP_2), size = 0.2, alpha = 0.8, color = "gray") +
   umap_theme() +
   #scale_colour_manual(name = "Timepoint", values = color_scanpy_default) +
   theme(aspect.ratio=1) +
   theme(legend.text=element_markdown(size=16),
         legend.title=element_text(size=16)) +
   guides(colour = guide_legend(override.aes = list(size=3))) +
   new_scale_color() +
   
   geom_point(data = all_fortify_can[(all_fortify_can$seurat_clusters_new1 %in% c("RT1","RT2","RT3","RT4","RT5")),], aes(x = UMAP_1, y = UMAP_2,color=seurat_clusters_new1),size = .2) +
   scale_colour_manual(name = "Only R2", values = color_scanpy_13[4:8], guide = guide_legend(override.aes = list(size=3), order = 1)) +
   theme(aspect.ratio=1,
         legend.text=element_markdown(size=14),
         legend.title=element_text(size=14)) +
   new_scale_color()+
   ggtitle("Only R2 cluster")
 
 pdf(file="Epi_onlyR2_highlight_UMAP.pdf",width = 5, height = 5)
 print(CC_OnlyR2_plot)
 dev.off()

CC_c1R2_plot <- ggplot() +
   geom_point(data = all_fortify_can[!(all_fortify_can$seurat_clusters_new1 %in% c("CRT1","CRT2","CRT3","Fibroblasts/stem_cells","cancer_stem_cells")),],
              aes(x = UMAP_1, y = UMAP_2), size = 0.2, alpha = 0.8, color = "gray") +
   umap_theme() +
   #scale_colour_manual(name = "Timepoint", values = color_scanpy_default) +
   theme(aspect.ratio=1) +
   theme(legend.text=element_markdown(size=16),
         legend.title=element_text(size=16)) +
   guides(colour = guide_legend(override.aes = list(size=3))) +
   new_scale_color() +
   
   geom_point(data = all_fortify_can[(all_fortify_can$seurat_clusters_new1 %in% c("CRT1","CRT2","CRT3","Fibroblasts/stem_cells","cancer_stem_cells")),], aes(x = UMAP_1, y = UMAP_2,color=seurat_clusters_new1),size = .2) +
   scale_colour_manual(name = "C1&R2", values = color_scanpy_13[9:13], guide = guide_legend(override.aes = list(size=3), order = 1)) +
   theme(aspect.ratio=1,
         legend.text=element_markdown(size=14),
         legend.title=element_text(size=14)) +
   new_scale_color()+
   ggtitle("C1&R2 cluster")
 
 pdf(file="Epi_c1R2_highlight_UMAP.pdf",width = 5, height = 5)
 print(CC_c1R2_plot)
 dev.off()

CC_c1R2_plot2 <- ggplot() +
   geom_point(data = all_fortify_can[!(all_fortify_can$seurat_clusters_new1 %in% c("CRT1","CRT2","CRT3","Fibroblasts/stem_cells","cancer_stem_cells")),],
              aes(x = UMAP_1, y = UMAP_2), size = 0.2, alpha = 0.8, color = "gray") +
   umap_theme() +
   #scale_colour_manual(name = "Timepoint", values = color_scanpy_default) +
   theme(aspect.ratio=1) +
   theme(legend.text=element_markdown(size=16),
         legend.title=element_text(size=16)) +
   guides(colour = guide_legend(override.aes = list(size=3))) +
   new_scale_color() +
   
   geom_point(data = all_fortify_can[(all_fortify_can$seurat_clusters_new1 %in% c(3,4,7,9,12)),], aes(x = UMAP_1, y = UMAP_2,color=samples),size = .2) +
   scale_colour_manual(name = "C1&R2", values =c("#DF711B","#64C9CF"), guide = guide_legend(override.aes = list(size=3), order = 1)) +
   theme(aspect.ratio=1,
         legend.text=element_markdown(size=14),
         legend.title=element_text(size=14)) +
   new_scale_color()+
   ggtitle("C1&R2 cluster")
 
 pdf(file="Epi_c1R2_highlight_UMAP2.pdf",width = 5, height = 5)
 print(CC_c1R2_plot2)
 dev.off()


cancer_scrna$seurat_clusters_new1<-factor(all_fortify_can$seurat_clusters_new1,levels=c("CT1","CT2","CT3","RT1","RT2","RT3","RT4","RT5","CRT1","CRT2","CRT3","Fibroblasts/stem_cells","cancer_stem_cells"))

 plots_clusters_2 <- DimPlot(cancer_scrna,group.by = "seurat_clusters_new1",pt.size=0.8,label = TRUE, repel=TRUE)+ umap_theme()+ 
 ggtitle("Epithelial and stem cells")+
  scale_colour_manual(name = "Cluster", values = color_scanpy_13) +
  theme(aspect.ratio=1)

pdf(file="Epi_stem_newcluster_umap_0.5.pdf",height=7,width=7)
print(plots_clusters_2)
dev.off()
}

 save(cancer_scrna,file="cancer_scrna.RData")

############################
###5.计算类别之间的marker
Idents(cancer_scrna)<-cancer_scrna$seurat_clusters_new1

plan("multiprocess", workers = N_WORKERS)

  cancer_markers <- parallelFindAllMarkers(cancer_scrna)
  names(cancer_markers)<-levels(Idents(cancer_scrna))

  save(cancer_markers,file="cancer_markers.RData")

##################################
#7.样本比例柱状图 _new
load("cancer_scrna.RData")
#cancer_scrna$seurat_clusters<- as.vector(paste0("T",cancer_scrna$seurat_clusters))

pdf("stack_plot_epicell_samples.pdf",width=3,height=5)
stack_plot(obj=cancer_scrna,options=c("seurat_clusters_tumor","samples"),colors= color_scanpy_viridis28[c(1,4,5,7,20,11,13,15)])
dev.off()

}

