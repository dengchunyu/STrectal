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
#"E:/OneDrive/colorectal_sc/CC_space/src/Functions.r"
#"E:/OneDrive/colorectal_sc/CC_space/3.Epi_cells_stage"
setwd("/share/pub/dengcy/Singlecell/CC_space/3.2.Immune_cells_stage")
load("/share/pub/dengcy/Singlecell/CC_space/1.scRNA_seq_prepare/CC_scrna.RData")

Immune_scrna<-subset(CC_scrna,annotation %in% c("B_cell","Macrophage","NK_cells","T_cells","Neutrophils"))
Immune_scrna<-cluster_pca_umap(Immune_scrna, assay = "SCT",reduction="pca",cluster_res = 0.5)

Immune_scrna$annotation
#cluster1<-table(Immune_scrna@meta.data$seurat_clusters)
Immune_scrna$annotation2<-paste0(Immune_scrna$annotation,"_",Immune_scrna$seurat_clusters)

###过滤低细胞群体
cluster_a1<- table(Immune_scrna$seurat_clusters)
Immune_scrna<- subset(Immune_scrna, subset = seurat_clusters %in% names(cluster_a1)[which(cluster_a1>50)])
#26496  9670
Immune_scrna$seurat_clusters <- as.numeric(Immune_scrna$seurat_clusters)

#CC_scrna$singler <- pred.brca$labels

library("celldex")
library("SingleR")
library("pracma")
ref <- HumanPrimaryCellAtlasData() 
#############
common <- intersect(rownames(ref), rownames(Immune_scrna))
Immune_scrna.singler <- Immune_scrna[common,]
ref <- ref[common,]
Immune_scrna.singler.sce <- SingleCellExperiment(assays = list(counts = Immune_scrna.singler))
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
write.csv(singler.final,file="singler.final.csv")

#####绘制注释后的细胞图
#
seurat_clusters<-as.vector(Immune_scrna$seurat_clusters)
annotation<-as.vector(seurat_clusters)

annotation[seurat_clusters == 2]<-"CD4.TCM"
annotation[seurat_clusters == 3]<-"Cytotoxic.CD8.T"
annotation[seurat_clusters == 4]<-"Terminally.Exhausted.T"
annotation[seurat_clusters == 6]<-"CD4.TCM"
annotation[seurat_clusters == 1]<-"Neutrophils"
annotation[seurat_clusters == 5]<-"NK.cell"
annotation[seurat_clusters == 7] <-"TAM"
annotation[seurat_clusters == 9]<-"TAM"
annotation[seurat_clusters == 8] <-"B.cell"
annotation[seurat_clusters == 10] <-"B.cell"
annotation[seurat_clusters == 11]<-"B.cell.activate.memory"
annotation[seurat_clusters == 12]<-"B.cell"

Immune_scrna$annotation<-annotation
############################
#####2.聚类
Immune_scrna$annotation<- factor(Immune_scrna$annotation,levels= c("CD4.TCM","Cytotoxic.CD8.T","Terminally.Exhausted CD8T","B.cell","B.cell.activite.memory","TAM","NK.cell","Neutrophils"))
Idents(Immune_scrna) <- Immune_scrna$annotation

pdf(file="Immune_cells_TSNE_0.5.pdf",height=8,width=8)
DimPlot(Immune_scrna,group.by = "annotation",pt.size=0.8,label = TRUE,reduction="tsne", repel=TRUE)+ 
umap_theme()+ ggtitle("Immune cells")+labs(x="TSNE",y="")+
        scale_colour_manual(name = "Cluster", values = color_scanpy_viridis28[18:28]) +
        theme(aspect.ratio=1)
dev.off()


#############
#差异基因热图
###################
markers= FindAllMarkers(Immune_scrna,logfc.threshold = 0)
markers<-markers[!grepl("MT-",markers$gene,ignore.case=F),]

save(markers,file="Immune_markers2.0.RData")
write.csv(markers,file="Immune_markers.csv")
#
save(Immune_scrna,file="Immune_scrna.RData")
top_markers = markers %>% group_by(cluster) %>% top_n(n=5,wt=avg_log2FC)  #For Seurat>4.0, set "wt=avg_logFC"

pdf("DotPlot_immuneMarkers.pdf",width=15,height=5)
DotPlot(Immune_scrna, features =top_markers$gene, cols = c("lightgrey", "#CD1818"))+ RotatedAxis()
dev.off()

Lineage<-c("CD3E","CD4","CD8A")
 Naive_memory<- c("SELL","CCR7","IL7R","FAS","CD27","CD28","ITGAE","KLRB1")
 Inhibitory<-c("PDCD1","TIGIT","HAVCR2","LAG3","CTLA4","TNFRSF14","BTLA","CD244","CD160")
Ectonucleases<-c("CD38","ENTPD1","NT5E")
 Activation<-c("CD69","IL2RA","ICOS","TNFRSF4","TNFRSF9","HLA-DRA","FASLG","CD40LG")
Cytolytic<-c("GZMA","GZMB","GZMH","GZMK","GZMM","PRF1","NKG7","GNLY","LAMP1")
 Cytokines<-c("IL2","IFNG","TNF")

pdf("FeaturePlot_Lineagegene.pdf",height=6)
FeaturePlot(Immune_scrna, features = Lineage,reduction="tsne",cols=c("#FAEEE0","#F14A16"))
dev.off()

pdf("FeaturePlot_Naive_memorygene.pdf",width=10,height=8)
FeaturePlot(Immune_scrna, features = Naive_memory,reduction="tsne",cols=c("#FAEEE0","#F14A16"))
dev.off()

pdf("FeaturePlot_Inhibitorygene.pdf",width=10,height=8)
FeaturePlot(Immune_scrna, features = Inhibitory,reduction="tsne",cols=c("#FAEEE0","#F14A16"))
dev.off()

pdf("FeaturePlot_Ectonucleasesgene.pdf")
FeaturePlot(Immune_scrna, features = Ectonucleases,reduction="tsne",cols=c("#FAEEE0","#F14A16"))
dev.off()

pdf("FeaturePlot_Activationgene.pdf",width=10,height=8)
FeaturePlot(Immune_scrna, features = Activation,reduction="tsne",cols=c("#FAEEE0","#F14A16"))
dev.off()

pdf("FeaturePlot_Cytolyticgene.pdf",width=10,height=8)
FeaturePlot(Immune_scrna, features = Cytolytic,reduction="tsne",cols=c("#FAEEE0","#F14A16"))
dev.off()

pdf("FeaturePlot_Cytokinesgene.pdf",height=6)
FeaturePlot(Immune_scrna, features = Cytokines,reduction="tsne",cols=c("#FAEEE0","#F14A16"))
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
save(Immune_scrna_pas,file="Immune_scrna_Immu.RData")

markers_Immu = FindAllMarkers(Immune_scrna_pas,logfc.threshold = 0)
save(markers_Immu,file="markers_Immu.RData")

pathways = markers_Immu %>% group_by(cluster) %>% top_n(n=3,wt=avg_log2FC)  #For Seurat>4.0, set "wt=avg_logFC"

Idents(Immune_scrna_pas) <- factor(as.vector(Idents(Immune_scrna_pas)),levels=sort(unique(as.vector(Idents(Immune_scrna_pas)))))

pdf("DoHeatmap_immunecells_immu.pdf",width=12,height=6)
DoHeatmap(Immune_scrna_pas,features=pathways$gene,,
  group.colors=color_scanpy_viridis28[18:28],size = 2)+ggtitle("ImmunePathways")
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
  group.colors=color_scanpy_viridis28[18:28],size = 2)+ggtitle("ReactomePathways")
dev.off()

#####################绘制免疫细胞比例图
load("Immune_scrna.RData")
load("Immune_markers.RData")

library(plyr)
library(ggplot2)


pdf("Immunecellproportion_in_samples2.pdf",width=8,height=3)
stack_plot(obj=Immune_scrna,options=c("annotation","samples"),colors= color_scanpy_viridis28[18:28],angle=FALSE)
dev.off()

save(Immune_scrna,file="Immune_scrna.RData")
##########################促炎基因和抗炎基因集合在不同样本，不同免疫细胞中的表达情况

#################两个得分的umap图

T_scrna<-subset(Immune_scrna,annotation %in% c("CD4.TCM","Cytotoxic.CD8.T","Terminally.Exhausted CD8T"))
Anti_inflammatory<-c("TIGIT","IDO1","LGALS3","PDCD1","FOXP3","ENTPD1","CD274","CSF2","CTLA4","CXCL12","CXCL5","IL8","MIF","PTGS2","VEGFA")
Pro_inflammatory<-c("IL1A","IL1B","TNF","IFNG","TBX21","CCL3","CCL4","PRF1","GZMA","GZMB","GZMK","GZMH","CD8A","FASLG","CCL2","CCL20","IL2","IL6","IL12a","IL17a","IL23a","PTGS2","TLR4","TNF")

Anti_inflammatory<-intersect(Anti_inflammatory,rownames(Immune_scrna))
Pro_inflammatory<-intersect(Pro_inflammatory,rownames(Immune_scrna))

T_scrna<-AddModuleScore(T_scrna,list(Anti_inflammatory,Pro_inflammatory),name=c("Anti_inflammatory","Pro_inflammatory"))

T_fortify_can<-fortify.Seurat(T_scrna)

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

