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
library(future)
options(future.globals.maxSize = 10000 * 1024^2)
source("/share/pub/dengcy/Singlecell/CC_space/src/Functions.r")
source("/share/pub/dengcy/Singlecell/CC_space/src/library.r")
setwd("/share/pub/dengcy/Singlecell/CC_space/1.scRNA_seq_prepare")

file_pre="/share/pub/dengcy/Singlecell/CC_space/gcy_analysis/Analysis/1.Basic_analysis/1.2.filtered_feature_bc_matrix"
output_pre="/share/pub/dengcy/Singlecell/CC_space/1.scRNA_seq_prepare/"

#tar czvf CRC_scrnaseq.tar.gz /share/pub/dengcy/Singlecell/CC_space/gcy_analysis/Analysis/1.Basic_analysis/1.2.filtered_feature_bc_matrix

if(F){
## human or mouse: cc_genes_hg.txt or cc_genes_hg.txt
cell_cycle_genes <- read.table(file = "/share/pub/xiongyc/project/scRNA/JiangFanChen/data/cc_genes_hg.txt",header = TRUE);
############################################
############ step1: data prepration
##########################################
  #### step1.1: create Seurat object
  data.10x<- Read10X(data.dir = paste(file_pre,sep="")); 
  #CRC_scrna = CreateSeuratObject(counts = data.10x, min.cells=0, min.features=0, project=projectName[file_pre]);
  CC_scrna = CreateSeuratObject(counts = data.10x,min.cells=3, min.features=200,plot_feature=F,
  cell_cycle_genes=cell_cycle_genes, project="CC_scRNA-seq");
  #28051 features across 9763 samples within 1 assay

  #### step1.2: mitochondrial percentage and rRNA percentage
  CC_scrna[["percent.mt"]] <- PercentageFeatureSet(CC_scrna, pattern = "^mt-")
  CC_scrna[["percent.ribo"]] <- PercentageFeatureSet(CC_scrna, pattern = "^Rp[sl][[:digit:]]")
  samples<-unlist(str_sub(rownames(CC_scrna@meta.data),18,18))
  samples[which(samples=="1")]<-"CC1"
  samples[which(samples=="2")]<-"CC2"

  CC_scrna@meta.data$samples <- samples
table( CC_scrna@meta.data$samples)
# CC1  CC2 
#3968 5795
  #### QC cutoff
  #nCount_RNA: the number of reads
  CC_scrna <- subset(CC_scrna, subset = percent.mt < 20 & nCount_RNA > 300 )
   s.genes = cc.genes$s.genes
   g2m.genes = cc.genes$g2m.genes
   CC_scrna <- NormalizeData(CC_scrna, normalization.method = "LogNormalize", scale.factor = 10000)
   CC_scrna <- CellCycleScoring(CC_scrna,  s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  
   #### normalization: 'nFeature_RNA', 'nCount_RNA',"percent.mt","percent.ribo","S.Score", "G2M.Score"
   CC_scrna=SCTransform(CC_scrna, vars.to.regress = c('nFeature_RNA', 'nCount_RNA',"percent.mt","percent.ribo","S.Score", "G2M.Score"), verbose = FALSE)
  #
#####1.4 pca
CC_scrna <- FindVariableFeatures(CC_scrna,nfeatures = 7000)
CC_scrna <- RunPCA(object = CC_scrna, assay = "SCT", npcs = 50)

############################################
############ step2: data cluster,并且初步展示样本和聚类类别
##########################################
#CC_scrna<-cluster_pca_umap(CC_scrna, assay = "SCT",reduction="pca",cluster_res = 0.5)
CC_scrna<-cluster_pca_umap(CC_scrna, assay = "SCT",reduction="pca",cluster_res = 0.3)

cluster_a1<- table(CC_scrna$seurat_clusters)

write.csv(cluster_a1,file="CC_scrna_raw_seurat_clusters.csv")
#26496  9763
 CC_scrna_filter <- subset(CC_scrna, subset = seurat_clusters %in% names(cluster_a1)[which(cluster_a1>100)])
#26496  9670
CC_scrna_filter$seurat_clusters <- as.numeric(CC_scrna_filter$seurat_clusters)
#cluster_a2<-table(CC_scrna$seurat_clusters)

#write.csv(cluster_a2,file="CC_scrna_filter_seurat_clusters.csv")
#########
#figure1.1展示聚类cluster
plots_clusters <- DimPlot(CC_scrna_filter,group.by = "seurat_clusters",pt.size=0.8,label = TRUE, repel=TRUE)+ 
umap_theme()+ ggtitle("Colon Cancer clusters")+labs(x="TSNE",y="")+
        scale_colour_manual(name = "Cluster", values = color_scanpy_viridis28[1:17]) +
        theme(aspect.ratio=1)

pdf(file=paste(output_pre,"allClusters_umap_0.3.pdf",sep=""),height=8,width=8)
print(plots_clusters)
dev.off()

plots_clusters2 <- DimPlot(CC_scrna_filter,group.by = "seurat_clusters", pt.size=0.8,reduction="tsne",label = TRUE, repel=TRUE)+ 
umap_theme()+ ggtitle("Colon Cancer clusters")+
        scale_colour_manual(name = "Cluster", values = color_scanpy_viridis28[1:17]) +labs(x="TSNE",y="")+
        theme(aspect.ratio=1)
pdf(file=paste(output_pre,"allClusters_tsne_0.3.pdf",sep=""),height=8,width=8)
print(plots_clusters2)
dev.off()

########
#figure1.2映射展示不同样本
plots_samples <- DimPlot(CC_scrna_filter,group.by = "samples",pt.size=0.5,reduction="tsne",label = F, repel=TRUE)+ 
umap_theme()+ ggtitle("Colon Cancer Samples")+labs(x="TSNE",y="")+
        scale_colour_manual(name = "Samples", values = c("#E05D5D","#00A19D")) +
        theme(aspect.ratio=1)

pdf(file=paste(output_pre,"allsamples_tsne_0.3.pdf",sep=""),height=8,width=8)
print(plots_samples)
dev.off()
#save(CC_scrna_filter,file="CC_scrna.RData")

############################################
############ step3: data samples vivalize，这里进一步凸显不同样本位置，去掉
##########################################
#########分别展示不同的样本位置：
all_fortify_cc <- fortify.Seurat(CC_scrna)

CC_c1_plot <- ggplot() +
        geom_point(data = all_fortify_crc[!(all_fortify_crc$samples %in% c("C1")),],
                   aes(x = UMAP_1, y = UMAP_2), size = 0.1, alpha = 0.8, color = "gray") +
        umap_theme() +
        #scale_colour_manual(name = "Timepoint", values = color_scanpy_default) +
        theme(aspect.ratio=1) +
        theme(legend.text=element_markdown(size=14),
              legend.title=element_text(size=14)) +
        guides(colour = guide_legend(override.aes = list(size=3))) +
        new_scale_color() +

        geom_point(data = all_fortify_crc[(all_fortify_crc$samples %in% c("C1")),], aes(x = UMAP_1, y = UMAP_2),color = "black", size = .1) +
        #scale_colour_manual(name = "Epithelium", values = epi_colors, guide = guide_legend(override.aes = list(size=3), order = 1)) +
        theme(aspect.ratio=1,
              legend.text=element_markdown(size=14),
              legend.title=element_text(size=14)) +
        new_scale_color()+
        ggtitle("CC C1")

        pdf(file="CC_c1_highlight_UMAP.pdf",width = 8, height = 8)
        print(CC_c1_plot)
        dev.off()

CC_R2_plot <- ggplot() +
        geom_point(data = all_fortify_crc[!(all_fortify_crc$samples %in% c("R2")),],
                   aes(x = UMAP_1, y = UMAP_2), size = 0.1, alpha = 0.8, color = "gray") +
        umap_theme() +
        #scale_colour_manual(name = "Timepoint", values = color_scanpy_default) +
        theme(aspect.ratio=1) +
        theme(legend.text=element_markdown(size=14),
              legend.title=element_text(size=14)) +
        guides(colour = guide_legend(override.aes = list(size=3))) +
        new_scale_color() +

        geom_point(data = all_fortify_crc[(all_fortify_crc$samples %in% c("R2")),], aes(x = UMAP_1, y = UMAP_2),color = "black", size = .1) +
        #scale_colour_manual(name = "Epithelium", values = epi_colors, guide = guide_legend(override.aes = list(size=3), order = 1)) +
        theme(aspect.ratio=1,
              legend.text=element_markdown(size=14),
              legend.title=element_text(size=14)) +
        new_scale_color()+
        ggtitle("CC R2")

        pdf(file="CC_R2_highlight_UMAP.pdf",width = 8, height = 8)
        print(CC_R2_plot)
        dev.off()


#########################################
#### step4: find cluster biomarkers
######################################
#利用MAST方法计算差异分析,并保存结果
#load("CC_scrna.RData")
CC_scrna<-CC_scrna_filter
rm(CC_scrna_filter)
plan("multiprocess", workers = N_WORKERS)
  Allmarkers <- parallelFindAllMarkers(CC_scrna)
 save(Allmarkers,file="All_cluster_markers.RData")
#Allmarkers[[0+15]]<-NULL
Allmarkers<-lapply(0:(length(Allmarkers)-1),function(x){
  df<-Allmarkers[[x+1]]
  df<-df[order(-df$avg_log2FC),]
  df$cluster<-rep(x,nrow(df))
  return(df)
  })
Allmarkers<-Reduce(function(dtf1, dtf2) rbind(dtf1, dtf2),Allmarkers)
write.csv(Allmarkers,file="17cluster_Allmarkers.csv")

#head(Allmarkers[[0 + 1]][order(-Allmarkers[[0 + 1]]$avg_log2FC),])
###################singleR
library("celldex")
library("SingleR")
library("pracma")
library(SingleCellExperiment)
ref <- HumanPrimaryCellAtlasData() 
#############er
common <- intersect(rownames(ref), rownames(CC_scrna))
CC_scrna.singler <- CC_scrna[common,]
ref <- ref[common,]
CC_scrna.singler.sce <- SingleCellExperiment(assays = list(counts = CC_scrna.singler))
CC_scrna.singler.sce <-SingleCellExperiment(assays = List(counts = GetAssayData(object =CC_scrna.singler, slot = "counts")))
CC_scrna.singler.sce <- logNormCounts(CC_scrna.singler.sce)

pred.brca <- SingleR(test = CC_scrna.singler.sce, 
  ref = ref, assay.type.test=1,labels = ref$label.main)
save(pred.brca,file="singler_pred.RData")


singler.results <- merge(data.frame(cell = rownames(pred.brca), singler = pred.brca$labels), 
                         data.frame(cell = rownames(CC_scrna@meta.data), 
                                    cluster = CC_scrna@meta.data$seurat_clusters), 
                         by = "cell", 
                         all.y = FALSE)
singler.results$cell <- NULL
singler.results$count <- 1
singler.results <- aggregate(count ~ ., singler.results, FUN = sum)
singler.final <- singler.results %>% group_by(cluster) %>% top_n(n = 1, wt = count)
write.csv(singler.final,file="singler.final.csv")
#CC_scrna$singler <- pred.brca$labels

#####绘制注释后的细胞图
#
seurat_clusters<-as.vector(CC_scrna$seurat_clusters)
annotation<-as.vector(seurat_clusters)

annotation[seurat_clusters %in% c(1 ,10)]<-"T_cells"
annotation[seurat_clusters == 8]<-"NK_cells"
annotation[seurat_clusters %in% c(2,3, 5,9,12,13,14,15)]<-"Epithelial_cells"
annotation[seurat_clusters == 6]<-"Macrophage"
annotation[seurat_clusters == 11]<-"Cancer_stem_cells"
annotation[seurat_clusters == 4]<-"Neutrophils"

annotation[seurat_clusters %in% c(7,16)]<-"B_cell"
annotation[seurat_clusters== 17]<-"Endothelial_cells"
#annotation[seurat_clusters== 15]<-"CMP"

annotation_col
pred.brca$labels<-annotation
pdf("singler_ScoreHeatmap_samples.pdf")
plotScoreHeatmap(pred.brca, clusters = CC_scrna@meta.data$samples)
dev.off()

pdf("singler_ScoreHeatmap_clusters.pdf")
plotScoreHeatmap(pred.brca, clusters = as.character(CC_scrna@meta.data$seurat_clusters))
dev.off()


CC_scrna@meta.data$annotation <- annotation
save(CC_scrna,file="CC_scrna.RData")

plots_annotation <- DimPlot(CC_scrna,group.by = "annotation",pt.size=0.5,reduction="tsne",label = TRUE, repel=TRUE)+ 
                    umap_theme()+ ggtitle("Annotation")+labs(x="TSNE",y="")+
        scale_colour_manual(name = "Annotation", values =color_scanpy_8 ) +
        theme(aspect.ratio=1)

pdf(file=paste(output_pre,"cluster_annotation_umap.pdf",sep=""),height=8,width=8)
print(plots_annotation)
dev.off()

##################输出结果
setwd("/share/pub/dengcy/Singlecell/CC_space/1.scRNA_seq_prepare")
load("CC_scrna.RData")
data_mat<- GetAssayData(CC_scrna, slot = "counts", assay ="RNA")
write.table(data_mat,file="SupplementaryData1.txt",quote=F)

meta_mat<- CC_scrna@meta.data
write.table(meta_mat,file="SupplementaryData2.txt",quote=F)
