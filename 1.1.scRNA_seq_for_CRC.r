library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(rtracklayer)
library(sctransform)
library(harmony)
library(cowplot)
library(stringr)
library(ggpubr)
source("/share/pub/dengcy/Singlecell/CC_space/src/Functions.r")
#source("/share/pub/dengcy/Singlecell/CC_space/src/library.r")
setwd("/share/pub/dengcy/Singlecell/CC_space/1.scRNA_seq_prepare")
###########
#1.E-MTAB-8107
##########
files<-c("scrEXT009.counts.csv","scrEXT013.counts.csv","scrEXT024.counts.csv","scrEXT028.counts.csv",
"scrEXT001.counts.csv","scrEXT010.counts.csv","scrEXT021.counts.csv","scrEXT025.counts.csv",
"scrEXT002.counts.csv","scrEXT018.counts.csv","scrEXT022.counts.csv",
"scrEXT012.counts.csv","scrEXT019.counts.csv","scrEXT027.counts.csv")
############################################
############ step1: data prepration
##########################################

cell_cycle_genes <- read.table(file = "/share/pub/xiongyc/project/scRNA/JiangFanChen/data/cc_genes_hg.txt",header = TRUE);
 #### step1.1: create Seurat object
CRC_Lambrechts_list<-lapply(files,function(x){
   count<-read.csv(paste0("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/colorectalcancer/E-MTAB-8107/",x))
   rownames(count)<-count[,1]
   count<-count[,-1]
   CRC_Lambrechts = CreateSeuratObject(counts = count,min.cells=3, min.features=200,plot_feature=F,
   cell_cycle_genes=cell_cycle_genes, project="CRC_Lambrechts");
   CRC_Lambrechts[["percent.mt"]] <- PercentageFeatureSet(CRC_Lambrechts, pattern = "^mt-")
   CRC_Lambrechts[["percent.ribo"]] <- PercentageFeatureSet(CRC_Lambrechts, pattern = "^Rp[sl][[:digit:]]")
   samples<-unlist(str_sub(rownames(CRC_Lambrechts@meta.data),18,18))
   return(CRC_Lambrechts)

 	})
 CRC_Lambrechts_merged <- merge(CRC_Lambrechts_list[[1]],y= CRC_Lambrechts_list[2:length(files)], project = "merged", merge.data = TRUE)
 CRC_Lambrechts_merged <- ScaleData(CRC_Lambrechts_merged)

 saveRDS(CRC_Lambrechts_merged,file="CRC_Lambrechts_merged.rds")


###########
#2.GSE164522
##########
#85052

#files<-c("GSE164522_CRLM_LN_expression.csv.gz",
#	"GSE164522_CRLM_PT_expression.csv.gz"
#)

 #### step1.1: create Seurat object
#  meta.data <- bigreadr::fread2("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/colorectalcancer/GSE164522/GSE164522_CRLM_metadata.csv.gz")
#rownames(meta.data)<- meta.data$V1
#CRC_Zhang_list<-lapply(files,function(x){
#   count <- bigreadr::fread2(paste0("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/colorectalcancer/GSE164522/",x))
#   rownames(count)<-count[,1]
#   count<-count[,-1]
#   CRC_Zhang = CreateSeuratObject(counts = count,min.cells=3, min.features=200,plot_feature=F,
#   cell_cycle_genes=cell_cycle_genes, project="CRC_Zhang");
#   return(CRC_Zhang)

# 	})
# CRC_Zhangs_merged <- merge(CRC_Zhang_list[[1]],y= CRC_Zhang_list[[2]], project = "merged", merge.data = TRUE)
 #rm(CRC_Zhang_list)
 
# meta.data <- bigreadr::fread2("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/colorectalcancer/GSE164522/GSE164522_CRLM_metadata.csv.gz")
 #meta.data<-  meta.data[ meta.data$sample %in% c("lymph node","primary tumor"),]
# rownames(meta.data)<-str_replace_all(meta.data$V1,"-",".")
# meta.data<-meta.data[rownames(meta.data) %in% colnames(CRC_Zhangs_merged),]
# meta.data<-meta.data[colnames(CRC_Zhangs_merged),]
# table(rownames(meta.data)==colnames(CRC_Zhangs_merged))
# CRC_Zhangs_merged@meta.data<-cbind(CRC_Zhangs_merged@meta.data,meta.data)
#CRC_Zhangs_merged <- ScaleData(CRC_Zhangs_merged)

# saveRDS(CRC_Zhangs_merged,file="CRC_Zhang_merged.rds")
################
##3.GSE144735
#############
#17678
counts<-bigreadr::fread2("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/colorectalcancer/GSE144735/GSE144735_processed_KUL3_CRC_10X_raw_UMI_count_matrix.txt.gz")
annotation<-bigreadr::fread2("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/colorectalcancer/GSE144735/GSE144735_processed_KUL3_CRC_10X_annotation.txt.gz")
  rownames(counts)<-counts[,1]
   counts<-counts[,-1]
   rownames(annotation)<-colnames(counts)
CRC_Park_ng1 = CreateSeuratObject(counts = counts,min.cells=3, 
	min.features=200,plot_feature=F, 
	project="CRC_Park",meta.data = annotation);
CRC_Park_ng1<-CRC_Park_ng1[,CRC_Park_ng1$Class!="Normal"]
CRC_Park_ng1 <- ScaleData(CRC_Park_ng1)

saveRDS(CRC_Park_ng1,file="CRC_Park_ng1.rds")

################
##4. GSE132257
#############
#7494
counts<-bigreadr::fread2("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/colorectalcancer/GSE132257/GSE132257_GEO_processed_protocol_and_fresh_frozen_raw_UMI_count_matrix.txt.gz")
annotation<-bigreadr::fread2("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/colorectalcancer/GSE132257/GSE132257_processed_protocol_and_fresh_frozen_cell_annotation.txt.gz")
  rownames(counts)<-counts[,1]
   counts<-counts[,-1]
   rownames(annotation)<-colnames(counts)
CRC_Park_ng2 = CreateSeuratObject(counts = counts,min.cells=3, 
	min.features=200,plot_feature=F, 
	project="CRC_Park2",meta.data = annotation);
CRC_Park_ng2<-CRC_Park_ng2[,CRC_Park_ng2$Class=="Tumor"]
CRC_Park_ng2<-CRC_Park_ng2[,CRC_Park_ng2$Status!="Frozen"]
CRC_Park_ng2 <- ScaleData(CRC_Park_ng2)

saveRDS(CRC_Park_ng2,file="CRC_Park_ng2.rds")

############
##5.GSE132465
############
#47285
counts<-bigreadr::fread2("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/colorectalcancer/GSE132465/GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz")
annotation<-bigreadr::fread2("/share/pub/dengcy/Cancer_Gwas/CollectedData/SingleCelldata/colorectalcancer/GSE132465/GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz")
  rownames(counts)<-counts[,1]
   counts<-counts[,-1]
   rownames(annotation)<-colnames(counts)
   counts<-counts[,annotation$Class=="Tumor"]
   annotation<-annotation[annotation$Class=="Tumor",]
CRC_Park_ng3 = CreateSeuratObject(counts = counts,min.cells=3, 
	min.features=200,plot_feature=F, 
	project="CRC_Park3",meta.data = annotation);
CRC_Park_ng3 <- ScaleData(CRC_Park_ng3)
saveRDS(CRC_Park_ng3,file="CRC_Park_ng3.rds")


##########
#二、整合
##########
setwd("/share/pub/dengcy/Singlecell/CC_space/1.scRNA_seq_prepare")
CRC_Lambrechts<-readRDS("CRC_Lambrechts_merged.rds")
#CRC_Zhang<-readRDS("CRC_Zhang_merged.rds")
CRC_Park1<-readRDS("CRC_Park_ng1.rds")
CRC_Park2<-readRDS("CRC_Park_ng2.rds")
CRC_Park3<-readRDS("CRC_Park_ng3.rds")
CRC_own <- readRDS("CRC_scRNAseq.rds")
#DefaultAssay(CRC_own)<-"RNA"
#CRC_own <- ScaleData(CRC_own)
#saveRDS(CRC_own,file="CRC_scRNAseq.rds")
crc_features <- SelectIntegrationFeatures(object.list = list(CRC_own,CRC_Lambrechts,CRC_Park1,CRC_Park2,CRC_Park3), nfeatures = 2000)
CRC_merged <- merge(CRC_own, y = list(CRC_Lambrechts,CRC_Park1,CRC_Park2,CRC_Park3), project = "merged", merge.data = TRUE)
rm(CRC_own,CRC_Lambrechts,CRC_Park1,CRC_Park2,CRC_Park3)
CRC_merged <- ScaleData(CRC_merged)
VariableFeatures(CRC_merged) <- crc_features
CRC_merged <- RunPCA(object = CRC_merged, assay = "RNA", features = crc_features, npcs = 50)
###
fortify_df<-fortify_df[fortify_df$orig.ident==,]

fortify_df$Samples<-fortify_df$orig.ident
fortify_df$Samples[!is.na(fortify_df$samples)]<-fortify_df$samples[!is.na(fortify_df$samples)]
#fortify_df$Patients<- fortify_df$Samples

a<-unlist(lapply(fortify_df$Samples,function(x) strsplit(x, "-")[[1]][1]))
a[a=='SMC171.T.BEL']<-'SMC171'
a[a=='SMC171.T.SGI']<-'SMC171'
a[a=='SMC171.T.SING']<-'SMC171'

fortify_df$Patients<-a
############################
####Harmony
############################
print("step3: Hamony: remove batch effect ##################################################################")
###there are little  difference between two samples
CRC_merged_batched<- RunHarmony(object = CRC_merged,
                              assay.use = "RNA",
                              reduction = "pca",
                              dims.use = 1:50,
                              group.by.vars = "orig.ident",
                              plot_convergence = TRUE)

saveRDS(CRC_merged_batched,file="CRC_merged_all.rds")

############
#step1:预处理,聚类可视化
#############
setwd("/share/pub/dengcy/Singlecell/CC_space/1.scRNA_seq_prepare")
source("/share/pub/dengcy/Singlecell/CC_space/src/Functions.r")

CRC_merged_all<-readRDS("CRC_merged_all.rds")
GetAssayData(CRC_merged_all,slot="data",assays="SCI")
CRC_merged_all@assays$SCT<-NULL
saveRDS(CRC_merged_all,file="CRC_merged_all.rds")

fortify_df<-fortify.Seurat(CRC_merged_all)
fortify_df$stage<-"StageIII"
fortify_df$stage[fortify_df$samples=="CC1"]<-"StageI"
fortify_df<-fortify_df[fortify_df$orig.ident!="KUL21-T",]
fortify_df<-fortify_df[fortify_df$orig.ident!="KUL21-B",]
fortify_df<-fortify_df[fortify_df$orig.ident!="SMC25-T",]
fortify_df<-fortify_df[fortify_df$orig.ident!="SMC21-T",]
fortify_df<-fortify_df[fortify_df$orig.ident!="scrEXT012",]
fortify_df<-fortify_df[fortify_df$orig.ident!="scrEXT013",]
save(fortify_df,file="CRC_fortify_df.RData")

load("CRC_fortify_df.RData")
CRC_merged_all<-CRC_merged_all[,rownames(fortify_df)]

CRC_merged_all<-cluster_pca_umap(CRC_merged_all, assay = "RNA",reduction="harmony",cluster_res = 0.2)
saveRDS(CRC_merged_all,file="CRC_merged_all.rds")
"#E1FFEE",
color_CRC<-c("#7FBCD2","#A5F1E9","#FFEEAF","#7FB77E","#F7F6DC","#B1D7B4","#FFC090",
	"#F3E0B5","#FFAE6D","#FECD70","#E3C770","#94B49F","#CEE5D0","#ECB390","#F2D7D9",
	"#D3CEDF","#B2C8DF","#C4D7E0","#F8F9D7","#76BA99","#ADCF9F",
	"#CED89E","#FFDCAE","#FFE3A9","#FFC3C3","#FF8C8C","#97C4B8","#F1F0C0","#B1BCE6")

plots_tsne <- DimPlot(CRC_merged_all,group.by = "seurat_clusters",reduction ="tsne",
	label = TRUE, pt.size =0.5,
  repel=TRUE,
	label.size = 4,
  label.color = "black",
  label.box = TRUE) + umap_theme()+ ggtitle("CRC merged")+
        scale_colour_manual(name = "Cluster", values = color_CRC)
pdf(file="step1_allClusters_tsne.pdf",height=10,width=10)
print(plots_tsne)
dev.off()

#######
#设置stage注释
fortify_df<-fortify.Seurat(CRC_merged_all)
fortify_df$stage<-"StageIII"
fortify_df$stage[fortify_df$samples=="CC1"]<-"StageI"
fortify_df$stage[fortify_df$samples=="CC2"]<-"StageIII"
fortify_df$stage[fortify_df$orig.ident %in% c("KUL01-B","KUL01-T","KUL28-B","KUL28-T","KUL30-B","KUL30-T",
  "SMC01-T","SMC05-T","SMC09-T","SMC10-T","SMC15-T","SMC18-T","scrEXT001","scrEXT002","scrEXT021",
  "scrEXT022","scrEXT027","scrEXT028")]<-"StageII"
fortify_df$stage[fortify_df$orig.ident %in% c("KUL31-B","KUL31-T","SMC07-T","SMC24-T","scrEXT018","scrEXT019")]<-"StageI"
table(fortify_df[,c("stage","orig.ident")])
save(fortify_df,file="CRC_fortify_df.RData")

CRC_merged_all$stage<-fortify_df$stage
########
#设置数据来源注释
CRC_merged_all$datasource<-"GSE132465"
CRC_merged_all$datasource[fortify_df$orig.ident %in% c("scrEXT001","scrEXT002","scrEXT009","scrEXT010","scrEXT018",
  "scrEXT019","scrEXT021","scrEXT022","scrEXT024","scrEXT025","scrEXT027","scrEXT028")]<-"E-MTAB-8107"
CRC_merged_all$datasource[fortify_df$orig.ident %in% c("KUL01-B","KUL01-T ","KUL19-B","KUL19-T","KUL28-B","KUL28-T",
  "KUL30-B","KUL30-T","KUL31-B","KUL31-T")]<-"GSE144735"
CRC_merged_all$datasource[fortify_df$orig.ident=="CC_scRNA-seq"]<-"twoCRCsamples"
CRC_merged_all$datasource[fortify_df$orig.ident %in% c("SMC13T.A1","SMC17-T","SMC171.T.BEL","SMC171.T.SGI","SMC171.T.SING")]<-"GSE132257"

fortify_df$datasource<-CRC_merged_all$datasource
save(fortify_df,file="CRC_fortify_df.RData")

plots<-lapply(unique(CRC_merged_all$datasource),function(x){
 p1 <- ggplot() +
        geom_point(data = fortify_df[!(fortify_df$datasource %in% x),],
                   aes(x = TSNE_1, y = TSNE_2), size = 0.2, alpha = 0.2, color = "gray") +
        umap_theme() +
        geom_point(data = fortify_df[(fortify_df$datasource %in% x),], 
          aes(x = TSNE_1, y = TSNE_2),color = "black", size = .2) +
        ggtitle(x)
  return(p1) 
  })

pdf(file="step1_datasource_split_stne.pdf",height=11,width=15)
CombinePlots(plots)
dev.off()

colors<-c("#89C9B8","#D9ADAD","#84A9AC","#DDDDDD")
names(colors)<-c("StageI","StageIII" ,"StageII","default")
plots<-lapply(c("StageI","StageII","StageIII"),function(x){
 p1 <- ggplot() +
        geom_point(data = fortify_df[!(fortify_df$stage %in% x),],
                   aes(x = TSNE_1, y = TSNE_2), size = 0.2, alpha = 0.2, color = "gray") +
        umap_theme() +
        geom_point(data = fortify_df[(fortify_df$stage %in% x),], 
          aes(x = TSNE_1, y = TSNE_2),color = colors[x], size = .2) +
        ggtitle(x)
  return(p1) 
  })

pdf(file="step1_stage_split_stne.pdf",height=6,width=17)
CombinePlots(plots, ncol =3)
dev.off()
library(ggtext)
Samples_plot <- ggplot() +
        geom_point(data = fortify_df,
                   aes(x = TSNE_1, y = TSNE_2,color=Samples), size = 0.03, alpha = 0.8) +
        umap_theme() +
        #scale_colour_manual(name = "Samples", values = color_CRC[]) +
        theme(aspect.ratio=1) +
        theme(legend.text=element_markdown(size=14),
              legend.title=element_text(size=14)) +
        guides(colour = guide_legend(override.aes = list(size=3))) +
        ggtitle("Samples")

        pdf(file="Samples_tsne.pdf",width = 12, height = 12)
        print(Samples_plot)
        dev.off()
p_plot <- ggplot() +
        geom_point(data = fortify_df,
                   aes(x = TSNE_1, y = TSNE_2,color=Patients), size = 0.03, alpha = 0.8) +
        umap_theme() +
        #scale_colour_manual(name = "Samples", values = color_CRC[]) +
        theme(aspect.ratio=1) +
        theme(legend.text=element_markdown(size=14),
              legend.title=element_text(size=14)) +
        guides(colour = guide_legend(override.aes = list(size=3))) +
        ggtitle("Patients")

        pdf(file="Patients_tsne.pdf",width = 12, height =12)
        print(p_plot)
        dev.off()
######
#细胞量统计
cluster_c<-table(fortify_df[,c("seurat_clusters")])
cluster_d<-table(fortify_df[,c("datasource")])
stage_s<-table(fortify_df[,c("stage")])
samples_s<-table(fortify_df[,c("Samples")])
Patients_s<-table(fortify_df[,c("Patients")])

table(fortify_df[,c("stage","Samples")])
length(unique(fortify_df$Samples[fortify_df$stage =="StageI"]))
#7
length(unique(fortify_df$Samples[fortify_df$stage =="StageII"]))
#18
length(unique(fortify_df$Samples[fortify_df$stage =="StageIII"]))
#24

library(dplyr)
library(xtable)
library(flextable)
library(officer)
m1 = xtable_to_flextable(xtable(cluster_c))
m2 = xtable_to_flextable(xtable(cluster_d))
m3 = xtable_to_flextable(xtable(stage_s))
m4 = xtable_to_flextable(xtable(samples_s))
m5 = xtable_to_flextable(xtable(Patients_s))

doc = read_docx()
doc = body_add_flextable(doc,m1)
print(doc,"/share/pub/dengcy/Singlecell/CC_space/1.scRNA_seq_prepare/cluster_stat.docx")
doc = read_docx()
doc = body_add_flextable(doc,m2)
print(doc,"/share/pub/dengcy/Singlecell/CC_space/1.scRNA_seq_prepare/datasource_stat.docx")
doc = read_docx()
doc = body_add_flextable(doc,m3)
print(doc,"/share/pub/dengcy/Singlecell/CC_space/1.scRNA_seq_prepare/stage_stat.docx")
doc = read_docx()
doc = body_add_flextable(doc,m4)
print(doc,"/share/pub/dengcy/Singlecell/CC_space/1.scRNA_seq_prepare/samples_stat.docx")
doc = read_docx()
doc = body_add_flextable(doc,m5)
print(doc,"/share/pub/dengcy/Singlecell/CC_space/1.scRNA_seq_prepare/Patients_stat.docx")


##############
#step2: 计算差异基因，细胞注释
################
setwd("/share/pub/dengcy/Singlecell/CC_space/1.scRNA_seq_prepare")
source("/share/pub/dengcy/Singlecell/CC_space/src/Functions.r")

load("CRC_merged_all.RData")
meta.data<-CRC_merged_all@meta.data
save(meta.data,file="meta.data_all.RData")
#plan("multiprocess", workers = 20)
CRC.markers <- FindAllMarkers(object = CRC_merged_all, only.pos = TRUE, 
                               min.pct = 0.25, 
                               thresh.use = 0.25)
save(CRC.markers,file="CRC.markers.RData")
#CRC.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
load("CRC.markers.RData")
top5 <- CRC.markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)

pdf("step2_CRC.markergenes.HEATMAP.pdf",height=12,width=10)
DoHeatmap(CRC_merged_all,features=top5$gene,group.colors=color_CRC,size = 3)
dev.off()

top10 <- CRC.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
write.csv(top10,file="/share/pub/dengcy/Singlecell/CC_space/1.scRNA_seq_prepare/firstcluster_markers.csv")

###singleR
library("celldex")
library("SingleR")
library("pracma")
library(scater)
library(SingleCellExperiment)
ref <- HumanPrimaryCellAtlasData() 
common <- intersect(rownames(ref), rownames(CRC_merged_all))
CRC_merged_all.singler <- CRC_merged_all[common,]

meta.data<-CRC_merged_all@meta.data
rm(CRC_merged_all)
gc()
ref <- ref[common,]

CRC_merged_all.singler.sce <-SingleCellExperiment(assays = list(counts = GetAssayData(object =CC_scrna.singler, slot = "counts")))
CRC_merged_all.singler.sce <- logNormCounts(CRC_merged_all.singler.sce)

pred.brca <- SingleR(test = CRC_merged_all.singler.sce, 
  ref = ref, assay.type.test=1,labels = ref$label.main)
save(pred.brca,file="singler_pred.RData")
singler.results <- merge(data.frame(cell = rownames(pred.brca), singler = pred.brca$labels), 
                         data.frame(cell = rownames(meta.data), 
                                    cluster = meta.data$seurat_clusters), 
                         by = "cell", 
                         all.y = FALSE)
singler.results$cell <- NULL
singler.results$count <- 1
singler.results <- aggregate(count ~ ., singler.results, FUN = sum)
singler.final <- singler.results %>% group_by(cluster) %>% top_n(n = 1, wt = count)
write.csv(singler.final,file="singler.final.csv")

#/share/pub/dengcy/Singlecell/CC_space/SupplementFiles10.26

#####step2:绘制注释后的细胞图
#
singler.final<-read.csv("/share/pub/dengcy/Singlecell/CC_space/1.scRNA_seq_prepare/singler.final.csv")
write.csv(singler.final,file="/share/pub/dengcy/Singlecell/CC_space/SupplementFiles10.26/TableS4.csv")

seurat_clusters<-as.vector(CRC_merged_all$seurat_clusters)
annotation<-as.vector(seurat_clusters)
load("singler_pred.RData")
annotation[seurat_clusters %in% singler.final$cluster[which(singler.final$singler=="Epithelial_cells")]]<-"Epithelial_cells"
annotation[seurat_clusters %in% singler.final$cluster[which(singler.final$singler=="T_cells")]]<-"T_cells"
annotation[seurat_clusters %in% singler.final$cluster[which(singler.final$singler=="B_cell")]]<-"B_cell"
annotation[seurat_clusters == 3]<-"Macrophage"
annotation[seurat_clusters == 6 ]<-"Fibroblasts"
annotation[seurat_clusters ==15]<-"Monocyte"
annotation[seurat_clusters== 8]<-"Endothelial_cells"
annotation[seurat_clusters== 14]<-"CMP"
annotation[seurat_clusters== 12]<-"Tissue_stem_cells"
CRC_merged_all$Initial_annotation<-annotation
save(CRC_merged_all,file="CRC_merged_all.RData")

pred.brca$labels<-annotation
pdf("step2_singler_ScoreHeatmap_clusters.pdf")
plotScoreHeatmap(pred.brca, clusters = as.character(CRC_merged_all@meta.data$seurat_clusters))
dev.off()


plots_annotation <- DimPlot(CRC_merged_all,group.by = "Initial_annotation",reduction ="tsne",
  label = TRUE, pt.size =0.5,
  repel=TRUE,
  label.size = 4,
  label.color = "black",
  label.box = F) + umap_theme()+ ggtitle("CRC merged")+
        scale_colour_manual(name = "Annotation", values = color_CRC)
pdf(file="step2_cluster_annotation_tsne.pdf",height=10,width=10)
print(plots_annotation)
dev.off()

CRC_merged_all$stage<-"StageIII"
CRC_merged_all$stage[CRC_merged_all$samples=="CC1"]<-"StageI/II"
CRC_merged_all$stage[CRC_merged_all$orig.ident %in% c("KUL01-B","KUL01-T","KUL28-B","KUL28-T","KUL30-B","KUL30-T",
  "SMC01-T","SMC05-T","SMC09-T","SMC10-T","SMC15-T","SMC18-T","scrEXT001","scrEXT002","scrEXT021",
  "scrEXT022","scrEXT027","scrEXT028")]<-"StageI/II"
CRC_merged_all$stage[CRC_merged_all$orig.ident %in% c("KUL31-B","KUL31-T","SMC07-T","SMC24-T","scrEXT018","scrEXT019")]<-"StageI/II"
CRC_merged_all$orig.ident[CRC_merged_all$samples=="CC1"]<-"CC1"
CRC_merged_all$orig.ident[CRC_merged_all$samples=="CC2"]<-"CC2"

pdf("allcellproportion_in_stage.pdf",width=8,height=3)
stack_plot(obj=CRC_merged_all,options=c("Initial_annotation","stage"),colors= rev(color_CRC[1:9]),angle=FALSE)
dev.off()
all_metadata<- CRC_merged_all@meta.data
save(all_metadata,file="all_metadata.RData")

##################输出结果
#setwd("/share/pub/dengcy/Singlecell/CC_space/1.scRNA_seq_prepare")
#load("CC_scrna.RData")
#data_mat<- GetAssayData(CC_scrna, slot = "counts", assay ="RNA")
#write.table(data_mat,file="SupplementaryData1.txt",quote=F)

#meta_mat<- CC_scrna@meta.data
#write.table(meta_mat,file="SupplementaryData2.txt",quote=F)

#####marker gene
pdf("epi_marker.pdf",width=9)
FeaturePlot(CRC_merged_all, reduction ="tsne",features = c('EPCAM', 'SFN','KRT19','KRT18'))
dev.off()