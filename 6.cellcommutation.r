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
setwd("/share/pub/dengcy/Singlecell/CC_space/6.cellcom")
load("/share/pub/dengcy/Singlecell/CC_space/3.Epi_cells_stage/cancer_scrna_cc1.RData")
load("/share/pub/dengcy/Singlecell/CC_space/3.Epi_cells_stage/cancer_scrna_cc2.RData")

load("/share/pub/dengcy/Singlecell/CC_space/3.2.Immune_cells_stage/Immune_scrna.RData")
########################CC1
Immune_scrna_cc1<- subset(Immune_scrna,samples == "CC1" )
scrna_cc1<-merge(Immune_scrna_cc1,cancer_scrna_cc1)
Idents(scrna_cc1)<-scrna_cc1$annotation
save(scrna_cc1,file="scrna_cc1.RData")

Immune_scrna_cc2<- subset(Immune_scrna,samples == "CC2" )
scrna_cc2<-merge(Immune_scrna_cc2,cancer_scrna_cc2)
Idents(scrna_cc2)<-scrna_cc2$annotation
save(scrna_cc2,file="scrna_cc2.RData")


library("SingleCellSignalR")
###1.singlecellsignalR
test1<-subset(scrna_cc1, annotation
 %in% c("CT1","CT2","Cytotoxic.CD8.T","Neutrophils","B.cell","Terminally.Exhausted CD8T","CD4.TCM"))
test1_data<-GetAssayData(object =test1, slot = "data")

cluster<-test1$annotation

cluster[cluster=="CT1"]<-1
cluster[cluster=="CT2"]<-2
cluster[cluster=="CT3"]<-3
cluster[cluster=="TAM"]<-4
cluster[cluster=="NK.cell"]<-5
cluster[cluster=="B.cell.activite.memory"]<-6
cluster[cluster=="Cytotoxic.CD8.T"]<-7
cluster[cluster=="CD4.TCM"]<-8
cluster[cluster=="Terminally.Exhausted CD8T"]<-9
cluster[cluster=="Neutrophils"]<-10
cluster[cluster=="B.cell"]<-11

signal_test1 <- cell_signaling(data = test1_data, genes = rownames(test1_data), cluster = cluster, write = FALSE)
inter.net <- inter_network(data = data, signal = signal, genes = rownames(test1_data), cluster = cluster, write = FALSE)
pdf("signal_test1.pdf")
visualize_interactions(signal = signal_test1)
dev.off()

#####################

load("scrna_cc2.RData")

test2<-subset(scrna_cc2, annotation %in% c("CT4","CT5","CT6","B.cell.activite.memory","TAM","Cytotoxic.CD8.T","NK.cell","Neutrophils","B.cell","Terminally.Exhausted CD8T","CD4.TCM"))
test2_data<-GetAssayData(object =test2, slot = "data")

cluster<-test2$annotation

cluster[cluster=="CT4"]<-1
cluster[cluster=="CT5"]<-2
cluster[cluster=="CT6"]<-3
cluster[cluster=="TAM"]<-4
cluster[cluster=="NK.cell"]<-5
cluster[cluster=="B.cell.activite.memory"]<-6
cluster[cluster=="Cytotoxic.CD8.T"]<-7
cluster[cluster=="CD4.TCM"]<-8
cluster[cluster=="Terminally.Exhausted CD8T"]<-9
cluster[cluster=="Neutrophils"]<-10
cluster[cluster=="B.cell"]<-11

signal_test2 <- cell_signaling(data = test2_data, genes = rownames(test2_data), cluster = cluster, write = FALSE)
inter.net <- inter_network(data = data, signal = signal, genes = rownames(test2_data), cluster = cluster, write = FALSE)
pdf("signal_test2.pdf")
visualize_interactions(signal = signal_test2)
dev.off()


########################
#cellchat

library(CellChat)
library(ggplot2)                  
library(patchwork)
library(igraph)
options(stringsAsFactors = FALSE)
#data_cc1  <- scrna_cc1@assays$RNA@data
#identity_cc1 = data.frame(group =scrna_cc1$annotation, row.names = names(scrna_cc1$annotation)) # create a dataframe consisting of the cell labels
#unique(identity_cc1$group) # check the cell labels
scrna_cc1_ct2<-subset(scrna_cc1,subset= annotation %in% 
	c("B.cell","CD4.TCM","CT2","Cytotoxic.CD8.T","Neutrophils","Terminally.Exhausted CD8T","NK.cell","TAM"))
cellchat <- createCellChat(object = scrna_cc1_ct2,group.by="annotation")
cellchat <- setIdent(cellchat, ident.use = "annotation")

cellchat_single(cellchat=cellchat,filename="cc1_ct2")

scrna_cc1_ct1<-subset(scrna_cc1,subset= annotation %in% 
	c("B.cell","CD4.TCM","CT1","Cytotoxic.CD8.T","Neutrophils","Terminally.Exhausted CD8T","NK.cell","TAM"))
cellchat <- createCellChat(object = scrna_cc1_ct1,group.by="annotation")
cellchat <- setIdent(cellchat, ident.use = "annotation")

cellchat_single(cellchat=cellchat,filename="cc1_ct1")


scrna_cc2_ct4<-subset(scrna_cc2,subset= annotation %in% c("B.cell","CD4.TCM","CT4","Cytotoxic.CD8.T","Neutrophils","Terminally.Exhausted CD8T","TAM","NK.cell","B.cell.activite.memory"))
cellchat <- createCellChat(object = scrna_cc2_ct4,group.by="annotation")
cellchat <- setIdent(cellchat, ident.use = "annotation")

cellchat_single(cellchat=cellchat,filename="cc2_ct4")


scrna_cc2_ct5<-subset(scrna_cc2,subset= annotation %in% c("B.cell","CD4.TCM","CT5","Cytotoxic.CD8.T","Neutrophils","Terminally.Exhausted CD8T","TAM","NK.cell","B.cell.activite.memory"))
cellchat <- createCellChat(object = scrna_cc2_ct5,group.by="annotation")
cellchat <- setIdent(cellchat, ident.use = "annotation")

cellchat_single(cellchat=cellchat,filename="cc2_ct5")


scrna_cc2_ct6<-subset(scrna_cc2,subset= annotation %in% c("B.cell","CD4.TCM","CT6","Cytotoxic.CD8.T","Neutrophils","Terminally.Exhausted CD8T","TAM","NK.cell","B.cell.activite.memory"))
cellchat <- createCellChat(object = scrna_cc2_ct6,group.by="annotation")
cellchat <- setIdent(cellchat, ident.use = "annotation")

cellchat_single(cellchat=cellchat,filename="cc2_ct6")


scrna_cc1<-subset(scrna_cc1,subset= annotation %in% 
	c("B.cell","CD4.TCM","CT1","CT2","Cytotoxic.CD8.T","Neutrophils","Terminally.Exhausted CD8T","NK.cell","TAM"))
cellchat <- createCellChat(object = scrna_cc1,group.by="annotation")
cellchat <- setIdent(cellchat, ident.use = "annotation")

cellchat<-cellchat_single(cellchat=cellchat,filename="cc1")

pdf("cc1_CT2_othercells_netVisual_bubble.pdf",width=4,height=5)
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(1,2,5,6,7,8,9), remove.isolate = FALSE)
#> Comparing communications on a single object
dev.off()

pdf("cc1_CT1_othercells_netVisual_bubble.pdf",width=4,height=5)
netVisual_bubble(cellchat, sources.use = 3, targets.use = c(1,2,5,6,7,8,9), remove.isolate = FALSE)
#> Comparing communications on a single object
dev.off()


scrna_cc2<-subset(scrna_cc2,subset= annotation %in% c("B.cell","CD4.TCM","CT4","CT5","CT6","Cytotoxic.CD8.T","Neutrophils","Terminally.Exhausted CD8T","TAM","NK.cell","B.cell.activite.memory"))
cellchat <- createCellChat(object = scrna_cc2,group.by="annotation")
cellchat <- setIdent(cellchat, ident.use = "annotation")
cellchat<-cellchat_single(cellchat=cellchat,filename="cc2")

pdf("cc2_CT4_othercells_netVisual_bubble.pdf",width=5)
netVisual_bubble(cellchat, sources.use = 4, targets.use =c(1:3,7:11), remove.isolate = FALSE)
#> Comparing communications on a single object
dev.off()

pdf("cc2_CT5_othercells_netVisual_bubble.pdf",width=5,height=3.5)
netVisual_bubble(cellchat, sources.use = 5, targets.use =c(1:3,7:11), remove.isolate = FALSE)
#> Comparing communications on a single object
dev.off()

pdf("cc2_CT6_othercells_netVisual_bubble.pdf",width=5,height=3)
netVisual_bubble(cellchat, sources.use = 6, targets.use =c(1:3,7:11), remove.isolate = FALSE)
#> Comparing communications on a single object
dev.off()


#############################
library(nichenetr)
library(Seurat) # please update to Seurat V4
library(tidyverse)
#Read in NicheNet’s ligand-target prior model, ligand-receptor network and weighted integrated networks:
ligand_target_matrix = readRDS("ligand_target_matrix.rds")
#ligand_target_matrix[1:5,1:5] # target genes in rows, ligands in columns
lr_network = readRDS("lr_network.rds")
weighted_networks = readRDS("weighted_networks.rds")
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))

 load("/share/pub/dengcy/Singlecell/CC_space/4.stage_visualize/cor_p.RData")
stage_dogenes
stage_upgenes

## receiver
scrna_cc1$celltype<-scrna_cc1$annotation
scrna_cc2$celltype<-scrna_cc2$annotation
Idents(scrna_cc1)<-scrna_cc1$annotation
Idents(scrna_cc2)<-scrna_cc2$annotation
 
scrna_cc1=SCTransform(scrna_cc1, vars.to.regress = c('nFeature_RNA', 'nCount_RNA',"percent.mt","percent.ribo","S.Score", "G2M.Score"), verbose = FALSE)
scrna_cc1 <- FindVariableFeatures(scrna_cc1,nfeatures = 7000)
scrna_cc1 <- RunPCA(object = scrna_cc1, assay = "SCT", npcs = 50)
scrna_cc1<- cluster_pca_umap(scrna_cc1, assay = "SCT",reduction="pca",cluster_res = 0.3)

scrna_cc2=SCTransform(scrna_cc2, vars.to.regress = c('nFeature_RNA', 'nCount_RNA',"percent.mt","percent.ribo","S.Score", "G2M.Score"), verbose = FALSE)
scrna_cc2 <- FindVariableFeatures(scrna_cc2,nfeatures = 7000)
scrna_cc2 <- RunPCA(object = scrna_cc2, assay = "SCT", npcs = 50)
scrna_cc2<- cluster_pca_umap(scrna_cc2, assay = "SCT",reduction="pca",cluster_res = 0.3)

save(scrna_cc1,file="scrna_cc1.RData")
save(scrna_cc2,file="scrna_cc2.RData")
load("scrna_cc1.RData")
load("scrna_cc2.RData")
#########

  EMT_genes<-read.delim("/share/pub/dengcy/Singlecell/CC_space/4.stage_visualize/EMTgene.txt")

stage_dogenes =intersect(stage_dogenes,rownames(scrna_cc1))
stage_upgenes =intersect(stage_upgenes,rownames(scrna_cc1))
emt_genes<-intersect(EMT_genes$x,rownames(scrna_cc1))

ligand_target_pheatmap<-function(scrna_cc1=scrna_cc1,
	geneset_oi=stage_dogenes,
	celltype1="Terminally.Exhausted CD8T",
	celltype2="CT2",
	colors="#11468F",
	filenames="TE-ct2_down_ligand_target_network.pdf",
	width=4,
	height=7,
	title1="Prioritized Terminally.Exhausted.T-ligands",
	title2="invasion-down genes in CT2 cells"
	){

expression = t(scrna_cc1@assays$SCT@data)
sample_info = scrna_cc1@meta.data # contains meta-information about the cells
sample_info$cellid<-rownames(sample_info)

tumors_remove = c("HN10","HN","HN12", "HN13", "HN24", "HN7", "HN8","HN23")
ids1 = sample_info %>% filter(annotation == celltype1) %>% pull(cellid)
ids2 = sample_info %>% filter(annotation ==celltype2 ) %>% pull(cellid)

expressed_genes_sender = expression[ids1,] %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}) %>% .[. >= 2] %>% names()

expressed_genes_receiver = expression[ids2,] %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}) %>% .[. >= 2] %>% names()

background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
ligands = lr_network %>% pull(from) %>% unique()
expressed_ligands = intersect(ligands,expressed_genes_sender)
receptors = lr_network %>% pull(to) %>% unique()
expressed_receptors = intersect(receptors,expressed_genes_receiver)
lr_network_expressed = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) 
potential_ligands = lr_network_expressed %>% pull(from) %>% unique()

ligand_activities = predict_ligand_activities(geneset = geneset_oi, 
	background_expressed_genes = background_expressed_genes, 
	ligand_target_matrix = ligand_target_matrix, 
	potential_ligands = potential_ligands)

best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand)

active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,
	geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 1000) %>% bind_rows()

active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, 
	ligand_target_matrix = ligand_target_matrix, cutoff = 0.2)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets = active_ligand_target_links_df$target %>% unique()
vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
p_ligand_target_network = vis_ligand_target %>% 
make_heatmap_ggplot(title1,title2, 
	color = colors,legend_position = "top", x_axis_position = "top",
	legend_title = "Regulatory potential") + 
scale_fill_gradient2(low = "whitesmoke",  high = colors, breaks = c(0,0.005,0.01)) + 
theme(axis.text.x = element_text(face = "italic"))

pdf(filenames,width=width,height=height)
print(p_ligand_target_network)
dev.off()

}


ligand_target_pheatmap(scrna_cc1=scrna_cc1,
	geneset_oi=stage_dogenes,
	celltype1="Terminally.Exhausted CD8T",
	celltype2="CT2",colors="#11468F",
	filenames="TE-ct2_down_ligand_target_network.pdf",
	width=4,height=7,
	title1="Prioritized Terminally.Exhausted.T-ligands",
	title2="invasion-down genes in CT2 cells"
	)

ligand_target_pheatmap(scrna_cc1=scrna_cc1,
	geneset_oi=stage_upgenes,
	celltype1="Terminally.Exhausted CD8T",
	celltype2="CT2",colors="#DA1212",
	filenames="TE-ct2_up_ligand_target_network.pdf",
	width=3,height=7,
	title1="Prioritized Terminally.Exhausted.T-ligands",
	title2="invasion-up genes in CT2 cells"
	)

ligand_target_pheatmap(scrna_cc1=scrna_cc1,
	geneset_oi=stage_dogenes,
	celltype1="TAM",
	celltype2="CT2",colors="#11468F",
	filenames="TAM-ct2_down_ligand_target_network.pdf",
	width=4,height=7,
	title1="Prioritized TAM-ligands",
	title2="invasion-down genes in CT2 cells"
	)
ligand_target_pheatmap(scrna_cc1=scrna_cc1,
	geneset_oi=stage_upgenes,
	celltype1="TAM",
	celltype2="CT2",colors="#DA1212",
	filenames="TAM-ct2_up_ligand_target_network.pdf",
	width=3,height=7,
	title1="Prioritized TAM-ligands",
	title2="invasion-up genes in CT2 cells"
	)
ligand_target_pheatmap(scrna_cc1=scrna_cc1,
	geneset_oi=stage_dogenes,
	celltype1="Cytotoxic.CD8.T",
	celltype2="CT2",colors="#11468F",
	filenames="CD8T-ct2_down_ligand_target_network.pdf",
	width=4,height=7,
	title1="Prioritized Cytotoxic.CD8.T-ligands",
	title2="invasion-down genes in CT2 cells"
	)
ligand_target_pheatmap(scrna_cc1=scrna_cc1,
	geneset_oi=stage_upgenes,
	celltype1="Cytotoxic.CD8.T",
	celltype2="CT2",colors="#DA1212",
	filenames="CD8T-ct2_up_ligand_target_network.pdf",
	width=4,height=7,
	title1="Prioritized Cytotoxic.CD8.T-ligands",
	title2="invasion-up genes in CT2 cells"
	)

ligand_target_pheatmap(scrna_cc1=scrna_cc1,
	geneset_oi=stage_dogenes,
	celltype1="Neutrophils",
	celltype2="CT2",colors="#11468F",
	filenames="Neutrophils-ct2_down_ligand_target_network.pdf",
	width=4,height=7,
	title1="Prioritized Neutrophils-ligands",
	title2="invasion-down genes in CT2 cells"
	)

ligand_target_pheatmap(scrna_cc1=scrna_cc1,
	geneset_oi=stage_upgenes,
	celltype1="Neutrophils",
	celltype2="CT2",colors="#DA1212",
	filenames="Neutrophils-ct2_up_ligand_target_network.pdf",
	width=4,height=7,
	title1="Prioritized Neutrophils-ligands",
	title2="invasion-up genes in CT2 cells"
	)

ligand_target_pheatmap(scrna_cc1=scrna_cc2,
	geneset_oi=stage_dogenes,
	celltype1="Neutrophils",
	celltype2="CT4",colors="#11468F",
	filenames="Neutrophils-ct4_down_ligand_target_network.pdf",
	width=4,height=7,
	title1="Prioritized Neutrophils-ligands",
	title2="invasion-down genes in CT4 cells"
	)
ligand_target_pheatmap(scrna_cc1=scrna_cc2,
	geneset_oi=stage_upgenes,
	celltype1="Neutrophils",
	celltype2="CT4",colors="#DA1212",
	filenames="Neutrophils-ct4_up_ligand_target_network.pdf",
	width=4,height=7,
	title1="Prioritized Neutrophils-ligands",
	title2="invasion-up genes in CT4 cells"
	)


ligand_target_pheatmap(scrna_cc1=scrna_cc2,
	geneset_oi=stage_dogenes,
	celltype1="B.cell.activite.memory",
	celltype2="CT4",colors="#11468F",
	filenames="BM-ct4_down_ligand_target_network.pdf",
	width=4,height=7,
	title1="Prioritized B.cell.activite.memory-ligands",
	title2="invasion-down genes in CT4 cells"
	)
ligand_target_pheatmap(scrna_cc1=scrna_cc2,
	geneset_oi=stage_upgenes,
	celltype1="B.cell.activite.memory",
	celltype2="CT4",colors="#DA1212",
	filenames="BM-ct4_up_ligand_target_network.pdf",
	width=4,height=7,
	title1="Prioritized B.cell.activite.memory-ligands",
	title2="invasion-up genes in CT4 cells"
	)

ligand_target_pheatmap(scrna_cc1=scrna_cc2,
	geneset_oi=stage_dogenes,
	celltype1="TAM",
	celltype2="CT4",colors="#11468F",
	filenames="TAM-ct4_down_ligand_target_network.pdf",
	width=4,height=7,
	title1="Prioritized TAM-ligands",
	title2="invasion-down genes in CT4 cells"
	)
ligand_target_pheatmap(scrna_cc1=scrna_cc2,
	geneset_oi=stage_upgenes,
	celltype1="TAM",
	celltype2="CT4",colors="#DA1212",
	filenames="TAM-ct4_up_ligand_target_network.pdf",
	width=4,height=7,
	title1="Prioritized TAM-ligands",
	title2="invasion-up genes in CT4 cells"
	)

ligand_target_pheatmap(scrna_cc1=scrna_cc2,
	geneset_oi=stage_dogenes,
	celltype1="Terminally.Exhausted CD8T",
	celltype2="CT4",colors="#11468F",
	filenames="TE-ct4_down_ligand_target_network.pdf",
	width=5,height=7,
	title1="Prioritized Terminally.Exhausted CD8T-ligands",
	title2="invasion-down genes in CT4 cells"
	)
ligand_target_pheatmap(scrna_cc1=scrna_cc2,
	geneset_oi=stage_upgenes,
	celltype1="Terminally.Exhausted CD8T",
	celltype2="CT4",colors="#DA1212",
	filenames="TE-ct4_up_ligand_target_network.pdf",
	width=7,height=7,
	title1="Prioritized Terminally.Exhausted CD8T-ligands",
	title2="invasion-up genes in CT4 cells"
	)


ligand_target_pheatmap(scrna_cc1=scrna_cc2,
	geneset_oi=stage_dogenes,
	celltype1="NK.cell",
	celltype2="CT4",colors="#11468F",
	filenames="NK-ct4_down_ligand_target_network.pdf",
	width=5,height=7,
	title1="Prioritized NK-ligands",
	title2="invasion-down genes in CT4 cells"
	)
ligand_target_pheatmap(scrna_cc1=scrna_cc2,
	geneset_oi=stage_upgenes,
	celltype1="NK.cell",
	celltype2="CT4",colors="#DA1212",
	filenames="NK-ct4_up_ligand_target_network.pdf",
	width=7,height=7,
	title1="Prioritized NK-ligands",
	title2="invasion-up genes in CT4 cells"
	)

ligand_target_pheatmap(scrna_cc1=scrna_cc2,
	geneset_oi=stage_dogenes,
	celltype1="Neutrophils",
	celltype2="CT6",colors="#11468F",
	filenames="Neutrophils-ct6_down_ligand_target_network.pdf",
	width=4,height=4,
	title1="Prioritized Neutrophils-ligands",
	title2="invasion-down genes in CT6 cells"
	)
ligand_target_pheatmap(scrna_cc1=scrna_cc2,
	geneset_oi=stage_upgenes,
	celltype1="Neutrophils",
	celltype2="CT6",colors="#DA1212",
	filenames="Neutrophils-ct6_up_ligand_target_network.pdf",
	width=4,height=4,
	title1="Prioritized Neutrophils-ligands",
	title2="invasion-up genes in CT6 cells"
	)

ligand_target_pheatmap(scrna_cc1=scrna_cc2,
	geneset_oi=stage_dogenes,
	celltype1="Neutrophils",
	celltype2="CT6",colors="#11468F",
	filenames="Neutrophils-ct6_down_ligand_target_network.pdf",
	width=4,height=4,
	title1="Prioritized Neutrophils-ligands",
	title2="invasion-down genes in CT6 cells"
	)
ligand_target_pheatmap(scrna_cc1=scrna_cc2,
	geneset_oi=stage_upgenes,
	celltype1="Neutrophils",
	celltype2="CT6",colors="#DA1212",
	filenames="Neutrophils-ct6_up_ligand_target_network.pdf",
	width=4,height=4,
	title1="Prioritized Neutrophils-ligands",
	title2="invasion-up genes in CT6 cells"
	)


ligand_target_pheatmap(scrna_cc1=scrna_cc2,
	geneset_oi=stage_dogenes,
	celltype1="Terminally.Exhausted CD8T",
	celltype2="CT6",colors="#11468F",
	filenames="TE-ct6_down_ligand_target_network.pdf",
	width=4,height=7,
	title1="Prioritized Terminally.Exhausted CD8T-ligands",
	title2="invasion-down genes in CT6 cells"
	)
ligand_target_pheatmap(scrna_cc1=scrna_cc2,
	geneset_oi=stage_upgenes,
	celltype1="Terminally.Exhausted CD8T",
	celltype2="CT6",colors="#DA1212",
	filenames="TE-ct6_up_ligand_target_network.pdf",
	width=4,height=7,
	title1="Prioritized Terminally.Exhausted CD8T-ligands",
	title2="invasion-up genes in CT6 cells"
	)
ligand_target_pheatmap(scrna_cc1=scrna_cc2,
	geneset_oi=stage_dogenes,
	celltype1="TAM",
	celltype2="CT6",colors="#11468F",
	filenames="TAM-ct6_down_ligand_target_network.pdf",
	width=4,height=7,
	title1="Prioritized TAM-ligands",
	title2="invasion-down genes in CT6 cells"
	)
ligand_target_pheatmap(scrna_cc1=scrna_cc2,
	geneset_oi=stage_upgenes,
	celltype1="TAM",
	celltype2="CT6",colors="#DA1212",
	filenames="TAM-ct6_up_ligand_target_network.pdf",
	width=4,height=7,
	title1="Prioritized TAM-ligands",
	title2="invasion-up genes in CT6 cells"
	)

ligand_target_pheatmap(scrna_cc1=scrna_cc2,
	geneset_oi=stage_dogenes,
	celltype1="Cytotoxic.CD8.T",
	celltype2="CT6",colors="#11468F",
	filenames="Cytotoxic.CD8.T-ct6_down_ligand_target_network.pdf",
	width=4,height=7,
	title1="Prioritized Cytotoxic.CD8.T-ligands",
	title2="invasion-down genes in CT6 cells"
	)
ligand_target_pheatmap(scrna_cc1=scrna_cc2,
	geneset_oi=stage_upgenes,
	celltype1="Cytotoxic.CD8.T",
	celltype2="CT6",colors="#DA1212",
	filenames="Cytotoxic.CD8.T-ct6_up_ligand_target_network.pdf",
	width=4,height=7,
	title1="Prioritized Cytotoxic.CD8.T-ligands",
	title2="invasion-up genes in CT6 cells"
	)
#scrna_cancer<-merge(scrna_cc1,scrna_cc2)
#head(scrna_cancer@meta.data)

#nichenet_Neutrophils_CC1vsCC2 = nichenet_seuratobj_aggregate(
#  seurat_obj = scrna_cancer, 
#  receiver = "Neutrophils", 
#  condition_colname = "samples", condition_oi = "CC1", condition_reference = "CC2", 
#  sender = c("B.cell","B.cell.activite.memory","CD4.TCM","Cytotoxic.CD8.T","NK.cell","TAM","Terminally.Exhausted CD8T"), 
#  ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, 
#  weighted_networks = weighted_networks, organism = "human")

#############"Cytotoxic.CD8.T"
# receiver
receiver = "Cytotoxic.CD8.T"
expressed_genes_receiver = get_expressed_genes(ident=receiver,scrna_cancer, pct = 0.10)
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
## sender
sender_celltypes = c("B.cell","B.cell.activite.memory","CD4.TCM",
	"Neutrophils","NK.cell","TAM","Terminally.Exhausted CD8T")
list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, scrna_cc1, 0.10) # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()

seurat_obj_receiver= subset(scrna_cc1, idents = receiver)
seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver$samples)

condition_oi = "CC1"
condition_reference = "CC2" 

DE_table_receiver = FindMarkers(object = seurat_obj_receiver, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = 0.10) %>% rownames_to_column("gene")

geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

