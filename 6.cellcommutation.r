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
setwd("/share/pub/dengcy/Singlecell/CC_space/6.cellcom")
load("/share/pub/dengcy/Singlecell/CC_space/3.Epi_cells_stage/cancer_scrna_stage1.RData")
load("/share/pub/dengcy/Singlecell/CC_space/3.Epi_cells_stage/cancer_scrna_stage2.RData")

load("/share/pub/dengcy/Singlecell/CC_space/3.2.Immune_cells_stage/Immune_scrna.RData")
load("scrna_stage1.RData")
load("scrna_stage1.4.RData")
load("scrna_stage2.RData")

########################
Immune_scrna_cc1<- subset(Immune_scrna,stage %in% c("StageI","StageII") )
Immune_scrna_cc2<- subset(Immune_scrna,stage == "StageIII" )

sc5<- subset(cancer_scrna_stage1,State == 5)
sc5$annotation<-"StageI/II-State5"
scrna_stage1<-merge(Immune_scrna_cc1,sc5)
Idents(scrna_stage1)<-scrna_stage1$annotation
save(scrna_stage1,file="scrna_stage1.RData")

sc4<- subset(cancer_scrna_stage1,State == 4)
sc4$annotation<-"StageI/II-State4"
scrna_stage1.4<-merge(Immune_scrna_cc1,sc4)
Idents(scrna_stage1.4)<-scrna_stage1$annotation
save(scrna_stage1.4,file="scrna_stage1.4.RData")

Immune_scrna_cc2<- subset(Immune_scrna,stage == "StageIII" )
sc3<- subset(cancer_scrna_stage2,State ==3)
sc3$annotation<-"StageIII-State3"
scrna_stage2<-merge(Immune_scrna_cc2,sc3)

Idents(scrna_stage2)<-scrna_stage2$annotation
save(scrna_stage2,file="scrna_stage2.RData")

##############
Immune_scrna_cc1<- subset(Immune_scrna,stage %in% c("StageI","StageII") )
sc2<- subset(cancer_scrna_stage1,Cluster == 2)
sc2$annotation<-"StageI/II-Cluster2"
scrna_stage1_c2<-merge(Immune_scrna_cc1,sc2)
Idents(scrna_stage1_c2)<-scrna_stage1_c2$annotation
save(scrna_stage1_c2,file="scrna_stage1_c2.RData")

sc7<- subset(cancer_scrna_stage1,Cluster == 7)
sc7$annotation<-"StageI/II-Cluster7"
scrna_stage1_c7<-merge(Immune_scrna_cc1,sc7)
Idents(scrna_stage1_c7)<-scrna_stage1_c7$annotation
save(scrna_stage1_c7,file="scrna_stage1_c7.RData")

Immune_scrna_cc2<- subset(Immune_scrna,stage == "StageIII" )
sc3<- subset(cancer_scrna_stage2,Cluster ==3)
sc3$annotation<-"StageIII-Cluster3"
scrna_stage2_c3<-merge(Immune_scrna2_cc2,sc3)
Idents(scrna_stage2_c3)<-scrna_stage2_c3$annotation
save(scrna_stage2_c3,file="scrna_stage2_c3.RData")

sc7<- subset(cancer_scrna_stage2,Cluster ==7)
sc7$annotation<-"StageIII-Cluster7"
scrna_stage2_c7<-merge(Immune_scrna_cc2,sc7)
Idents(scrna_stage2_c7)<-scrna_stage2_c7$annotation
save(scrna_stage2_c7,file="scrna_stage2_c7.RData")

#####################
########################
#cellchat

library(CellChat)
library(ggplot2)
library(patchwork)
library(igraph)
options(stringsAsFactors = FALSE)

cellchat <- createCellChat(object = scrna_stage1,group.by="annotation")
cellchat <- setIdent(cellchat, ident.use = "annotation")
cellchat <- cellchat_single(cellchat=cellchat,filename="stage1_immune")
pdf("stage1.5_netVisual_bubble.pdf",width=4,height=5)
netVisual_bubble(cellchat, sources.use =2, targets.use = c(1,3,4,5,6,7,8,9), remove.isolate = FALSE)
dev.off()
save(cellchat,file="stage1.5_cellchat.RData")

cellchat <- createCellChat(object = scrna_stage1.4,group.by="annotation")
cellchat <- setIdent(cellchat, ident.use = "annotation")
cellchat <- cellchat_single(cellchat=cellchat,filename="stage1.4_immune")
cellchat@net
pdf("stage1.4_netVisual_bubble.pdf",width=4,height=3)
netVisual_bubble(cellchat, sources.use =4, targets.use = c(1,2,3,5,6,7,8,9), remove.isolate = FALSE)
dev.off()
save(cellchat,file="stage1.4_cellchat.RData")

cellchat <- createCellChat(object = scrna_stage2,group.by="annotation")
cellchat <- setIdent(cellchat, ident.use = "annotation")
cellchat<-cellchat_single(cellchat=cellchat,filename="stage2_immune")
save(cellchat,file="stage2_cellchat.RData")
cellchat@net
pdf("stage2_netVisual_bubble.pdf",width=4,height=5)
netVisual_bubble(cellchat, sources.use = 2, targets.use = c(1,3,4,5,6,7,8,9), remove.isolate = FALSE)
dev.off()


cellchat <- createCellChat(object = scrna_stage1_c2,group.by="annotation")
cellchat <- setIdent(cellchat, ident.use = "annotation")
cellchat <- cellchat_single(cellchat=cellchat,filename="stage1_c2")
save(cellchat,file="stage1_c2_cellchat.RData")
pdf("stage1_c2_netVisual_bubble.pdf",width=4,height=5)
netVisual_bubble(cellchat, sources.use =4, targets.use = c(1,2,3,5,6,7,8,9), remove.isolate = FALSE)
dev.off()

cellchat <- createCellChat(object = scrna_stage1_c7,group.by="annotation")
cellchat <- setIdent(cellchat, ident.use = "annotation")
cellchat <- cellchat_single(cellchat=cellchat,filename="stage1_c7")
save(cellchat,file="stage1_c7_cellchat.RData")
pdf("stage1_c7_netVisual_bubble.pdf",width=4,height=5)
netVisual_bubble(cellchat, sources.use =4, targets.use = c(1,2,3,5,6,7,8,9), remove.isolate = FALSE)
dev.off()

cellchat <- createCellChat(object = scrna_stage2_c3,group.by="annotation")
cellchat <- setIdent(cellchat, ident.use = "annotation")
cellchat <- cellchat_single(cellchat=cellchat,filename="stage2_c3")
save(cellchat,file="stage2_c3_cellchat.RData")
pdf("stage2_c3_netVisual_bubble.pdf",width=4,height=5)
netVisual_bubble(cellchat, sources.use =4, targets.use = c(1,2,3,5,6,7,8,9), remove.isolate = FALSE)
dev.off()


cellchat <- createCellChat(object = scrna_stage2_c7,group.by="annotation")
cellchat <- setIdent(cellchat, ident.use = "annotation")
cellchat <- cellchat_single(cellchat=cellchat,filename="stage2_c7")
save(cellchat,file="stage2_c7_cellchat.RData")
pdf("stage2_c7_netVisual_bubble.pdf",width=4,height=5)
netVisual_bubble(cellchat, sources.use =4, targets.use = c(1,2,3,5,6,7,8,9), remove.isolate = FALSE)
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

load("/share/pub/dengcy/Singlecell/CC_space/4.stage_visualize/emt_gene.RData")

sc5<- subset(cancer_scrna_stage1,State == 5)
sc5$annotation<-"StageI/II-State5"
sc4<- subset(cancer_scrna_stage1,State == 4)
sc4$annotation<-"StageI/II-State4"
scrna_stage1<-merge(Immune_scrna_cc1,list(sc4,sc5))
Idents(scrna_stage1)<-scrna_stage1$annotation

sc3<- subset(cancer_scrna_stage2,Cluster == 3)
sc3$annotation<-"StageIII-Cluster3"
sc7<- subset(cancer_scrna_stage2,Cluster == 7)
sc7$annotation<-"StageIII-Cluster7"
scrna_stage2<-merge(Immune_scrna_cc2,list(sc3,sc7))
Idents(scrna_stage2)<-scrna_stage2$annotation

save(scrna_stage1,file="scrna_stage1_45.RData")
save(scrna_stage2,file="scrna_stage2_37.RData")

load("scrna_stage1_45.RData")
load("scrna_stage2_37.RData")

#########

EMT_genes<-read.delim("/share/pub/dengcy/Singlecell/CC_space/4.stage_visualize/EMTgene.txt")

stage_invasion_genes =intersect(stage_upgenes,rownames(scrna_stage1))
#stage_upgenes =intersect(stage_upgenes,rownames(scrna_cc1))
emt_genes<-intersect(EMT_genes$x,rownames(scrna_stage1))

ligand_target_pheatmap<-function(scrna_cc1=scrna_stage1,
	geneset_oi=stage_invasion_genes,
	celltype1="TAM",
	celltype2="StageI/II-State5",
	colors="#11468F",
	filenames="TAM-Stage1-state5_ligand_target_network.pdf",
	width=4,
	height=4,
	title1="Prioritized TAM",
	title2="Invasion genes in StageI/II-State5"
	){

expression = t(scrna_cc1@assays$RNA@data)
sample_info = scrna_cc1@meta.data # contains meta-information about the cells
sample_info$cellid<-rownames(sample_info)

ids1 = sample_info %>% filter(annotation == celltype1) %>% pull(cellid)
ids2 = sample_info %>% filter(annotation ==celltype2 ) %>% pull(cellid)

expressed_genes_sender = expression[ids1,] %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x,na.rm=T) + 1)}) %>% .[. >= 2] %>% names()

expressed_genes_receiver = expression[ids2,] %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x,na.rm=T) + 1)}) %>% .[. >= 2] %>% names()

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
active_ligand_target_links_df<-na.omit(active_ligand_target_links_df)
ligand_target_matrix<-na.omit(ligand_target_matrix)
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


ligand_target_pheatmap(scrna_cc1=scrna_stage1,
	geneset_oi=stage_invasion_genes,
	celltype1="TAM",
	celltype2="StageI/II-State5",colors="#ed6663",
	filenames="tam_s5_ligand_target_network.pdf",
	width=3,height=4,
	title1="Prioritized TAM-ligands",
	title2="invasion genes in stageI/II"
	)

ligand_target_pheatmap(scrna_cc1=scrna_stage1,
	geneset_oi=stage_invasion_genes,
	celltype1="Neutrophils",
	celltype2="StageI/II-State5",colors="#ed6663",
	filenames="Neu_s5_ligand_target_network.pdf",
	width=3,height=4,
	title1="Prioritized Neutrophils-ligands",
	title2="invasion genes in stageI/II"
	)

ligand_target_pheatmap(scrna_cc1=scrna_stage1,
	geneset_oi=stage_invasion_genes,
	celltype1="Treg",
	celltype2="StageI/II-State5",colors="#ed6663",
	filenames="Treg-s5_down_ligand_target_network.pdf",
	width=3,height=4,
	title1="Prioritized Treg-ligands",
	title2="invasion genes in stageI/II"
	)
ligand_target_pheatmap(scrna_cc1=scrna_stage1,
	geneset_oi=stage_invasion_genes,
	celltype1="Tex",
	celltype2="StageI/II-State5",colors="#ed6663",
	filenames="Tex-s5_ligand_target_network.pdf",
	width=3,height=4,
	title1="Prioritized Tex-ligands",
	title2="invasion-up genes in stageI/II"
	)

ligand_target_pheatmap(scrna_cc1=scrna_stage1,
	geneset_oi=stage_invasion_genes,
	celltype1="TAM",
	celltype2="StageI/II-State4",colors="#ffa372",
	filenames="tam_s4_ligand_target_network.pdf",
	width=3,height=4,
	title1="Prioritized TAM-ligands",
	title2="invasion genes in stageI/II"
	)

ligand_target_pheatmap(scrna_cc1=scrna_stage1,
	geneset_oi=stage_invasion_genes,
	celltype1="Neutrophils",
	celltype2="StageI/II-State4",colors="#ffa372",
	filenames="Neu_s4_ligand_target_network.pdf",
	width=3,height=4,
	title1="Prioritized Neutrophils-ligands",
	title2="invasion genes in stageI/II"
	)

ligand_target_pheatmap(scrna_cc1=scrna_stage1,
	geneset_oi=stage_invasion_genes,
	celltype1="Treg",
	celltype2="StageI/II-State4",colors="#ffa372",
	filenames="Treg-s4_down_ligand_target_network.pdf",
	width=3,height=4,
	title1="Prioritized Treg-ligands",
	title2="invasion genes in stageI/II"
	)
ligand_target_pheatmap(scrna_cc1=scrna_stage1,
	geneset_oi=stage_invasion_genes,
	celltype1="Tex",
	celltype2="StageI/II-State4",colors="#ffa372",
	filenames="Tex-s4_ligand_target_network.pdf",
	width=3,height=4,
	title1="Prioritized Tex-ligands",
	title2="invasion-up genes in stageI/II"
	)
####################

ligand_target_pheatmap(scrna_cc1=scrna_stage2,
	geneset_oi=stage_invasion_genes,
	celltype1="TAM",
	celltype2="StageIII-Cluster3",colors="#4e89ae",
	filenames="tam_c3_ligand_target_network.pdf",
	width=3,height=4,
	title1="Prioritized TAM-ligands",
	title2="invasion genes in stageIII"
	)

ligand_target_pheatmap(scrna_cc1=scrna_stage2,
	geneset_oi=stage_invasion_genes,
	celltype1="Neutrophils",
	celltype2="StageIII-Cluster3",colors="#4e89ae",
	filenames="Neu_c3_ligand_target_network.pdf",
	width=3,height=4,
	title1="Prioritized Neutrophils-ligands",
	title2="invasion genes in stageIII"
	)

ligand_target_pheatmap(scrna_cc1=scrna_stage2,
	geneset_oi=stage_invasion_genes,
	celltype1="Treg",
	celltype2="StageIII-Cluster3",colors="#4e89ae",
	filenames="Treg-c3_down_ligand_target_network.pdf",
	width=3,height=4,
	title1="Prioritized Treg-ligands",
	title2="invasion genes in stageIII"
	)
ligand_target_pheatmap(scrna_cc1=scrna_stage2,
	geneset_oi=stage_invasion_genes,
	celltype1="Tex",
	celltype2="StageIII-Cluster3",colors="#4e89ae",
	filenames="Tex-c3_ligand_target_network.pdf",
	width=3,height=4,
	title1="Prioritized Tex-ligands",
	title2="invasion-up genes in stageIII"
	)


ligand_target_pheatmap(scrna_cc1=scrna_stage2,
	geneset_oi=stage_invasion_genes,
	celltype1="TAM",
	celltype2="StageIII-Cluster7",colors="#3e978b",
	filenames="tam_c7_ligand_target_network.pdf",
	width=3,height=4,
	title1="Prioritized TAM-ligands",
	title2="invasion genes in stageIII"
	)

ligand_target_pheatmap(scrna_cc1=scrna_stage2,
	geneset_oi=stage_invasion_genes,
	celltype1="Neutrophils",
	celltype2="StageIII-Cluster7",colors="#3e978b",
	filenames="Neu_c7_ligand_target_network.pdf",
	width=3,height=4,
	title1="Prioritized Neutrophils-ligands",
	title2="invasion genes in stageIII"
	)

ligand_target_pheatmap(scrna_cc1=scrna_stage2,
	geneset_oi=stage_invasion_genes,
	celltype1="Treg",
	celltype2="StageIII-Cluster7",colors="#3e978b",
	filenames="Treg-c7_down_ligand_target_network.pdf",
	width=3,height=4,
	title1="Prioritized Treg-ligands",
	title2="invasion genes in stageIII"
	)
ligand_target_pheatmap(scrna_cc1=scrna_stage2,
	geneset_oi=stage_invasion_genes,
	celltype1="Tex",
	celltype2="StageIII-Cluster7",colors="#3e978b",
	filenames="Tex-c7_ligand_target_network.pdf",
	width=3,height=4,
	title1="Prioritized Tex-ligands",
	title2="invasion-up genes in stageIII"
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

