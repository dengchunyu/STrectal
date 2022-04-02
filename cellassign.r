setwd("/share/pub/dengcy/BRCA_organoid/runtime/1.dataProgress")

source("/share/pub/dengcy/BRCA_organoid/runtime/src/Functions.r")
source("/share/pub/dengcy/BRCA_organoid/runtime/src/library.r")
setwd("/share/pub/dengcy/BRCA_organoid/runtime/1.dataProgress")
load("BRCA_merged.RData")

library(SingleCellExperiment)
library(cellassign)
CancerDiffGene<-read.table("BrcaCancerDiffGeneTop50.txt",header=T,sep="\t")
immune_marker<-read.table("Cancer_Res_immune_marker_30773503.txt",header=T,sep="\t")
cancerstemcell<-read.csv("CellMarker_cancerstemcell.csv")
epicell<-read.csv("CellMarker_epicell.csv")
progenitor<-read.csv("CellMarker_progenitor.csv")

rownames(immune_marker)<-immune_marker$gene
immune_marker<-immune_marker[intersect(rownames(BRCA_merged),immune_marker$gene),]
immune_marker<-tapply(immune_marker$gene,factor(immune_marker$cell_type),function(cell){
   return(cell)
  })

immune_marker<-immune_marker[c("Activated_B_cells","Central_memory_CD4","Central_memory_CD8",
"Cytotoxic_cells","DC","Effector_memory_CD8","Eosinophil","iDC","Immature_B_cells",
"Macrophages","Mast_cells","MDSC","Memory_B_cells","Monocytes","Neutrophils",
"NK","NKT" ,"T_cells","TGD","Th1","Th17","Th2","Treg")]


BRCATME_markers<-list(CancerDiffGene=intersect(rownames(BRCA_merged),CancerDiffGene$Gene.Symbol),
  cancerstemcell=intersect(rownames(BRCA_merged),cancerstemcell$Cell.Marker),
  epicell=intersect(rownames(BRCA_merged),epicell$Cell.Marker),progenitor=intersect(rownames(BRCA_merged),progenitor$Cell.Marker))
BRCATME_markers<-c(BRCATME_markers,immune_marker)
marker_mat <- marker_list_to_mat(BRCATME_markers)
BRCA_merged.sce <-SingleCellExperiment(assays = List(counts = GetAssayData(object =BRCA_merged, slot = "counts")))
library(scran)
s_factor <-calculateSumFactors(BRCA_merged.sce, assay.type = "counts")

fit <- cellassign(exprs_obj = BRCA_merged.sce[rownames(marker_mat),],
  marker_gene_info = marker_mat,s = s_factor,
  learning_rate = 1e-2, shrinkage = TRUE,verbose = FALSE)
save(fit,file="cellassign_fit.RData")
load("cellassign_fit.RData")
assign_celltype<-celltypes(fit)
names(assign_celltype)<-colnames(BRCA_merged)

all_fortify$cellassign_celltype<-assign_celltype
all_fortify[,c("seurat_clusters","cellassign_celltype")]
tapply(all_fortify[,"cellassign_celltype"],factor(all_fortify[,"seurat_clusters"]),function(assign) sort(table(assign),decreasing=T))
