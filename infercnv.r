#options(future.globals.maxSize = 10000 * 1024^2)
#file_pre="/share/pub/dengcy/BRCA_organoid/RawData/Matrix/"
setwd("/share/pub/dengcy/BRCA_organoid/runtime/1.dataProgress")
source("/share/pub/dengcy/BRCA_organoid/runtime/src/Data_input_progress.r")
source("/share/pub/dengcy/BRCA_organoid/runtime/src/Functions.r")
source("/share/pub/dengcy/BRCA_organoid/runtime/src/library.r")
output_pre="/share/pub/dengcy/BRCA_organoid/runtime/1.dataProgress/"
setwd("/share/pub/dengcy/BRCA_organoid/runtime/1.dataProgress")
load("BRCA_merged.RData")
library(infercnv)
all_fortify <- fortify.Seurat(BRCA_merged)
counts_matrix<-GetAssayData(object =BRCA_merged, slot = "counts")
celltype<-rep("PTPRC-",ncol(counts_matrix))
celltype[counts_matrix["PTPRC",]!=0]<-"PTPRC+"

cellAnnota <- data.frame(celltype=celltype,seurat_clusters=all_fortify$seurat_clusters,Type=all_fortify$Type)

rownames(cellAnnota)<-colnames(counts_matrix)
load("/share/pub/dengcy/Singlecell/COVID19/1.rolypoly_result/geneid_df1.RData")
geneid_df1$chrom<-paste0("chr",geneid_df1$chrom)
geneLocate<-geneid_df1[,c(1,2,3)]
rownames(geneLocate)<-geneid_df1$label

infercnv_obj = infercnv::CreateInfercnvObject(
                  raw_counts_matrix = counts_matrix,
                  annotations_file = cellAnnota,
                  gene_order_file = geneLocate,
                  ref_group_names = c("PTPRC+"))

infercnv_obj<-infercnv::run(infercnv_obj,cutoff=0.1,cluster_by_groups=T,denoise=TRUE,HMM=TRUE,
  out_dir="/share/pub/dengcy/BRCA_organoid/runtime/1.dataProgress/infercnv",no_prelim_plot = T,png_res = 300,num_threads=40)
save(infercnv_obj,file="infercnv_obj2.0.RData")
