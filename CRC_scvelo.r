#########################
#velocyto.R读取文件
#BiocManager::install("hdf5r")    #安装包
#library(hdf5r)
library(Seurat)
library(loomR)
#data_sample <- Read10X_h5("Women/GSE118127_RAW/GSM3319032_sample_1-1_filtered_gene_bc_matrices_h5.h5")  #导入数据
#data_seurat <- CreateSeuratObject(data_sample,project = "data_sample") #后面就可以单细胞处理的标准流程啦
################
#读取原始文件转化为loom文件
################
setwd("/share/pub/dengcy/Singlecell/CRC_space/6.scvelo")
file_pre="/share/pub/dengcy/Singlecell/CRC_space/gcy_analysis/Analysis/1.Basic_analysis/1.1.raw_feature_bc_matrix"
cell_cycle_genes <- read.table(file = "/share/pub/xiongyc/project/scRNA/JiangFanChen/data/cc_genes_hg.txt",header = TRUE);
  #### step1.1: create Seurat object
  data.10x<- Read10X(data.dir = paste(file_pre,sep="")); 
  #CRC_scrna = CreateSeuratObject(counts = data.10x, min.cells=0, min.features=0, project=projectName[file_pre]);
  CRC_scrna_raw = CreateSeuratObject(counts = data.10x,min.cells=3, min.features=200,plot_feature=F,cell_cycle_genes=cell_cycle_genes, project="CRC_scRNA-seq");

  sdata.loom <- as.loom(x = CRC_scrna_raw, filename = "CRC_scrna_raw.loom", verbose = FALSE)# Always remember to close loom files when done
  sdata.loom$close_all()