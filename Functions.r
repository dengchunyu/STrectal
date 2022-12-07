#########################################################
#FUNCTION
cellchat_single<-function(cellchat=cellchat,filename="cc1"){

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
#future::plan("multiprocess", workers = 4) # do parallel

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# project gene expression data onto PPI network (optional)
cellchat <- projectData(cellchat, PPI.human)

# Inference of cell-cell communication network
cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
#df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5))
#df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb"))
cellchat <- computeCommunProbPathway(cellchat,thresh = 0.05)
cellchat <- aggregateNet(cellchat)

#paste0("netVisual_circle_",filename,".pdf")
pdf(paste0("netVisual_circle_",filename,".pdf"))
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

mat <- cellchat@net$weight

pdf(paste0("netVisual_eachgroup_circle_",filename,".pdf"),width=10,height=10)
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()
#Visualization of cell-cell communication network
#Visualize each signaling pathway using Hierarchy plot, Circle plot or Chord diagram

# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')

return(cellchat)
}


#' Run the common Seurat clustring, PCA, and UMAP
#'
#' @param dims_umap dimensions to run UMAP
#' @param dims_neighbors How many dims for neighbor finding
#' @param k_param k for finding neighbors
#' @param cluster_res Resolution for clustering
#' @return Seurat object
cluster_pca_umap <- function(obj,assay=NULL, reduction,cluster_res = 0.3){
  #obj2 <- RunPCA(obj, assay = "SCT", reduction = "harmony",verbose = F)
  obj2 <- RunTSNE(object = obj,assay = assay, reduction = reduction, dims = 1:50,check_duplicates = FALSE)
  obj2 <- RunUMAP(object = obj2, assay =assay, reduction = reduction, dims = 1:50,check_duplicates = FALSE)
  obj2 <- FindNeighbors(object=obj2, assay = assay, reduction = reduction, dims = 1:50)
  obj2 <- FindClusters(object=obj2, resolution = cluster_res)
  return(obj2)
}
####################
##

pred_final<-function(pred.brca=pred.brca,CRC_scrna=CRC_scrna){
  singler.results <- merge(data.frame(cell = rownames(pred.brca), singler = pred.brca$labels), 
                         data.frame(cell = rownames(CRC_scrna@meta.data), 
                                    cluster = CRC_scrna@meta.data$seurat_clusters), 
                         by = "cell", 
                         all.y = FALSE)
singler.results$cell <- NULL
singler.results$count <- 1
singler.results <- aggregate(count ~ ., singler.results, FUN = sum)
singler.final <- singler.results %>% group_by(cluster) %>% top_n(n = 1, wt = count)
return(singler.final)
}

######vinplot 画图代码
##感谢https://www.jianshu.com/p/db0d0927e613
modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.text.y = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.ticks.y = element_line(),
          axis.title.y = element_text(size = rel(1), angle = 0, vjust = 0.5),
          plot.margin = plot.margin )
  return(p)
}

## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0,
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.ticks.x = element_line())
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}
#' This can use a ton of memory, but it will save a lot of time
#' Important: The output is the same order as levels(Idents(obj))
#'
#' This needs to export the Seurat object to all workers - this takes a lot of RAM.
#' Default
#'       logfc.threshold = 0.25,
#'       test.use = "wilcox",
#'       min.pct = 0.1,
#'       min.diff.pct = -Inf,
#'       verbose = TRUE,
#'       only.pos = FALSE,
#'       max.cells.per.ident = Inf,
#'       random.seed = 1,
#'      latent.vars = NULL,
#'      min.cells.feature = 3,
#'       min.cells.group = 3,

#' @param obj Seurat object
#' @param n_cor number of CPU cores
#' @return List of data tables, one for each numeric cluster
parallelFindAllMarkers <- function(obj){

  all_markers <- lapply(levels(Idents(obj)), function(x){ # Five expression patterns
    FindMarkers(obj, ident.1 = x, ident.2 = NULL, test.use = "wilcox")
  })

  return(value(all_markers))
}
#'
fortify.Seurat <- function(x){
  xy <- as.data.frame(Embeddings(x, reduction = "umap"))
  colnames(xy) <- c("UMAP_1", "UMAP_2")
  xy$UMAP_1 <- as.numeric(xy$UMAP_1)
  xy$UMAP_2 <- as.numeric(xy$UMAP_2)

  xy2 <- as.data.frame(Embeddings(x, reduction = "tsne"))
  colnames(xy2) <- c("TSNE_1", "TSNE_2")
  xy2$TSNE_1 <- as.numeric(xy2$TSNE_1)
  xy2$TSNE_2 <- as.numeric(xy2$TSNE_2)
  xy<-cbind(xy,xy2)
  return(cbind(xy, as.data.frame(x@meta.data)))
}
#'
umap_theme <- function(){
  theme_grey() %+replace%
    theme(panel.background = element_rect(fill = "white", colour = "white", size = 0.1),
          axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y  = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.key=element_blank())
}

featureplot_theme <- function(){
  theme_grey() %+replace%
    theme(panel.background = element_rect(fill = "white", colour = "black", size = 1),
          axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y  = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
}


prism_theme <- function(){
  theme_grey() %+replace%
    theme(panel.border = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_line(),
          panel.grid.minor = element_blank(),
          legend.key=element_blank(),
          axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
          axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
          axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 0.95),
          axis.text.y = element_text(size = 14),
          axis.ticks = element_line(size = 0.5),
          axis.ticks.length = unit(0.25, "cm"),
          axis.title.y = element_text(size = 16, vjust = 4, angle = 90),
          axis.title.x = element_blank(),
          plot.margin = margin(t = 5, r = 5, b = 5, l = 20),
          legend.text = element_text(size = 14),
          legend.key.size = unit(2,"line"),
          legend.title = element_blank())
}

singler_func<-function(CRC_scrna,ref){
common <- intersect(rownames(ref), rownames(CRC_scrna))
CRC_scrna.singler <- CRC_scrna[common,]
ref <- ref[common,]
CRC_scrna.singler.sce <- SingleCellExperiment(assays = list(counts = CRC_scrna.singler))
CRC_scrna.singler.sce <-SingleCellExperiment(assays = List(counts = GetAssayData(object =CRC_scrna.singler, slot = "counts")))
CRC_scrna.singler.sce <- logNormCounts(CRC_scrna.singler.sce)

pred.brca <- SingleR(test = CRC_scrna.singler.sce, 
  ref = ref, assay.type.test=1,labels = ref$label.main)
CRC_scrna$singler <- pred.brca$labels
return(CRC_scrna)
}

color_scanpy_patient <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22")

color_scanpy_viridis28 <- c("#D9DD6B","#ECEFA4","#D54C4C","#8D2828","#FDD2BF","#E98580","#DF5E5E","#492F10","#334257","#476072","#548CA8",
"#00A19D","#ECD662","#5D8233","#284E78","#3E215D","#835151","#F08FC0","#C6B4CE","#BB8760","#FFDADA","#3C5186",
"#558776","#E99497","#FFBD9B","#0A1D37","#01937C","#464660","#368B85")
color_scanpy_13 <- c("#FDD2BF","#DF5E5E","#492F10","#753422","#628395","#262A53","#ECD662","#5D8233","#284E78","#3E215D","#DF711B","#FFB740","#64C9CF")

color_scanpy_type3 <-c("#D9DD6B","#D54C4C","#548CA8")
color_scanpy_8 <- c("#FDD2BF","#753422","#297F87","#DF5E5E","#628395","#262A53","#ECD662","#5D8233")

set.seed(42) # For reproducability

N_WORKERS <- 2
n_print <- 1:20
options(future.globals.maxSize=15*1024*1024^2)


saveTiff <- function(path, image, width = 5, height = 5, dpi = 300, units = "in"){
  if (!file.exists(path)){
    dir.create(dirname(path), showWarnings = FALSE)
  }

  if (Sys.info()["sysname"]=='Darwin'){
    # lzw doesn't work on mac with quartz
    tmp_path <- suppressWarnings(normalizePath(paste0(path, "_tmp.tiff")))
    out_path <- suppressWarnings(normalizePath(path))

    tiff(tmp_path, width = width, height = height, units = units, res = dpi, compression = "lzw", type = "quartz")
    print(image)
    dev.off()
    # requires imagemagick
    Sys.sleep(0.5)
    system(paste0("convert ", tmp_path, " -density ", dpi,  " -resize ", width * dpi, "x", height * dpi, "\\> -compress lzw ", out_path), ignore.stdout = TRUE)

    if (file.exists(out_path)) {
      #Delete file if it exists
      file.remove(tmp_path)
    }

  } else {
    tiff(path, width = width, height = height, units = units, res = dpi, compression = "lzw")
    print(image)
    dev.off()
  }

}


renameFactorIdent <- function(obj, start, end){
  levels(obj)[levels(obj) == start] <- end
  return(obj)
}

#' Make a wide DF that has some summary data about the clusters, useful for heatmaps
#'
#' @param obj seurat object to run on, needs to have a $timepoint variable
#' @param marker_genes List of marker genes to include
#' @param include_zeros Include zero expression clusters
#' @param sort_genes Use hclust to sort the order of the genes
#' @param sort_clusters Use hclust to sort the order of the clusters
#' @return Returns Wide dataframe
#' @examples
#'
#'heatmap_df <- make_heatmap_df(meso, marker_genes)
#'
#'
#'heatmap_df <- heatmap_df %>% group_by(gene,cluster) %>% summarise(expression_mean=mean(expression))
#'ggplot(heatmap_df, aes(x = gene, y = cluster, fill = expression_mean)) +
#'        geom_tile(color = "white", size = 0.1) +
#'        scale_fill_distiller(palette = "Blues", direction = 1, trans = "sqrt") +
#'        ggtitle("Expression") +
#'        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#'        coord_equal() +


#' plot the plot for proportion of cell types
#'
#' @param obj seurat object to run on
#' @param options 
#' @param 
#' @param colors

stack_plot<-function(obj,options,colors,coord=T,angle=T,label.size=4){
    ggplot_df<-obj@meta.data[,options]
    colnames(ggplot_df)<-c("index1","index2")

a<-tapply(ggplot_df$index1,factor(ggplot_df$index2),function(index2){table(index2) })
dfll <- do.call(rbind,lapply(a, data.frame))
dfll$index1 <- rep(names(a),sapply(a,length))
# Calculate the percentages
dfll = plyr::ddply(dfll, .(index1), transform, percent = Freq/sum(Freq) * 100)

# Format the labels and calculate their positions
dfll = plyr::ddply(dfll, .(index1), transform, pos = cumsum(percent)/100)
dfll$label = paste0(sprintf("%.0f", dfll$percent), "%")
dfll$index2<-factor(dfll$index2,levels = rev(unique(as.vector(dfll$index2))))
p<-ggplot(dfll, aes(x =index1, y = percent, fill =index2)) +
  geom_bar(stat = "identity", position= 'fill',size=0.3) +
  geom_text(aes(y = pos, label = label), size = label.size)+
  theme_classic()+
        scale_fill_manual(values=colors)+
        #
        labs(x = options[1], title = options[2])
  if(angle){
    p<-p+theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1))
  }
   if(coord){
    p<-p+coord_flip()
   }     
    #p<-ggplot(data = dfll, mapping = aes(x=index1,y=Freq,fill =factor(index2))) + 
    #geom_bar(stat = 'identity',position = 'fill',color="grey",size=0)+
return(p)
}

###空间转录组得标准化流程
#
normalize_spacial<-function(Spatial_data,lowqspot=0.02,mitper=25,geneExprMin=15,spot_meta="RST2bei"){
    ###标准化
    ##Quantification of spots and genes were carried out by using subset function. 
    #We filtered out ~2% of low-quality spots. Spots with over 25% mitochondrial gene expression were also discarded
    feature_nums <- colSums(as.matrix(Spatial_data@assays$Spatial@counts) >0)
    mycut_feature <- as.numeric(quantile(feature_nums, lowqspot))
    count_nums <- colSums(as.matrix(Spatial_data@assays$Spatial@counts))
    mycut_count <- as.numeric(quantile(count_nums, lowqspot))
    Spatial_data_f1 <- subset(Spatial_data, subset = nFeature_Spatial>=mycut_feature | nCount_Spatial>=mycut_count)# & percent.mt < mitper )
    message("cut the count is success!")

    #Genes expressed in fewer than 15 spots were excluded.
    min.spots <- geneExprMin
    num.spots <- rowSums(as.matrix(Spatial_data@assays$Spatial@counts) > 0)
    genes.use <- names(num.spots[which(num.spots >= min.spots)])
    mykeepgene <- c(1:nrow(Spatial_data_f1))[rownames(Spatial_data_f1)%in%as.character(genes.use)]
    Spatial_data_f2 <- subset(Spatial_data_f1,features=mykeepgene)
    message("cut the spot is success!")
    # Filter out contaminated genes with specified names.
    #Genes related to hemoglobin (considerable variation from blood contents) and Y-chromosome linked genes were removed.
    filter_genes <-c("MALAT1", "SLC4A1", "KDM5D", "ANK1", "DDX3Y", "EIF2AK1", "HBQ1", "FTL", "GATA1", "KLF1", "USP9Y", "NFE2", "MT1G", "RPS4Y1", "HBZ", "GYPC", "HEMGN", "SLC25A37", "ALAS2", "EPB41", "AHSP", "GYPA", "UTY", "HBA2", "HBG2", "EIF1AY", "HBA1", "HBM", "HBE1", "HBG1", "MTRNR2L4", "HBB", "MTRNR2L5", "MTRNR2L8", "MTRNR2L10", "MTRNR2L3", "MTRNR2L1", "MTRNR2L7", "MTRNR2L12", "MTRNR2L11", "MTRNR2L13", "MTRNR2L6")
    keepgenes <- c(1:nrow(Spatial_data_f2))[!(rownames(Spatial_data_f2)%in%as.character(filter_genes))]
    Spatial_data_f3 <- subset(Spatial_data_f2,features=keepgenes)

    write.csv(Spatial_data_f3@meta.data, file = paste0(spot_meta,"-spots-metadata-clean.csv")) # Export the clean spots metadata
    message("cut the gene is success!")
    # SCT normalization
    #The clean expression matrix data were normalized using regularized negative binomial regression.
    Spatial_data_f3  <- SCTransform(Spatial_data_f3 , assay = "Spatial", return.only.var.genes = FALSE)
    DefaultAssay(Spatial_data_f3) <- "SCT"
    write.csv(Spatial_data_f3@meta.data, file = paste0(spot_meta,"-spots-metadata-clean-SCT.csv"))
    message("SCT normalization is success!")
    # Dimensionality reduction (PCA)
    #PCA was performed and the 10 most significant components was determined by the DimHeatmap and ElbowPlot function. 
    Spatial_data_f3 <- RunPCA(Spatial_data_f3)
    message("runpca is success!")
    pdf(paste0(spot_meta,"_dimheatmap.pdf"))
    DimHeatmap(Spatial_data_f3, dims = 1:15, cells = 2000, balanced = TRUE)
    dev.off()
    message("DimHeatmap is success!")
    p1 <- ElbowPlot(Spatial_data_f3)
    ggsave(paste0(spot_meta,"-PCA-ElbowPlot.pdf"), plot=p1, width = 5, height =5)
    message("ElbowPlot is success!")
    return(Spatial_data_f3)
}



ggplot_ScatterPlot1<-function(data_m,x_index="",y_index="",size_index="",color_index="",Theme="horizon",title="",shape_index=NULL){
  ##
  library(ggthemes)
  
    if(is.null(size_index)){
      p1<-ggplot(data_m, aes(x=data_m[,x_index], y=data_m[,y_index], color=data_m[,color_index]))
      if(is.null(shape_index)){
        p2<-geom_point(shape=20,alpha = 0.7,size=3)
      }else{
        p2<-geom_point(aes(shape=data_m[,shape_index]),alpha = 0.7,size=3)
      }
      
    }else{
      p1<-ggplot(data_m, aes(x=data_m[,x_index], y=data_m[,y_index], size=data_m[,size_index], color=data_m[,color_index]))
      if(is.null(shape_index)){
      p2<-geom_point(shape=20,alpha = 0.7)
      }else{
        p2<-geom_point(aes(shape=data_m[,shape_index]),alpha = 0.7) 
      }
    }
  
  p<-p1+p2
  # if(!isTRUE(legend)){
  #   p<-p+ theme(legend.position="none")
  # }
  ##是否加入x轴虚线
   library(ggthemes)
  if(Theme=="classic"){
    p<-p+theme_classic()##无网格主图
  }else if(Theme=="horizon"){
    p<-p+theme_calc()##横线主题
  }else if(Theme=="grid"){
    p<-p+theme_gdocs()
  }
  p<-p+labs(x = x_index,y=y_index,title = title)
  
  # p4<-geom_hline(aes(yintercept=), colour="#990000", linetype="dashed",alpha = 0.6)
 
  return(p)
}

