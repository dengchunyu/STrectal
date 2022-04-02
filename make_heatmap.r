

make_heatmap_df <- function(obj, marker_genes, include_zeros = FALSE, sort_genes = FALSE, sort_clusters = TRUE) {
  marker_genes <- marker_genes
  obj <- obj
  cluster <- factor(BRCA_merged@meta.data$seurat_clusters)
  #timepoint_ordered <- ordered(timepoint_factor, levels = unique(obj$timepoint))
  cell_ident_df <- data.frame(#timepoint = timepoint_ordered,
                              cluster = Idents(obj))
  count_table <- table(#cell_ident_df$timepoint, 
    cell_ident_df$cluster)
  prop_celltypes <- round(prop.table(count_table, 1) * 100, 1)

  counts_markers <- t(FetchData(object = obj, slot = "data", vars = marker_genes))
  #rownames(counts_markers) <- gsub("^rna_","",rownames(counts_markers))

  heatmap_df <- data.frame()
  # Todo: Replace this with an apply
  for (cluster in levels(Idents(obj))) {
    
      if (as.data.frame.matrix(count_table)[cluster] > 3) {

        cell_expression <- melt(t(counts_markers[, Idents(obj) == cluster]))
        colnames(cell_expression) <- c("cell", "gene", "expression")
        n_rep <- nrow(cell_expression)
        cell_expression[,'cluster'] <- rep(cluster, n_rep)
        
        cell_expression[,'prop'] <- rep(as.data.frame.matrix(prop_celltypes)[cluster] + 1, n_rep)

        heatmap_df <- rbind(heatmap_df, cell_expression)
      } else if (include_zeros){
        gene_rep <- length(marker_genes)
        zero_df <- data.frame(cell = rep("NA", gene_rep), gene = marker_genes, expression = rep(0, gene_rep),
                              cluster = rep(cluster, gene_rep), time = rep(time, gene_rep), prop = rep(0, gene_rep))
        heatmap_df <- rbind(heatmap_df, zero_df)
      }
    }
  }


  # Use hclust to sort the genes by complete linkage
  if (sort_genes) {
    expression_summary_wide <- acast(heatmap_df, cluster+time~gene, value.var = "expression")
    hc <- hclust(dist(t(expression_summary_wide)))
    heatmap_df$gene <- ordered(as.factor(heatmap_df$gene), colnames(expression_summary_wide)[hc$order])
  } else {
    #Order by defined order..
    heatmap_df$gene <- ordered(as.factor(heatmap_df$gene), marker_genes)
  }

  if (sort_clusters) {
    # Use hclust to sort the clusters by complete linkage
    expression_summary_wide_cluster <- acast(heatmap_df, gene ~ cluster, value.var = "expression", fun.aggregate = mean)
    expression_summary_wide_cluster[is.na(expression_summary_wide_cluster)] <- 0
    hc_cluster <- hclust(dist(t(expression_summary_wide_cluster)))
    heatmap_df$cluster <- ordered(as.factor(heatmap_df$cluster), colnames(expression_summary_wide_cluster)[hc_cluster$order])
  } else {
    heatmap_df$cluster <- ordered(as.factor(heatmap_df$cluster), unique(heatmap_df$cluster))
  }

  return(heatmap_df)
}



heatmapRearrange <- function(data_m, x_index = "", x_sequence, y_index = "", y_sequence,
                             fill_index = "", is.text, title = "", xlab = NULL, ylab = NULL, pdffile = "", pdfwidth, pdfheight){
  #1.进行行列重排
  if(!is.null(x_sequence)){
    data_m <- data_m[which(data_m[, x_index] %in% x_sequence),]
    data_m[, x_index] <- factor(data_m[, x_index], levels = x_sequence)
  }
  if(!is.null(y_sequence)){
    data_m <- data_m[which(data_m[, y_index] %in% y_sequence),]
    data_m[, y_index] <- factor(data_m[, y_index], levels = y_sequence)
  }
  
  #2.是否在热图色块上面显示数值
  p <- ggplot(data = data_m, aes(x = data_m[, x_index], y = data_m[, y_index], fill = data_m[, fill_index])) + 
    geom_tile(color="white") + 
    scale_fill_gradient(low = "white", high = "blue", limit = c(0.5,1), na.value = "white") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size=12), axis.text.y = element_text(size=12))+ 
    labs(x = xlab, y = ylab, title = title)
  
  if(is.text){
    p <- p + geom_text(aes(x = data_m[, x_index], y = data_m[, y_index], label = data_m[, fill_index]), color = "black", size = 5)
  }
  
  pdf(file = pdffile, width = pdfwidth, height = pdfheight)
  print(p)
  dev.off()
  
  return(p)
}
