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
library(slingshot)
source("/share/pub/dengcy/Singlecell/CC_space/src/Functions.r")
source("/share/pub/dengcy/Singlecell/CC_space/src/library.r")
setwd("/share/pub/dengcy/Singlecell/CC_space/4.stage_visualize")
load("/share/pub/dengcy/Singlecell/CC_space/3.Epi_cells_stage/cancer_scrna.RData")
load("/share/pub/dengcy/Singlecell/CC_space/3.Epi_cells_stage/Slingshot_scrna1.RData")
load("/share/pub/dengcy/Singlecell/CC_space/3.Epi_cells_stage/Slingshot_scrna2.RData")

load("/share/pub/dengcy/Singlecell/CC_space/3.Epi_cells_stage/cancer_scrna_cc1.RData")
load("/share/pub/dengcy/Singlecell/CC_space/3.Epi_cells_stage/cancer_scrna_cc2.RData")

##############
##1.基于不同的基因集合获得stage基因
##################
###1）.pseud time

    color_ptime1<- cancer_scrna_cc1$Pseudotime
    color_ptime2<- cancer_scrna_cc2$Pseudotime

#2）.hallmark genes
  x <- readLines("/share/pub/dengcy/Singlecell/CC_space/3.Epi_cells_stage/h.all.v7.4.symbols.gmt")
  res <- strsplit(x, "\t")
  names(res) <- vapply(res, function(y) y[1], character(1))
  genes.by.pathway.h <- lapply(res, "[", -c(1:2))
  save(genes.by.pathway.h,file="genes.by.pathway.h.RData")
  genes.pathway.h <- unique(unlist(genes.by.pathway.h))

#3）差异表达基因
   ##肿瘤正常差异基因输入p<0.01
  degene<-read.delim("/share/pub/dengcy/Singlecell/CC_space/4.stage_visualize/table_degenes.txt")
  #degene_up<-degene[degene$Log2.Fold.Change.>0,]
  #degene_down<-degene[degene$Log2.Fold.Change.<0,]
#4）EMT基因
  EMT_genes<-read.delim("/share/pub/dengcy/Singlecell/CC_space/4.stage_visualize/EMTgene.txt")
#将以上基因取交集
#5) ct1/ct8 的差异基因
#ct1_markers<-markers[markers$avg_log2FC>0 & markers$cluster %in% "CT1","gene"]
#ct8_markers<-markers[markers$avg_log2FC>0 & markers$cluster %in% "CT8","gene"]

  dehallgene<-intersect(degene$Gene.Symbol,EMT_genes$x)#,genes.pathway.h)

#####数据矩阵
  cancer_cc1mat<-GetAssayData(object =cancer_scrna_cc1, slot = "data")
  cancer_cc2mat<-GetAssayData(object =cancer_scrna_cc2, slot = "data")

  dehallgene<-intersect(rownames(cancer_cc1mat),dehallgene)
  #do_dehallgene<-intersect(rownames(cancer_cc1mat),do_dehallgene)

  scrna_cor1<- cancer_cc1mat[dehallgene,]
  #scrna_cor_do1<- cancer_cc1mat[do_dehallgene,]
  scrna_cor2<- cancer_cc2mat[dehallgene,]
  #scrna_cor_do2<- cancer_cc2mat[do_dehallgene,]

##相关性
##up
cor_p<- function(genelist=up_dehallgene,mat,Pseudotime,cor_f=0.1,p_f=0.001,type="up"){

  ps1<-unlist(lapply(genelist,function(gene){
  p<-cor.test(Pseudotime,mat[gene,])$p.value
  return(p)
  	}))
  cors1<-unlist(lapply(genelist,function(gene){
  c<- cor(Pseudotime,mat[gene,])
  return(c)
  	}))
  df<-data.frame(gene=genelist,cor=cor_f,p=p_f)
  rownames(df)<-genelist
  if(type=="up"){
   genes<-as.vector(na.omit(genelist[cors1>cor_f & ps1<p_f]))
    }else{
   genes<-as.vector(na.omit(genelist[cors1<cor_f & ps1<p_f]))
    }
  
  return(df[genes,])
}

up_cc1_genes<-cor_p(genelist=dehallgene,mat=scrna_cor1,Pseudotime=cancer_scrna_cc1$Pseudotime,cor_f=0.1,p_f=0.01,type="up")
down_cc1_genes<-cor_p(genelist=dehallgene,mat=scrna_cor1,Pseudotime=cancer_scrna_cc1$Pseudotime,cor_f=-0.1,p_f=0.01,type="down")
write.csv(rbind(up_cc1_genes,down_cc1_genes),file="cc1_cor_invariongenes.csv")

up_cc2_genes<-cor_p(genelist=dehallgene,mat=scrna_cor2,Pseudotime=cancer_scrna_cc2$Pseudotime,cor_f=0.1,p_f=0.01,type="up")
down_cc2_genes<-cor_p(genelist=dehallgene,mat=scrna_cor2,Pseudotime=cancer_scrna_cc2$Pseudotime,cor_f=-0.1,p_f=0.01,type="down")
write.csv(rbind(up_cc2_genes,down_cc2_genes),file="cc2_cor_invariongenes.csv")

stage_upgenes<-intersect(up_cc1_genes$gene,up_cc2_genes$gene)
stage_dogenes<-intersect(down_cc1_genes$gen,down_cc2_genes$gen)

write.csv(stage_upgenes,file="stage_upgenes.csv")
write.csv(stage_dogenes,file="stage_downgenes.csv")

  save(stage_dogenes,stage_upgenes,file="cor_p.RData")

load("cor_p.RData")

#####观察基因和hallmark基因的重合数量
library(ggplot2) 
library(cowplot) 
library(ggthemes) 

a<-list(down_genes=as.vector(stage_dogenes),up_genes=as.vector(stage_upgenes))

p_list<-lapply(names(a),function(x){
  stage_genes<-a[[x]]
  stage_pathway<-lapply(genes.by.pathway.h,function(x){
   a<-intersect(stage_genes,x)
   if(length(a)>=1){return(a)}
  })
  stage_pathway<-stage_pathway[sapply(stage_pathway,function(x) !is.null(x))]

num_df<- data.frame(num=unlist(lapply(stage_pathway,length)),
  pathway=names(stage_pathway),
  genes=unlist(lapply(stage_pathway,function(x) paste0(x,collapse=",")))
  )

#num_df<-num_df[num_df$num!=0,]
######### 
p1<-ggplot(num_df, aes(x=reorder(pathway, X = num), y=num))+   
geom_bar(position = "dodge",stat = "identity",fill="#22577E")+theme_classic()+   
labs(x = "pathway", y = "number of gene", title = x)+   
  geom_text(aes(label=genes),stat = "identity",vjust=0.5,hjust=1, 
            color="white", size=3)+
#scale_x_discrete(limits=factor(gg_gp$Description)) 
coord_flip()+theme_minimal()

  })

pdf("down_genes_number_pathway_stage.pdf",width=6,height=4)
print(p_list[[1]])
dev.off()
pdf("up_genes_number_pathway_stage.pdf",height=7,width=7)
print(p_list[[2]])
dev.off()

#######################
#2.绘制不同基因在不同肿瘤细胞中的点图
###################
#Idents(cancer_scrna)<- as.vector(cancer_scrna$seurat_clusters_tumor)
#save(cancer_scrna,file="/share/pub/dengcy/Singlecell/CC_space/3.Epi_cells_stage/cancer_scrna.RData")

#pdf("StackedVlnPlot_downgenes.pdf",height=20)
#StackedVlnPlot(cancer_scrna, features=as.vector(stage_dogenes),cols =color_scanpy_viridis28[c(1,4,5,7,20,11,13,15)])
#dev.off()

#pdf("StackedVlnPlot_upgenes.pdf")
#StackedVlnPlot(cancer_scrna, features=as.vector(stage_upgenes),cols =color_scanpy_viridis28[c(1,4,5,7,20,11,13,15)])
#dev.off()

pdf("DotPlot_cc1_stage_dogenes.pdf",width=8,height=3)
DotPlot(cancer_scrna_cc1, features =as.vector(stage_dogenes), cols = c("lightgrey", "#0F2C67"))+ RotatedAxis()
dev.off()

pdf("DotPlot_cc1_stage_upgenes.pdf",width=12,height=3)
DotPlot(cancer_scrna_cc1, features =as.vector(stage_upgenes), cols = c("lightgrey", "#CD1818"))+ RotatedAxis()
dev.off()

Idents(cancer_scrna_cc2)<- factor(cancer_scrna_cc2$annotation,levels=c("CT5","CT4","CT6"))

pdf("DotPlot_cc2_stage_dogenes.pdf",width=8,height=3)
DotPlot(cancer_scrna_cc2, features =as.vector(stage_dogenes), cols = c("lightgrey", "#0F2C67"))+ RotatedAxis()
dev.off()

pdf("DotPlot_cc2_stage_upgenes.pdf",width=12,height=3)
DotPlot(cancer_scrna_cc2, features =as.vector(stage_upgenes), cols = c("lightgrey", "#CD1818"))+ RotatedAxis()
dev.off()
#############################
#3.绘制侵袭得分的伪时间相关性
#######################
#Invasion
cancer_scrna <- AddModuleScore(cancer_scrna,list(stage_upgenes,stage_dogenes),name=c("Invasion_up_score","Invasion_down_score"))
cancer_scrna$Invasion_up_score<-cancer_scrna$Invasion_up_score1
cancer_scrna$Invasion_down_score<-cancer_scrna$Invasion_down_score2
all_fortify_can <- fortify.Seurat(cancer_scrna)

cancer_scrna_cc1$Invasion_up_score<- all_fortify_can[colnames(cancer_scrna_cc1),]$Invasion_up_score
cancer_scrna_cc1$Invasion_down_score<- all_fortify_can[colnames(cancer_scrna_cc1),]$Invasion_down_score
cancer_scrna_cc2$Invasion_up_score<- all_fortify_can[colnames(cancer_scrna_cc2),]$Invasion_up_score
cancer_scrna_cc2$Invasion_down_score<- all_fortify_can[colnames(cancer_scrna_cc2),]$Invasion_down_score

save(cancer_scrna_cc1,file="/share/pub/dengcy/Singlecell/CC_space/3.Epi_cells_stage/cancer_scrna_cc1.RData")
save(cancer_scrna_cc2,file="/share/pub/dengcy/Singlecell/CC_space/3.Epi_cells_stage/cancer_scrna_cc2.RData")

###1.相关性绘图
p_cor1 <- ggplot(cancer_scrna_cc1@meta.data,aes(x =Pseudotime ,
                y =Invasion_up_score ))+
          geom_point(color="#676FA3",size=0.5)+
          geom_smooth(se = F,
                      color = "#8E0505")+
          labs(x = "Pseudotime",y="Invasion_up_score",title="CC1")+
          theme_classic()+
          stat_cor(data =cancer_scrna_cc1@meta.data,
           method = "pearson")

p_cor2 <- ggplot(cancer_scrna_cc1@meta.data,aes(x = Pseudotime,
                y =Invasion_down_score ))+
          geom_point(color="#676FA3",size=0.5)+
          geom_smooth(se = F,
                     color = "#8E0505")+
          labs(x = "Pseudotime",y="Invasion_down_score",title="CC1")+
          theme_classic()+
          stat_cor(data = cancer_scrna_cc1,
                   method = "pearson")
p_cor3 <- ggplot(cancer_scrna_cc2@meta.data,aes(x =Pseudotime ,
                y =Invasion_up_score ))+
          geom_point(color="#676FA3",size=0.5)+
          geom_smooth(se = F,
                      color = "#8E0505")+
          labs(x = "Pseudotime",y="Invasion_up_score",title="CC2")+
          theme_classic()+
          stat_cor(data =cancer_scrna_cc2@meta.data,
           method = "pearson")

p_cor4 <- ggplot(cancer_scrna_cc2@meta.data,aes(x = Pseudotime,
                y =Invasion_down_score ))+
          geom_point(color="#676FA3",size=0.5)+
          geom_smooth(se = F,
                     color = "#8E0505")+
          labs(x = "Pseudotime",y="Invasion_down_score",title="CC2")+
          theme_classic()+
          stat_cor(data = cancer_scrna_cc2,
                   method = "pearson")
pdf("Ptime_score.pdf",width=8,height=4)
ggarrange(p_cor1,p_cor2,p_cor3,p_cor4, ncol = 2)
dev.off()
#########为时间在不同组中的比较

pdf("VlnPlot_cc1_ptimeincluster.pdf",width=6,height=5)
VlnPlot(object = cancer_scrna_cc1, features =c("Pseudotime"),cols =color_scanpy_viridis28[c(1,4)])+
ggtitle("CC1 Pseudotime")
dev.off()


pdf("VlnPlot_cc2_ptimeincluster.pdf",width=6,height=5)
VlnPlot(object = cancer_scrna_cc2, features =c("Pseudotime"),cols =color_scanpy_viridis28[c(20,7,11)])+
ggtitle("CC2 Pseudotime")
dev.off()

pdf("VlnPlot_cc1_Invasion_up_score.pdf",width=6,height=5)
VlnPlot(object = cancer_scrna_cc1, features =c("Invasion_up_score"),cols =color_scanpy_viridis28[c(1,4)])+
   stat_compare_means(method ="wilcox.test",
                       hide.ns=T,label="p.signif",
                       color="black",size=4,
                       symnum.args= list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                         symbols = c("****", "***", "**", "*", "ns")))+

ggtitle("CC1 Invasion_up_score")
dev.off()

pdf("VlnPlot_cc1_Invasion_down_score.pdf",width=6,height=5)
VlnPlot(object = cancer_scrna_cc1, features =c("Invasion_down_score"),cols =color_scanpy_viridis28[c(1,4)])+
   stat_compare_means(method ="wilcox.test",
                       hide.ns=T,label="p.signif",
                       color="black",size=4,
                       symnum.args= list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                         symbols = c("****", "***", "**", "*", "ns")))+

ggtitle("CC1 Invasion_down_score")
dev.off()


pdf("VlnPlot_cc2_Invasion_up_score.pdf",width=6,height=5)
VlnPlot(object = cancer_scrna_cc2, features =c("Invasion_up_score"),cols =color_scanpy_viridis28[c(20,7,11)])+
   stat_compare_means(method ="kruskal.test",
        #comparisons=list(c('CT4','CT5'),c('CT4','CT6'),c('CT5','CT6')),
                       hide.ns=T,label="p.signif",
                       color="black",size=4,
                       symnum.args= list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                         symbols = c("****", "***", "**", "*", "ns")))+

ggtitle("CC2 Invasion_up_score")
dev.off()

pdf("VlnPlot_cc2_Invasion_down_score.pdf",width=6,height=5)
VlnPlot(object = cancer_scrna_cc2, features =c("Invasion_down_score"),cols =color_scanpy_viridis28[c(20,7,11)])+
   stat_compare_means(method ="kruskal.test",
                       hide.ns=T,label="p.signif",
                       color="black",size=4,
                       symnum.args= list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                         symbols = c("****", "***", "**", "*", "ns")))+

ggtitle("CC2 Invasion_down_score")
dev.off()

########2.整体排序分布得分
#stage_df <- cancer_scrna@meta.data$stage_score1#xxx存储的是每个细胞的恶性程度的打分
#stage_df <- as.matrix(stage_df)
#rownames(stage_df)<-rownames(cancer_scrna@meta.data)#行名是细胞，列是每个细胞的打分
#colnames(stage_df)<-"stage_score"
#stage_df<-as.data.frame(stage_df)
all_fortify_can <- fortify.Seurat(cancer_scrna)
save(all_fortify_can,file="all_fortify_can.RData")
if(F){
stage_df1<-all_fortify_can[,c("Invasion_up_score","annotation")]
stage_df2<-all_fortify_can[,c("Invasion_down_score","annotation")]

#stage_df$stage_score <- exp(stage_df$stage_score1)

stage_df1 <- stage_df1[order(stage_df1$Invasion_up_score,decreasing=FALSE),]#打分从小到大排一个序
stage_df2 <- stage_df2[order(stage_df2$Invasion_down_score,decreasing=FALSE),]#打分从小到大排一个序

save(stage_df1,stage_df2,file="stage_df.RData")
stage_df1$index<-1:nrow(stage_df1)
stage_df2$index<-1:nrow(stage_df2)
#stage_df<-as.matrix(stage_df)

####排序点图
pdf("/share/pub/dengcy/Singlecell/CC_space/4.stage_visualize/pot_line_stagescore_up.pdf")
ggplot(data = stage_df1) + 
  geom_point(mapping = aes(x = index, y = Invasion_up_score, color =annotation))+
   #scale_color_manual(values=color_scanpy_viridis28[c(1,4,5,7,20,11,13,15)])+
  theme_classic()
dev.off()

pdf("/share/pub/dengcy/Singlecell/CC_space/4.stage_visualize/pot_line_stagescore_down.pdf")
ggplot(data = stage_df2) + 
  geom_point(mapping = aes(x = index, y = Invasion_down_score, color =annotation))+
  #scale_color_manual(values=color_scanpy_viridis28[c(1,4,5,7,20,11,13,15)])+
  theme_classic()
dev.off()

###########
##计算排序

median_stage1<-unlist(lapply(unique(stage_df1$annotation),function(x){
  a<-median(stage_df1[stage_df1$annotation %in% x,3])
  return(a)
  }))
names(median_stage1)<-unique(stage_df1$annotation)

#Idents(cancer_scrna)<-cancer_scrna$annotation
col_index<-c(1,4,5,7,20,11,13,15)[c(7,8,6,4,1,5,3,2)]

  plots <- VlnPlot(cancer_scrna,features ="Invasion_up_score",sort="increasing")+ggtitle("Invasion_up_score")+
    stat_compare_means(method ="kruskal.test",
                       hide.ns=T,label="p.signif",
                       color="red",size=4,
                       symnum.args= list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                         symbols = c("****", "***", "**", "*", "ns")))+
    scale_fill_manual(name = "T", values =color_scanpy_viridis28[col_index]) +
    theme(legend.position="none",aspect.ratio=1)
pdf("Invasion_up_score_vinplot.pdf")
print(plots)
dev.off()
}
###########stage umap图

load("/share/pub/dengcy/Singlecell/CC_space/3.Epi_cells_stage/Slingshot_embedded.RData")
CC1_fortify_can <- fortify.Seurat(cancer_scrna_cc1)
CC1_stage_plot1 <- ggplot() +
        geom_point(data = CC1_fortify_can,
                   aes(x = UMAP_1, y = UMAP_2,color =Invasion_up_score), size = 0.2, alpha = 1) +
        umap_theme() +
        scale_colour_gradient(low="#ECDBBA",high="#C84B31")+
        theme(aspect.ratio=1) +
        theme(legend.text=element_markdown(size=14),
              legend.title=element_text(size=14)) +
        guides(colour = guide_legend(override.aes = list(size=3)))  +
       geom_path(data=embedded1.1, aes(x=UMAP_1, y=UMAP_2), size=0.3,color="black",alpha=1)+
                     geom_path(data=embedded1.2, aes(x=UMAP_1, y=UMAP_2), size=0.3,color="black",alpha=1)+

        ggtitle("CC1 Invasion_up_score")

        pdf(file="CC1_Invasion_up_score_plot.pdf",width = 5, height = 5)
        print(CC1_stage_plot1)
        dev.off()
CC1_stage_plot2 <- ggplot() +
        geom_point(data = CC1_fortify_can,
                   aes(x = UMAP_1, y = UMAP_2,color =Invasion_down_score), size = 0.2, alpha = 1) +
        umap_theme() +
        scale_colour_gradient(low="#ECDBBA",high="#161853")+
        theme(aspect.ratio=1) +
        theme(legend.text=element_markdown(size=14),
              legend.title=element_text(size=14)) +
        guides(colour = guide_legend(override.aes = list(size=3)))  +
       geom_path(data=embedded1.1, aes(x=UMAP_1, y=UMAP_2), size=0.3,color="black",alpha=1)+
              geom_path(data=embedded1.2, aes(x=UMAP_1, y=UMAP_2), size=0.3,color="black",alpha=1)+

        ggtitle("CC1 Invasion_down_score")

        pdf(file="CC1_Invasion_down_score_plot.pdf",width = 5, height = 5)
        print(CC1_stage_plot2)
        dev.off()

#######################
CC2_fortify_can <- fortify.Seurat(cancer_scrna_cc2)
CC2_stage_plot1 <- ggplot() +
        geom_point(data = CC2_fortify_can,
                   aes(x = UMAP_1, y = UMAP_2,color =Invasion_up_score), size = 0.2, alpha = 1) +
        umap_theme() +
        scale_colour_gradient(low="#ECDBBA",high="#C84B31")+
        theme(aspect.ratio=1) +
        theme(legend.text=element_markdown(size=14),
              legend.title=element_text(size=14)) +
        guides(colour = guide_legend(override.aes = list(size=3)))  +
       geom_path(data=embedded2.1, aes(x=UMAP_1, y=UMAP_2), size=0.3,color="black",alpha=1)+
       geom_path(data=embedded2.2, aes(x=UMAP_1, y=UMAP_2), size=0.3,color="black",alpha=1)+
        ggtitle("CC2 Invasion_up_score")

        pdf(file="CC2_Invasion_up_score_plot.pdf",width = 5, height = 5)
        print(CC2_stage_plot1)
        dev.off()
CC2_stage_plot2 <- ggplot() +
        geom_point(data = CC2_fortify_can,
                   aes(x = UMAP_1, y = UMAP_2,color =Invasion_down_score), size = 0.2, alpha = 1) +
        umap_theme() +
        scale_colour_gradient(low="#ECDBBA",high="#161853")+
        theme(aspect.ratio=1) +
        theme(legend.text=element_markdown(size=14),
              legend.title=element_text(size=14)) +
        guides(colour = guide_legend(override.aes = list(size=3)))  +
       geom_path(data=embedded2.1, aes(x=UMAP_1, y=UMAP_2), size=0.3,color="black",alpha=1)+
              geom_path(data=embedded2.2, aes(x=UMAP_1, y=UMAP_2), size=0.3,color="black",alpha=1)+

        ggtitle("CC2 Invasion_down_score")

        pdf(file="CC2_Invasion_down_score_plot.pdf",width = 5, height = 5)
        print(CC2_stage_plot2)
        dev.off()
