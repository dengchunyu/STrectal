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

library(survival)
#library(dplyr)
library(survminer)
library(ggplot2, quietly = TRUE)
library(ggpubr)

library(monocle)
library(Seurat)
source("/share/pub/dengcy/Singlecell/CC_space/src/Functions.r")
#source("/share/pub/dengcy/Singlecell/CC_space/src/library.r")
setwd("/share/pub/dengcy/Singlecell/CC_space/4.stage_visualize")
#load("/share/pub/dengcy/Singlecell/CC_space/3.Epi_cells_stage/cancer_scrna.RData")
load("/share/pub/dengcy/Singlecell/CC_space/3.Epi_cells_stage/emt_pseudotim_gene.RData")
load("/share/pub/dengcy/Singlecell/CC_space/3.Epi_cells_stage/monocle/cds_cancer_scrna_stage2.RData")
load("/share/pub/dengcy/Singlecell/CC_space/3.Epi_cells_stage/monocle/cds_cancer_scrna_stage1.RData")
load("/share/pub/dengcy/Singlecell/CC_space/7.TCGA-READ/TCGA_COADREAD_data.RData")
load("/share/pub/dengcy/Singlecell/CC_space/4.stage_visualize/genes.by.pathway.h.RData")
##############
##1.get the emt stage for different sample
##################
#预后
ct_os2 <- list()
for (Gene in intersect(emt_pseudotim_gene,rownames(expr))) {
  pdata$gene_level <- unlist(expr[Gene,])
  fit <- coxph(Surv(OS.time, OS) ~ gene_level, data = pdata)
  s <- summary(fit)
  ct_os2[[Gene]] = c(
    s$conf.int[1],
    s$conf.int[3],
    s$conf.int[4],
    s$logtest[3]
  )
}

ct_os2 <- t(as.data.frame(ct_os2, check.names = FALSE))
colnames(ct_os2) <- c("hr", "lower.95", "upper.95", "pvalue")
ct_os2 <- as.data.frame(ct_os2)
ct_os2<-ct_os2[ct_os2$pvalue<0.05,]
save(ct_os2,file="stage_genes_os.RData")
stage_genes_os<-ct_os2
write.csv(ct_os2,file="stage_genes_os.csv")
emt_pseudotim_gene<-rownames(stage_genes_os)




##相关性
##up
#pData(cds1)
cor_p<- function(genelist=up_dehallgene,mat,Pseudotime,cor_f=0.1,p_f=0.001,type="up"){

  ps1<-unlist(lapply(genelist,function(gene){
  p<-cor.test(Pseudotime,mat[gene,])$p.value
  return(p)
  	}))
  cors1<-unlist(lapply(genelist,function(gene){
  c<- cor(Pseudotime,mat[gene,])
  return(c)
  	}))
  df<-data.frame(gene=genelist,cor=cors1,p=ps1)
  rownames(df)<-genelist
 # if(type=="up"){
 #  genes<-as.vector(na.omit(genelist[cors1>cor_f & ps1<p_f]))
  #  }else{
  # genes<-as.vector(na.omit(genelist[cors1<cor_f & ps1<p_f]))
  #  }
  
  return(df)
}
cds1pdata<-pData(cds1)
cc1_genes<-cor_p(genelist=emt_pseudotim_gene,mat=exprs(cds1),Pseudotime=cds1pdata$Pseudotime)
up_cc1_genes<-cc1_genes[cc1_genes$cor > 0.1 & cc1_genes$p<0.05,]
down_cc1_genes<-cc1_genes[cc1_genes$cor < -0.1 & cc1_genes$p<0.05,]
write.csv(cc1_genes,file="cc1_cor_invariongenes.csv")

cds2data<-pData(cds2)
cc2_genes<-cor_p(genelist=emt_pseudotim_gene,mat=exprs(cds2),Pseudotime=cds2data$Pseudotime)
write.csv(cc2_genes,file="cc2_cor_invariongenes.csv")
cc2_genes<-read.csv("cc2_cor_invariongenes.csv")
up_cc2_genes<-cc2_genes[cc2_genes$cor > 0.1 & cc1_genes$p<0.05,]
down_cc2_genes<-cc2_genes[cc2_genes$cor < -0.1 & cc1_genes$p<0.05,]

stage_upgenes<-intersect(up_cc1_genes$gene,up_cc2_genes$gene)
stage_dogenes<-intersect(down_cc1_genes$gene,down_cc2_genes$gene)

write.csv(stage_upgenes,file="stage_upgenes.csv")
#write.csv(stage_dogenes,file="stage_downgenes.csv")

save(stage_dogenes,stage_upgenes,file="emt_gene.RData")

for(i in stage_upgenes){
  pdata$level<-"High"
  a<-unlist(expr[i,])
  b<-median(a)
  pdata$level[a<b]<-"Low"
  
  pdf(paste0(i,"_read_TCGA_os.pdf"), width =10, height = 8)
  print(ggsurvplot(
    survfit(Surv(OS.time, OS) ~ level, data = pdata),
    font.legend=14,xlab="OS Time (days)",legend.title="cluster",pval=TRUE,
    font.y=c(20),font.x=c(20),font.tickslab =c(18),pval.size=7,
    palette = c("#297F87", "#F6D167"),risk.table = TRUE,
    ggtheme = theme_survminer()
  ))
  dev.off()
  
}

#####hallmark genes
library(ggplot2) 
library(cowplot) 
library(ggthemes) 
load("emt_gene.RData")
a<-list(#down_genes=as.vector(stage_dogenes),
  up_genes=as.vector(stage_upgenes))

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

#pdf("down_genes_number_pathway_stage.pdf",width=6,height=4)
#print(p_list[[1]])
#dev.off()
pdf("up_genes_number_pathway_stage.pdf",height=5,width=7)
print(p_list[[1]])
dev.off()

#######################
#2.绘制不同基因在不同肿瘤细胞中的点图
###################
library(monocle)
load("/share/pub/dengcy/Singlecell/CC_space/3.Epi_cells_stage/cancer_scrna_stage1.RData")
load("/share/pub/dengcy/Singlecell/CC_space/3.Epi_cells_stage/cancer_scrna_stage2.RData")

cancer_scrna_stage1@meta.data<-pData(cds1)
cancer_scrna_stage2@meta.data<-pData(cds2)
save(cancer_scrna_stage1,file="/share/pub/dengcy/Singlecell/CC_space/3.Epi_cells_stage/cancer_scrna_stage1.RData")
save(cancer_scrna_stage2,file="/share/pub/dengcy/Singlecell/CC_space/3.Epi_cells_stage/cancer_scrna_stage2.RData")



Idents(cancer_scrna_stage1)<-cancer_scrna_stage1$State
#pdf("DotPlot_cc1_stage_dogenes.pdf",width=5,height=3)
#DotPlot(cancer_scrna_stage1, features =as.vector(stage_dogenes), cols = c("lightgrey", "#0F2C67"))+ RotatedAxis()
#dev.off()

pdf("DotPlot_cc1_stage_upgenes.pdf",width=5,height=3)
DotPlot(cancer_scrna_stage1, features =as.vector(stage_upgenes), cols = c("lightgrey", "#CD1818"))+ RotatedAxis()
dev.off()

Idents(cancer_scrna_stage2)<- cancer_scrna_stage2$State

#pdf("DotPlot_cc2_stage_dogenes.pdf",width=5,height=3)
#DotPlot(cancer_scrna_stage2, features =as.vector(stage_dogenes), cols = c("lightgrey", "#0F2C67"))+ RotatedAxis()
#dev.off()

pdf("DotPlot_cc2_stage_upgenes.pdf",width=5,height=3)
DotPlot(cancer_scrna_stage2, features =as.vector(stage_upgenes), cols = c("lightgrey", "#CD1818"))+ RotatedAxis()
dev.off()

p1<-plot_genes_in_pseudotime(cds1[stage_upgenes,],color_by="State")+scale_colour_manual(name = "State", values = c("#f38181","#fce38a","#eaffd0","#95e1d3","#64958f") )
p2<-plot_genes_in_pseudotime(cds2[stage_upgenes,],color_by="State")+scale_colour_manual(name = "State", values = c("#f38181","#fce38a","#eaffd0","#95e1d3","#64958f") )
ggsave("Up_gene_in_pseudotime_1.pdf",plot= p1|p2,width=8,height=10)


#############################
#3.绘制侵袭得分的伪时间相关性
#######################
#Invasion
cancer_scrna_stage1 <- AddModuleScore(cancer_scrna_stage1,list(stage_upgenes),name=c("Invasion_score"))
cancer_scrna_stage2 <- AddModuleScore(cancer_scrna_stage2,list(stage_upgenes),name=c("Invasion_score"))

cancer_scrna_stage1$Invasion_score<-cancer_scrna_stage1$Invasion_score1
cancer_scrna_stage1$Invasion_score[cancer_scrna_stage1$Invasion_score>10]<-10
cancer_scrna_stage1$Invasion_score[cancer_scrna_stage1$Invasion_score< -10]<- -10

cancer_scrna_stage2$Invasion_score<-cancer_scrna_stage2$Invasion_score1
cancer_scrna_stage2$Invasion_score[cancer_scrna_stage2$Invasion_score>10]<-10
cancer_scrna_stage2$Invasion_score[cancer_scrna_stage2$Invasion_score< -10]<- -10
#cancer_scrna_stage2$Invasion_down_score<-cancer_scrna_stage2$Invasion_down_score2
pData(cds1)$Invasion_score<-cancer_scrna_stage1$Invasion_score
#pData(cds1)$Invasion_down_score<-cancer_scrna_stage1$Invasion_down_score
pData(cds2)$Invasion_score<-cancer_scrna_stage2$Invasion_score
#pData(cds2)$Invasion_down_score<-cancer_scrna_stage2$Invasion_down_score

exprData1 <- cds1@assayData$exprs
exprData1 <- LogNormalize(exprData1) #
exprData2 <- cds2@assayData$exprs
exprData2 <- LogNormalize(exprData2) 
#对数据normaliza一下，会比较显著，记得引用seurat包

for (i in stage_upgenes) {
cds1$emtgenes <- exprData1[i,] #输入想看的基因，cds是monocle的对象
cds2$emtgenes <- exprData2[i,] #输入想看的基因，cds是monocle的对象
p1<-plot_cell_trajectory(cds1, color_by ="emtgenes",cell_size=1,show_backbone=TRUE)+scale_color_gradient(low = "#EBEBEB",high = "#F5A25D",limits=c(0,2.5))
p2<-plot_cell_trajectory(cds2, color_by ="emtgenes",cell_size=1,show_backbone=TRUE)+scale_color_gradient(low = "#EBEBEB",high = "#F5A25D",limits=c(0,2.5))
ggsave(paste0(i,"_for_cell_trajectory.pdf"),plot= p1|p2,width=14,height=7)
}

p1<-plot_cell_trajectory(cds1,cell_size=1,color_by="Invasion_score")+ scale_color_gradient(low = "#f5efef",high = "#cf7500")

p2<-plot_cell_trajectory(cds2,cell_size=1,color_by="Invasion_score")+ scale_color_gradient(low = "#f5efef",high = "#cf7500")

ggsave("Invasion_score_for_cell_trajectory2.pdf",plot= p1|p2,width=14,height=7)

all_fortify_can1 <- fortify.Seurat(cancer_scrna_stage1)
all_fortify_can2 <- fortify.Seurat(cancer_scrna_stage2)

save(cancer_scrna_stage1,file="/share/pub/dengcy/Singlecell/CC_space/3.Epi_cells_stage/cancer_scrna_stage1.RData")
save(cancer_scrna_stage2,file="/share/pub/dengcy/Singlecell/CC_space/3.Epi_cells_stage/cancer_scrna_stage2.RData")
save(cds1,file="cds_cancer_scrna_stage1.RData")
save(cds2,file="cds_cancer_scrna_stage2.RData")

###1.相关性绘图
library(ggpubr) 
all_fortify_can1 <- fortify.Seurat(cancer_scrna_stage1)
all_fortify_can2 <- fortify.Seurat(cancer_scrna_stage2)


p_cor1 <- ggplot(all_fortify_can1,aes(x =Pseudotime ,
                y =Invasion_score))+
          geom_point(aes(colour=State),size=0.5)+
          geom_smooth(se = F,
                      color = "#8E0505")+
          labs(x = "Pseudotime",y="Invasion_score",title="StageI/II")+
          theme_classic()+
          stat_cor(data =all_fortify_can1,
           method = "pearson")+
          scale_colour_manual(name = "State", values = c("#f38181","#fce38a","#eaffd0","#95e1d3","#64958f") )

p_cor3 <- ggplot(all_fortify_can2,aes(x =Pseudotime ,
                y =Invasion_score ))+
          geom_point(aes(colour=State),size=0.5)+
          geom_smooth(se = F,
                      color = "#8E0505")+
          labs(x = "Pseudotime",y="Invasion_score",title="StageIII")+
          theme_classic()+
          stat_cor(data =all_fortify_can2,
           method = "pearson")+
          scale_colour_manual(name = "State", values = c("#f38181","#fce38a","#eaffd0","#95e1d3","#64958f") )

pdf("Ptime_score.pdf",width=8,height=7)
ggarrange(p_cor1,p_cor3, ncol = 1)
dev.off()


a1<-all_fortify_can1[all_fortify_can1$State==5,]
a2<-all_fortify_can2[all_fortify_can2$State ==3,]
a3<-all_fortify_can1[all_fortify_can1$State ==4,]

p_cor2 <- ggplot(a1,aes(x =Pseudotime ,
                y =Invasion_score))+
          geom_point(colour="#64958f",size=0.5)+
          geom_smooth(se = F,
                      color = "#8E0505")+
          labs(x = "Pseudotime",y="Invasion_score",title="StageI/II")+
          theme_classic()+
          stat_cor(data =a1,
           method = "pearson")
         # scale_colour_manual(name = "State", values = c("#f38181","#fce38a","#eaffd0","#95e1d3","#64958f") )

p_cor4 <- ggplot(a2,aes(x =Pseudotime ,
                y =Invasion_score ))+
          geom_point(colour="#eaffd0",size=0.5)+
          geom_smooth(se = F,
                      color = "#8E0505")+
          labs(x = "Pseudotime",y="Invasion_score",title="StageIII")+
          theme_classic()+
          stat_cor(data =a2,
           method = "pearson")
          #scale_colour_manual(name = "State", values = c("#f38181","#fce38a","#eaffd0","#95e1d3","#64958f") )
p_cor3 <- ggplot(a3,aes(x =Pseudotime ,
                y =Invasion_score))+
          geom_point(colour="#95e1d3",size=0.5)+
          geom_smooth(se = F,
                      color = "#8E0505")+
          labs(x = "Pseudotime",y="Invasion_score",title="StageI/II")+
          theme_classic()+
          stat_cor(data =a3,
                       method = "pearson")
         # scale_colour_manual(name = "State", values =            method = "pearson")
#c("#f38181","#fce38a","#eaffd0","#95e1d3","#64958f") )
pdf("Ptime_score_sub.pdf",width=5,height=12)
ggarrange(p_cor2,p_cor3,p_cor4, ncol = 1)
dev.off()
#########为时间在不同组中的比较

pdf("VlnPlot_StageI_Invasion_score.pdf",width=6,height=5)
VlnPlot(object = cancer_scrna_stage1,group.by ="State", features =c("Invasion_score"),
  cols =c("#f38181","#fce38a","#eaffd0","#95e1d3","#64958f"))+
   stat_compare_means(method ="kruskal.test",
                       hide.ns=T,label="p.signif",
                       color="black",size=4,
                       symnum.args= list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                         symbols = c("****", "***", "**", "*", "ns")))+

ggtitle("StageI/II Invasion_score")
dev.off()

pdf("VlnPlot_StageIII_Invasion_score.pdf",width=6,height=5)
VlnPlot(object = cancer_scrna_stage2,group.by ="State", features =c("Invasion_score"),
  cols =c("#f38181","#fce38a","#eaffd0","#95e1d3","#64958f"))+
   stat_compare_means(method ="kruskal.test",
        #comparisons=list(c('CT4','CT5'),c('CT4','CT6'),c('CT5','CT6')),
                       hide.ns=T,label="p.signif",
                       color="black",size=4,
                       symnum.args= list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                         symbols = c("****", "***", "**", "*", "ns")))+

ggtitle("StageIII Invasion_score")
dev.off()

library(ggpubr)

color_CRC<-c("#7FBCD2","#A5F1E9","#FFEEAF","#7FB77E","#F7F6DC","#B1D7B4","#FFC090",
  "#F3E0B5","#FFAE6D","#FECD70","#E3C770","#94B49F","#CEE5D0","#ECB390","#F2D7D9",
  "#D3CEDF","#B2C8DF","#C4D7E0","#F8F9D7","#76BA99","#ADCF9F",
  "#CED89E","#FFDCAE","#FFE3A9","#FFC3C3","#FF8C8C","#97C4B8","#F1F0C0","#B1BCE6")

pdf("VlnPlot_StageI_clusters_score.pdf",width=6,height=5)
ggplot(all_fortify_can1,aes(x=Cluster,y=Invasion_score,fill=Cluster))+
  geom_boxplot(outlier.size = 0.6,alpha=0.6)+ #箱式图异常值大小调整
    scale_fill_manual(values=color_CRC[1:7])+ #填充颜色调整
    stat_compare_means(method ="kruskal.test",
                       hide.ns=T,label="p.signif",
                       color="black",size=4,
                       symnum.args= list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                         symbols = c("****", "***", "**", "*", "ns")))+
     theme_classic()+
    ggtitle("StageI/II Invasion_score")
dev.off()

pdf("VlnPlot_StageIII_clusters_score.pdf",width=6,height=5)
ggplot(all_fortify_can1,aes(x=Cluster,y=Invasion_score,fill=Cluster))+
  geom_boxplot(outlier.size = 0.6,alpha=0.6)+ #箱式图异常值大小调整
    scale_fill_manual(values=color_CRC[1:7])+ #填充颜色调整
    stat_compare_means(method ="kruskal.test",
                       hide.ns=T,label="p.signif",
                       color="black",size=4,
                       symnum.args= list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                         symbols = c("****", "***", "**", "*", "ns")))+
     theme_classic()+
    ggtitle("StageI/II Invasion_score")
dev.off()


save.image(".RData")
###############################################################删除######################################################################
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

###############CXCL1,2,3差异分析
a1<-cancer_scrna_stage1@meta.data$stage
a2<-cancer_scrna_stage2@meta.data$stage
s1<-cancer_scrna_stage1$State
s2<-cancer_scrna_stage2$Cluster
cancer_scrna_stage1<-NormalizeData(
  cancer_scrna_stage1,
  normalization.method = "LogNormalize"
)
cancer_scrna_stage2<-NormalizeData(
  cancer_scrna_stage2,
  normalization.method = "LogNormalize"
)
b1<-GetAssayData(cancer_scrna_stage1, assay="RNA")
b2<-GetAssayData (cancer_scrna_stage2, assay="RNA")
b1<-t(b1[c('CXCL1','CXCL2','CXCL3','CXCL8','SRF','IFNG','IGF2','SEMA3C'),])
b2<-t(b2[c('CXCL1','CXCL2','CXCL3','CXCL8','SRF','IFNG','IGF2','SEMA3C'),])

df<-as.data.frame(rbind(b1,b2))

df$stage<-c(a1,a2)
df$state<-c(paste0("State",s1),paste0("Cluster",s2))
df$stage[df$stage != "StageIII"]<-"StageI/II"
df<-df[df$state %in% c("State4","State5","Cluster3","Cluster7"),]

pdf("CXCL1_for_Stage_clusters.pdf",width=3.5,height=5)
ggplot(df,aes(x=stage,y=CXCL1,fill=stage))+
  geom_boxplot(outlier.size = 0.4,alpha=0.6)+ #箱式图异常值大小调整
    scale_fill_manual(values=c("#4e89ae","#ed6663"))+ #填充颜色调整
    stat_compare_means(method ="wilcox.test",
                       hide.ns=T,label="p.signif",
                       color="black",size=4,
                       symnum.args= list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                         symbols = c("****", "***", "**", "*", "ns")))+
     theme_classic()+
    ggtitle("CXCL1")
dev.off()

pdf("SRF_for_Stage_clusters.pdf",width=3.5,height=5)
ggplot(df,aes(x=stage,y=SRF,fill=stage))+
  geom_boxplot(outlier.size = 0.6,alpha=0.6)+ #箱式图异常值大小调整
    scale_fill_manual(values=c("#4e89ae","#ed6663"))+ #填充颜色调整
    stat_compare_means(method ="wilcox.test",
                       hide.ns=T,label="p.signif",
                       color="black",size=4,
                       symnum.args= list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                         symbols = c("****", "***", "**", "*", "ns")))+
     theme_classic()+
    ggtitle("SRF")
dev.off()

pdf("CXCL8_for_Stage_clusters.pdf",width=3.5,height=5)
ggplot(df,aes(x=stage,y=CXCL8,fill=stage))+
  geom_boxplot(outlier.size = 0.6,alpha=0.6)+ #箱式图异常值大小调整
    scale_fill_manual(values=c("#4e89ae","#ed6663"))+ #填充颜色调整
    stat_compare_means(method ="wilcox.test",
                       hide.ns=T,label="p.signif",
                       color="black",size=4,
                       symnum.args= list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                         symbols = c("****", "***", "**", "*", "ns")))+
     theme_classic()+
    ggtitle("CXCL8")
dev.off()


pdf("CXCL2_for_Stage_clusters.pdf",width=6,height=5)
ggplot(df,aes(x=stage,y=CXCL2,fill=stage))+
  geom_boxplot(outlier.size = 0.6,alpha=0.6)+ #箱式图异常值大小调整
    scale_fill_manual(values=c("#4e89ae","#ed6663"))+ #填充颜色调整
    stat_compare_means(method ="wilcox.test",
                       hide.ns=T,label="p.signif",
                       color="black",size=4,
                       symnum.args= list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                         symbols = c("****", "***", "**", "*", "ns")))+
     theme_classic()+
    ggtitle("CXCL2")
dev.off()

pdf("CXCL3_for_Stage_clusters.pdf",width=6,height=5)
ggplot(df,aes(x=stage,y=CXCL3,fill=stage))+
  geom_boxplot(outlier.size = 0.6,alpha=0.6)+ #箱式图异常值大小调整
    scale_fill_manual(values=c("#4e89ae","#ed6663"))+ #填充颜色调整
    stat_compare_means(method ="wilcox.test",
                       hide.ns=T,label="p.signif",
                       color="black",size=4,
                       symnum.args= list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                         symbols = c("****", "***", "**", "*", "ns")))+
     theme_classic()+
    ggtitle("CXCL3")
dev.off()

pdf("IFNG_for_Stage_clusters.pdf",width=6,height=5)
ggplot(df,aes(x=stage,y=IFNG,fill=stage))+
  geom_boxplot(outlier.size = 0.6,alpha=0.6)+ #箱式图异常值大小调整
    scale_fill_manual(values=c("#4e89ae","#ed6663"))+ #填充颜色调整
    stat_compare_means(method ="wilcox.test",
                       hide.ns=T,label="p.signif",
                       color="black",size=4,
                       symnum.args= list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                         symbols = c("****", "***", "**", "*", "ns")))+
     theme_classic()+
    ggtitle("IFNG")
dev.off()

###############通路差异分析
x <- readLines("/share/pub/dengcy/Singlecell/CC_space/4.stage_visualize/h.all.v7.5.1.symbols.gmt")
res <- strsplit(x, "\t")
names(res) <- vapply(res, function(y) y[1], character(1))
genes.by.hallmark.pathway <- lapply(res, "[", -c(1:2))

cancer_scrna_stage<-merge(cancer_scrna_stage1,cancer_scrna_stage2)
cancer_scrna_stage<-NormalizeData(
  cancer_scrna_stage,
  normalization.method = "LogNormalize"
)

cancer_scrna_stage <- AddModuleScore(cancer_scrna_stage,genes.by.hallmark.pathway)
s1<-cancer_scrna_stage1$State
s2<-cancer_scrna_stage2$Cluster
df2<-cancer_scrna_stage@meta.data
df2$states<-c(paste0("State",s1),paste0("Cluster",s2))
df2$stage[df2$stage != "StageIII"]<-"StageI/II"
df2<-df2[df2$states %in% c("State4","State5","Cluster3","Cluster7"),]
df2<-df2[,c(25,45:94)]
colnames(df2)[2:51]<-names(genes.by.hallmark.pathway)

df_gg<-reshape2::melt(df2,id.vars='stage')

p <- ggplot(df_gg, aes(x = variable, y =value,fill=stage)) + 
  geom_boxplot(outlier.size = 0.3,alpha=0.8)+theme_classic()+
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),legend.position = "top")+
  scale_fill_manual(values=c("#4e89ae","#ed6663"))+
  stat_compare_means(method ="wilcox.test",
                             hide.ns=T,label="p.signif",
                             color="black",size=4,
                             symnum.args= list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                               symbols = c("****", "***", "**", "*", "ns")))


#ggsave('Figure_tcga_SRF.pdf',width=4,height=4)
library(ggpubr) 
pdf('Figure_stage_hallmark.pdf',width=15,height=5)
p
dev.off()

pa<-c("stage","HALLMARK_TNFA_SIGNALING_VIA_NFKB","HALLMARK_HYPOXIA","HALLMARK_INTERFERON_ALPHA_RESPONSE","HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")
df3<-df2[,pa]
df_gg2<-reshape2::melt(df3,id.vars='stage')

p <- ggplot(df_gg2, aes(x = variable, y =value,fill=stage)) + 
  geom_boxplot(outlier.size = 0.3,alpha=0.8)+theme_classic()+
  theme(axis.text.x = element_text(angle = 90,vjust = 1,hjust = 1),legend.position = "top")+
  scale_fill_manual(values=c("#4e89ae","#ed6663"))+
  stat_compare_means(method ="wilcox.test",
                             hide.ns=T,label="p.signif",
                             color="black",size=4,
                             symnum.args= list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                               symbols = c("****", "***", "**", "*", "ns")))


#ggsave('Figure_tcga_SRF.pdf',width=4,height=4)
library(ggpubr) 
pdf('Figure_stage_hallmark_select.pdf',width=4,height=7)
p
dev.off()