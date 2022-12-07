
tcga_gdc_expr_process<-function(expr,pdata,sdata,output=".RData"){
  expr$Ensembl_ID<-unlist(str_sub(expr$Ensembl_ID,1,15))
  expr<-expr[!duplicated(expr$Ensembl_ID),]
  rownames(expr)<-expr$Ensembl_ID
  expr<-expr[,-1]
  #colnames(expr) <- unlist(str_sub(colnames(expr),1,15))
  ###
  a<-data.frame(ensembl_id=rownames(expr))
  g2s=toTable(org.Hs.egSYMBOL)
  g2e=toTable(org.Hs.egENSEMBL)
  b=merge(a,g2e,by="ensembl_id",all.x=T)
  d=merge(b,g2s,by="gene_id",all.x=T)
  expr<-expr[d$ensembl_id,]
  expr$gene<-d$symbol
  expr<-expr[!duplicated(expr$gene),]
  expr<-expr[!is.na(expr$gene),]
  rownames(expr)<-expr$gene
  expr<-expr[,1:(ncol(expr)-1)]
  
  pdata$tumor_stage.diagnoses[pdata$tumor_stage.diagnoses=="not reported"]<-NA
  pdata$tumor_stage.diagnoses[pdata$tumor_stage.diagnoses==""]<-NA
  pdata$tumor_stage.diagnoses[pdata$tumor_stage.diagnoses %in% c("stage iia","stage ii","stage iib","stage iic")]<-'stageII'
  pdata$tumor_stage.diagnoses[pdata$tumor_stage.diagnoses %in% c("stage i","stage ia")]<-'stageI'
  pdata$tumor_stage.diagnoses[pdata$tumor_stage.diagnoses %in% c("stage iii","stage iiia","stage iiib","stage iiic")]<-'stageIII'
  pdata$tumor_stage.diagnoses[pdata$tumor_stage.diagnoses %in% c("stage iv","stage iva","stage ivb")]<-'stageIV'
  rownames(pdata)<-pdata$submitter_id.samples
  rownames(sdata)<-sdata$sample
  
  rownames(sdata)<-str_replace_all(sdata$sample,"-",".")
  rownames(sdata)<-unlist(str_sub(rownames(sdata),1,nchar(colnames(expr)[1])))
  
  rownames(pdata)<-str_replace_all(pdata$submitter_id.samples,"-",".")
  rownames(pdata)<-unlist(str_sub(rownames(pdata),1,nchar(colnames(expr)[1])))
  
  ss<-intersect(intersect(colnames(expr),rownames(sdata)),rownames(pdata))
  expr<-expr[,ss]
  pdata<-pdata[ss,]
  sdata<-sdata[ss,]
  pdata<-cbind(pdata,sdata)
  pdata<-pdata[pdata$sample_type.samples=="Primary Tumor",]
  expr<-expr[,rownames(pdata)]
  save(expr,pdata,file = output)
}

library(readr)
library(data.table)
library(stringr)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)
library(enrichplot)
f_data<-list(
setwd='D:/OneDrive/colorectal_sc/CC_space/7.TCGA-READ',
tcga_expr='TCGA-READ.htseq_fpkm.tsv',
tcga_expr2='TCGA-COAD.htseq_fpkm.tsv',

tcga_p='TCGA-READ.GDC_phenotype.tsv',
tcga_p2='TCGA-COAD.GDC_phenotype.tsv',

tcga_s='TCGA-READ.survival.tsv',
tcga_s2='TCGA-COAD.survival.tsv',

emt_gene="D:/OneDrive/colorectal_sc/CC_space/4.stage_visualize/emt_gene.RData")

invasion_emtgenes<-c('CXCL1','CXCL8','SRF','MIF','IL10','IL6','TGFB1','IL4','TNF','IFNB1','HMGB1','TNFSF10','APOE','IL12B','IFNG','NAMPT','EDN1','MMP9','CXCL9','CXCL3','TRH','CXCL11','TGFB1','MIF','IL4','TNF','TNFSF10','APOE','IFNG','EDN1','NAMPT','MMP9','HBEGF','MIF','IL10','TGFB1','IL6','TNF','IL24','TNFSF10','HMGB1','IFNB1','APOE','IL1B','NAMPT','MMP9','IL1RN','IL1A','MIF','IL10','TGFB1','IL6','IL24','TNF','IFNB1','TNFSF10','HMGB1','APOE','IL1B','NAMPT','EDN1','MMP9','IL1RN','IL1A','ALCAM','MIF','TGFB1','HMGB1','TNFSF10','APOE','NAMPT','EDN1','TGFB1','MIF','TNFSF10','APOE','EDN1','NAMPT','TNFSF13B','ICAM1','MIF','TGFB1','HMGB1','IL1B','NAMPT','IL1RN','IL1A','TNFSF13B','ICAM1','MIF','TGFB1','HMGB1','IL1B','NAMPT','IL1RN','IL1A','ANXA1','PMCH','XCL1','MIF','IL10','TGFB1','TNF','HMGB1','APOE','IFNG','NAMPT','ANXA1','PMCH','XCL1','TGFB1','MIF','TNF','APOE','IFNG','NAMPT','ANXA1','XCL1','MIF','TGFB1','TNF','TNFSF10','HMGB1','IL1B','IFNG','NAMPT','IL17A','IL13','MIF','IL10','TGFB1','IL6','IL24','TNF','IFNB1','TNFSF10','HMGB1','APOE','IL1B','NAMPT','EDN1','MMP9','IL1RN','IL1A','TNFSF13B','ICAM1','MIF','TGFB1','HMGB1','IL1B','NAMPT','IL1RN','IL1A','ANXA1','ANXA1','PMCH','XCL1','MIF','TGFB1','TNF','TNFSF10','HMGB1','APOE','IFNG','NAMPT','IL17A','ANXA1','PMCH','XCL1','MIF','TGFB1','TNF','TNFSF10','HMGB1','APOE','IFNG','NAMPT','IL17A','ANXA1','MIF','TGFB1','TNF','HMGB1','TNFSF10','IFNG','NAMPT','IL13','ANXA1','XCL1','TGFB1','MIF','TNF','TNFSF10','IFNG','NAMPT','IL13','ANXA1','XCL1','MIF','TGFB1','TNF','TNFSF10','HMGB1','IL1B','IFNG','NAMPT','IL17A','IL13','ANXA1','XCL1','MIF','TGFB1','TNF','TNFSF10','HMGB1','IL1B','IFNG','NAMPT','IL17A','IL13')

#############
#读取数据
setwd(f_data$setwd)
expr<-read.table(f_data$tcga_expr,header = T,stringsAsFactors = T)
pdata<-read.delim(f_data$tcga_p,header = T,sep = "\t",)
sdata<-read.delim(f_data$tcga_s,header = T,sep = "\t",)

tcga_gdc_expr_process(expr,pdata,sdata,output="TCGA_READ_data.RData")

pdata2<-read.delim(f_data$tcga_p2,header = T,sep = "\t",)
expr2<-read.table(f_data$tcga_expr2,header = T,stringsAsFactors = T)
sdata2<-read.delim(f_data$tcga_s2,header = T,sep = "\t",)
tcga_gdc_expr_process(expr=expr2,pdata=pdata2,sdata=sdata2,output="TCGA_COAD_data.RData")
expr2<-expr
pdata2<-pdata

expr<-cbind(expr1,expr2)
aa<-c("age_at_initial_pathologic_diagnosis","gender.demographic","tissue_or_organ_of_origin.diagnoses","prior_treatment.diagnoses","disease_code","OS","OS.time","tumor_stage.diagnoses","days_to_last_follow_up.diagnoses")
pdata<-rbind(pdata1[,aa],pdata2[,aa])
##
save(expr,pdata,file = "TCGA_COADREAD_data.RData")

load("TCGA_COADREAD_data.RData")
#######
#tcga数据验证emt得分的预后功能
load(f_data$emt_gene)
#构造得分
library(survival)
#library(dplyr)
library(survminer)
library(ggplot2, quietly = TRUE)
library(ggpubr)
#library(GSVA)
library(ggalluvial, quietly = TRUE)

# ct<-ct[ct$pvalue<0.01,]
# goodGenes<-rownames(ct[ct$hr<1,])
# badGenes<-rownames(ct[ct$hr>1,])

##cox回归

ct_os2 <- list()
for (Gene in intersect(c(stage_dogenes,stage_upgenes),rownames(expr))) {
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
#ct_os2$Necroptosis_signature<-rownames(ct_os2)
stage_genes_os<-ct_os2
write.csv(stage_genes_os,file = "stage_genes_os.csv")
ct_os2<-ct_os2[ct_os2$pvalue<0.05,]

for(i in c(stage_dogenes,stage_upgenes)){
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


invasion_emtgenes<-unique(invasion_emtgenes)
invasion_emtgenes<-intersect(rownames(expr),invasion_emtgenes)
data1<-expr[c(invasion_emtgenes),]
#data2<-expr[c(stage_dogenes),]

com1 <- prcomp(t(data1), center = TRUE,scale. = TRUE)
#com2 <- prcomp(t(data2), center = TRUE,scale. = TRUE)

pdata$emt_Pseudotime_score<- com1$x[,1]

pdata$emt_Pseudotime_score_level<-rep("High",nrow(pdata))
pdata$emt_Pseudotime_score_level[which(pdata$emt_Pseudotime_score <= 
                                    median(pdata$emt_Pseudotime_score))]<-"Low"
survNecr <- ggsurvplot(
  survfit(Surv(OS.time, OS) ~ emt_Pseudotime_score_level, data = pdata),
  font.legend=14,xlab="OS Time (days)",legend.title="emt_Pseudotime_score_level",pval=TRUE,
  font.y=c(20),font.x=c(20),font.tickslab =c(18),
  pval.size=7,palette = c("#FF7F00", "#377EB8"),
  risk.table = T,
  ggtheme = theme_survminer(),
  title="TCGA Data OS survival analysis"
)
pdf("TCGAscore_level_OS.pdf",width = 7, onefile=T)
print(survNecr)
dev.off()

#######
#tcga数据验证
#http://timer.comp-genomics.org/
timer_data<-read.csv("D:/OneDrive/PTjob/肺腺癌铜死亡/F220401001-肺腺癌铜死亡模式解析/runtime/1.mian_process/infiltration_estimation_for_tcga.csv")

pdata$CXCL1 <- unlist(expr['CXCL1',])
pdata$CXCL8 <- unlist(expr['CXCL8',])
pdata$SRF <- unlist(expr['SRF',])

pdata$tumor_stage.diagnoses<-factor(pdata$tumor_stage.diagnoses,levels=c("stageI","stageII","stageIII","stageIV"))

pdata<-pdata[!is.na(pdata$tumor_stage.diagnoses)]

p1<-ggboxplot(pdata, x = "tumor_stage.diagnoses", y = "CXCL1",
                 fill = 'tumor_stage.diagnoses', 
                 palette = c('#f6eec9','#fed39f','#fe8761','#af460f'))+
  stat_compare_means(method = "anova",
                     label = "p.signif", 
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                                            symbols = c("****", "***", "**", "*", "ns")))+
  labs(x='Stage',y='Expression')

print(p1)
ggsave('Figure_tcga_CXCL1.pdf',width=4,height=4)

p2<-ggboxplot(pdata, x = "tumor_stage.diagnoses", y = "CXCL8",
              fill = 'tumor_stage.diagnoses', 
              palette = c('#f6eec9','#fed39f','#fe8761','#af460f'))+
  stat_compare_means(method = "anova",
                     label = "p.signif", 
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                        symbols = c("****", "***", "**", "*", "ns")))+
  labs(x='Stage',y='Expression')

print(p2)
ggsave('Figure_tcga_CXCL8.pdf',width=4,height=4)

p3<-ggboxplot(pdata, x = "tumor_stage.diagnoses", y = "SRF",
              fill = 'tumor_stage.diagnoses', 
              palette = c('#f6eec9','#fed39f','#fe8761','#af460f'))+
  stat_compare_means(method = "anova",
                     label = "p.signif", 
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                        symbols = c("****", "***", "**", "*", "ns")))+
  labs(x='Stage',y='Expression')

print(p3)
#ggsave('Figure_tcga_SRF.pdf',width=4,height=4)
library(ggpubr) 
pdf('Figure_tcga_stage.pdf',width=9,height=5)
ggarrange(p1, p2, p3,nrow =1)
dev.off()

##########################
#emtg高低得分组的免疫细胞浸润差异
library(stringr) 

timer_data$cell_type<-str_replace_all(timer_data$cell_type,"-",".")
rownames(timer_data)<-timer_data$cell_type
ss<-intersect(rownames(pdata),rownames(timer_data))
pdata<-cbind(pdata[ss,],timer_data[ss,])
xcelldata<-pdata[,90:128]

xcelldata$emt_Pseudotime_score_level <- pdata$emt_Pseudotime_score_level
xcelldata2<-xcelldata[,-c(1,37,38,39,22,23,13,14,16,19,20)]
xcelldatagg<-reshape2::melt(xcelldata2,id.vars='emt_Pseudotime_score_level')


xcelldatagg$variable<-unlist(lapply(as.vector(xcelldatagg$variable),function(x){
  return(strsplit(x, split = "_",fixed=T)[[1]][1])
}))

p <- ggplot(xcelldatagg, aes(x = variable, y =value,fill=emt_Pseudotime_score_level)) + 
  geom_boxplot(outlier.size = 0.3,alpha=0.8)+theme_classic()+
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),legend.position = "top")+
  scale_fill_manual(values=c('#ff9a76','#679b9b'))+
  stat_compare_means(method ="wilcox.test",
                             hide.ns=T,label="p.signif",
                             color="black",size=4,
                             symnum.args= list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                               symbols = c("****", "***", "**", "*", "ns")))


#ggsave('Figure_tcga_SRF.pdf',width=4,height=4)
library(ggpubr) 
pdf('Figure_tcga_scorelevel.pdf',width=8,height=5)
p
dev.off()

