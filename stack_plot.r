################
#绘制stack plot
library(plyr)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(ggtext)
color_CRC<-c("#7FBCD2","#A5F1E9","#FFEEAF","#7FB77E","#F7F6DC","#B1D7B4","#FFC090",
             "#F3E0B5","#FFAE6D","#FECD70","#E3C770","#94B49F","#CEE5D0","#ECB390","#F2D7D9",
             "#D3CEDF","#B2C8DF","#C4D7E0","#F8F9D7","#76BA99","#ADCF9F",
             "#CED89E","#FFDCAE","#FFE3A9","#FFC3C3","#FF8C8C","#97C4B8","#F1F0C0","#B1BCE6")


setwd("D:/OneDrive/colorectal_sc/CC_space/3.2.Immune_cells_stage/")
load("metadata.RData")
options<-c("annotation","stage")
label.size<-4
colors<-rev(color_CRC[14:21])
metadata$stage<-"StageIII"
metadata$stage[metadata$samples=="CC1"]<-"StageI/II"
metadata$stage[metadata$orig.ident %in% c("KUL01-B","KUL01-T","KUL28-B","KUL28-T","KUL30-B","KUL30-T",
                                                      "SMC01-T","SMC05-T","SMC09-T","SMC10-T","SMC15-T","SMC18-T","scrEXT001","scrEXT002","scrEXT021",
                                                      "scrEXT022","scrEXT027","scrEXT028")]<-"StageI/II"
metadata$stage[metadata$orig.ident %in% c("KUL31-B","KUL31-T","SMC07-T","SMC24-T","scrEXT018","scrEXT019")]<-"StageI/II"
metadata$orig.ident[metadata$samples=="CC1"]<-"CC1"
metadata$orig.ident[metadata$samples=="CC2"]<-"CC2"

save(metadata,file="metadata.RData")

ggplot_df<-metadata[,options]
colnames(ggplot_df)<-c("index1","index2")

a<-tapply(ggplot_df$index1,factor(ggplot_df$index2),function(index2){table(index2) })
dfll <- do.call(rbind,lapply(a, data.frame))
dfll$index1 <- rep(names(a),sapply(a,length))
# Calculate the percentages
dfll = plyr::ddply(dfll, .(index1), transform, percent = Freq/sum(Freq) * 100)

# Format the labels and calculate their positions
dfll = plyr::ddply(dfll, .(index1), transform, pos = cumsum(percent)/100)
dfll$label = paste0(sprintf("%.0f", dfll$percent), "%")
dfll$index2<-factor(dfll$index2,levels = rev(as.vector(dfll$index2[1:8])))
p<-ggplot(dfll, aes(x =index1, y = percent, fill =index2)) +
  geom_bar(stat = "identity", position= 'fill',size=0.3) +
  geom_text(aes(y = pos, label = label), size = label.size)+
  theme_classic()+
  scale_fill_manual(values=colors)+
  #
  labs(x = options[1], title = options[2])

  p<-p+theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1))

  p<-p+coord_flip()
  
pdf("Immuneproportion_in_stage.pdf",width = 8,height = 3)
print(p)
dev.off()

##################
#样本的差异分析
