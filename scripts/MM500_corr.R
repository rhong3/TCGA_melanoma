# MM500 New Prot-RNA correlations
library(readxl)
ProteinAbundanceRank <- read_excel("Data/ProteinAbundanceRank.xlsx")
Protein = ProteinAbundanceRank[, 3:4]
s <- strsplit(Protein$Gene, split = ";")
Protein = data.frame(PROT_Median = rep(Protein$MedianIntensity, sapply(s, length)), `Gene name` = unlist(s))
RNAseq <- read.csv("Data/log2_RNAseq.csv")

# CORUM protein complex
library(MoonFinder)
data("corumModuleList")
corum_gene=c()
for (corum in corumModuleList){
  corum_gene=c(corum_gene, corum)
}

# Take mean
library(matrixStats)
library(dplyr)
RNAseq$RNA_Mean = rowMeans(RNAseq[,3:445], na.rm = TRUE)
RNAseq$RNA_Median = rowMedians(as.matrix(RNAseq[,3:445]), na.rm = TRUE)
RNAseq$RNA_STD = rowSds(as.matrix(RNAseq[,3:445]), na.rm = TRUE)
RNAseq$RNA_absVMR = abs(RNAseq$RNA_STD*RNAseq$RNA_STD/RNAseq$RNA_Mean)
RNAseq = RNAseq[,-c(2:445)]
RNAseq = na.omit(RNAseq)
colnames(RNAseq) = gsub('Hugo_Symbol', 'Gene.name', colnames(RNAseq))

RNA_PROT = merge(RNAseq, Protein, by = 'Gene.name')
write.csv(RNA_PROT, "Data/RNA_PROT_sum_new.csv", row.names=FALSE)

# Correlation plots
library(ggplot2)
library(ggpubr)
library(ggrepel)
RNA_PROT_sum = read.csv("Data/RNA_PROT_sum_new.csv")
RNA_PROT_sum$Mean_outlier = "NA"
RNA_PROT_sum$Median_outlier = "NA"
RNA_PROT_sum$CORUM = "No"
for (row in 1:nrow(RNA_PROT_sum)){
  if (RNA_PROT_sum[row, "PROT_Median"] < 10 | RNA_PROT_sum[row, "PROT_Median"] > 30){
    RNA_PROT_sum[row, "Median_outlier"] = as.character(RNA_PROT_sum[row, "Gene.name"])
  }
  if (RNA_PROT_sum[row, "Gene.name"] %in% corum_gene){
    RNA_PROT_sum[row, "CORUM"] = "Yes"
  }
}
RNA_PROT_sum[RNA_PROT_sum=="NA"]=NA
median_cor = cor(RNA_PROT_sum$RNA_Median, RNA_PROT_sum$PROT_Median, method='pearson')

RNA_PROT_sum = RNA_PROT_sum[RNA_PROT_sum$RNA_Median>0,]

pdf(file='Results/MM500/Median_cor_all.pdf',
    width=15,height=10,useDingbats=FALSE)
ggscatter(RNA_PROT_sum, x = "RNA_Median", y = "PROT_Median", 
          add = "reg.line", conf.int = TRUE, 
          xlab = "log2 RNA Median", ylab = "Protein Median", size=0.05, title="Correlation of RNA and Protein Median Expression")+
  stat_cor(method = "pearson", label.x = 12, label.y = 30) + 
  geom_point(color='#21908dff', size=0.05)+ 
  geom_text_repel(label=RNA_PROT_sum$Median_outlier, nudge_x = 0.35, size = 3)+
  theme(plot.title = element_text(hjust=0.5))
dev.off()


RNA_PROT_sum_corum = RNA_PROT_sum[RNA_PROT_sum$CORUM == "Yes",]
RNA_PROT_sum_non_corum= RNA_PROT_sum[RNA_PROT_sum$CORUM == "No",]

median_cor = cor(RNA_PROT_sum_corum$RNA_Median, RNA_PROT_sum_corum$PROT_Median, method='pearson')

pdf(file='Results/MM500/Median_cor_CORUM.pdf',
    width=15,height=10,useDingbats=FALSE)
ggscatter(RNA_PROT_sum_corum, x = "RNA_Median", y = "PROT_Median", 
          add = "reg.line", conf.int = TRUE, 
          xlab = "log2 RNA Median", ylab = "Protein Median", size=0.05, title="Correlation of RNA and Protein Median Expressionin in Complex")+
  stat_cor(method = "pearson", label.x = 12, label.y = 30) + 
  geom_point(color='#21908dff', size=0.05)+ 
  geom_text_repel(label=RNA_PROT_sum_corum$Median_outlier, nudge_x = 0.35, size = 3)+
  theme(plot.title = element_text(hjust=0.5))
dev.off()

median_cor = cor(RNA_PROT_sum_non_corum$RNA_Median, RNA_PROT_sum_non_corum$PROT_Median, method='pearson')

pdf(file='Results/MM500/Median_cor_nonCORUM.pdf',
    width=15,height=10,useDingbats=FALSE)
ggscatter(RNA_PROT_sum_non_corum, x = "RNA_Median", y = "PROT_Median", 
          add = "reg.line", conf.int = TRUE, 
          xlab = "log2 RNA Median", ylab = "Protein Median", size=0.05, title="Correlation of RNA and Protein Median Expressionin not in Complex")+
  stat_cor(method = "pearson", label.x = 12, label.y = 30) + 
  geom_point(color='#21908dff', size=0.05)+ 
  geom_text_repel(label=RNA_PROT_sum_non_corum$Median_outlier, nudge_x = 0.35, size = 3)+
  theme(plot.title = element_text(hjust=0.5))
dev.off()



