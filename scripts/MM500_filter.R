# filter log2RNA < 5
# Correlation plots
library(ggplot2)
library(ggpubr)
library(readxl)
library(biomaRt)
library(org.Hs.eg.db)
library(stringr)
library(readr)
library(plyr)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(robustHD)
# CORUM protein complex
library(MoonFinder)
data("corumModuleList")
corum_gene=c()
for (corum in corumModuleList){
  corum_gene=c(corum_gene, corum)
}

RNA_PROT_sum = read.csv("Data/RNA_PROT_sum.csv")
RNA_PROT_sum$Mean_outlier = "NA"
RNA_PROT_sum$Median_outlier = "NA"
RNA_PROT_sum$CORUM = "No"
for (row in 1:nrow(RNA_PROT_sum)){
  if (RNA_PROT_sum[row, "PROT_Mean"] < -10 | RNA_PROT_sum[row, "PROT_Mean"] > 10){
    RNA_PROT_sum[row, "Mean_outlier"] = RNA_PROT_sum[row, "Gene.name"]
  }
  if (RNA_PROT_sum[row, "PROT_Median"] < -10 | RNA_PROT_sum[row, "PROT_Mean"] > 10){
    RNA_PROT_sum[row, "Median_outlier"] = RNA_PROT_sum[row, "Gene.name"]
  }
  if (RNA_PROT_sum[row, "Gene.name"] %in% corum_gene){
    RNA_PROT_sum[row, "CORUM"] = "Yes"
  }
}
RNA_PROT_sum[RNA_PROT_sum=="NA"]=NA
mean_cor = cor(RNA_PROT_sum$RNA_Mean, RNA_PROT_sum$PROT_Mean, method='pearson')
median_cor = cor(RNA_PROT_sum$RNA_Median, RNA_PROT_sum$PROT_Median, method='pearson')

pdf(file='Results/MM500/Mean_cor_all_filter.pdf',
    width=15,height=10,useDingbats=FALSE)
ggscatter(RNA_PROT_sum[RNA_PROT_sum$RNA_Mean > 5, ], x = "RNA_Mean", y = "PROT_Mean", 
          add = "reg.line", conf.int = TRUE, 
          xlab = "log2 RNA Mean", ylab = "Protein Mean", size=0.05, title="Correlation of RNA and Protein Mean Expression")+
  stat_cor(method = "pearson", label.x = 3, label.y = 15) + 
  geom_point(color='#21908dff', size=0.05)+ 
  geom_text(label=RNA_PROT_sum[RNA_PROT_sum$RNA_Mean > 5, ]$Mean_outlier, nudge_x = 0.35, size = 3)+
  theme(plot.title = element_text(hjust=0.5))
dev.off()

pdf(file='Results/MM500/Median_cor_all_filter.pdf',
    width=15,height=10,useDingbats=FALSE)
ggscatter(RNA_PROT_sum[RNA_PROT_sum$RNA_Median > 5, ], x = "RNA_Median", y = "PROT_Median", 
          add = "reg.line", conf.int = TRUE, 
          xlab = "log2 RNA Median", ylab = "Protein Median", size=0.05, title="Correlation of RNA and Protein Median Expression")+
  stat_cor(method = "pearson", label.x = 3, label.y = 15) + 
  geom_point(color='#21908dff', size=0.05)+ 
  geom_text(label=RNA_PROT_sum[RNA_PROT_sum$RNA_Median > 5, ]$Median_outlier, nudge_x = 0.35, size = 3)+
  theme(plot.title = element_text(hjust=0.5))
dev.off()


RNA_PROT_sum_corum = RNA_PROT_sum[RNA_PROT_sum$CORUM == "Yes",]
RNA_PROT_sum_non_corum= RNA_PROT_sum[RNA_PROT_sum$CORUM == "No",]

mean_cor = cor(RNA_PROT_sum_corum$RNA_Mean, RNA_PROT_sum_corum$PROT_Mean, method='pearson')
median_cor = cor(RNA_PROT_sum_corum$RNA_Median, RNA_PROT_sum_corum$PROT_Median, method='pearson')

pdf(file='Results/MM500/Mean_cor_CORUM_filter.pdf',
    width=15,height=10,useDingbats=FALSE)
ggscatter(RNA_PROT_sum_corum[RNA_PROT_sum_corum$RNA_Mean > 5, ], x = "RNA_Mean", y = "PROT_Mean", 
          add = "reg.line", conf.int = TRUE, 
          xlab = "log2 RNA Mean", ylab = "Protein Mean", size=0.05, title="Correlation of RNA and Protein Mean Expression in Complex")+
  stat_cor(method = "pearson", label.x = 3, label.y = 10) + 
  geom_point(color='#21908dff', size=0.05)+ 
  geom_text(label=RNA_PROT_sum_corum[RNA_PROT_sum_corum$RNA_Mean > 5, ]$Mean_outlier, nudge_x = 0.35, size = 3)+
  theme(plot.title = element_text(hjust=0.5))
dev.off()

pdf(file='Results/MM500/Median_cor_CORUM_filter.pdf',
    width=15,height=10,useDingbats=FALSE)
ggscatter(RNA_PROT_sum_corum[RNA_PROT_sum_corum$RNA_Median > 5, ], x = "RNA_Median", y = "PROT_Median", 
          add = "reg.line", conf.int = TRUE, 
          xlab = "log2 RNA Median", ylab = "Protein Median", size=0.05, title="Correlation of RNA and Protein Median Expressionin in Complex")+
  stat_cor(method = "pearson", label.x = 3, label.y = 10) + 
  geom_point(color='#21908dff', size=0.05)+ 
  geom_text(label=RNA_PROT_sum_corum[RNA_PROT_sum_corum$RNA_Median > 5, ]$Median_outlier, nudge_x = 0.35, size = 3)+
  theme(plot.title = element_text(hjust=0.5))
dev.off()


mean_cor = cor(RNA_PROT_sum_non_corum$RNA_Mean, RNA_PROT_sum_non_corum$PROT_Mean, method='pearson')
median_cor = cor(RNA_PROT_sum_non_corum$RNA_Median, RNA_PROT_sum_non_corum$PROT_Median, method='pearson')

pdf(file='Results/MM500/Mean_cor_nonCORUM_filter.pdf',
    width=15,height=10,useDingbats=FALSE)
ggscatter(RNA_PROT_sum_non_corum[RNA_PROT_sum_non_corum$RNA_Mean > 5, ], x = "RNA_Mean", y = "PROT_Mean", 
          add = "reg.line", conf.int = TRUE, 
          xlab = "log2 RNA Mean", ylab = "Protein Mean", size=0.05, title="Correlation of RNA and Protein Mean Expression not in Complex")+
  stat_cor(method = "pearson", label.x = 3, label.y = 15) + 
  geom_point(color='#21908dff', size=0.05)+ 
  geom_text(label=RNA_PROT_sum_non_corum[RNA_PROT_sum_non_corum$RNA_Mean > 5, ]$Mean_outlier, nudge_x = 0.35, size = 3)+
  theme(plot.title = element_text(hjust=0.5))
dev.off()

pdf(file='Results/MM500/Median_cor_nonCORUM_filter.pdf',
    width=15,height=10,useDingbats=FALSE)
ggscatter(RNA_PROT_sum_non_corum[RNA_PROT_sum_non_corum$RNA_Median > 5, ], x = "RNA_Median", y = "PROT_Median", 
          add = "reg.line", conf.int = TRUE, 
          xlab = "log2 RNA Median", ylab = "Protein Median", size=0.05, title="Correlation of RNA and Protein Median Expressionin not in Complex")+
  stat_cor(method = "pearson", label.x = 3, label.y = 15) + 
  geom_point(color='#21908dff', size=0.05)+ 
  geom_text(label=RNA_PROT_sum_non_corum[RNA_PROT_sum_non_corum$RNA_Median > 5, ]$Median_outlier, nudge_x = 0.35, size = 3)+
  theme(plot.title = element_text(hjust=0.5))
dev.off()

# VMR corrrelations
pdf(file='Results/MM500/VMR_cor.pdf',
    width=15,height=10,useDingbats=FALSE)
ggscatter(RNA_PROT_sum[RNA_PROT_sum$PROT_absVMR < 100, ], x = "RNA_absVMR", y = "PROT_absVMR", 
          add = "reg.line", conf.int = TRUE, 
          xlab = "log2 RNA abs of Variance to Mean Ratio", ylab = "Protein abs of Variance to Mean Ratio", size=0.05, title="Correlation of RNA and Protein Expression Variance to Mean Ratio")+
  stat_cor(method = "pearson", label.x = 3, label.y = 15) + 
  geom_point(color='#21908dff', size=0.05)+ 
  theme(plot.title = element_text(hjust=0.5))
dev.off()

# Top 500 VMR
RNA_PROT_sum = RNA_PROT_sum[RNA_PROT_sum$RNA_Mean > 5, ]
RNA_PROT_sum = RNA_PROT_sum[order(-RNA_PROT_sum['RNA_absVMR']), ]
colnames(RNA_PROT_sum) = gsub('Gene name', 'Gene.name', colnames(RNA_PROT_sum))

RNA_clinical = read.csv("Data/clinical_RNAseq.csv")
colnames(RNA_clinical)[1] = 'ENTREZID'
clinical = RNA_clinical[1:8,]
RNA_clinical = RNA_clinical[-c(1:8), ]
annotation <- select(org.Hs.eg.db, keys=as.character(unique(RNA_clinical$ENTREZID)), columns=c('ENTREZID', 'SYMBOL'), keytype="ENTREZID")
RNA_clinical = merge(RNA_clinical, annotation, by='ENTREZID')
RNA_clinical = RNA_clinical[,-c(1,445)]
RNA_clinical = cbind(RNA_clinical[,444], RNA_clinical[,1:443])
colnames(RNA_clinical)[1] = 'Gene.name'
colnames(clinical)[1]  = 'Gene.name'
clinical = clinical[,-445]

RNA_PROT_sum_clinical = RNA_clinical[RNA_clinical$Gene.name %in% RNA_PROT_sum$Gene.name[1:500], ]
row.names(RNA_PROT_sum_clinical) = RNA_PROT_sum_clinical$Gene.name
RNA_PROT_sum_clinical = RNA_PROT_sum_clinical[,-1]
RNA_PROT_sum_clinical <- sapply(RNA_PROT_sum_clinical, as.numeric)


anno = HeatmapAnnotation(age = as.numeric(as.matrix(clinical[1,2:444])), 
                         gender = as.character(as.matrix(clinical[2,2:444])), 
                         AJCC_stage = as.character(as.matrix(clinical[3,2:444])),
                         weight = as.numeric(as.matrix(clinical[4,2:444])), 
                         survival_months = as.numeric(as.matrix(clinical[5,2:444])), 
                         cancer_type = as.character(as.matrix(clinical[6,2:444])),
                         aneuploidy_score = as.numeric(as.matrix(clinical[7,2:444])), 
                         sample_type = as.character(as.matrix(clinical[8,2:444])) )
# col = list(days = col_fun2, stage = col_fun3, NRAS = col_fun4))

breaksList = seq(min(as.matrix(RNA_PROT_sum_clinical)), max(as.matrix(RNA_PROT_sum_clinical)), by=2)
col = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(length(breaksList))
pdf("Results/MM500/HM_RNA_500.pdf", height = 20, width = 30)
hp = Heatmap(as.matrix(RNA_PROT_sum_clinical), col = col, column_title = paste("RNAseq log2-transformed Median Expression"), top_annotation = anno,
             cluster_rows = TRUE, cluster_columns = TRUE, show_row_names = FALSE, show_column_names = FALSE, name = "value", heatmap_legend_param = list(direction = "vertical"))
draw(hp, heatmap_legend_side = "bottom", 
     annotation_legend_side = "right", merge_legend = TRUE,)
dev.off()




