# filter log2RNA < 5
RNA_Prospective_sum = read.csv('Results/Mitochondria/RNA_Prospective.csv')
RNA_Postmortem_sum = read.csv('Results/Mitochondria/RNA_Postmortem.csv')

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

pdf(file='Results/Mitochondria/Mean_cor_prospective_filter.pdf',
    width=15,height=10,useDingbats=FALSE)
ggscatter(RNA_Prospective_sum[RNA_Prospective_sum$RNA_Mean > 5, ], x = "RNA_Mean", y = "PROT_Mean", 
          add = "reg.line", conf.int = TRUE, 
          xlab = "log2 RNA Mean", ylab = "Protein Mean", size=0.1, title="Correlation of RNA and Prospective Protein Mean Expression")+
  stat_cor(method = "pearson", label.x = 3, label.y = 30) + 
  geom_point(color='#21908dff', size=0.1)+ 
  
  theme(plot.title = element_text(hjust=0.5))
dev.off()

pdf(file='Results/Mitochondria/Median_cor_prospective_filter.pdf',
    width=15,height=10,useDingbats=FALSE)
ggscatter(RNA_Prospective_sum[RNA_Prospective_sum$RNA_Median > 5, ], x = "RNA_Median", y = "PROT_Median", 
          add = "reg.line", conf.int = TRUE, 
          xlab = "log2 RNA Median", ylab = "Protein Median", size=0.1, title="Correlation of RNA and Prospective Protein Median Expression")+
  stat_cor(method = "pearson", label.x = 3, label.y = 30) + 
  geom_point(color='#21908dff', size=0.1)+ 
  
  theme(plot.title = element_text(hjust=0.5))
dev.off()


pdf(file='Results/Mitochondria/Mean_cor_postmortem_filter.pdf',
    width=15,height=10,useDingbats=FALSE)
ggscatter(RNA_Postmortem_sum[RNA_Postmortem_sum$RNA_Mean > 5, ], x = "RNA_Mean", y = "PROT_Mean", 
          add = "reg.line", conf.int = TRUE, 
          xlab = "log2 RNA Mean", ylab = "Protein Mean", size=0.1, title="Correlation of RNA and Postmortem Protein Mean Expression")+
  stat_cor(method = "pearson", label.x = 3, label.y = 30) + 
  geom_point(color='#21908dff', size=0.1)+ 
  
  theme(plot.title = element_text(hjust=0.5))
dev.off()

pdf(file='Results/Mitochondria/Median_cor_postmortem_filter.pdf',
    width=15,height=10,useDingbats=FALSE)
ggscatter(RNA_Postmortem_sum[RNA_Postmortem_sum$RNA_Median > 5, ], x = "RNA_Median", y = "PROT_Median", 
          add = "reg.line", conf.int = TRUE, 
          xlab = "log2 RNA Median", ylab = "Protein Median", size=0.1, title="Correlation of RNA and Postmortem Protein Median Expression")+
  stat_cor(method = "pearson", label.x = 3, label.y = 30) + 
  geom_point(color='#21908dff', size=0.1)+ 
  
  theme(plot.title = element_text(hjust=0.5))
dev.off()


# Top 500 VMR
RNA_Prospective_sum = RNA_Prospective_sum[RNA_Prospective_sum$RNA_Mean > 5, ]
RNA_Postmortem_sum = RNA_Postmortem_sum[RNA_Postmortem_sum$RNA_Mean > 5, ]
RNA_Prospective_sum = RNA_Prospective_sum[order(-RNA_Prospective_sum['RNA_absVMR']), ]
RNA_Postmortem_sum = RNA_Postmortem_sum[order(-RNA_Postmortem_sum['RNA_absVMR']), ]

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

RNA_clinical_prospective = RNA_clinical[RNA_clinical$Gene.name %in% RNA_Prospective_sum$Gene.name[1:500], ]
row.names(RNA_clinical_prospective) = RNA_clinical_prospective$Gene.name
RNA_clinical_prospective = RNA_clinical_prospective[,-1]
RNA_clinical_prospective <- sapply(RNA_clinical_prospective, as.numeric)
RNA_clinical_postmortem = RNA_clinical[RNA_clinical$Gene.name %in% RNA_Postmortem_sum$Gene.name[1:500], ]
row.names(RNA_clinical_postmortem) = RNA_clinical_postmortem$Gene.name
RNA_clinical_postmortem = RNA_clinical_postmortem[,-1]
RNA_clinical_postmortem <- sapply(RNA_clinical_postmortem, as.numeric)


anno = HeatmapAnnotation(age = as.numeric(as.matrix(clinical[1,2:444])), 
                         gender = as.character(as.matrix(clinical[2,2:444])), 
                         AJCC_stage = as.character(as.matrix(clinical[3,2:444])),
                         weight = as.numeric(as.matrix(clinical[4,2:444])), 
                         survival_months = as.numeric(as.matrix(clinical[5,2:444])), 
                         cancer_type = as.character(as.matrix(clinical[6,2:444])),
                         aneuploidy_score = as.numeric(as.matrix(clinical[7,2:444])), 
                         sample_type = as.character(as.matrix(clinical[8,2:444])) )
# col = list(days = col_fun2, stage = col_fun3, NRAS = col_fun4))

breaksList = seq(min(as.matrix(RNA_clinical_prospective)), max(as.matrix(RNA_clinical_prospective)), by=2)
col = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(length(breaksList))
pdf("Results/Mitochondria/HM_RNA_clinical_prospective_500.pdf", height = 20, width = 30)
hp = Heatmap(as.matrix(RNA_clinical_prospective), col = col, column_title = paste("RNAseq log2-transformed Median Expression"), top_annotation = anno,
             cluster_rows = TRUE, cluster_columns = TRUE, show_row_names = FALSE, show_column_names = FALSE, name = "value", heatmap_legend_param = list(direction = "vertical"))
draw(hp, heatmap_legend_side = "bottom", 
     annotation_legend_side = "right", merge_legend = TRUE,)
dev.off()

breaksList = seq(min(as.matrix(RNA_clinical_postmortem)), max(as.matrix(RNA_clinical_postmortem)), by=2)
col = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(length(breaksList))
pdf("Results/Mitochondria/HM_RNA_clinical_postmortem_500.pdf", height = 20, width = 30)
hp = Heatmap(as.matrix(RNA_clinical_postmortem), col = col, column_title = paste("RNAseq log2-transformed Median Expression"), top_annotation = anno,
             cluster_rows = TRUE, cluster_columns = TRUE, show_row_names = FALSE, show_column_names = FALSE, name = "value", heatmap_legend_param = list(direction = "vertical"))
draw(hp, heatmap_legend_side = "bottom", 
     annotation_legend_side = "right", merge_legend = TRUE,)
dev.off()

