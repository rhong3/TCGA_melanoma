# Mitochondria
# Data preparation
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

Prospective <- read_excel("Data/Melanoma proteome/MitochondrialProteinListProspectivePostmortem.xlsx", sheet = "Prospective")
Postmortem <- read_excel("Data/Melanoma proteome/MitochondrialProteinListProspectivePostmortem.xlsx", sheet = "Postmortem")
Prospective = Prospective[-c(1:11), -c(78:84, 86:91)]
Postmortem = Postmortem[-c(1:13), -c(75:81, 83:87)]
Prospective[,1:77] <- lapply(Prospective[,1:77], function(x) as.numeric(gsub("NaN", NA, x)))
Postmortem[,1:74] <- lapply(Postmortem[,1:74], function(x) as.numeric(gsub("NaN", NA, x)))
temp = strsplit(Prospective$PG.ProteinAccession, ';')
for (i in 1:length(temp)){Prospective[i,78] = temp[[i]][1]}
temp = strsplit(Postmortem$`0PG.ProteinAccessions`, ';')
for (i in 1:length(temp)){Postmortem[i,75] = temp[[i]][1]}

annotation <- select(org.Hs.eg.db, keys=unique(Prospective$PG.ProteinAccession), columns=c('UNIPROT', 'SYMBOL'), keytype="UNIPROT")
colnames(Prospective) = gsub('PG.ProteinAccession', "UNIPROT", colnames(Prospective))
Prospective = merge(Prospective, annotation, by="UNIPROT")

annotation <- select(org.Hs.eg.db, keys=unique(Postmortem$`0PG.ProteinAccessions`), columns=c('UNIPROT', 'SYMBOL'), keytype="UNIPROT")
colnames(Postmortem) = gsub('0PG.ProteinAccessions', "UNIPROT", colnames(Postmortem))
Postmortem = merge(Postmortem, annotation, by="UNIPROT")


RNAseq <- read_csv("Data/log2_RNAseq.csv")
RNAseq = na.omit(RNAseq)

# Venn plot of RNA, PROT
library(VennDiagram)
Prospective_unq = unique(na.omit(Prospective$SYMBOL))
RNA_unq = unique(na.omit(RNAseq$Hugo_Symbol))
Postmortem_unq = unique(na.omit(Postmortem$SYMBOL))
venn.diagram(
  x = list(Prospective_unq, Postmortem_unq, RNA_unq),
  category.names = c("Genes in prospective (1684)" , "Genes in postmortem (1800)", "Genes in transcriptomics (20501)"),
  filename = 'Results/Mitochondria/PROT_RNA_venn.png',
  output = TRUE ,
  imagetype="png" ,
  height = 700 , 
  width = 700 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  col=c("#440154ff", '#21908dff', '#e6550d'),
  fill = c(alpha("#440154ff",0.5), alpha('#21908dff',0.5), alpha('#e6550d',0.5)),
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.3,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  cat.col = c("#440154ff", '#21908dff', '#e6550d')
)

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

colnames(Prospective) = gsub('SYMBOL', 'Gene.name', colnames(Prospective))
colnames(Postmortem) = gsub('SYMBOL', 'Gene.name', colnames(Postmortem))
Prospective = Prospective[,-1]
Postmortem = Postmortem[,-1]

Prospective$PROT_Mean = rowMeans(Prospective[,1:77], na.rm = TRUE)
Prospective$PROT_Median = rowMedians(as.matrix(Prospective[,1:77]), na.rm = TRUE)
Prospective$PROT_STD = rowSds(as.matrix(Prospective[,1:77]), na.rm = TRUE)
Prospective$PROT_absVMR = abs(Prospective$PROT_STD*Prospective$PROT_STD/Prospective$PROT_Mean)
Prospective = Prospective[,-c(1:77)]
Prospective = na.omit(Prospective)

Postmortem$PROT_Mean = rowMeans(Postmortem[,1:74], na.rm = TRUE)
Postmortem$PROT_Median = rowMedians(as.matrix(Postmortem[,1:74]), na.rm = TRUE)
Postmortem$PROT_STD = rowSds(as.matrix(Postmortem[,1:74]), na.rm = TRUE)
Postmortem$PROT_absVMR = abs(Postmortem$PROT_STD*Postmortem$PROT_STD/Postmortem$PROT_Mean)
Postmortem = Postmortem[,-c(1:74)]
Postmortem = na.omit(Postmortem)

RNA_Prospective = merge(RNAseq, Prospective, by = 'Gene.name')

RNA_Prospective_sum = RNA_Prospective  %>%
  group_by(`Gene.name`) %>%
  summarize(
    RNA_Mean = mean(RNA_Mean)
    , RNA_Median = mean(RNA_Median)
    , RNA_STD = mean(RNA_STD)
    , RNA_absVMR = mean(RNA_absVMR)
    , PROT_Mean = mean(PROT_Mean)
    , PROT_Median = mean(PROT_Median)
    , PROT_STD = mean(PROT_STD)
    , PROT_absVMR = mean(PROT_absVMR)
  )


RNA_Postmortem = merge(RNAseq, Postmortem, by = 'Gene.name')

RNA_Postmortem_sum = RNA_Postmortem  %>%
  group_by(`Gene.name`) %>%
  summarize(
    RNA_Mean = mean(RNA_Mean)
    , RNA_Median = mean(RNA_Median)
    , RNA_STD = mean(RNA_STD)
    , RNA_absVMR = mean(RNA_absVMR)
    , PROT_Mean = mean(PROT_Mean)
    , PROT_Median = mean(PROT_Median)
    , PROT_STD = mean(PROT_STD)
    , PROT_absVMR = mean(PROT_absVMR)
  )

ENT <- select(org.Hs.eg.db, keys=as.character(unique(RNA_Prospective_sum$Gene.name)), columns=c("ENTREZID", "CHR",  'MAP'), "ALIAS")
colnames(ENT) = gsub('ALIAS', 'Gene.name', colnames(ENT))
RNA_Prospective_sum = merge(ENT, RNA_Prospective_sum, by='Gene.name')
write_csv(RNA_Prospective_sum, 'Results/Mitochondria/RNA_Prospective.csv')

ENT <- select(org.Hs.eg.db, keys=as.character(unique(RNA_Postmortem_sum$Gene.name)), columns=c("ENTREZID", "CHR",  'MAP'), "ALIAS")
colnames(ENT) = gsub('ALIAS', 'Gene.name', colnames(ENT))
RNA_Postmortem_sum = merge(ENT, RNA_Postmortem_sum, by='Gene.name')
write_csv(RNA_Postmortem_sum, 'Results/Mitochondria/RNA_Postmortem.csv')

# Correlation plots
library(ggplot2)
library(ggpubr)

pdf(file='Results/Mitochondria/Mean_cor_prospective.pdf',
    width=15,height=10,useDingbats=FALSE)
ggscatter(RNA_Prospective_sum, x = "RNA_Mean", y = "PROT_Mean", 
          add = "reg.line", conf.int = TRUE, 
          xlab = "log2 RNA Mean", ylab = "Protein Mean", size=0.1, title="Correlation of RNA and Prospective Protein Mean Expression")+
  stat_cor(method = "pearson", label.x = 3, label.y = 30) + 
  geom_point(color='#21908dff', size=0.1)+ 
  
  theme(plot.title = element_text(hjust=0.5))
dev.off()

pdf(file='Results/Mitochondria/Median_cor_prospective.pdf',
    width=15,height=10,useDingbats=FALSE)
ggscatter(RNA_Prospective_sum, x = "RNA_Median", y = "PROT_Median", 
          add = "reg.line", conf.int = TRUE, 
          xlab = "log2 RNA Median", ylab = "Protein Median", size=0.1, title="Correlation of RNA and Prospective Protein Median Expression")+
  stat_cor(method = "pearson", label.x = 3, label.y = 30) + 
  geom_point(color='#21908dff', size=0.1)+ 
  
  theme(plot.title = element_text(hjust=0.5))
dev.off()


pdf(file='Results/Mitochondria/Mean_cor_postmortem.pdf',
    width=15,height=10,useDingbats=FALSE)
ggscatter(RNA_Postmortem_sum, x = "RNA_Mean", y = "PROT_Mean", 
          add = "reg.line", conf.int = TRUE, 
          xlab = "log2 RNA Mean", ylab = "Protein Mean", size=0.1, title="Correlation of RNA and Postmortem Protein Mean Expression")+
  stat_cor(method = "pearson", label.x = 3, label.y = 30) + 
  geom_point(color='#21908dff', size=0.1)+ 
  
  theme(plot.title = element_text(hjust=0.5))
dev.off()

pdf(file='Results/Mitochondria/Median_cor_postmortem.pdf',
    width=15,height=10,useDingbats=FALSE)
ggscatter(RNA_Postmortem_sum, x = "RNA_Median", y = "PROT_Median", 
          add = "reg.line", conf.int = TRUE, 
          xlab = "log2 RNA Median", ylab = "Protein Median", size=0.1, title="Correlation of RNA and Postmortem Protein Median Expression")+
  stat_cor(method = "pearson", label.x = 3, label.y = 30) + 
  geom_point(color='#21908dff', size=0.1)+ 
  
  theme(plot.title = element_text(hjust=0.5))
dev.off()


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

RNA_clinical_prospective = RNA_clinical[RNA_clinical$Gene.name %in% unique(RNA_Prospective_sum$Gene.name), ]
row.names(RNA_clinical_prospective) = RNA_clinical_prospective$Gene.name
RNA_clinical_prospective = RNA_clinical_prospective[,-1]
RNA_clinical_prospective <- sapply(RNA_clinical_prospective, as.numeric)
RNA_clinical_postmortem = RNA_clinical[RNA_clinical$Gene.name %in% unique(RNA_Postmortem_sum$Gene.name), ]
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
pdf("Results/Mitochondria/HM_RNA_clinical_prospective.pdf", height = 20, width = 30)
hp = Heatmap(as.matrix(RNA_clinical_prospective), col = col, column_title = paste("RNAseq log2-transformed Median Expression"), top_annotation = anno,
             cluster_rows = TRUE, cluster_columns = TRUE, show_row_names = FALSE, show_column_names = FALSE, name = "value", heatmap_legend_param = list(direction = "vertical"))
draw(hp, heatmap_legend_side = "bottom", 
     annotation_legend_side = "right", merge_legend = TRUE,)
dev.off()

breaksList = seq(min(as.matrix(RNA_clinical_postmortem)), max(as.matrix(RNA_clinical_postmortem)), by=2)
col = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(length(breaksList))
pdf("Results/Mitochondria/HM_RNA_clinical_postmortem.pdf", height = 20, width = 30)
hp = Heatmap(as.matrix(RNA_clinical_postmortem), col = col, column_title = paste("RNAseq log2-transformed Median Expression"), top_annotation = anno,
             cluster_rows = TRUE, cluster_columns = TRUE, show_row_names = FALSE, show_column_names = FALSE, name = "value", heatmap_legend_param = list(direction = "vertical"))
draw(hp, heatmap_legend_side = "bottom", 
     annotation_legend_side = "right", merge_legend = TRUE,)
dev.off()


